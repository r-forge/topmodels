hxlr <- function(formula, data, subset = NULL, na.action = NULL, weights, 
  thresholds, control = hxlr.control(...), ...)
{
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  

  ## extract terms, response, model matrix
  mt <- terms(formula, data = data)
  mtX <- delete.response(terms(formula, data = data, rhs = 1L))
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf)
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
  
  


  ## obtain correct subset of predvars/dataClasses to terms
  .add_predvars_and_dataClasses <- function(terms, model.frame) {
    ## original terms
    rval <- terms
    ## terms from model.frame
    nval <- if(inherits(model.frame, "terms")) model.frame else terms(model.frame)

    ## associated variable labels
    ovar <- sapply(as.list(attr(rval, "variables")), deparse)[-1]
    nvar <- sapply(as.list(attr(nval, "variables")), deparse)[-1]
    if(!all(ovar %in% nvar)) stop(paste("The following terms variables are not part of the model.frame:",
      paste(ovar[!(ovar %in% nvar)], collapse = ", ")))
    ix <- match(ovar, nvar)
  
    ## subset predvars
    if(!is.null(attr(rval, "predvars"))) warning("terms already had 'predvars' attribute, now replaced")
    attr(rval, "predvars") <- attr(nval, "predvars")[1L + c(0L, ix)]

    ## subset dataClasses
    if(!is.null(attr(rval, "dataClasses"))) warning("terms already had 'dataClasses' attribute, now replaced")
    attr(rval, "dataClasses") <- attr(nval, "dataClasses")[ix]
  
    return(rval)
  }
  mt  <- .add_predvars_and_dataClasses(mt,  mf)
  mtX <- .add_predvars_and_dataClasses(mtX, mf)
  mtZ <- .add_predvars_and_dataClasses(mtZ, mf)


  ## factorize response if numeric 
  if(!is.factor(Y)) Y <- cut(Y, c(-Inf, thresholds, Inf))
  ## new formula with factorized response
  mformula <- reformulate(attr(mtX, "term.labels"), response = "Y")
  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(any(summary(Y)==0)) {
    warning("Intervals with no data! Corresponding threshold ignored")
    thresholds <- thresholds[-which(summary(data$ff.cat)==0)]
    Y <- droplevels(Y)
  }
    


  ## convenience variables
  n <- length(Y)

  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  ## weights must not be 0 for clm
  data <- data[weights > 0,]
  Y <- Y[weights > 0]
  weights <- weights[weights > 0]

  ## prepare data
  data <- get_all_vars(mt, data)
  if(!is.null(subset)) data <- data[subset,]
  if(!is.null(na.action)) data <- na.action(data)
  else data <- na.omit(data)

  stopifnot(requireNamespace("ordinal"))
  ## get environment from clm
  env <- ordinal::clm(formula = mformula, scale = mtZ, data = data, weights = weights, doFit = FALSE, ...)
  
  ## thresholds can also be data.frame with several columns (predictor variables for intercept model)
  q <- model.matrix(~ thresholds)
  p <- ncol(q)  # number of columns (predictor variables)

  ## new objective function
  ## intercept are replaced with a[1] + a[2] * thresholds
  nll2 <- function(par, envir) {
    a <- par[1:p] ## parameters for the thresholds
    ## resulting intercepts:
    theta <- q %*% a ## = a[1] + a[2] * thresholds
    ## set parameters in envir:
    envir$par <- c(theta, par[-(1:p)])
    ## evaluate negative log-likelihood:
    envir$clm.nll(envir)
  }
    

  ## control parameters
  ocontrol <- control
  method <- control$method
  hessian <- control$hessian
  start <- control$start
  control$method <- control$hessian <- control$start <- NULL

  ## starting values
  if(is.null(start)) {
    ## starting values from clm fit, threshold coefficients set to 0 and 1
    strt <- ordinal::clm(formula = mformula, scale = mtZ, data = data, weights = weights, ...)
    strt <- c(0, rep(1, p-1)/(p-1), strt$beta, strt$zeta)
  }
  if(is.list(start)) start <- do.call("c", start) 
  if(length(start) > p + length(attr(mt, "term.labels"))) {
    warning(paste("too many entries in start! only first", length(attr(mt, "term.labels")) + p, "entries are considered"))
    start <- start[1: (length(attr(mt, "term.labels")) + p)]
  }

  ## estimation
  opt <- optim(par = strt, fn = nll2, envir = env, method = method, hessian = hessian, control = control)

  if(opt$convergence > 0) {
    converged <- FALSE
    warning("optimization failed to converge")
  } else {
    converged <- TRUE
  }

  ## compute Hessian
  vcov <- if (hessian) solve(as.matrix(opt$hessian)) else NULL

  if (hessian) {
    colnames(vcov) <- rownames(vcov) <- c(
      colnames(q), tail(colnames(X),-1),
      tail(colnames(Z),-1)
      )
  }
  
  ## coefficients
  par <- opt$par
  nbeta <- length(par) - p - length(attr(mtZ, "term.labels"))
  ## intercept coefficients
  alpha <- par[1:p]
  names(alpha) <- colnames(q)
  ## logit coefficients
  beta <- par[seq.int(length.out = nbeta) + p]
  ## scale coefficients
  delta2 <- tail(par, length(attr(mtZ, "term.labels")))



  ## coefficients for location and scale (gamma, delta)
  if(p ==2) {
    gamma <- c(-alpha[1], beta)/alpha[2]
    delta <- c(-log(alpha[2]), delta2)
    names(delta)[1] <- "(Intercept)"
  } else {gamma <- delta <- NULL}

  ## output
  rval <- list(
    coefficients = list(location = gamma, scale = delta),
    coefficients.CLM = list(intercept = alpha, location = beta, scale = delta2),
    fitted.values = list(location = drop(X %*% gamma),  scale = exp(drop(Z %*% delta))),
    optim = opt, 
    method = method,
    control = control,
    start = strt,
    weights = weights, 
    nobs = n,
    loglik = - opt$value,
    vcov = vcov,
    converged = converged, 
    iterations = as.vector(tail(na.omit(opt$counts), 1)),
    call = cl,
    formula = formula,
    terms = list(location = mtX, scale = mtZ, full = mt),
    levels = list(location = .getXlevels(mtX, mf), scale = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf)),
    thresholds = thresholds
  )
  
  class(rval) <- "hxlr"
  return(rval)
}

hxlr.control <- function(method = "BFGS", maxit = 5000, hessian = TRUE, trace = FALSE, start = NULL, ...)
{
  rval <- list(method = method, maxit = maxit, hessian = hessian, trace = trace, start = start)
  rval <- c(rval, list(...))
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  rval$fnscale <- 1
  if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.2)
  rval
}


logLik.hxlr <- function(object, ...)
  structure(object$loglik, df = sum(sapply(object$coefficients, length)), class = "logLik")

predict.hxlr <- function(object, newdata = NULL,
  type = c("class", "probability", "cumprob", "location", "scale"),
  thresholds = object$thresholds, na.action = na.pass, ...)
{
  type <- match.arg(type)
  
  if(missing(newdata)) {
    mu <- object$fitted.values$location
    sigma <- object$fitted.values$scale
  } else {
    ## get coefficients
    gamma <- object$coefficients$location
    delta <- object$coefficients$scale

    mf <- model.frame(delete.response(object$terms$full), newdata, na.action = na.action, xlev = object$levels$full)

    ## matrices for location and scale model from newdata
    x <- model.matrix(delete.response(object$terms$location), mf)
    z <- model.matrix(object$terms$scale, mf)
    
    ## location and scale 
    mu <- drop(x %*% gamma)
    sigma <- exp(drop(z %*% delta))
  }

  ## cumulative probabilities P(y<threshold[i]|x)
  cumprob <- NULL
  for(i in 1:length(thresholds)) { 
    cumprob <- cbind(cumprob, plogis((thresholds[i] - mu)/sigma))
  }
  
  ## category probabilities P(threshold[i-1]<=y<threshold[i]|x)
  prob <- t(apply(cbind(0, cumprob, 1), 1, diff))

  rval <- switch(type,
      "location" = mu,
      "scale" = sigma,
      "cumprob" = cumprob,
      "probability" = prob,
      "class" = apply(prob, 1, which.max)
    )

  return(rval)
}

print.hxlr <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if(length(x$coefficients.CLM$location) | length(x$coefficients.CLM$intercept)) {
      cat(paste("Coefficients (location model):\n", sep = ""))
       print.default(format(c(x$coefficients.CLM$intercept, x$coefficients.CLM$location), digits = digits), print.gap = 2, quote = FALSE)
       cat("\n")
    } else cat("No coefficients\n\n")
    if(length(x$coefficients.CLM$scale)) {
      cat(paste("Coefficients (scale model with log link):\n", sep = ""))
      print.default(format(x$coefficients.CLM$scale, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in scale model)\n\n")
    cat("---\n")
    cat("converted coefficients for location and scale of latent distribution:\n")
    if(length(x$coefficients$location)) {
      cat(paste("Coefficients (location model):\n", sep = ""))
      print.default(format(x$coefficients$location, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in location model)\n\n")
    if(length(x$coefficients$scale)) {
      cat(paste("Coefficients (scale model with log link):\n", sep = ""))
      print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in scale model)\n\n")
    cat("\n")
  }
  invisible(x)
}

summary.hxlr <- function(object, ...)
{
  ## extend coefficient table
  k <- length(object$coefficients.CLM$intercept)
  l <- length(object$coefficients.CLM$location)
  m <- length(object$coefficients.CLM$scale)
  cf <- as.vector(do.call("c", object$coefficients.CLM))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if(length(object$coefficients.CLM$scale)) {
    cf <- list(
      intercept = cf[seq.int(length.out = k), , drop = FALSE],
      location = cf[seq.int(length.out = l) + k, , drop = FALSE],
      scale = cf[seq.int(length.out = m) + l + k, , drop = FALSE])
    rownames(cf$scale) <- names(object$coefficients.CLM$scale)
  } else {
    cf <- list(intercept = cf[seq.int(length.out = k), , drop = FALSE], location = cf[seq.int(length.out = l) + k, , drop = FALSE])
  }
  rownames(cf$intercept) <- names(object$coefficients.CLM$intercept)
  rownames(cf$location) <- names(object$coefficients.CLM$location)
  object$coefficients.CLM <- cf

  ## delete some slots
  object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL

  ## return
  class(object) <- "summary.hxlr"
  object
}

print.summary.hxlr <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if(NROW(x$coefficients.CLM$location) | NROW(x$coefficients.CLM$intercept)) {
      cat(paste("\nCoefficients:\n", sep = ""))
      printCoefmat(rbind(x$coefficients.CLM$intercept, x$coefficients.CLM$location), digits = digits, signif.legend = FALSE)
    } 

    if(NROW(x$coefficients.CLM$scale)) {
      cat(paste("\nlog-scale coefficients:\n", sep = ""))
      printCoefmat(x$coefficients.CLM$scale, digits = digits, signif.legend = FALSE)
    } 

    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients.CLM)[, 4L] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    cat("Log-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df\n")
    cat(paste("Number of iterations in", x$method, "optimization:", x$iterations[1L], "\n"))
  }
  invisible(x)
}

terms.hxlr <- function(x, model = c("location", "scale", "full"), ...) x$terms[[match.arg(model)]]

fitted.hxlr <- function(object, type = c("location", "scale"), ...) object$fitted.values[[match.arg(type)]]

coef.hxlr <- function(object, model = c("full", "intercept", "location", "scale"), type = c("CLM", "latent"), ...) {
  type<- match.arg(type)
  model <- match.arg(model)
  cf <- switch(type, 
    "CLM" = {
      object$coefficients.CLM
    },
    "latent" = {
      object$coefficients
    }
  )
  switch(model,
    "location" = {
      c(cf$location, cf$intercept)
    },
    "scale" = {
      cf$scale
    },
    "full" = {
      cf <- c(cf$intercept, cf$location, cf$scale)
      names(cf) <- colnames(object$vcov)
      cf
    }
  )
}


vcov.hxlr <- function(object, model = c("full", "intercept", "location", "scale"), ...) {
  vc <- object$vcov
  k <- length(object$coefficients.CLM$intercept)
  l <- length(object$coefficients.CLM$location)
  m <- length(object$coefficients.CLM$scale)

  model <-  match.arg(model)

  switch(model,
    "intercept" = {
      vc[seq.int(length.out = k) , seq.int(length.out = k), drop = FALSE]
    },
    "location" = {
      vc[seq.int(length.out = l) + k, seq.int(length.out = l) + k, drop = FALSE]
    },
    "scale" = {
      vc <- vc[seq.int(length.out = m) + k + l, seq.int(length.out = m) + k + l, drop = FALSE]
      colnames(vc) <- rownames(vc) <- names(object$coefficients.CLM$scale)
      vc
    },
    "full" = {
      vc
    }
  )
}


