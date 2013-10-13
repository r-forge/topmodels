crch <- function(formula, data, subset, na.action, weights, offset,
  dist = c("student", "gaussian", "logistic"), df = NULL,
  left = -Inf, right = Inf, control = crch.control(...),
  model = TRUE, x = FALSE, y = FALSE, ...)
{
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
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

  ## distribution
  dist <- match.arg(dist)

  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(dist == "student") {
    if(!is.null(df) && df <= 0) stop("'df' must be positive")
    if(!is.null(df) && !is.finite(df)) dist <- "gaussian"
  }

  ## convenience variables
  n <- length(Y)


  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)

  ## offsets
  expand_offset <- function(offset) {
    if(is.null(offset)) offset <- 0
    if(length(offset) == 1) offset <- rep.int(offset, n)
    as.vector(offset)
  }
  ## in location part of formula
  offsetX <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 1L, terms = TRUE)))
  ## in scale part of formula
  offsetZ <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 2L, terms = TRUE)))
  ## in offset argument (used for location)
  if(!is.null(cl$offset)) offsetX <- offsetX + expand_offset(mf[, "(offset)"])
  ## collect
  offset <- list(location = offsetX, scale = offsetZ)
 

  ## call the actual workhorse: crch.fit() 
  rval <- crch.fit(x = X, y = Y, z = Z, left = left, right = right, dist = dist, df = df, weights = weights, offset = offset, control = control)

  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(location = mtX, scale = mtZ, full = mt)
  rval$levels <- list(location = .getXlevels(mtX, mf), scale = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(location = attr(X, "contrasts"), scale = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(location = X, scale = Z)

  class(rval) <- "crch"
  return(rval)
}


crch.control <- function(method = "BFGS", maxit = 5000, hessian = TRUE, trace = FALSE, start = NULL, ...)
{
  rval <- list(method = method, maxit = maxit, hessian = hessian, trace = trace, start = start)
  rval <- c(rval, list(...))
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  rval$fnscale <- 1
  if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.2)
  rval
}


crch.fit <- function(x, z, y, left, right, dist = "student", df = NULL,  
  weights = NULL, offset = NULL, control = crch.control()) 
{
  ## response and regressor matrix
  n <- NROW(x)  
  k <- NCOL(x)
  if(is.null(weights)) weights <- rep.int(1, n)
  nobs <- sum(weights > 0)
  dfest <- dist == "student" & is.null(df)
  if(is.null(offset)) offset <- rep.int(0, n)
  if(!is.list(offset)) offset <- list(location = offset, scale = rep.int(0, n))
  if(is.null(z)) {
    q <- 1L
    z <- matrix(1, ncol = q, nrow = n)
    colnames(z) <- "(Intercept)"
    rownames(z) <- rownames(x)
  } else {
    q <- NCOL(z)
    if(q < 1L) stop("scale regression needs to have at least one parameter")
  }

  dist <- match.arg(dist, c("student", "gaussian", "logistic"))
  ddist <- switch(dist,
    "student"  = function(x, location, scale, df) dt((x - location)/scale, df = df, log = TRUE) - log(scale),
    "gaussian" = function(x, location, scale, df) dnorm((x - location)/scale, log = TRUE) - log(scale),
    "logistic" = function(x, location, scale, df) dlogis((x - location)/scale, log = TRUE) - log(scale)
  )
  pdist <- switch(dist,
    "student"  = function(x, location, scale, df, lower.tail = TRUE) pt((x - location)/scale, df = df, lower.tail = lower.tail, log.p = TRUE),
    "gaussian" = function(x, location, scale, df, lower.tail = TRUE) pnorm((x - location)/scale, lower.tail = lower.tail, log.p = TRUE),
    "logistic" = function(x, location, scale, df, lower.tail = TRUE) plogis((x - location)/scale, lower.tail = lower.tail, log.p = TRUE)
  )

  ## control parameters
  ocontrol <- control
  method <- control$method
  hessian <- control$hessian
  start <- control$start
  control$method <- control$hessian <- control$start <- NULL

  ## starting values
  if(is.null(start)) {
    auxreg <- lm.wfit(x, y, w = weights, offset = offset[[1L]])
    beta <- auxreg$coefficients
    gamma <- c(log(sqrt(sum(weights * auxreg$residuals^2)/auxreg$df.residual)), rep(0, q-1))
    start <- if(dfest) c(beta, gamma, log(10)) else c(beta, gamma)
  }
  if(is.list(start)) start <- do.call("c", start) 
  if(length(start) > k + q + dfest) {
    warning(paste("too many entries in start! only first", k + q + dfest, "entries are considered"))
    start <- start[1: (k + q + dfest)]
  }

  ## various fitted quantities (parameters, linear predictors, etc.)
  fitfun <- function(par) {
    beta <- par[seq.int(length.out = k)]
    gamma <- par[seq.int(length.out = q) + k]
    delta <- if(dfest) tail(par, 1) else NULL
    mu <- drop(x %*% beta) + offset[[1L]]
    sigma <- exp(drop(z %*% gamma)) + offset[[2L]]
    df <- if(dfest) exp(delta) else df
    list(
      beta = beta,
      gamma = gamma,
      delta = delta,
      mu = mu,
      sigma = sigma,
      df = df
    )
  }

  ## objective function
  loglikfun <- function(par) {
    fit <- fitfun(par)
    ll <- with(fit, ifelse(y <= left, pdist(left, mu, sigma, df, lower.tail = TRUE),
      ifelse(y >= right, pdist(right, mu, sigma, df, lower.tail = FALSE),
      ddist(y, mu, sigma, df))))

    return(-sum(weights * ll))
  }

  opt <- suppressWarnings(optim(par = start, fn = loglikfun, method = method, hessian = hessian, control = control))
  if(opt$convergence > 0) {
    converged <- FALSE
    warning("optimization failed to converge")
  } else {
    converged <- TRUE
  }

  par <- opt$par
  fit <- fitfun(par)
  beta <- fit$beta
  gamma <- fit$gamma
  delta <- fit$delta
  mu <- fit$mu
  sigma <- fit$sigma
  vcov <- if (hessian) solve(as.matrix(opt$hessian)) else NULL
  ll <- -opt$value
  df <- if(dfest) exp(delta) else df

  names(beta) <- colnames(x)
  names(gamma) <- colnames(z)
  if (hessian) {
    colnames(vcov) <- rownames(vcov) <- c(
      colnames(x),
      paste("(scale)", colnames(z), sep = "_"),
      if(dfest) "(Log(df))")
  }

  rval <- list(
    coefficients = list(location = beta, scale = gamma, df = delta),
    df = df,
    residuals = y - mu,
    fitted.values = list(location = mu, scale = sigma),
    dist = dist,
    cens = list(left = left, right = right),
    optim = opt,  
    method = method,
    control = ocontrol,
    start = start,  
    weights = if(identical(as.vector(weights), rep.int(1, n))) NULL else weights,
    offset = list(location = if(identical(offset[[1L]], rep.int(0, n))) NULL else offset[[1L]],
    n = n,
    nobs = nobs,
    scale = if(identical(offset[[2L]], rep.int(0, n))) NULL else offset[[2L]]),
    loglik = ll,
    vcov = vcov,
    converged = converged,
    iterations = as.vector(tail(na.omit(opt$counts), 1))
  )
  
  class(rval) <- "crch"
  return(rval)
}

print.crch <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
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
    cat(paste("Distribution: ", x$dist, "\n", sep = ""))
    if(length(x$df)) {
      cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
    }
    cat("\n")
  }

  invisible(x)
}


summary.crch <- function(object, ...)
{
  ## residuals
  object$residuals <- object$residuals/object$fitted.values$scale

  ## extend coefficient table
  k <- length(object$coefficients$location)
  m <- length(object$coefficients$scale)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if(length(object$coefficients$df)) {
    cf <- list(location = cf[seq.int(length.out = k), , drop = FALSE], scale = cf[seq.int(length.out = m) + k, , drop = FALSE], 
      df = cf[nrow(cf), , drop = FALSE])
    rownames(cf$df) <- names(object$coefficients$df)
  } else {
    cf <- list(location = cf[seq.int(length.out = k), , drop = FALSE], scale = cf[seq.int(length.out = m) + k, , drop = FALSE])
  }
  rownames(cf$location) <- names(object$coefficients$location)
  rownames(cf$scale) <- names(object$coefficients$scale)
  object$coefficients <- cf

  ## delete some slots
  object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL

  ## return
  class(object) <- "summary.crch"
  object
}


print.summary.crch <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    cat(paste("Standardized residuals:\n", sep = ""))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")))

    if(NROW(x$coefficients$location)) {
      cat(paste("\nCoefficients (location model):\n", sep = ""))
      printCoefmat(x$coefficients$location, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients (in location model)\n")

    if(NROW(x$coefficients$scale)) {
      cat(paste("\nCoefficients (log(scale) model):\n", sep = ""))
      printCoefmat(x$coefficients$scale, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients (in log(scale) model)\n")

    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

    cat(paste("\nDistribution: ", x$dist, "\n", sep = ""))
    if(length(x$df)) {
      cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
    }
    cat("Log-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df\n")
    cat(paste("Number of iterations in", x$method, "optimization:", x$iterations[1L], "\n"))
  }

  invisible(x)
}

terms.crch <- function(x, model = c("location", "scale", "full"), ...) x$terms[[match.arg(model)]]

fitted.crch <- function(object, type = c("location", "scale"), ...) object$fitted.values[[match.arg(type)]]

predict.crch <- function(object, newdata = NULL,
  type = c("response", "location", "scale", "quantile"), na.action = na.pass, at = 0.5, ...)
{
  type <- match.arg(type)
  ## response/location are synonymous
  if(type == "location") type <- "response"
  ## currently all distributions supported are symmetric, hence:
  if(type == "quantile" & identical(at, 0.5)) type <- "response"
  
  if(type == "quantile") {
    qdist <- switch(object$dist,
    "student"  = function(at, location, scale, df) {
      rval <- sapply(at, function(p) qt(p, df = df) * scale + location)
      if(length(at) > 1L) {
        if(NCOL(rval) == 1L) rval <- matrix(rval, ncol = length(at),
	  dimnames = list(unique(names(rval)), NULL))
        colnames(rval) <- paste("q_", at, sep = "")
      } else {
        rval <- drop(rval)
      }
      rval
    },
    "gaussian" = function(at, location, scale, df) {
      rval <- sapply(at, function(p) qnorm(p) * scale + location)
      if(length(at) > 1L) {
        if(NCOL(rval) == 1L) rval <- matrix(rval, ncol = length(at),
	  dimnames = list(unique(names(rval)), NULL))
        colnames(rval) <- paste("q_", at, sep = "")
      } else {
        rval <- drop(rval)
      }
      rval 
    },
    "logistic" = function(at, location, scale, df) {
      rval <- sapply(at, function(p) qlogis(p) * scale + location)
      if(length(at) > 1L) {
        if(NCOL(rval) == 1L) rval <- matrix(rval, ncol = length(at),
	  dimnames = list(unique(names(rval)), NULL))
        colnames(rval) <- paste("q_", at, sep = "")
      } else {
        rval <- drop(rval)
      }
      rval
    })
  }
  
  if(missing(newdata)) {

    rval <- switch(type,
      "response" = {
        object$fitted.values$location
      },
      "scale" = {
        object$fitted.values$scale
      },
      "quantile" = {
        qdist(at, object$fitted.values$location, object$fitted.values$scale, object$df)
      },
    )
    return(rval)

  } else {

    tnam <- switch(type,
      "response" = "location",
      "scale" = "scale",
      "quantile" = "full"
    )

    mf <- model.frame(delete.response(object$terms[[tnam]]), newdata, na.action = na.action, xlev = object$levels[[tnam]])
    newdata <- newdata[rownames(mf), , drop = FALSE]
    offset <- list(location = rep.int(0, nrow(mf)), scale = rep.int(0, nrow(mf)))

    if(type %in% c("response", "quantile")) {
      X <- model.matrix(delete.response(object$terms$location), mf, contrasts = object$contrasts$location)
      if(!is.null(object$call$offset)) offset[[1L]] <- offset[[1L]] + eval(object$call$offset, newdata)
      if(!is.null(off.num <- attr(object$terms$location, "offset"))) {
        for(j in off.num) offset[[1L]] <- offset[[1L]] + eval(attr(object$terms$location, "variables")[[j + 1L]], newdata)
      }
    }
    if(type %in% c("scale", "quantile")) {
      Z <- model.matrix(object$terms$scale, mf, contrasts = object$contrasts$scale)
      if(!is.null(off.num <- attr(object$terms$scale, "offset"))) {
        for(j in off.num) offset[[2L]] <- offset[[2L]] + eval(attr(object$terms$scale, "variables")[[j + 1L]], newdata)
      }
    }
   
    rval <- switch(type,
      "response" = {
        drop(X %*% object$coefficients$location + offset[[1L]])
      },
      "scale" = {
        exp(drop(Z %*% object$coefficients$scale + offset[[2L]]))
      },
      "quantile" = {
        mu <- drop(X %*% object$coefficients$location + offset[[1L]])
        sigma <- exp(drop(Z %*% object$coefficients$scale + offset[[2L]]))
        df <- object$df
        qdist(at, mu, sigma, df)       
      }
    )
    return(rval)

  }
}


coef.crch <- function(object, model = c("full", "location", "scale", "df"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
    "location" = {
      cf$location
    },
    "scale" = {
      cf$scale
    },
    "df" = {
      cf$df
    },
    "full" = {
      nam1 <- names(cf$location)
      nam2 <- names(cf$scale)
      cf <- c(cf$location, cf$scale, cf$df)
      names(cf) <- colnames(object$vcov)
      cf
    }
  )
}


vcov.crch <- function(object, model = c("full", "location", "scale", "df"), ...) {
  vc <- object$vcov
  k <- length(object$coefficients$location)
  m <- length(object$coefficients$scale)
  l <- length(object$coefficients$df)

  model <-  match.arg(model)

  switch(model,
    "location" = {
      vc[seq.int(length.out = k), seq.int(length.out = k), drop = FALSE]
    },
    "scale" = {
      vc <- vc[seq.int(length.out = m) + k, seq.int(length.out = m) + k, drop = FALSE]
      colnames(vc) <- rownames(vc) <- names(object$coefficients$scale)
      vc
    },
    "df" = {
      vc <- vc[seq.int(length.out = l) + k + m, seq.int(length.out = l) + k + m, drop = FALSE]
      colnames(vc) <- rownames(vc) <- names(object$coefficients$df)
      vc
    },
    "full" = {
      vc
    }
  )
}

logLik.crch <- function(object, ...) structure(object$loglik, df = sum(sapply(object$coefficients, length)), class = "logLik")

residuals.crch <- function(object, type = c("standardized", "pearson", "response"), ...) {
  if(match.arg(type) == "response") object$residuals else object$residuals/object$fitted.values$scale
}

if(FALSE) {
  ## load code and data
  library("zoo")
  library("AER")
  load("Neuhof_24.RData")
  # load("pottenbrunn_24.RData")
  source("transform.R")
  source("crch.R")

  ## transform data
  d <- as.data.frame(turbine.data)
  p2w <- make_p2w(pc)
  w2p <- make_w2p(pc)
  d$wnow <- p2w(d$pnow)
  d <- d[order(d$fflin),]

  ## fit models
  mg <- crch(wnow ~ fflin | 1, data = d, dist = "gaussian", left = 3, right = 17)
  ml <- crch(wnow ~ fflin | 1, data = d, dist = "logistic", left = 3, right = 17)
  ms <- crch(wnow ~ fflin | 1, data = d, dist = "student",  left = 3, right = 17)

  ## visualize fitted means (virtually identical)
  plot(wnow ~ fflin, data = d)
  abline(coef(mg, "location"))
  abline(coef(ml, "location"), col = 2)
  abline(coef(ms, "location"), col = 4)

  ## visualize fitted densities (gaussian differen, student and logistic similar)
  x <- seq(-6, 6, by = 0.1)
  plot(x, dt(x/exp(coef(ms, "scale")), df = exp(ms$coefficients$df))/exp(coef(ms, "scale")), type = "l", col = 4)
  lines(x, dnorm(x/exp(coef(mg, "scale")))/exp(coef(mg, "scale")), col = 1)
  lines(x, dlogis(x/exp(coef(ml, "scale")))/exp(coef(ml, "scale")), col = 2)

  ## likelihood ratio test between student and gaussian model
  lrtest(mg, ms)
}
