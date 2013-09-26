hxlr <- function(formula, scale = ~ NULL, data, weights = NULL, thresholds, ...)
{
  cl <- match.call()
  if(is.null(data)) data <- cbind(model.frame(formula), model.frame(scale))
  nobs <- nrow(data)
  
  ## weights must not be 0 for clm
  if(!is.null(weights)) {
    data <- data[weights>0,]
    weights <- weights[weights>0]
  }

  stopifnot(requireNamespace("ordinal"))

  ## get environment from clm
  env <- ordinal::clm(as.character(formula), scale = scale, data = data, weights = weights, doFit = FALSE, ...)
  
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
    
  ## starting values from clm fit, threshold coefficients set to 0 and 1
  strt <- c(0, rep(1, p-1)/(p-1), ordinal::clm(formula, scale = scale, data = data, ...)$coefficients[-(1:nrow(q))])

  ## estimation
  opt <- nlminb(start=strt, objective=nll2, envir=env)

  ## converged?
  if(opt$convergence > 0) {
    converged <- FALSE
    warning("optimization failed to converge")
  } else {
    converged <- TRUE
  }

  ## compute Hessian
  hessian <- numDeriv::hessian(nll2, opt$par, envir = env)
  vcov <- solve(as.matrix(hessian)) 
  
  ## coefficients
  par <- opt$par
  nbeta <- length(par) - p - length(all.vars(scale))
  
  ## intercept coefficients
  alpha <- par[1:p]
  names(alpha) <- colnames(q)
  ## logit coefficients
  beta <- par[seq.int(length.out = nbeta) + p]
  ## scale coefficients
  delta2 <- tail(par, length(all.vars(scale)))

  ## likelihood
  ll <- -opt$objective


  ## coefficients for location and scale (gamma, delta)
  if(p ==2) {
    gamma <- c(-alpha[1], beta)/alpha[2]
    delta <- c(-log(alpha[2]), delta2)
    names(delta)[1] <- "(Intercept)"
  } else {gamma <- delta <- NULL}
  
  ## output
  rval <- list(
    coefficients = list(intercept = alpha, location = beta, scale = delta2),
    coefficients2 = list(location = gamma, scale = delta),
    optim = opt,  
    loglik = - opt$objective,
    start = start, 
    thresholds = thresholds,
    data = data,
    formula = formula,
    scale = scale,
    nobs = nobs,
    call = cl,
    converged = converged,
    vcov = vcov
  )
  
  class(rval) <- "hxlr"
  return(rval)
}

logLik.hxlr <- function(object, ...)
  structure(object$loglik, df = sum(sapply(object$coefficients, length)), class = "logLik")

predict.hxlr <- function(object, newdata = object$data,
  type = c("class", "probability", "cumprob", "location", "scale"),
  thresholds = object$thresholds, ...)
{
  type <- match.arg(type)
  
  ## get coefficients
  gamma <- object$coefficients2$location
  delta <- object$coefficients2$scale

  ## matrices for location and scale model from newdata
  x <- model.matrix(delete.response(terms(object$formula)), newdata)
  z <- model.matrix(object$scale, newdata)
  
  ## location and scale 
  mu <- drop(x %*% gamma)
  sigma <- exp(drop(z %*% delta))

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
    if(length(x$coefficients$location) | length(x$coefficients$intercept)) {
      cat(paste("Coefficients:\n", sep = ""))
       print.default(format(c(x$coefficients$intercept, x$coefficients$location), digits = digits), print.gap = 2, quote = FALSE)
       cat("\n")
    } else cat("No coefficients\n\n")
    if(length(x$coefficients$scale)) {
      cat(paste("Coefficients (scale model with log link):\n", sep = ""))
      print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in scale model)\n\n")
    cat("---\n")
    cat("Converted coefficients for location and scale:\n")
    if(length(x$coefficients2$location)) {
      cat(paste("Coefficients (location model):\n", sep = ""))
       print.default(format(x$coefficients2$location, digits = digits), print.gap = 2, quote = FALSE)
       cat("\n")
    } else cat("No coefficients (in location model)\n\n")
    if(length(x$coefficients2$scale)) {
      cat(paste("Coefficients (scale model with log link):\n", sep = ""))
      print.default(format(x$coefficients2$scale, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in scale model)\n\n")
    cat("\n")
  }
  invisible(x)
}

summary.hxlr <- function(object, ...)
{
  ## extend coefficient table
  k <- length(object$coefficients$intercept)
  l <- length(object$coefficients$location)
  m <- length(object$coefficients$scale)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if(length(object$coefficients$scale)) {
    cf <- list(
      intercept = cf[seq.int(length.out = k), , drop = FALSE],
      location = cf[seq.int(length.out = l) + k, , drop = FALSE],
      scale = cf[seq.int(length.out = m) + k + l, , drop = FALSE])
    rownames(cf$scale) <- names(object$coefficients$scale)
  } else {
    cf <- list(intercept = cf[seq.int(length.out = k), , drop = FALSE], location = cf[seq.int(length.out = l) + k, , drop = FALSE])
  }
  rownames(cf$location) <- names(object$coefficients$location)
  rownames(cf$intercept) <- names(object$coefficients$intercept)
  object$coefficients <- cf

  ## info
  object$info <- cbind(nobs = object$nobs, logLik = object$loglik, AIC = AIC(object), niter = object$optim$iterations)
  rownames(object$info) = ""

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
    print(x$info)
    if(NROW(x$coefficients$location) | NROW(x$coefficients$intercept)) {
      cat(paste("\nCoefficients:\n", sep = ""))
      printCoefmat(rbind(x$coefficients$intercept, x$coefficients$location), digits = digits, signif.legend = FALSE)
    } 

    if(NROW(x$coefficients$scale)) {
      cat(paste("\nlog-scale coefficients:\n", sep = ""))
      printCoefmat(x$coefficients$scale, digits = digits, signif.legend = FALSE)
    } 

    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
  }
  invisible(x)
}

