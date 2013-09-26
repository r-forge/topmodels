hxlr <-
function(formula, scale=~NULL, data, weights = NULL, thresholds, ...) {
  cl <- match.call()
  if(is.null(data)) data <- cbind(model.frame(formula), model.frame(scale))
  nobs <- nrow(data)
  
  ## weights must not be 0 for clm
  if(!is.null(weights)) {
    data <- data[weights>0,]
    weights <- weights[weights>0]
  }

  ## get environment from clm
  env <- clm(as.character(formula), scale = scale, data = data, weights = weights, doFit=FALSE, ...)
  
  ## thresholds can also be data.frame with several columns (predictor variables for intercept model)
  q <- model.matrix(~thresholds)
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
  strt <- c(0, rep(1, p-1)/(p-1), clm(formula, scale = scale, data = data, ...)$coefficients[-(1:nrow(q))])

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
    coefficients = list(intercept = alpha, predictor = beta, scale = delta2),
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
