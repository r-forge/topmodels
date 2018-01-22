## function to set some control parameters for boosting (replaces crch.control() for boosting)
crch.glmnet <- function(maxit = 100000, nlambda = NULL, 
  lambda.min.ratio = 0.0001, lambda= NULL,
  dot = "separate", mstop = c("max", "aic", "bic", "cv"),  nfolds = 10, 
  foldid = NULL, reltol = 1E-7, ...)
{ 
  if(is.null(nlambda)) {
    if(!is.null(lambda)) {
      nlambda <- length(lambda)
      warning("nlambda is ignored if lambda is supplied")
    } else nlambda <- 100
  }
  lambda <- if(is.null(lambda)) -1
  if(is.numeric(mstop)) {
    maxit <- mstop
    mstop <- "max"
  }
  rval <- list(maxit = maxit, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
    lambda = lambda, nfolds = nfolds, foldid = foldid, reltol = reltol,
    start = start, dot = dot, mstop = match.arg(mstop), fit = "crch.glmnet.fit", ...)
  rval
}


## Function to (re)standardize vectors and model matrices:
standardize.matrix <- function(x, restandardize = FALSE, center = NULL, scale = NULL,
  weights = rep(1, NROW(x))) {
  if(is.null(center)) center <- attr(x, "standardize:center")
  if(is.null(scale)) scale <- attr(x, "standardize:scale")
  x <- as.matrix(x)
  if(is.null(center))
    center <- colSums(x*weights)/sum(weights)

  center[grep("(Intercept)", colnames(x))] <- 0
  mcenter <- matrix(rep(center, nrow(x)), nrow(x), ncol(x), byrow = TRUE)

  if(is.null(scale))
    scale <- sqrt(colSums(weights*(x-mcenter)^2)/(sum(weights)))

  scale[grep("(Intercept)", colnames(x))] <- 1
  mscale  <- matrix(rep(scale , nrow(x)), nrow(x), ncol(x), byrow = TRUE)
 
  if(!restandardize) { 
    x <- (x - mcenter)/mscale
  } else {
    x <- x*mscale + mcenter
  }
  attr(x, "standardize:center") <- center
  attr(x, "standardize:scale") <- scale
  x  
}

## Function to (re)standardize coefficients
standardize.coefficients <- function(coef, restandardize = FALSE, center, scale) {
  interceptind <- grep("(Intercept)", names(coef))
  if(length(interceptind > 0)) {
    intercept <- coef[interceptind]
    center <- center[-interceptind]
    scale <- scale[-interceptind]
    coef <- coef[-interceptind]
  } else {
    intercept <- 0
    names(intercept) <- "(Intercept)"
  }
  if(restandardize) {
    intercept <- intercept - sum(coef*center/scale)
    coef <- coef/scale
  } else {
    coef <- coef*scale
    intercept <- intercept + sum(coef*center/scale)
  }
  c(intercept, coef)  
}


## actual fitting function
crch.glmnet.fit <- function(x, z, y, left, right, truncated = FALSE, 
  dist = "gaussian", df = NULL, link.scale = "log", type = "ml",
  weights = NULL, offset = NULL, control = crch.glmnet()) 
{
  ## type = "crps" currently not supported
  if(type == "crps") stop("type = 'crps' currently not supported for boosting")
  ## response and regressor matrix
  n <- NROW(x)  
  k <- NCOL(x)
  if(is.null(weights)) weights <- rep.int(1, n)
  nobs <- sum(weights > 0)
  dfest <- identical(dist, "student") & is.null(df)
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

  ## control parameters
  maxit <- control$maxit
  nlambda <- control$nlambda
  lambda.min.ratio <- control$lambda.min.ratio
  lambda <- control$lambda
  start <- control$start
  mstop <- control$mstop
  reltol <- control$reltol

  ## extend left and right to vectors of length n
  left2 <- left
  right2 <- right
  if(length(left) == 1) {
    left <- rep(left, n)
  }
  if(length(right) == 1) {
    right <- rep(right, n)
  }

  ## link
linkfun2 <- switch(link.scale, "log" = 1L, "identity" = 2L, "quadratic" = 3L) 
  if(is.character(link.scale)) {
    linkstr <- link.scale
    if(linkstr != "quadratic") {
      linkobj <- make.link(linkstr)
      linkobj$dmu.deta <- switch(linkstr, 
        "identity" = function(eta) rep.int(0, length(eta)), 
        "log" = function(eta) pmax(exp(eta), .Machine$double.eps))
    } else {
      linkobj <- structure(list(
        linkfun = function(mu) mu^2,
        linkinv = function(eta) sqrt(eta),
        mu.eta = function(eta) 1/2/sqrt(eta),
        dmu.deta = function(eta) -1/4/sqrt(eta^3),
        valideta = function(eta) TRUE,
        name = "quadratic"
      ), class = "link-glm")
    }
  } else {
    linkobj <- link.scale
    linkstr <- link.scale$name
    if(is.null(linkobj$dmu.deta) & !hessian) {
      warning("link.scale needs to provide dmu.deta component for analytical Hessian. Numerical Hessian is employed.")
      hessian <- TRUE
    }
  }
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta
  dmu.deta <- linkobj$dmu.deta



  ## scale response and regressors
  x <- standardize.matrix(x, weights = weights)
  z <- standardize.matrix(z, weights = weights)
  basefit <- crch.fit(x[,1, drop = FALSE], z[,1, drop = FALSE], y, left, right, 
      truncated, dist, df, link.scale, type, weights, offset)
  y <- standardize.matrix(y, center = basefit$coefficients$location, scale = linkinv(basefit$coefficients$scale))
    
  standardize <- list(
    x = list(center = attr(x, "standardize:center"), scale = attr(x, "standardize:scale")),
    z = list(center = attr(z, "standardize:center"), scale = attr(z, "standardize:scale")),
    y = list(center = attr(y, "standardize:center"), scale = attr(y, "standardize:scale")))

  ## standardize left, right, and offset
  left <- (left - standardize$y$center)/standardize$y$scale
  right <- (right - standardize$y$center)/standardize$y$scale
  offset[[1L]]  <- offset[[1L]]/standardize$y$scale
  
 if(is.character(dist)){
    ## distribution functions
    if(truncated) {
      ddist2 <- switch(dist, 
        "student"  = dtt, "gaussian" = dtnorm, "logistic" = dtlogis)
      sdist2 <- switch(dist, 
        "student"  = stt, "gaussian" = stnorm, "logistic" = stlogis)
      hdist2 <- switch(dist, 
        "student"  = htt, "gaussian" = htnorm, "logistic" = htlogis)
    } else {
      ddist2 <- switch(dist, 
        "student"  = dct, "gaussian" = dcnorm, "logistic" = dclogis)
      sdist2 <- switch(dist, 
        "student"  = sct, "gaussian" = scnorm, "logistic" = sclogis)
      hdist2 <- switch(dist, 
        "student"  = hct, "gaussian" = hcnorm, "logistic" = hclogis)
    }
    ddist <- if(dist == "student") ddist2 else function(..., df) ddist2(...)
    sdist <- if(dist == "student") sdist2 else function(..., df) sdist2(...)
    hdist <- if(dist == "student") hdist2 else function(..., df) hdist2(...)


  } else { 
    ## for user defined distribution (requires list with ddist, sdist (optional)
    ## and hdist (optional), ddist, sdist, and hdist must be functions with
    ## arguments x, mean, sd, df, left, right, and log)
    ddist <- dist$ddist
    sdist <- if(is.null(dist$sdist)) NULL else  dist$sdist
    if(is.null(dist$hdist)) {
      if(hessian == FALSE) warning("no analytic hessian available. Hessian is set to TRUE and numerical Hessian from optim is employed")
      hessian <- TRUE     
    } else hdist <- dist$hdist 
    dist <- "user defined"
  }




 
  ## various fitted quantities (parameters, linear predictors, etc.)
  fitfun <- function(par, subset = 1:n) {
    beta <- par[seq.int(length.out = k)]
    gamma <- par[seq.int(length.out = q) + k]
    mu <- (drop(x %*% beta) + offset[[1L]])[subset]
    zgamma <- (drop(z %*% gamma) + offset[[2L]])[subset]
    sigma <- linkinv(zgamma)
    list(
      beta = beta,
      gamma = gamma,
      mu = mu,
      zgamma = zgamma,
      sigma = sigma
    )
  }

  ## objective function
  loglikfun <- function(par, subset = 1:n, sum = TRUE) {
    fit <- fitfun(par, subset)
    ll <- with(fit,  
        ddist(y[subset], mu, sigma, df = df, left = left[subset], right = right[subset], log = TRUE))
    if(sum) if(any(!is.finite(ll))) NaN else -sum(weights[subset] * ll)  
    else weights[subset] * ll
  }

  ## functions to evaluate gradient
  if(dfest | is.null(sdist)) {
    stop("score function required for boosting not available")
  } else { 

    wobs <- function(coef, which = "mu") {
      fit <- fitfun(par)
      grad <- with(fit, 
        sdist(y, mu, sigma, df = df, left = left, right = right, which = which))
      
      hess <- with(fit, hdist(y, mu, sigma, left = left, right = right,
        df = df, which = which))      

      if(which == "sigma") {
        grad <- grad * mu.eta(fit$zgamma)
        hess <- hess*mu.eta(fit$zgamma)^2 + grad*dmu.deta(fit$zgamma)
        which <- "zgamma"
      }


      wo <- fit[[which]]  - 1 / hess * grad
      rval <- data.frame(wo, -hess)
      names(rval) <- c("wo", "weights")
      return(rval)
    }
  }
  
  ## initialize intercepts
  par <- rep(0, k+q)
  if(link.scale %in% c("identity", "quadratic")) par[k+1] <- 1 
  coefpath <- par
  loglikpath <- -loglikfun(par)


  ## penalized IWLS
  fit <- fitfun(par)
  wobs2 <- wobs(par, which = "mu")
  w <- wobs2$weights  #weights
  wo <- wobs2$wo
  lambdamax <- max(apply(x, 2, function(x) abs(x%*%wo)/n))

  wobs2 <- wobs(par, which = "sigma")
  w <- wobs2$weights  #weights
  wo <- wobs2$wo
#  lambdamax <- max(c(lambdamax, max(apply(z, 2, function(z) abs(sum(w*z*wo))))))
  lambdaseq <- exp(seq(log(lambdamax), log(lambdamax*lambda.min.ratio), length.out = nlambda))
iter <- 0


browser()
a <- .Call("crchglmnet", x, z, y, left, right, lambdaseq, as.integer(maxit), reltol)
coefpath <- a[[1]]
loglikpath <- a[[2]]


  ## cross validation to find optimum mstop
  if(mstop == "cv") {
    nfolds <- control$nfolds
    foldid <- control$foldid
    control$standardize <- FALSE
    control$mstop <- "max"
    control$start <- NULL
    if(is.null(foldid)) foldid <- (sample(1:nfolds, size = length(y), replace = TRUE))
  
    lltestall <- matrix(0, length(y), nlambda)
    for(i in 1:nfolds) {
      train <- foldid != i
      test <- foldid == i
      glmnet <- crch.glmnet.fit(x = x[train, , drop = FALSE], y = y[train], z = z[train, , drop = FALSE], left = left[train], 
          right = right[train], link.scale = link.scale, dist = dist, df = df, weights = weights[train], 
          offset = list(offset[[1L]][train], offset[[2L]][train]), 
          control = control, truncated = truncated)
      lltest <- apply(glmnet$coefpath, 1, loglikfun, subset = test, sum = FALSE)
#      lltest <- cbind(lltest, matrix(rep(lltest[,NCOL(lltest)], maxit + 1 - NCOL(lltest)), nrow = NROW(lltest)))
      lltestall[test,] <- lltest
    }

    mstopopt<- which.max(colMeans(lltestall))
  } else {
    mstopopt <- NULL
  }

#browser()     
  colnames(coefpath) <- c(colnames(x), paste("scale_", colnames(z), sep = ""))
  ## restandardize
  left  <-  left*standardize$y$scale + standardize$y$center
  right <- right*standardize$y$scale + standardize$y$center
  offset[[1L]] <- offset[[1L]]*standardize$y$scale
  beta <- coefpath[,seq.int(length.out = k), drop = FALSE]
  gamma <- coefpath[,seq.int(length.out = q) + k, drop = FALSE]
  beta[,] <- t(apply(beta, 1, standardize.coefficients, 
    center = standardize$x$center, scale = standardize$x$scale, restandardize = TRUE))*standardize$y$scale
  beta[, grep("(Intercept)", colnames(beta))] <- 
    beta[, grep("(Intercept)", colnames(beta)), drop = FALSE] + standardize$y$center
  gamma[,] <- t(apply(gamma, 1, standardize.coefficients, 
    center = standardize$z$center, scale = standardize$z$scale, restandardize = TRUE))
  gamma[, grep("(Intercept)", colnames(gamma))] <- 
    gamma[, grep("(Intercept)", colnames(gamma)), drop = FALSE] + linkfun(standardize$y$scale)
  coefpath <- cbind(beta, gamma)
  y <- standardize.matrix(y, restandardize = TRUE)
  x <- standardize.matrix(x, restandardize = TRUE)
  z <- standardize.matrix(z, restandardize = TRUE)

  ## optimum stopping iterations
  mstopopt <- c("cv" = mstopopt, "max" = NROW(coefpath), 
      "aic" = which.min(2*rowSums(coefpath!=0) - 2* loglikpath),
      "bic" = which.min(log(nobs)*rowSums(coefpath!=0) - 2* loglikpath)) - 1

  ## coefficients, fitted values and likelihood for optimum stopping iteration
  par <- coefpath[mstopopt[[mstop]] + 1,]
  fit <- fitfun(par)
  beta <- fit$beta
  gamma <- fit$gamma
  mu <- fit$mu
  sigma <- fit$sigma
  ll <- - loglikfun(par)

  names(beta) <- colnames(x)
  names(gamma) <- colnames(z)
  rownames(coefpath) <- c(0:(NROW(coefpath) - 1))

  rval <- list(
    coefficients = list(location = beta, scale = gamma),
    df = df,
    residuals = y - mu,
    fitted.values = list(location = mu, scale = sigma),
    dist = dist,
    cens = list(left = left2, right = right2), 
    control = control,
    weights = if(identical(as.vector(weights), rep.int(1, n))) NULL else weights,
    offset = list(location = if(identical(offset[[1L]], rep.int(0, n))) NULL else 
      offset[[1L]],
      scale = if(identical(offset[[2L]], rep.int(0, n))) NULL else offset[[2L]]),
    n = n,
    nobs = nobs,
    loglik = ll,
    link = list(scale = linkobj),
    truncated = truncated,
#    converged = converged,
    coefpath = coefpath,
    loglikpath = loglikpath,
    iterations = NROW(coefpath) - 1,
    lambda = lambdaseq,
    mstop = mstop,
    mstopopt = mstopopt,
    standardize = standardize
  )
  
  class(rval) <- c("crch.glmnet", "crch.boost", "crch")
  return(rval)
}





plot.crch.glmnet <- function(object, 
  standardize = TRUE, which = c("both", "location", "scale"), mstop = NULL,
  coef.label = TRUE, col = NULL, ylim = NULL, ...) 
{
  which <- match.arg(which)
  k <- length(object$coefficient$location)
  q <- length(object$coefficient$scale)
  if(is.null(col)) col <- if(which == "both") c(1,2) else c(1,1)

  if(is.null(mstop)) {
    mstop <- object$mstop
  } else {
    mstop <- match.arg(mstop, c("max", "aic", "bic", "cv", "all", "no"))
  }
  if(mstop %in% c("no", "max")) {
    mstop <- NULL 
  } else {
    if(mstop == "all") mstop <- c("aic", "bic", "cv")
    mstop <- object$mstopopt[mstop] 
  }

  
  beta <- object$coefpath[,seq.int(length.out = k), drop = FALSE]
  gamma <- object$coefpath[,seq.int(length.out = q) + k, drop = FALSE]
  if(standardize) {
    beta[,] <- t(apply(beta, 1, standardize.coefficients, 
      center = object$standardize$x$center, scale = object$standardize$x$scale))
    beta[, grep("(Intercept)", colnames(beta))] <- 
      beta[, grep("(Intercept)", colnames(beta))] - 
      object$standardize$y$center 
    beta <- beta/object$standardize$y$scale  
    gamma[,] <- t(apply(gamma, 1, standardize.coefficients, 
      center = object$standardize$z$center, scale = object$standardize$z$scale))
    gamma[, grep("(Intercept)", colnames(gamma))] <- 
      gamma[, grep("(Intercept)", colnames(gamma))] - 
      object$link$scale$linkfun(object$standardize$y$scale)  
    colnames(gamma) <- names(object$coefficients$scale)
  } 
  path <- cbind(beta, gamma)
  ylab <- if(standardize) "standardized coefficients" else "coefficients"
 

  if(which == "location") {path <- path[,1:k]; q <- 0}
  if(which == "scale") {path <- path[,(k+1):(k+q)]; k <- 0}
  
  if(is.logical(coef.label)) {
    coef.label <- if(coef.label) {
        list(location = colnames(path)[1:k], scale = colnames(path)[(k+1):(q+k)])
    } else NULL
  }
  ## adjust figure margins
  rmar <- if(length(coef.label)) max(strwidth(colnames(object$coefpath), units = "in"))+0.5 else 0.42
  umar <- if(is.null(mstop)) 0.82 else 1.02
  par(mai = c(1.02, 0.82, umar, rmar))
  if(is.null(ylim)) ylim <- c(min(path), max(path))

  plot(-log(object$lambda), path[,1], ylab = ylab, xlab = "log-lambda", col = col[1], 
    type = "l", xaxt = "n", ylim = ylim, ...)
  for(i in 1:k) lines(-log(object$lambda), path[,i], col = col[1])
  for(i in 1:q+k) lines(-log(object$lambda), path[,i], col = col[2])
  axis(1, axTicks(1), -axTicks(1))
  ## label paths
  if(length(coef.label)) {
    if(k!=0) axis(4, at = as.data.frame(tail(path, 1))[1:k][coef.label$location], 
      labels = coef.label$location, las = 1, col.axis = col[1], col.ticks = col[1])
    if(q!=0) axis(4, at = as.data.frame(tail(path, 1))[(k+1):(q+k)][coef.label$scale], 
      labels = coef.label$scale, las = 1, col.axis = tail(col, 1), col.ticks = tail(col,1))
  }

  ## vertical line at mstop
  if(!is.null(mstop)) {
    abline(v = -log(object$lambda[mstop]), lty = 2)
    axis(3, at = -log(object$lambda[mstop]), labels = names(mstop))
  }
}
