## function to set some control parameters for boosting (replaces crch.control() for boosting)
crch.boost <- function(method = "boosting", mstop = 100, nu = 0.1, standardize = TRUE, 
  nfolds = 10, foldid = NULL, dot = "separate",  ...)
{
  rval <- list(method = method, mstop = mstop, nu = nu, standardize = standardize, 
    nfolds = nfolds, foldid = foldid, dot = dot, fit = "crch.boost.fit")
  rval <- c(rval, list(...))
  rval
}

## actual fitting function
crch.boost.fit <- function(x, z, y, left, right, truncated = FALSE, 
  dist = "gaussian", df = NULL, link.scale = "log",
  weights = NULL, offset = NULL, control = crch.boost()) 
{
  ## standardize response and regressors
  standardize <- control$standardize
  if(standardize) {
    centerx <- matrix(rep(c(0, colMeans(x)[-1]), nrow(x)), nrow(x), ncol(x), byrow = TRUE)
    scalex <- matrix(rep(c(1, apply(x, 2, sd)[-1]), nrow(x)), nrow(x), ncol(x), byrow = TRUE)
    x <- (x - centerx)/scalex
    
    if(!is.null(z)) {
      centerz <- matrix(rep(c(0, colMeans(z)[-1]), nrow(z)), nrow(z), ncol(z), byrow = TRUE)
      scalez <- matrix(rep(c(1, apply(z, 2, sd)[-1]), nrow(z)), nrow(z), ncol(z), byrow = TRUE)
      z <- (z - centerz)/scalez
    }

    centery <- mean(y)
    scaley <- sd(y)
    y <- (y - centery)/scaley
    left <- (left - centery)/scaley
    right <- (right - centery)/scaley
  }
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
  ocontrol <- control
  method <- control$method
  mstop <- control$mstop
  nu <- control$nu

  if(is.character(dist)){
    ## distribution functions
    if(truncated) {
      ddist2 <- switch(dist, 
        "student"  = dtt, "gaussian" = dtnorm, "logistic" = dtlogis)
      sdist2 <- switch(dist, 
        "student"  = stt, "gaussian" = stnorm, "logistic" = stlogis)
    } else {
      ddist2 <- switch(dist, 
        "student"  = dct, "gaussian" = dcnorm, "logistic" = dclogis)
      sdist2 <- switch(dist, 
        "student"  = sct, "gaussian" = scnorm, "logistic" = sclogis)
    }
    ddist <- if(dist == "student") ddist else function(..., df) ddist2(...)
    sdist <- if(dist == "student") sdist else function(..., df) sdist2(...)


  } else { 
    ## for user defined distribution (requires list with ddist, sdist (optional)
    ## and hdist (optional), ddist, sdist, and hdist must be functions with
    ## arguments x, mean, sd, df, left, right, and log)
    ddist <- dist$ddist
    sdist <- if(is.null(dist$sdist)) NULL else  dist$sdist
    dist <- "user defined"
  }


  ## link
  if(is.character(link.scale)) {
    linkstr <- link.scale
    if(linkstr != "quadratic") {
      linkobj <- make.link(linkstr)
    } else {
      linkobj <- structure(list(
        linkfun = function(mu) mu^2,
        linkinv = function(eta) sqrt(eta),
        mu.eta = function(eta) 1/2/sqrt(eta),
        valideta = function(eta) TRUE,
        name = "quadratic"
      ), class = "link-glm")
    }
  } else {
    linkobj <- link.scale
    linkstr <- link.scale$name
  }
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta

 
  ## various fitted quantities (parameters, linear predictors, etc.)
  fitfun <- function(par, subset = 1:n) {
    beta <- par[seq.int(length.out = k)]
    gamma <- par[seq.int(length.out = q) + k]
    mu <- drop(x[subset,] %*% beta) + offset[[1L]][subset]
    zgamma <- drop(z[subset,] %*% gamma) + offset[[2L]][subset]
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
        ddist(y[subset], mu, sigma, df = df, left = left, right = right, log = TRUE))
    if(sum) if(any(!is.finite(ll))) NaN else -sum(weights[subset] * ll)  
    else weights[subset] * ll
  }

  ## functions to evaluate gradients and hessian
  if(dfest | is.null(sdist)) {
    stop("score function required for boosting not available")
  } else { 
    gradfun <- function(par) {
      fit <- fitfun(par)
      grad <- with(fit, 
        sdist(y, mu, sigma, df = df, left = left, right = right))
      grad <- cbind(grad[,1], grad[,2] * mu.eta(fit$zgamma))
      return(-weights * grad)
    }
  }


  ## initialize intercepts
  par <- rep(0, k+q)
  par[1] <- mean(y)
  par[k+1] <- log(sd(y))
  ## actual boosting
  coefpath <- NULL
  xtx <- apply(x, 2, function(x) sum(x^2))
  ztz <- apply(z, 2, function(x) sum(x^2))
  for(i in 1:mstop) {
    ## gradient of negative likelihood
    grad <- gradfun(par)

    ## location
    basefits <- apply(x, 2, function(x) t(x) %*% -grad[,1])/xtx
    rss <- colSums((matrix(basefits, n, k, byrow = TRUE)*x + matrix(rep(grad[,1], k), n, k))^2)
    par2 <- par
    minind <- which.min(rss)
    par2[minind] <- par2[minind] + nu*basefits[minind]
    ## scale
    basefits <- apply(z, 2, function(z) t(z) %*% -grad[,2])/ztz
    rss <- colSums((matrix(basefits, n, q, byrow = TRUE)*z + matrix(rep(grad[,2], q), n, q))^2)
    par3 <- par
    minind <- which.min(rss)
    par3[k + minind] <- par3[k + minind] + nu*basefits[minind]


    ## to compare mu and sigma improvements loglik must be used instead of rss
    par <- if(-loglikfun(par3) > -loglikfun(par2)) par3 else par2
    coefpath <- rbind(coefpath, par)
  }

  ## cross validation to find optimum mstop
  cv <- NULL
  if(method == "cvboosting") {
    nfolds <- control$nfolds
    foldid <- control$foldid
    control$standardize <- FALSE
    control$method <- "boosting"
    if(is.null(foldid)) foldid <- (sample(1:nfolds, size = length(y), replace = TRUE))
  
    lltestall <- matrix(0, length(y), mstop)
    for(i in 1:nfolds) {
      train <- foldid != i
      test <- foldid == i
      boost <- crch.boost.fit(x = x[train, ], y = y[train], z = z[train, ], left = left, 
          right = right, link.scale = link.scale, dist = dist, df = df, weights = weights[train], 
          offset = list(offset[[1L]][train], offset[[2L]][train]), 
          control = control, truncated = truncated)
      lltest <- apply(boost$coefpath, 1, loglikfun, subset = test, sum = FALSE)
      lltestall[test,] <- lltest
    }

    bootll <- NULL
    for(i in 1:100) {
      bootindex <- sample(1:length(y), replace = TRUE)
      bootll <- rbind(bootll, colMeans(lltestall[bootindex,]))
    }
  
    mstopopt <- which.max(colMeans(bootll))

    opt1sd <- max(colMeans(bootll)) - sd(bootll[,mstopopt])
    mstopopt1se <- which(colMeans(bootll)>opt1sd)[1]
    cv <- list(bootll = bootll, mstopopt = mstopopt, mstopopt1se = mstopopt1se)
  }


  ## destandardize
  if(standardize) {
    beta <- coefpath[,seq.int(length.out = k)]
    gamma <- coefpath[,seq.int(length.out = q) + k]
    centerpathx <- matrix(rep(centerx[1,], mstop), mstop, k, byrow = TRUE)
    centerpathz <- matrix(rep(centerz[1,], mstop), mstop, q, byrow = TRUE)
    scalepathx  <- matrix(rep( scalex[1,], mstop), mstop, k, byrow = TRUE)
    scalepathz  <- matrix(rep( scalez[1,], mstop), mstop, q, byrow = TRUE)
    beta  <- cbind( beta[,1]*scaley - rowSums(as.matrix(beta[,-1]*centerpathx[,-1]/scalepathx[,-1]))*
      scaley + centery, beta[,-1]*scaley/scalepathx[,-1])
    if(NCOL(z) > 1) {
      gamma <- cbind(gamma[,1] - rowSums(as.matrix(gamma[,-1]*centerpathz[,-1]/scalepathz[,-1])) 
        + linkfun(scaley),  gamma[,-1]/scalepathz[,-1])
    } else {
      gamma <- gamma + linkfun(scaley)
    }
    coefpath <- cbind(beta, gamma)
    par <- if(method == "cvboosting") coefpath[mstopopt,] else tail(coefpath, 1)
    y <- y*scaley + centery
    x <- x*scalex + centerx
    z <- z*scalez + centerz
    left <- left*scaley + centery
    right <- right*scaley + centery
  }
  fit <- fitfun(par)
  beta <- fit$beta
  gamma <- fit$gamma
  mu <- fit$mu
  sigma <- fit$sigma
  ll <- -loglikfun(par)

  names(beta) <- colnames(x)
  names(gamma) <- colnames(z)
  colnames(coefpath) <- c(colnames(x), paste("scale_", colnames(z), sep = ""))
  rownames(coefpath) <- c(1:mstop)

  rval <- list(
    coefficients = list(location = beta, scale = gamma),
    df = df,
    residuals = y - mu,
    fitted.values = list(location = mu, scale = sigma),
    dist = dist,
    cens = list(left = left, right = right), 
    method = method, 
    control = ocontrol,
    weights = if(identical(as.vector(weights), rep.int(1, n))) NULL else weights,
    offset = list(location = if(identical(offset[[1L]], rep.int(0, n))) NULL else 
      offset[[1L]],
      scale = if(identical(offset[[2L]], rep.int(0, n))) NULL else offset[[2L]]),
    n = n,
    nobs = nobs,
    loglik = ll,
    vcov = matrix(NA, k+q, n+k),
    link = list(scale = linkobj),
    truncated = truncated,
    converged = TRUE,
    coefpath = coefpath,
    loglikpath = - apply(coefpath, 1, loglikfun),
    iterations = mstop,
    stepsize = nu,
    cv = cv
  )
  
  class(rval) <- c("crch.boost", "crch")
  return(rval)
}

