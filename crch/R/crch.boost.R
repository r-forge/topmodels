## function to set some control parameters for boosting (replaces crch.control() for boosting)
crch.boost <- function(method = "boosting", mstop = 100, nu = 0.1, scale = TRUE, 
  nfolds = 10, foldid = NULL, start = NULL, dot = "separate",  ...)
{
  rval <- list(method = method, mstop = mstop, nu = nu, scale = scale, 
    nfolds = nfolds, foldid = foldid, start = start, dot = dot, fit = "crch.boost.fit")
  rval <- c(rval, list(...))
  rval
}

## actual fitting function
crch.boost.fit <- function(x, z, y, left, right, truncated = FALSE, 
  dist = "gaussian", df = NULL, link.scale = "log",
  weights = NULL, offset = NULL, control = crch.boost()) 
{
  
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
  start <- control$start


  ## scale response and regressors
  scale <- control$scale
  if(scale) {
    x <- scale.matrix(x)
    z <- scale.matrix(z)
    y <- scale.matrix(y)
    
    scale <- list(
      x = list(center = attr(x, "scale:center"), scale = attr(x, "scale:scale")),
      z = list(center = attr(z, "scale:center"), scale = attr(z, "scale:scale")),
      y = list(center = attr(y, "scale:center"), scale = attr(y, "scale:scale")))

    left <- (left - scale$y$center)/scale$y$scale
    right <- (right - scale$y$center)/scale$y$scale
    
    ## rescale coefpath if start is supplied
    if(!is.null(start)) {
      beta <- start$coefpath[,seq.int(length.out = k)]
      beta <- t(apply(beta, 1, scale.coefficients, center = scale$x$center, scale = scale$x$scale))
      beta[, grep("(Intercept)", colnames(beta))] <- beta[, grep("(Intercept)", colnames(beta))] - scale$y$center 
      beta <- beta/scale$y$scale  
      gamma <- start$coefpath[,seq.int(length.out = q) + k]
      gamma <- t(apply(gamma, 1, scale.coefficients, center = scale$z$center, scale = scale$z$scale))
      gamma[, grep("(Intercept)", colnames(gamma))] <- gamma[, grep("(Intercept)", colnames(gamma))] - start$link$scale$linkfun(scale$y$scale)  
      start$coefpath <- cbind(beta, gamma)
    } 
  }

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
  if(is.null(start)) {
    par <- rep(0, k+q)
    par[1] <- mean(y)
    par[k+1] <- log(sd(y))
    coefpath <- NULL
    startit <- 1
  } else {
    par <- tail(start$coefpath, 1)
    coefpath <- start$coefpath
    startit <- start$iterations + 1 
  }
  

  ## actual boosting
  xtx <- apply(x, 2, function(x) sum(x^2))
  ztz <- apply(z, 2, function(x) sum(x^2))
  for(i in startit:mstop) {
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
    control$scale <- FALSE
    control$method <- "boosting"
    control$start <- NULL
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
  
    mstopopt <- which.max(colMeans(lltestall))

    opt1sd <- max(colMeans(bootll)) - sd(bootll[,mstopopt])
    mstopopt1se <- which(colMeans(bootll)>opt1sd)[1]
    cv <- list(bootll = bootll, mstopopt = mstopopt, mstopopt1se = mstopopt1se)
  }

  colnames(coefpath) <- c(colnames(x), paste("scale_", colnames(z), sep = ""))
  ## rescale
  if(control$scale) {
    beta <- coefpath[,seq.int(length.out = k), drop = FALSE]
    gamma <- coefpath[,seq.int(length.out = q) + k, drop = FALSE]
    beta[,] <- t(apply(beta, 1, scale.coefficients, 
      center = scale$x$center, scale = scale$x$scale, rescale = TRUE))*scale$y$scale
    beta[, grep("(Intercept)", colnames(beta))] <- 
      beta[, grep("(Intercept)", colnames(beta)), drop = FALSE] + scale$y$center

    gamma[,] <- t(apply(gamma, 1, scale.coefficients, 
      center = scale$z$center, scale = scale$z$scale, rescale = TRUE))
    gamma[, grep("(Intercept)", colnames(gamma))] <- 
      gamma[, grep("(Intercept)", colnames(gamma)), drop = FALSE] + linkfun(scale$y$scale)
    coefpath <- cbind(beta, gamma)




#    centerpathx <- matrix(rep(centerx[1,], mstop), mstop, k, byrow = TRUE)
#    centerpathz <- matrix(rep(centerz[1,], mstop), mstop, q, byrow = TRUE)
#    scalepathx  <- matrix(rep( scalex[1,], mstop), mstop, k, byrow = TRUE)
#    scalepathz  <- matrix(rep( scalez[1,], mstop), mstop, q, byrow = TRUE)
#    beta  <- cbind( beta[,1]*scaley - rowSums(as.matrix(beta[,-1]*centerpathx[,-1]/scalepathx[,-1]))*
#      scaley + centery, beta[,-1]*scaley/scalepathx[,-1])
#    if(NCOL(z) > 1) {
#      gamma <- cbind(gamma[,1] - rowSums(as.matrix(gamma[,-1]*centerpathz[,-1]/scalepathz[,-1])) 
#        + linkfun(scaley),  gamma[,-1]/scalepathz[,-1])
#    } else {
#      gamma <- gamma + linkfun(scaley)
#    }
#    coefpath <- cbind(beta, gamma)
    par <- if(method == "cvboosting") coefpath[mstopopt,] else tail(coefpath, 1)
    y <- scale.matrix(y, rescale = TRUE)
    x <- scale.matrix(x, rescale = TRUE)
    z <- scale.matrix(z, rescale = TRUE)
    left  <-  left*scale$y$scale + scale$y$center
    right <- right*scale$y$scale + scale$y$center
  }
  fit <- fitfun(par)
  beta <- fit$beta
  gamma <- fit$gamma
  mu <- fit$mu
  sigma <- fit$sigma
  ll <- -loglikfun(par)

  names(beta) <- colnames(x)
  names(gamma) <- colnames(z)
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

continue <- function(object, ...) UseMethod("continue")
continue.crch.boost<- function(object, mstop) {
  call <- object$call
  if(is.null(call$control)) {
    call$mstop <- mstop 
    call$start <- object
  } else {
    call$control$mstop <- mstop
    call$control$start <- object
  }
  rval <- eval(call, parent.frame())
  rval$call <- call
  rval
}

scale.matrix <- function(x, rescale = FALSE, center = NULL, scale = NULL) {
  center <- attr(x, "scale:center")
  scale <- attr(x, "scale:scale")
  x <- as.matrix(x)
  if(is.null(center))
    center <- colMeans(x)
  if(is.null(scale))
    scale <- apply(x, 2, sd)

  
  center[grep("(Intercept)", colnames(x))] <- 0
  scale[grep("(Intercept)", colnames(x))] <- 1

  
  mcenter <- matrix(rep(center, nrow(x)), nrow(x), ncol(x), byrow = TRUE)
  mscale  <- matrix(rep(scale , nrow(x)), nrow(x), ncol(x), byrow = TRUE)
  

  if(!rescale) { 
    x <- (x - mcenter)/mscale
  } else {
    x <- x*mscale + mcenter
  }
  attr(x, "scale:center") <- center
  attr(x, "scale:scale") <- scale
  x  
}

scale.coefficients <- function(coef, rescale = FALSE, center, scale) {
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
  if(rescale) {
    intercept <- intercept - sum(coef*center/scale)
    coef <- coef/scale
  } else {
    coef <- coef*scale
    intercept <- intercept + sum(coef*center/scale)
  }
  c(intercept, coef)  
}
  


plot.crch.boost <- function(object, quantity = c("coefficient", "likelihood"), scale = TRUE) {
  quantity <- match.arg(quantity)
  par(mar = c(5.1, 4.1, 2.1, 10.1))

  if(quantity == "coefficient") {
    if(scale) {
      mstopprev <- object$iterations
      beta <- object$coefpath[,seq.int(length.out = k)]
      gamma <- object$coefpath[,seq.int(length.out = q) + k]
      centerpathx <- matrix(rep(centerx[1,], mstopprev), mstopprev, k, byrow = TRUE)
      centerpathz <- matrix(rep(centerz[1,], mstopprev), mstopprev, q, byrow = TRUE)
      scalepathx  <- matrix(rep( scalex[1,], mstopprev), mstopprev, k, byrow = TRUE)
      scalepathz  <- matrix(rep( scalez[1,], mstopprev), mstopprev, q, byrow = TRUE)
                     
      beta  <- cbind(beta[,1]/scaley + 
        rowSums(as.matrix(beta[,-1]/scaley*scalepathx[,-1]*centerpathx[,-1]/scalepathx[,-1]))
        - centery/scaley, beta[,-1]/scaley*scalepathx[,-1])
      if(NCOL(z) > 1) {
        gamma <- cbind(gamma[,1] + rowSums(as.matrix(gamma[,-1]*centerpathz[,-1])) 
          - start$link$scale$linkfun(scaley),  gamma[,-1]*scalepathz[,-1])
      } else {
        gamma <- gamma - linkfun(scaley)
      }
      start$coefpath <- cbind(beta, gamma)
    } 

    plot.ts(object$coefpath, plot.type = "single", xlab = "boosting iteration", ylab = "coefficients", type = "s")
    axis(4, at = tail(object$coefpath, 1), labels = colnames(object$coefpath), las = 1)
  } else {
    coefupdate <- as.data.frame(apply(object$coefpath, 2, diff)!=0)
    coefupdate[["(Intercept)"]][rowSums(coefupdate) == 2] <- coefupdate[["scale_(Intercept)"]][rowSums(coefupdate) == 2] <- FALSE
    loglikpath <- object$loglikpath
    llcontrib <- diff(loglikpath)*coefupdate
    llcontrib <- apply(llcontrib, 2, cumsum)
    plot.ts(llcontrib, plot.type = "single", xlab = "boosting iteration", ylab = "likelihood contribution", type = "s", col = rep(c(1,2), c(k, q))
    axis(4, at = tail(llcontrib, 1), labels = colnames(llcontrib), las = 1)
  }

  if(!is.null(object$cv)) {
    abline(v = object$cv$mstopopt, col = 2)
    axis(1, at = object$cv$mstopopt, labels = "mstopopt")
  }
}


logLik.crch.boost <- function(object, ...) structure(object$loglik, df = NA, class = "logLik")

