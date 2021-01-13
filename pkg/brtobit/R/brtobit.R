brtobit <- function(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE,
  control = brtobit_control(...), ...)
{
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf)

  ## sanity check
  if(length(Y) < 1) stop("empty model")
  n <- length(Y)

  ## call the actual workhorse: brtobit_fit()
  rval <- brtobit_fit(X, Y, control)

  ## further model information
  rval$call <- cl
  rval$terms <- mt
  rval$levels <- .getXlevels(mt, mf)
  rval$contrasts <- attr(X, "contrasts")
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- X
  class(rval) <- "brtobit"
  return(rval)
}

brtobit_control <- function(fsmaxit = 100, start = NULL, epsilon = 1e-08, type = "BR", ...)
{
  type <- match.arg(type, c("BR", "ML"))
  list(fsmaxit = fsmaxit, start = start, epsilon = epsilon, type = type, ...)
}

## bias reduction
brtobit_fit <- function(x, y, control = brtobit_control())
{
  ## basic data properties
  n <- NROW(x)
  zero <- y <= 0
  
  ## starting values by default via OLS
  if (is.null(control$start)) {
    m <- lm.fit(x, y)
    par <- c(m$coefficients, "(Variance)" = sum(m$residuals^2)/m$df.residual)
  } else {
    par <- control$start
  }
  
  ## extract standard quantities for parameter vector
  fitfun <- function(par) {
    beta <- par[1L:(length(par) - 1L)]
    sigma <- if(par[length(par)] >= 0) sqrt(par[length(par)]) else NaN
    mu <- drop(x %*% beta)
    F <- pnorm(mu / sigma)
    f <- dnorm(mu / sigma)
    list(beta = beta, sigma = sigma, mu = mu, F = F, f = f)
  }
  
  ## (negative) log-likelihood for given parameter vector
  loglikfun <- function(par, fit = NULL){
    if(is.null(fit)) fit <- fitfun(par)
    ll <- crch::dcnorm(y, mean = fit$mu, sd = fit$sigma, left = 0, log = TRUE)
    if(any(!is.finite(ll))) NaN else -sum(ll)
  }
  
  ## corresponding gradient (contributions)
  gradfun <- function(par, sum = TRUE, fit = NULL){
    if(is.null(fit)) fit <- fitfun(par)
    sc <- crch:::scnorm(y, mean = fit$mu, sd = fit$sigma, left = 0)
    sc <- -cbind(sc[, 1] * x, sc[, 2]/(2 * fit$sigma))
    colnames(sc) <- names(par)
    if (sum) colSums(sc) else sc
  }
  
  ## corresponding Hessian
  hessfun <- function(par, fit = NULL) {
    if(is.null(fit))
      fit <- fitfun(par)
    mu <- fit$mu
    sigma <- fit$sigma
    f <- fit$f
    F <- fit$F
    
    K <- f / (1 - F)^2
    sinv <- 1 / sigma
    z <- mu / sigma
    A <- f^2 / (1 - F)
    
    ## Observed information
    o1 <- rep(NA, n)
    o1[zero] <- (- sinv * (K * (f * sinv - (sinv^2) * (1 - F) * mu)))[zero]
    o1[!zero] <- - sinv^2
    o1 <- t(x) %*% (as.vector(o1) * x)
    
    o2 <- rep(NA, n)
    o2[zero] <- ((-1 / (2 * sigma^3)) * K * (sinv^2 * (1 - F) * mu^2 - (1 - F) - z * f))[zero]
    o2[!zero] <- (- sinv^4 * (y - mu))[!zero]
    o2 <- t(x) %*% as.vector(o2)
    
    o3 <- rep(NA, n)
    o3[zero] <- ((1 / (4 * sigma^5)) * K * ( sinv^2 * (1 - F) * mu^3 - 3 * (1 - F) * mu - (mu^2) * f / sigma))[zero]
    o3[!zero] <- (1 / (2 * sigma^4) - (sinv^6) * (y - mu)^2)[!zero]
    o3 <- sum(o3)
    
    hessian <- rbind(cbind(o1, o2), c(t(o2), o3))
    obsinfo <- -hessian
    
    ## Expected information
    a <- (-1 / sigma^2) * (z * f - A - F)
    b <- (1 / (2 * sigma^3)) * (f * z^2 + f - A * z)
    c <- (-1 / (4 * sigma^4)) * (f * z^3 + f * z - A * z^2 - 2 * F)
    
    e1 <- t(x) %*% (as.vector(a) * x)
    e2 <- t(as.vector(b) %*% x)
    e3 <- sum(c)
    e <- cbind(e1, e2)
    expinfo <- rbind(e, c(t(e2), e3))
    
    dimnames(hessian) <- dimnames(obsinfo) <- dimnames(expinfo) <- list(names(par), names(par))
    list(hessian = hessian, obsinfo = obsinfo, expinfo = expinfo)
  }
  
  biasfun <- function(par, fit = NULL){
    if(is.null(fit)) {
      fit <- fitfun(par)
    }
    estcoef <- fit$beta
    estsigma <- fit$sigma
    mu <- fit$mu
    phi <- fit$f
    Phi <- fit$F
    expinfo <- hessfun(par, fit = fit)$expinfo
    
    ## Term Q
    qa <- (-phi / estsigma^3) * (phi^2 / (1 - Phi)^2 - (1 / estsigma) * (phi / (1 - Phi)) * mu - 1 )
    QA <- lapply(1:length(estcoef), function(i) t(x) %*% ((as.vector(qa) * x[, i]) * x))
    qb <- ( (phi^2 * mu) / (2 * estsigma^5 * (1 - Phi)) * (phi / (1 - Phi) - mu / estsigma) - (mu * phi) / (2 * estsigma^5) )
    QB <-  t(x) %*% (as.vector(qb) * x)
    qc <- (1 / (2 * estsigma^4)) * (phi^2 / (1 - Phi)) * ((- 1 / estsigma^2) * mu^2 + 1 + (1 / estsigma) * ((mu * phi) / (1 - Phi))) + (1 / estsigma^4) * (Phi - (mu * phi) / estsigma)
    QC <- lapply(1:length(estcoef), function(i) t(x) %*% ((as.vector(qc) * x[, i])))
    qd <- (phi^2 * mu) / (4 * estsigma^6 * (1 - Phi)) * (mu^2 / estsigma^2 - 1 - (phi * mu) / (estsigma * (1 - Phi))) + (phi / (2 * estsigma^5)) * (1 + mu^2 / estsigma^2)
    QD <- t(x) %*% (as.vector(qd))
    qe <- (1 / estsigma^5) * ((mu / (4 * estsigma)) * (phi^2 / (1 - Phi)) * (mu^2 / estsigma^2 - 3 - (mu / estsigma) * phi / (1 - Phi)) + (phi * mu^2) / estsigma^2 + (3 * phi) / 2 )
    QE <-  lapply(1:length(estcoef), function(i) (t(x[, i]) %*% as.vector(qe)))
    qf <- (mu^2 / (8 * estsigma^8)) * (phi / (1 - Phi)) * (- (mu^2 * phi) / estsigma^2 + 3 * phi + (mu / estsigma) * (phi^2 / (1 - Phi)) ) + Phi / estsigma^6 - (3 * mu * phi) / (4 * estsigma^7) - (phi * mu^3) / (2 * estsigma^9)
    QF <- sum(qf)
    
    QT <- lapply(1:length(estcoef), function(i) -rbind(cbind(QA[[i]], QC[[i]]), c(t(QC[[i]]), QE[[i]])))
    QTplus1 <- cbind(QB, QD)
    QTplus1 <- -rbind(QTplus1, c(t(QD), QF))
    
    ## Term P 
    pg <- (- 1 / estsigma^3) * (phi^3 / (1 - Phi)^2) + (phi / estsigma^5) * (mu^2 + 2 * estsigma^2)
    PG <- lapply(1:length(estcoef), function(i) t(x) %*% ((as.vector(pg) * x[, i]) * x))
    ph <- ((mu * phi) / (2 * estsigma^5)) * (phi^2 / (1 - Phi)^2 - 2 - mu^2 / estsigma^2) + Phi / estsigma^4
    PH <-  t(x) %*% (as.vector(ph) * x)
    pi <- ((phi * mu) / (2 * estsigma^5)) * (phi^2 / (1 - Phi)^2 - 2) + Phi / estsigma^4 - (phi * mu^3) / (2 * estsigma^7)
    PI <-  lapply(1:length(estcoef), function(i) t(x) %*% ((as.vector(pi) * x[, i])))
    pj <- ((phi * mu^2) / (2 * estsigma^7)) * (- phi^2 / (2 * (1 - Phi)^2) + 1) + phi / (4 * estsigma^5) * (5 + mu^4 / estsigma^4)
    PJ <- t(x) %*% (as.vector(pj))
    pk <- (- mu^2 * phi) / (4 * estsigma^7) * (phi^2 / (1 - Phi)^2 - 2 - mu^2 / estsigma^2) + (5 * phi) / (4 * estsigma^5)
    PK <-  lapply(1:length(estcoef), function(i) (t(x[, i]) %*% as.vector(pk)))
    pl <- (1 / (8 * estsigma^9)) * mu^3 * phi * ( phi^2 / (1 - Phi)^2 - 2 - mu^2 / estsigma^2) + Phi / estsigma^6 - (9 * phi * mu) / (8 * estsigma^7)
    PL <- sum(pl)
    
    PT <- lapply(1:length(estcoef), function(i) rbind(cbind(PG[[i]], PI[[i]]), c(t(PI[[i]]), PK[[i]])))
    PTplus1 <- cbind(PH, PJ)
    PTplus1 <- rbind(PTplus1, c(t(PJ), PL))
    
    ## Bias correction factor
    Acoef <- lapply(1:length(estcoef), function(i) 0.5 * sum(diag(solve(expinfo) %*% (PT[[i]] + QT[[i]]))))
    Avar <- 0.5 * sum(diag(solve(expinfo) %*% (PTplus1 + QTplus1)))
    A <- cbind(c(unlist(Acoef), Avar))
    rownames(A) <- rownames(expinfo)
    BC <- -solve(expinfo) %*% A
    
    list(adjustment = A, bias = BC)
  }
    
  ## opt <- optim(par = par, fn = loglikfun, gr = gradfun, method = "BFGS",
  ##   control = list(maxit = 5000, reltol = .Machine$double.eps^(1/1.2)))
  ## par <- opt$par
  step <- .Machine$integer.max
  iter <- 0
  if (control$fsmaxit > 0) {
    for (iter in 1L:control$fsmaxit) {
      stepPrev <- step
      stepFactor <- 0
      testhalf <- TRUE
      while (testhalf & stepFactor < 11) {
        fiti <- fitfun(par)
        scores <- -gradfun(par, fit = fiti)
        expinfo <- hessfun(par, fit = fiti)$expinfo
        bias <- if(control$type == "BR") biasfun(par, fit = fiti)$bias else 0
        par <- par + 2^(-stepFactor) * (step <- solve(expinfo) %*% scores - bias)
        stepFactor <- stepFactor + 1
        testhalf <- drop(crossprod(stepPrev) < crossprod(step))
      }
      if (all(abs(step) < control$epsilon)) {
        break
      }
    }
    fit <- fitfun(par)
    scores <- gradfun(par, fit = fit)
    expinfo <- hessfun(par, fit = fit)$expinfo
    
    if(control$type == "BR") {
      bias <- biasfun(par)$bias
      scores <- scores - expinfo %*% bias
    } else {
      bias <- biasfun(par)$bias
    }
  }
  
  ## nice formatting
  par <- drop(par)
  bias <- drop(bias)
  grad <- drop(scores)
  vc <- solve(expinfo)
  rownames(vc) <- colnames(vc) <- names(bias) <- names(grad) <- names(par) <- c(colnames(x), "(Variance)")
  
  rval <- list(
    coefficients = par,
    bias = bias,
    vcov = vc,
    loglik = -loglikfun(par),
    df = length(par),
    nobs = n,
    grad = grad,
    control = control,
    iterations = iter,
    converged = iter <= control$fsmaxit
  )
  return(rval)
}


logLik.brtobit <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

coef.brtobit <- function(object, model = c("full", "location", "scale"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
    "location" = {
      cf[-length(cf)]
    },
    "scale" = {
      cf[length(cf)]
    },
    "full" = {
      cf
    }
  )
}

vcov.brtobit <- function(object, model = c("full", "location", "scale"), ...) {
  model <- match.arg(model)
  vc <- object$vcov
  switch(model,
    "location" = {
      vc[-nrow(vc), -nrow(vc), drop = FALSE]
    },
    "scale" = {
      vc[nrow(vc), nrow(vc), drop = FALSE]
    },
    "full" = {
      vc
    }
  )
}

print.brtobit <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat(sprintf("%s tobit model\n\n", if(x$control$type == "BR") "Bias-reduced" else "Maximum likelihood"))
  if(!x$converged) {
    cat("Model did not converge\n")
  } else {
    if(length(x$coefficients) > 1L) {
      cat("Coefficients:\n")
      print.default(format(x$coefficients[-length(x$coefficients)], digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients\n\n")
    }
    if(length(x$coefficients) > 0L) {
      cat("Variance: ")
      cat(format(as.numeric(x$coefficients[length(x$coefficients)]), digits = digits))
      cat("\n\n")
    } else {
      cat("No variance\n\n")
    }
    cat(paste("Log-likelihood: ", format(x$loglik, digits = digits), "\n", sep = ""))
    if(length(x$df)) {
      cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
    }
  }

  invisible(x)
}

model.frame.brtobit <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  NextMethod()
} 

model.matrix.brtobit <- function(object, ...) {
  rval <- if(!is.null(object$x)) object$x
    else model.matrix(object$terms, model.frame(object), contrasts = object$contrasts)
  return(rval)
}

fitted.brtobit <- function(object, type = c("location", "scale"), ...) {
  predict(object, type = match.arg(type, c("location", "scale")), ...)
}

predict.brtobit <- function(object, newdata = NULL,
  type = c("response", "location", "scale", "parameter", "probability", "quantile"),
  na.action = na.pass, at = 0.5, ...)
{
  ## types of prediction
  ## response/location are synonymous
  type <- match.arg(type, c("response", "location", "scale", "parameter", "probability", "quantile"))
  if(type == "location") type <- "response"

  ## obtain model.frame/model.matrix
  if(type != "scale") {
    if(is.null(newdata)) {
      X <- model.matrix(object)
    } else {
      mf <- model.frame(delete.response(object$terms), newdata, na.action = na.action, xlev = object$levels)
      X <- model.matrix(delete.response(object$terms), mf, contrasts = object$contrasts)
    }
    n <- NROW(X)
    location <- drop(X %*% object$coefficients[-length(object$coefficients)])
  } else {
    n <- if(is.null(newdata)) object$nobs else NROW(newdata)
  }
  if(type != "response") {
    scale <- object$coefficients[length(object$coefficients)]
    scale <- rep.int(scale, n) 
  }

  ## compute result
  rval <- switch(type,
    "response" = location,
    "scale" = scale,
    "parameter" = data.frame(location, scale),
    "probability" = pnorm(at, mean = location, sd = scale),
    "quantile" = pmax(0, qnorm(at, mean = location, sd = scale))
  )
  return(rval)
}

bread.brtobit <- function(x, ...) x$vcov * x$nobs

summary.brtobit <- function(object, ...)
{
  ## extend coefficient table
  cf <- object$coefficients
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients <- cf

  ## delete some slots
  object$terms <- object$levels <- object$contrasts <- NULL

  ## return
  class(object) <- "summary.brtobit"
  object
}

print.summary.brtobit <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if(NROW(x$coefficients) > 0L) {
      cat(paste("Coefficients:\n", sep = ""))
      printCoefmat(x$coefficients, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients\n")
    if(getOption("show.signif.stars") & any(x$coefficients[, 4L] < 0.1, na.rm = TRUE))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df\n")
    cat(paste("Number of iterations in", x$control$type, "optimization:", x$iterations, "\n"))
  }

  invisible(x)
}

getSummary.brtobit <- function(obj, alpha = 0.05, ...) {
  ## coefficient matrix and confidence interval
  cf <- coef(obj)
  se <- sqrt(diag(vcov(obj)))
  cf <- cbind(
    "est" = cf,
    "se" = se,
    "stat" = cf/se,
    "p" = 2 * pnorm(abs(cf/se), lower.tail = FALSE),
    "lwr" = cf + qnorm(alpha/2) * se,
    "upr" = cf + qnorm(1 - alpha/2) * se
  )

  ## further summary statistics
  sstat <- c(
    "numdf" = obj$df,
    "N" = obj$nobs,
    "logLik" = obj$loglik,
    "AIC" = AIC(obj),
    "BIC" = BIC(obj)
  )

  ## return everything
  return(list(
    coef = cf,
    sumstat = sstat,
    contrasts = obj$contrasts,
    xlevels = obj$xlevels,
    call = obj$call
  ))
}

## setSummaryTemplate("brtobit" = c(
##   "Log-likelihood" = "($logLik:f#)",
##   "AIC" = "($AIC:f#)",
##   "BIC" = "($BIC:f#)",
##   "N" = "($N:d)"
## ))

