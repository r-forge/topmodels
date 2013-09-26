summary.hxlr <-
function(object, ...)
{
  ## extend coefficient table
  k <- length(object$coefficients$intercept)
  l <- length(object$coefficients$predictor)
  m <- length(object$coefficients$scale)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if(length(object$coefficients$scale)) {
    cf <- list(intercept = cf[seq.int(length.out = k), , drop = FALSE], predictor = cf[seq.int(length.out = l) + k, , drop = FALSE], scale = cf[seq.int(length.out = m) + k + l, , drop = FALSE])
    rownames(cf$scale) <- names(object$coefficients$scale)
  } else {
    cf <- list(intercept = cf[seq.int(length.out = k), , drop = FALSE], predictor = cf[seq.int(length.out = l) + k, , drop = FALSE])
  }
  rownames(cf$predictor) <- names(object$coefficients$predictor)
  rownames(cf$intercept) <- names(object$coefficients$intercept)
  object$coefficients <- cf

  ## info
  object$info <- cbind(nobs = object$nobs, logLik = object$loglik, AIC = AIC(object), niter = object$optim$iterations)
  rownames(object$info) = ""

  ## return
  class(object) <- "summary.hxlr"
  object
}
