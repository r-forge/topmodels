newresponse <- function(object, ...) {
  UseMethod("newresponse")
}

newresponse.default <- function(object, newdata, na.action = na.pass, ...)
{
  ## FIXME: use expand.model.frame() instead of this hand-crafted code
  if(missing(newdata) || is.null(newdata)) {
    mf <- model.frame(object)
  } else {
    mf <- model.frame(update(terms(object), . ~ 1), newdata, na.action = na.action, ...)
  }
  
  ## TODO: (ML) Is this so correct and should we do that (needed for rootogram)
  w <- model.weights(mf)
  if(is.null(w)) w <- rep(1, NROW(mf))

  y <- model.response(mf)

  ## TODO: (ML) Should we really catch response values above/below censor points
  ## FIXME: (ML) Solve quick workaround, must work with with vector of new censor points
  ##   might work employing `procast()` which uses `predict()`
  tmp <- procast(object)  
  left <- attr(tmp, "cens")$left
  right <- attr(tmp, "cens")$right
  if(!is.null(left)) y[y < left] <- left
  if(!is.null(right)) y[y > right] <- right

  attr(y, "weights") <- w
  return(y)
}


newresponse.glm <- function(object, 
                            newdata, 
                            na.action = na.pass, 
                            ...) {

  ## FIXME: use expand.model.frame() instead of this hand-crafted code
  if(missing(newdata) || is.null(newdata)) {
    mf <- model.frame(object)
  } else {
    mf <- model.frame(update(terms(object), . ~ 1), newdata, na.action = na.action, ...)
  }

  y <- model.response(mf)

  ## TODO: (ML) What is the difference to `weigths()`
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- rep(1, NROW(mf))

  nobs <- nobs(object)
  etastart <- NULL
  mustart <- NULL
  n <- NULL
  eval(object$family$initialize)

  attr(y, "nobs") <- nobs
  attr(y, "etastart") <- etastart
  attr(y, "mustart") <- mustart
  attr(y, "n") <- n
  attr(y, "weights") <- weights
  return(y)
}
