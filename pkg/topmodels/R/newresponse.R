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
  tmp <- procast(object)  ## TODO: Quick workaround
  left <- attr(tmp, "cens")$left
  right <- attr(tmp, "cens")$right
  if(!is.null(left)) y[y < left] <- left
  if(!is.null(right)) y[y > right] <- right

  attr(y, "weights") <- w
  return(y)
}
