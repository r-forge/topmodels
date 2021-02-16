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

  rval <- model.response(mf)
  attr(rval, "weights") <- w
  return(rval)

}
