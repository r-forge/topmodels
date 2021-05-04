newresponse <- function(object, ...) {
  UseMethod("newresponse")
}

newresponse.default <- function(object, newdata, na.action = na.pass, ...)
{

  stopifnot(
    "`object` must be a (model) OO object:\n  * you have supplied a base object" = 
      is.object(object),
      !all(class(object) == "list"),
    "`object` must be a (model) OO object:\n  * you have supplied an object of the single class `data.frame`" = 
      !all(class(object) == "data.frame"),
     "`object` must be a (model) OO object:\n  * you have supplied a numeric" =  
      !is.numeric(object) # TODO: (ML) is this case already caught by `is.object()`?
  )

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

  ## TODO: (ML) What is the difference to `weigths()` (seems to work correct)
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- rep(1, NROW(mf))

  nobs <- nobs(object)
  etastart <- NULL
  mustart <- NULL
  n <- NULL
  eval(object$family$initialize)

  attr(y, "nobs") <- nobs
  attr(y, "n") <- n
  attr(y, "weights") <- weights
  return(y)
}
