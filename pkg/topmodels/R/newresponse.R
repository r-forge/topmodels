newresponse <- function(object, ...) {
  UseMethod("newresponse")
}

newresponse.default <- function(object, 
                                newdata, 
                                ## FIXME: (ML) 
                                ## * Does it work as it should? 
                                ## * Should response have no NAs when covariates have NAs. 
                                ## * With na.pass, response and forecast won't match anymore
                                ## * Compare rootogram()
                                na.action = na.pass, 
                                ...) {

  ## Sanity checks
  stopifnot(
    "`object` must be a (model) OO object:\n  * you have supplied a base object" = 
      is.object(object),
      !all(class(object) == "list"),
    "`object` must be a (model) OO object:\n  * you have supplied an object of the single class `data.frame`" = 
      !all(class(object) == "data.frame"),
     "`object` must be a (model) OO object:\n  * you have supplied a numeric" =  
      !is.numeric(object) # TODO: (ML) is this case already caught by `is.object()`?
  )

  ## Get model.frame
  ## FIXME: (Z) Use expand.model.frame() instead of this hand-crafted code
  if(missing(newdata) || is.null(newdata)) {
    mf <- model.frame(object)
  } else {
    ## FIXME: (ML) See FIXME for `na.action`: 
    ## * Commented code of Z: this produces full response also when NAs in covariates (for na.action = na.omit)
    ## * Uncommented code is similar to countreg
    #mf <- model.frame(update(terms(object), . ~ 1), newdata, na.action = na.action, ...)
    mf <- model.frame(terms(object), newdata, na.action = na.action)
  }
  
  ## Get weights
  ## FIXME: (ML) Is this so correct and should we do that (needed for rootogram)
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- rep(1, NROW(mf))

  ## Get response
  y <- model.response(mf)

  ## FIXME: (ML)
  ## * Not very nice to return attributes here
  ## * How to guess properly if response is continous?!
  attr(y, "response_type") <- "continuous"
  attr(y, "weights") <- weights
  return(y)
}


newresponse.glm <- function(object, 
                            newdata, 
                            ## FIXME: (ML) 
                            ## * Does it work as it should? 
                            ## * Should response have no NAs when covariates have NAs. 
                            ## * With na.pass, response and forecast won't match anymore
                            ## * Compare rootogram()
                            na.action = na.pass, 
                            ...) {

  ## Get model.frame
  ## FIXME: (Z) use expand.model.frame() instead of this hand-crafted code
  if(missing(newdata) || is.null(newdata)) {
    mf <- model.frame(object)
  } else {
    ## FIXME: (ML) See FIXME for `na.action`: 
    ## * Commented code of Z: this produces full response also when NAs in covariates (for na.action = na.omit)
    ## * Uncommented code is similar to countreg
    #mf <- model.frame(update(terms(object), . ~ 1), newdata, na.action = na.action, ...)
    mf <- model.frame(terms(object), newdata, na.action = na.action)
  }

  ## Get weights
  ## FIXME: (ML) Is this so correct and should we do that (needed for rootogram)
  ## FIXME: (ML) What is the difference to `weights()` (seems to work correct)
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- rep(1, NROW(mf))

  ## Get response
  y <- model.response(mf)

  ## Initialize family to get nobs and n
  nobs <- nobs(object)
  etastart <- NULL
  mustart <- NULL
  n <- NULL
  eval(object$family$initialize)

  ## Get response type
  if (object$family$family %in% "gaussian") {
    response_type <- "continuous"
  } else {
    response_type <- "discrete"
  }

  attr(y, "nobs") <- nobs
  attr(y, "n") <- n
  attr(y, "weights") <- weights
  attr(y, "response_type") <- response_type
  return(y)
}
