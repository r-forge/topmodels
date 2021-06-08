#' Extract Observed Responses from New Data
#' 
#' Generic function and methods for computing (randomized) quantile residuals.
#' 
#' This is a convenience function to accompany \code{\link{procast}}, e.g., for
#' computing predicted probabilities at the observed responses (also known as
#' the probability integral transform, PIT).
#' 
#' @aliases newresponse newresponse.default newresponse.glm
#' @param object a fitted model object. For the \code{default} method this
#' needs to needs to be \code{formula}-based so that
#' \code{\link[stats]{model.frame}} can be used to extract the response from
#' the original data the model was fitted to or \code{\link[stats]{terms}} can
#' be used to set up the response on \code{newdata}.
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param na.action function determining what should be done with missing
#' values in \code{newdata}.  The default is to employ \code{NA}.
#' @param \dots further parameters passed to methods.
#' @return A vector of new responses.
#' @seealso \code{\link[stats]{terms}}, \code{\link[stats]{model.frame}}
#' @keywords regression
#' @examples
#' 
#' ## linear regression models (homoscedastic Gaussian response)
#' m <- lm(dist ~ speed, data = cars)
#' newresponse(m)
#' newresponse(m, newdata = cars[1:3, ])
#' 
#' @export
newresponse <- function(object, ...) {
  UseMethod("newresponse")
}

#' @rdname newresponse
#' @method newresponse default
#' @export
newresponse.default <- function(object,
                                newdata,
                                na.action = na.pass,
                                ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## Sanity checks
  ## * `newaction` w/i `model.frame()`
  stopifnot(
    "`object` must be a (model) OO object:\n  * you have supplied a base object" =
      is.object(object),
    !all(class(object) == "list"),
    "`object` must be a (model) OO object:\n  * you have supplied an object of the single class `data.frame`" =
      !all(class(object) == "data.frame"),
    "`object` must be a (model) OO object:\n  * you have supplied a numeric" =
      !is.numeric(object) # TODO: (ML) is this case already caught by `is.object()`?
  )

  # -------------------------------------------------------------------
  # PREPARE RESPONSE
  # -------------------------------------------------------------------
  ## Get model.frame
  ## FIXME: (Z) Use expand.model.frame() instead of this hand-crafted code
  if (missing(newdata) || is.null(newdata)) {
    mf <- model.frame(object)
  } else {
    mf <- model.frame(update(terms(object), . ~ 1), newdata, na.action = na.action, ...)
  }

  ## Get weights
  ## FIXME: (ML)
  ## * Is this so correct and should we do that (needed for rootogram)?
  ## * What about na.action? What should the weights be for NAs (compare implementation below)?
  weights <- model.weights(mf)
  if (is.null(weights)) weights <- rep(1, NROW(mf))

  ## Get response
  y <- model.response(mf)

  ## Set weights to NA for NAs in response
  weights[is.na(y)] <- NA

  # -------------------------------------------------------------------
  # RETURN
  # -------------------------------------------------------------------
  ## FIXME: (ML)
  ## * Not very nice to return attributes here
  ## * How to guess properly if response is continous?!
  attr(y, "response_type") <- "continuous"
  attr(y, "weights") <- weights
  return(y)
}


#' @rdname newresponse
#' @method newresponse glm
#' @export
newresponse.glm <- function(object,
                            newdata,
                            na.action = na.pass,
                            ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## Sanity checks
  ## * `newaction` w/i `model.frame()`

  # -------------------------------------------------------------------
  # PREPARE RESPONSE
  # -------------------------------------------------------------------
  ## Get model.frame
  ## FIXME: (Z) use expand.model.frame() instead of this hand-crafted code
  if (missing(newdata) || is.null(newdata)) {
    mf <- model.frame(object)
  } else {
    mf <- model.frame(update(terms(object), . ~ 1), newdata, na.action = na.action, ...)
  }

  ## Get weights
  ## FIXME: (ML)
  ## * Is this so correct and should we do that (needed for rootogram)?
  ## * What about na.action? What should the weights be for NAs (compare implementation below)?
  weights <- model.weights(mf)
  if (is.null(weights)) weights <- rep(1, NROW(mf))

  ## Get response
  y <- model.response(mf)

  ## Set weights to NA for NAs in response
  weights[is.na(y)] <- NA

  ## Initialize family to get nobs and n
  nobs <- nobs(object)
  etastart <- NULL
  mustart <- NULL
  n <- NULL
  eval(object$family$initialize)

  ## Get response type
  ## FIXME: (ML) This could/should be improved
  if (object$family$family %in% "gaussian") {
    response_type <- "continuous"
  } else {
    response_type <- "discrete"
  }

  # -------------------------------------------------------------------
  # RETURN
  # -------------------------------------------------------------------
  attr(y, "nobs") <- nobs
  attr(y, "n") <- n
  attr(y, "weights") <- weights
  attr(y, "response_type") <- response_type
  return(y)
}
