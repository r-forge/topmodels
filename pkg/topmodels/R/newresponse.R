#' Extract Observed Responses from New Data
#' 
#' Generic function and methods for extracting response variables from new data
#' based on fitted model objects.
#' 
#' \code{newresponse} is a convenience function that supports functions like
#' \code{\link{proscore}} or \code{\link{proresiduals}} which assess
#' discrepancies between predictions/forecasts on new data and the corresponding
#' observed response variables.
#'
#' The default method takes an approach that is similar to many \code{\link[stats]{predict}}
#' methods which rebuild the \code{\link[stats]{model.frame}} after dropping the
#' response from the \code{\link[stats]{terms}} of a model object. However, here
#' only the response variable is preserved and all explanatory variables are dropped.
#' Missing values values are typically preserved (i.e., using \code{\link[stats]{na.pass}}).
#'
#' If the new \code{model.frame} contains a variable \code{"(weights)"},
#' it is preserved along with the response variable(s).
#' 
#' A method for \code{distribution} objects is provided which expects that
#' \code{newdata} is essentially already the corresponding new response.
#' Thus, it needs to be a vector (or data frame) of the same length as \code{distribution}.
#' If it is not a data frame, yet, it is transformed to one but no further
#' modifications are made.
#' 
#' @param object a fitted model object. For the \code{default} method this
#' needs to needs to be \code{formula}-based so that
#' \code{\link[stats]{model.frame}} can be used to extract the response from
#' the original data the model was fitted to or \code{\link[stats]{terms}} can
#' be used to set up the response on \code{newdata}.
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param na.action function determining how to handle missing
#' values in \code{newdata}, by default these are preserved.
#' @param \dots further arguments passed to methods.
#' @param initialize logical. Should the response variable from \code{\link[stats]{glm}}
#' objects be initialized using the corresponding expression from the \code{family}?
#' If \code{NULL} (the default), the initialization is only used for \code{\link[stats]{binomial}}
#' and \code{\link[stats]{quasibinomial}} families.
#' 
#' @return A \code{data.frame} (\code{model.frame}) containing the response variable
#' (and optionally a variable with \code{"(weights)"}).
#' 
#' @seealso \code{\link[stats]{terms}}, \code{\link[stats]{model.frame}}
#' 
#' @keywords regression
#' 
#' @examples
#' ## linear regression model
#' m <- lm(dist ~ speed, data = cars)
#' 
#' ## extract response variable on data used for model fitting 
#' newresponse(m)
#' 
#' ## extract response variable on "new" data
#' newresponse(m, newdata = cars[1:3, ])
#' 
#' @export
newresponse <- function(object, ...) {
  UseMethod("newresponse")
}

#' @rdname newresponse
#' @method newresponse default
#' @export
newresponse.default <- function(object, newdata, na.action = na.pass, ...) {
  ## extract model.frame (after dropping explanatory variables)
  if (missing(newdata) || is.null(newdata)) {
    mf <- model.frame(object, ...)
    mt <- try(terms(mf), silent = TRUE)
    if(inherits(mt, "try-error") || !inherits(mt, "terms")) mt <- try(terms(object), silent = TRUE)
    if(inherits(mt, "try-error") || !inherits(mt, "terms")) stop("terms() cannot be extracted from object")
  } else {
    mt <- try(terms(object), silent = TRUE)
    if(inherits(mt, "try-error") || !inherits(mt, "terms")) stop("terms() cannot be extracted from object")
    mt <- update(mt, . ~ 1)
    mf <- model.frame(mt, newdata, na.action = na.action, ...)
    mt <- terms(mf)
  }

  ## response column
  y <- names(mf)[attr(mt, "response")]
  
  ## weights (if any)
  if("(weights)" %in% names(mf)) y <- c(y, "(weights)")
  
  ## return model.frame subset
  y <- mf[, y, drop = FALSE]
  return(y)
}


#' @rdname newresponse
#' @method newresponse glm
#' @export
newresponse.glm <- function(object, newdata, na.action = na.pass, initialize = NULL, ...) {
  newy <- if (missing(newdata) || is.null(newdata)) {
    newresponse.default(object, ...)
  } else {
    newresponse.default(object, newdata = newdata, na.action = na.action, ...)
  }

  ## initialize response via family
  binom <- !is.null(object$family) && inherits(object$family, "family") && object$family$family %in% c("binomial", "quasibinomial")
  if (is.null(initialize)) initialize <- binom
  if (initialize) {
    y <- newy[[1L]]
    weights <- newy[["(weights)"]]
    nobs <- nobs(object)
    etastart <- NULL
    mustart <- NULL
    n <- NULL
    eval(object$family$initialize)
    if (binom) y <- round(y * n) ## FIXME: "correct" transformation?
    newy[[1L]] <- y
  }
  
  return(newy)
}

#' @rdname newresponse
#' @method newresponse distribution
#' @export
newresponse.distribution <- function(object, newdata, ...) {
  if (missing(newdata) || is.null(newdata) || NROW(newdata) != length(object)) {
    stop("for a distribution 'object' the 'newdata' needs to be a vector (or data frame) of the same length")
  }
  if (!is.data.frame(newdata)) {
    newy <- data.frame(y = seq_along(object))
    names(newy) <- deparse(substitute(newdata))
    newy[[1L]] <- newdata
    return(newy)
  } else {
    return(newdata)
  }
}
