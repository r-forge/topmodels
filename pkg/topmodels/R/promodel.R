#' Predictions and Residuals Dispatch for Probabilistic Models
#' 
#' The function \code{promodel} is a wrapper for dispatching the base
#' \code{\link[stats]{predict}} and \code{\link[stats]{residuals}} methods to
#' the \code{\link[topmodels]{procast}} and \code{\link[topmodels]{proresiduals}}
#' functions for probabilistic forecasts and probabilistic residuals, respectively.
#' 
#' The default methods for \code{\link[topmodels]{procast}} and
#' \code{\link[topmodels]{proresiduals}} in this package make a wide range of
#' different probabilistic forecasts and probabilistic residuals
#' available for many fitted model object classes. However, it may sometimes be
#' useful to call these flexible methods via the base
#' \code{\link[stats]{predict}} and \code{\link[stats]{residuals}} methods.
#' For example, this may be useful in combination with other packages that
#' rely on the base functions such as \pkg{marginaleffects}.
#'
#' Therefore, the \code{promodel} wrapper function simply adds an additional
#' class \code{"promodel"} (probabilistic model) to the original class of an
#' object. Then the methods for \code{predict} and \code{residuals} then
#' strip off this class again before calling \code{procast} and \code{proresiduals},
#' respectively.
#'
#' @aliases promodel
#' @param object a fitted model object for which \code{\link[topmodels]{procast}}
#'   and/or \code{\link[topmodels]{proresiduals}} work.
#' @param \dots further arguments passed on to \code{\link[topmodels]{procast}}
#'   or \code{\link[topmodels]{proresiduals}}, respectively.
#'
#' @keywords regression
#' @examples
#' ## Poisson regression model for FIFA 2018 data:
#' ## number of goals scored by each team in each game, explained by
#' ## predicted ability difference of the competing teams
#' data("FIFA2018", package = "distributions3")
#' m <- glm(goals ~ difference, data = FIFA2018, family = poisson)
#' 
#' ## prediction using a new data set (final of the tournament)
#' final <- tail(FIFA2018, 2)
#'
#' ## base predict method computes linear predictor on link scale (here in logs)
#' predict(m, newdata = final)
#'
#' ## procast-based method computes distribution object by default
#' pm <- promodel(m)
#' predict(pm, newdata = final)
#'
#' ## all other procast types are available as well
#' predict(pm, newdata = final, type = "density", at = 0:4)
#' predict(pm, newdata = final, type = "cdf", at = 0:4)
#'
#' ## the base residuals method defaults to deviance residuals
#' ## but the proresiduals-based method defaults to quantile residuals
#' head(residuals(m))
#' head(residuals(pm))

#' @export 
promodel <- function(object) {
  class(object) <- c("promodel", class(object))
  return(object)
}

#' @rdname promodel
#' @export 
residuals.promodel <- function(object, ...) {
  class(object) <- setdiff(class(object), "promodel")
  proresiduals(object, ...)
}

#' @rdname promodel
#' @export 
predict.promodel <- function(object, ...) {
  class(object) <- setdiff(class(object), "promodel")
  procast(object, ...)
}
