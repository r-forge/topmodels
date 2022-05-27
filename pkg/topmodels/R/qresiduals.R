#' (Randomized) Quantile Residuals
#' 
#' Generic function and methods for computing (randomized) quantile residuals.
#' 
#' (Randomized) quantile residuals are simply the theoretical standard normal
#' quantiles evaluated at the PIT residuals as suggested by Dunn and Smyth (1996).
#' For regression models with a continuous response distribution these are exact;
#' for discrete distributions, PIT residuals are drawn from the range of
#' probabilities corresponding to each observation, hence quantile residuals must
#' be random as well.
#' 
#' The default \code{qresiduals} method calls \code{\link{pitresiduals}} with
#' \code{scale} equal \code{"normal"}, as employed in normal Q-Q
#' plots (\code{\link{qqrplot}}). 
#' 
#' @aliases qresiduals qresiduals.default
#' @param object an object. For the \code{default} method this needs to be
#' either a specification of probabilities (vector or 2-dimensional matrix of
#' probabilities) or an object from which the these can be obtained with
#' \code{\link{procast}}.
#' @param scale On which scale should the PIT residuals be shown; on the probability scale 
#' (\code{"uniform"}) or on the normal scale (\code{"normal"}).
#' Here, for (randomized) quantile residuals, the
#' quantiles of the standard normal distribution are computed per default.
#' @param \dots further parameters passed to  \code{\link{pitresiduals}}.
#' @return A vector or matrix of quantile residuals.
#' @note Note that there is also a \code{\link[statmod]{qresiduals}} function
#' in the \pkg{statmod} package that is not generic and always returns a single
#' random quantile residual.
#' @seealso \code{\link{pitresiduals}}, \code{\link[stats]{qnorm}}, \code{\link{qqrplot}}
#' @references Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile
#' Residuals.} \emph{Journal of Computational and Graphical Statistics},
#' \bold{5}(3), 236--244. \doi{10.2307/1390802}
#' @keywords regression
#' @examples
#' 
#' ## linear regression models (homoscedastic Gaussian response)
#' m <- lm(dist ~ speed, data = cars)
#' qresiduals(m)
#' 
#' @export
qresiduals <- function(object, ...) {
  UseMethod("qresiduals")
}

#' @rdname qresiduals
#' @method qresiduals default
#' @export
qresiduals.default <- function(object, 
                               scale = "normal",
                               ...) {
  pitresiduals(object, scale = scale, ...)
}
