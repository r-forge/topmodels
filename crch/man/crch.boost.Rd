\name{crch.boost}
\alias{crch.boost}
\alias{crch.boost.fit}

\title{Auxiliary functions to fit \code{crch} models via boosting.}

\description{
  Auxiliary functions to fit \code{crch} models via boosting
}
\usage{
crch.boost(maxit = 100, nu = 0.1, start = NULL, dot = "separate", 
  mstopopt = c("max", "aic", "bic", "cv"),  nfolds = 10, foldid = NULL, 
  reltol = sqrt(.Machine$double.eps)

crch.boost.fit(x, z, y, left, right, truncated = FALSE, dist = "gaussian",
  df = NULL, link.scale = "log", weights = NULL, offset = NULL, 
  control = crch.boost())
}

\arguments{
  \item{maxit}{the maximum number of boosting iterations.}
  \item{nu}{boosting step size. Default is 0.1.}
  \item{start}{a previously boosted but not converged \code{"crch.boost"} 
    object to continue.}
  \item{dot}{character specifying how to process formula parts with a dot
    (\code{.}) on the right-hand side. This can either be \code{"separate"}
    so that each formula part is expanded separately or \code{"sequential"}
    so that the parts are expanded sequentially conditional on all prior parts.
    Default is \code{"separate"}}
  \item{mstop}{method to find optimum stopping iteration. Default is \code{"max"}
    which is \code{maxit} or the boosting iteration at which the algorithm
    converges. Alternatives are \code{"aic"} and \code{"bic"} for AIC and BIC
    optimization and \code{"cv"} for cross validation.
    \code{mstop} can also be an integer to set the number of boosting
    iterations. Then \code{maxit} is set to \code{mstop} and \code{mstop="max"}.}
  \item{nfolds}{if \code{mstopopt = "cv"}, number of folds in cross validation.}
  \item{foldid}{if \code{mstopopt = "cv"}, an optional vector of values between 1 
    and \code{nfold} identifying what fold each observation is in. If supplied,
    \code{nfold} can be missing.}
  \item{reltol}{Relative convergence tolerance. The boosting stops if
    it is unable to reduce the negative log-likelihood by a factor of 
    \code{reltol * (- loglik + reltol)} at a step.  Default is 
    \code{sqrt(.Machine$double.eps)}. }
  \item{x, z, y, left, right, truncated, dist, df, link.scale, weights, offset, control}{see \code{\link{crch.fit}} for details.}
}

\details{
\code{crch.boost.fit} extends \code{crch} to fit censored (tobit) or 
truncated regression models with conditional heteroscedasticy with
boosting. If \code{crch.boost} is supplied as \code{control} in 
\code{crch} then \code{crch.boost.fit} is used as lower level fitting
function.
}

\value{
  For \code{crch.boost}: A list with components named as the arguments.
  For \code{crch.boost.fit}: An object of class \code{"crch.boost"}, 
    i.e., a list with the following elements.
  \item{coefficients}{list of coefficients for location and scale. Scale
  coefficients are in log-scale. Coefficients are of optimum stopping
  stopping iteration specified by \code{mstop}.}
  \item{df}{if \code{dist = "student"}: degrees of freedom of student-t
  distribution. else \code{NULL}.}
  \item{residuals}{the residuals, that is response minus fitted values.}
  \item{fitted.values}{list of fitted location and scale parameters at
  optimum stopping iteration specified by \code{mstop}.}
  \item{dist}{assumed distribution for the dependent variable \code{y}.}
  \item{cens}{list of censoring points.}
  \item{control}{list of control parameters.}
  \item{weights}{case weights used for fitting.}
  \item{offset}{list of offsets for location and scale.}
  \item{n}{number of observations.}
  \item{nobs}{number of observations with non-zero weights.}
  \item{loglik}{log-likelihood.}
  \item{link}{a list with element \code{"scale"} containing
          the link objects for the scale model.}
  \item{truncated}{logical indicating wheter a truncated model has been 
  fitted.}
  \item{converged}{logical variable whether boosting has converged or not.}
  \item{iterations}{number of boosting iterations.}
  \item{stepsize}{boosting stepsize \code{nu}.}
  \item{mstop}{criterion used to find optimum stopping iteration.}
  \item{mstopopt}{optimum stopping iterations for different criteria.}
  \item{standardize}{list of center and scale values used to standardize 
  response and regressors.}
}

\seealso{\code{\link{crch}}, \code{\link{crch.control}}}

\examples{
data("RainIbk")
RainIbk <- sqrt(RainIbk)
RainIbk$ensmean <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, mean)
RainIbk$enssd <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, sd)
RainIbk <- subset(RainIbk, enssd > 0)
boost <- crch(rain ~ ensmean + enssd | ensmean + enssd, RainIbk,
  left = 0, dist = "logistic", method = "boosting", mstop = "aic", maxit = 1000)
summary(boost)
plot(boost)
}

\keyword{regression}
