\name{crch}
\alias{crch}
\alias{crch.fit}

\title{Censored Regression with Conditional Heteroscedasticy}

\description{
  Fitting censored regression (tobit) models with conditional heteroscedasticy. 
}
\usage{
crch(formula, data, subset, na.action, weights, offset, 
  dist = c("student", "gaussian", "logistic"), df = NULL, 
  left = -Inf, right = Inf, control = crch.control(\dots), 
  model = TRUE, x = FALSE, y = FALSE, \dots)

crch.fit(x, z, y, left, right, dist = "student", df = NULL, 
  weights = NULL, offset = NULL, control = crch.control()) 
}
\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x | z} where
  \code{y} is the response and \code{x} and \code{z} are regressor variables
  for the location and the scale of the fitted distribution respectively.}
  \item{data}{an optional data frame containing the variables occurring in the
  formulas.}
  \item{subset}{an optional vector specifying a subset of observations to be
  used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
  contain \code{NA}s.}
  \item{weights}{optional case weights in fitting.}
  \item{offset}{optional numeric vector with \emph{a priori} known component to
  be included in the linear predictor for the location.  For \code{crch.fit}, 
  \code{offset} can also be a list of 2 offsets used for the location and scale
  respectively.}
  \item{dist}{assumed distribution for the dependent variable \code{y}.}
  \item{df}{optional degrees of freedom for \code{dist="student"}. If omitted
  the degrees of freedom are estimated.}
  \item{left}{left limit for the censored dependent variable \code{y}.  If set
  to \code{-Inf}, \code{y} is assumed not to be left-censored.}
  \item{right}{right limit for the censored dependent variable \code{y}.  If
  set to \code{Inf}, \code{y} is assumed not to be right-censored.}
  \item{control}{a list of control parameters passed to \code{\link{optim}}.
  Default is \code{\link{crch.control}()}}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
  included as a component of the returned value.}
  \item{x, y}{for \code{crch}: logical. If \code{TRUE} the model matrix and
  response vector used for fitting are returned as components of the returned
  value.

  for \code{crch.fit}: \code{x} is a design matrix with regressors for the
  location and \code{y} is a vector of observations.}
  \item{z}{a design matrix with regressors for the scale.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
  if it is not supplied directly.}
}

\details{
\code{crch} fits censored regression (tobit) models with conditional
heteroscedasticy with maximum likelihood estimation. Student-t, Gaussian, and
logistic distributions can be fitted to left- and/or right-censored responses.
Different regressors can be used to model the location and the scale of this
distribution. \code{crch.fit} is the lower level function where the actual
fitting takes place.
}

\value{
  An object of class \code{"crch"}, i.e., a list with the following elements.
  \item{coefficients}{list of coefficients for location, scale, and df. Scale
  and df coefficients are in log-scale.}
  \item{df}{if \code{dist = "student"}: degrees of freedom of student-t
  distribution. else \code{NULL}.}
  \item{residuals}{the residuals, that is response minus fitted values.}
  \item{fitted.values}{list of fitted location and scale parameters.}
  \item{dist}{assumed distribution for the dependent variable \code{y}.}
  \item{cens}{list of censoring points.}
  \item{optim}{output from optimization from \code{\link{optim}}.}
  \item{method}{optimization method used for \code{\link{optim}}}
  \item{control}{list of control parameters passed to \code{\link{optim}}}
  \item{start}{starting values of coefficients used in the optimization.}
  \item{weights}{case weights used for fitting.}
  \item{offset}{list of offsets for location and scale.}
  \item{n}{number of observations.}
  \item{nobs}{number of observations with non-zero weights.}
  \item{loglik}{log-likelihood.}
  \item{vcov}{covariance matrix.}
  \item{converged}{logical variable whether optimization has converged or not.}
  \item{iterations}{number of iterations in optimization.}
  \item{call}{function call.}
  \item{formula}{the formula supplied.}
  \item{terms}{the \code{terms} objects used.}
  \item{levels}{list of levels of the factors used in fitting for location and
  scale respectively.}
  \item{contrasts}{(where relevant) the contrasts used.}
  \item{y}{if requested, the response used.}
  \item{x}{if requested, the model matrix used.}
  \item{model}{if requested, the model frame used.}
}

\references{
  Messner JW, Zeileis A, Broecker J, Mayr GJ (2013). Probabilistic Wind Power
    Forecasts with an Inverse Power Curve Transformation and Censored Regression.
    \emph{Wind Energy}, in press. Preprint version available from
    \url{http://EconPapers.RePEc.org/RePEc:inn:wpaper:2013-01}.
}

\seealso{\code{\link{predict.crch}}, \code{\link{crch.control}}}

\examples{
data("RainIbk")
## mean and standard deviation of square root transformed ensemble forecasts
RainIbk$sqrtensmean <- 
  apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]), 1, mean)
RainIbk$sqrtenssd <- 
  apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]), 1, sd)

## fit linear regression model with Gaussian distribution 
CRCH <- crch(sqrt(rain) ~ sqrtensmean, data = RainIbk, dist = "gaussian")
## same as lm?
all.equal(coef(lm(sqrt(rain) ~ sqrtensmean, data = RainIbk)),
  head(coef(CRCH), -1))

## print
CRCH
## summary
summary(CRCH)

## left censored regression model with censoring point 0:
CRCH2 <- crch(sqrt(rain) ~ sqrtensmean, data = RainIbk, 
  dist = "gaussian", left = 0)

## left censored regression model with censoring point 0 and 
## conditional heteroscedasticy:
CRCH3 <- crch(sqrt(rain) ~ sqrtensmean|sqrtenssd, data = RainIbk, 
  dist = "gaussian",  left = 0)

## left censored regression model with censoring point 0 and 
## conditional heteroscedasticy with logistic distribution:
CRCH4 <- crch(sqrt(rain) ~ sqrtensmean|sqrtenssd, data = RainIbk, 
  dist = "logistic", left = 0)

## compare AIC 
AIC(CRCH, CRCH2, CRCH3, CRCH4)
}
\keyword{regression}
