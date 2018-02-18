\name{p4}

\alias{p4}
\alias{p4.lm}
\alias{p4_setup}

\title{Probabilistic Forecasting}

\description{
  Generic function and methods for computing various kinds of forecasts
  from probabilistic (regression) models (probabilistic 4-casting).
}

\usage{
p4(object, newdata = NULL, na.action = na.pass, type = "quantile", at = 0.5, \dots)

\method{p4}{lm}(object, newdata = NULL, na.action = na.pass,
  type = c("quantile", "mean", "variance", "parameter", "density", "probability", "score"),
  at = 0.5, drop = FALSE, \dots)

p4_setup(pars, FUN, at = NULL, drop = FALSE, type = "p4", \dots)
}

\arguments{
  \item{object}{a fitted model object. For the \code{default} method this needs to
    needs to be \code{formula}-based so that \code{\link[stats]{model.frame}}
    can be used to extract the response from the original data the model
    was fitted to or \code{\link[stats]{terms}} can be used to set up
    the response on \code{newdata}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{na.action}{function determining what should be done with missing
    values in \code{newdata}.  The default is to employ \code{NA}.}

  \item{type}{character specifying the type of probabilistic forecast to compute.
    In \code{p4_setup} the \code{type} is only used for nice labels of the
    returned data frame.}
  \item{at}{specification of values at which the forecasts should be evaluated,
    typically a numeric vector but possibly also a matrix or data frame.
    Additionally, \code{at} can be the character string \code{"function"}
    or \code{"list"}, see details below.}
  \item{drop}{logical. Should forecasts be returned in a data frame (default) or
    (if possible) dropped to a vector, see details below.}
  \item{FUN}{function to be used for forecasts. Either of type \code{FUN(pars, \dots)}
    or \code{FUN(at, pars, \dots)}, see details below.}
  \item{pars}{a data frame of predicted distribution parameters.}
  \item{\dots}{further parameters passed to methods.}
}

\details{
The function \code{p4} provides a unified framework for probabilistic
4-casting based on probabilistic (regression) models, also known
as distributional regression approaches. Typical types of predictions
include quanitles, probabilities, (conditional) expectations, variances,
(log-)densities, or scores. Internally, \code{p4} methods typically
compute the predicted parameters for each observation and then transform
these to the desired outcome. Some quantities (e.g., expectations or
variances) can be computed directly from the predicted parameters of
the distribution while others require an additional argument \code{at}
which the distribution is evaluated (e.g., the probability of a quantile
or an observation of the response. The argument \code{at} can also be
the character \code{"function"} or \code{"list"} so that a single function
or list of functions is set up that can be evaluated \code{at} different
values later on.

The function \code{p4_setup} is a convenience wrapper that make setting
up \code{p4} methods easier for package developers. It takes a data frame
of predicted parameters \code{pars} and a function \code{FUN} which is to
be evaluated at the parameters. This can either have the interface
\code{FUN(pars, \dots)} when the desired quantity can be predicted directly
from the predicted parameters -- or the interface \code{FUN(at, pars, \dots)}
if an additional argument \code{at} is needed. \code{p4_setup} takes care
of suitable expanding \code{at} to the dimensions of \code{pars} and
optionally setting up a (list of) function(s) to be returned.
}

\value{
  Either a \code{data.frame} of predictions (in case of multivariate forecasts,
  or if \code{drop = FALSE}, default) or a vector (in case of a univariate forecast
  and additionally \code{drop = TRUE}). Unless \code{at} is the character string
  \code{"function"} or \code{"list"} in which case a (list of) function(s) is
  returned.
}

\examples{
## linear regression models (homoscedastic Gaussian response)
m <- lm(dist ~ speed, data = cars)

## medians on observed data
p4(m)
p4(m, drop = TRUE)

## probability integral transform (PIT) on observed data
p4(m, type = "probability", at = cars$dist)

## log-likelihood contributions
p4(m, type = "density", at = cars$dist, log = TRUE)

## log-likelihood sum
sum(p4(m, type = "density", at = cars$dist, log = TRUE))
logLik(m)


## medians on new data
nd <- data.frame(speed = c(10, 15, 20))
p4(m, newdata = nd)

## different quantile for each observation
p4(m, newdata = nd, at = c(0.25, 0.5, 0.75))

## all combinations of quantiles and observations
p4(m, newdata = nd, at = rbind(c(0.25, 0.5, 0.75)))

## function for computing quantiles (vectorized)
qnt1 <- p4(m, newdata = nd, at = "function")
## as before
qnt1(0.5)
qnt1(c(0.25, 0.5, 0.75))
qnt1(rbind(c(0.25, 0.5, 0.75)))

## list of functions
qnt2 <- p4(m, newdata = nd, at = "list")
qnt2[[1]]
qnt2[[1]](0.5)
qnt2[[1]](c(0.25, 0.5, 0.75))

}

\keyword{regression}
