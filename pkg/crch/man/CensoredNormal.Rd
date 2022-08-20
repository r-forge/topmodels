\name{CensoredNormal}
\alias{CensoredNormal}
\alias{mean.CensoredNormal}
\alias{variance.CensoredNormal}
\alias{skewness.CensoredNormal}
\alias{kurtosis.CensoredNormal}
\alias{pdf.CensoredNormal}
\alias{log_pdf.CensoredNormal}
\alias{cdf.CensoredNormal}
\alias{quantile.CensoredNormal}
\alias{random.CensoredNormal}
\alias{crps.CensoredNormal}
\alias{support.CensoredNormal}

\title{Create a Censored Normal Distribution}

\description{
Class and methods for left-, right-, and interval-censored normal
distributions using the workflow from the \pkg{distributions3} package.
}

\usage{
CensoredNormal(mu = 0, sigma = 1, left = -Inf, right = Inf)
}

\arguments{
\item{mu}{numeric. The location parameter of the underlying \emph{uncensored}
  normal distribution, typically written \eqn{\mu} in textbooks.
  Can be any real number, defaults to \code{0}.}
\item{sigma}{numeric. The scale parameter (standard deviation) of
  the underlying \emph{uncensored} normal distribution,
  typically written \eqn{\sigma} in textbooks.
  Can be any positive number, defaults to \code{1}.}
\item{left}{numeric. The left censoring point. Can be any real number,
  defaults to \code{-Inf} (uncensored). If set to a finite value, the
  distribution has a point mass at \code{left} whose probability corresponds
  to the cumulative probability function of the \emph{uncensored}
  normal distribution at this point.}
\item{right}{numeric. The right censoring point. Can be any real number,
  defaults to \code{Inf} (uncensored). If set to a finite value, the
  distribution has a point mass at \code{right} whose probability corresponds
  to 1 minus the cumulative probability function of the \emph{uncensored}
  normal distribution at this point.}
}

\details{
The constructor function \code{CensoredNormal} sets up a distribution
object, representing the censored normal probability distribution by the
corresponding parameters: the latent mean \code{mu} = \eqn{\mu} and
latent standard deviation \code{sigma} = \eqn{\sigma} (i.e., the parameters
of the underlying \emph{uncensored} normal variable), the \code{left} censoring
point (with \code{-Inf} corresponding to uncensored), and the
\code{right} censoring point (with \code{Inf} corresponding to uncensored). 

The censored normal distribution has probability density function (PDF) \eqn{f(x)}:

\tabular{ll}{
\eqn{\Phi((left - \mu)/\sigma)} \tab if \eqn{x \le left}\cr
\eqn{1 - \Phi((right - \mu)/\sigma)} \tab if \eqn{x \ge right}\cr
\eqn{\phi((x - \mu)/\sigma)/\sigma} \tab if \eqn{left < x < right}
}

where \eqn{\Phi} and \eqn{\phi} are the cumulative distribution function
and probability density function of the standard normal distribution
respectively.

All parameters can also be vectors, so that it is possible to define a vector
of censored normal distributions with potentially different parameters.
All parameters need to have the same length or must be scalars (i.e.,
of length 1) which are then recycled to the length of the other parameters.

For the \code{CensoredNormal} distribution objects there is a wide range
of standard methods available to the generics provided in the \pkg{distributions3}
package: \code{\link[distributions3]{pdf}} and \code{\link[distributions3]{log_pdf}}
for the (log-)density (PDF), \code{\link[distributions3]{cdf}} for the probability
from the cumulative distribution function (CDF), \code{quantile} for quantiles,
\code{\link[distributions3]{random}} for simulating random variables,
\code{\link[scoringRules]{crps}} for the continuous ranked probability score
(CRPS), and \code{\link[distributions3]{support}} for the support interval
(minimum and maximum). Internally, these methods rely on the usual d/p/q/r
functions provided for the censored normal distributions in the \pkg{crch}
package, see \code{\link{dcnorm}}, and the \code{\link[scoringRules]{crps_cnorm}}
function from the \pkg{scoringRules} package.

See the examples below for an illustration of the workflow for the class and methods.
}

\value{
A \code{CensoredNormal} distribution object.
}

\seealso{\code{\link{dcnorm}}, \code{\link[distributions3]{Normal}}, \code{\link{TruncatedNormal}},
  \code{\link{CensoredLogistic}}, \code{\link{CensoredStudentsT}}}

\examples{
\dontshow{ if(!requireNamespace("distributions3")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }
}
## package and random seed
library("distributions3")
set.seed(6020)

## three censored normal distributions:
## - uncensored standard normal
## - left-censored at zero (Tobit) with latent mu = 1 and sigma = 1
## - interval-censored in [0, 5] with latent mu = 1 and sigma = 2
X <- CensoredNormal(
  mu    = c(   0,   1, 1),
  sigma = c(   1,   1, 2),
  left  = c(-Inf,   0, 0),
  right = c( Inf, Inf, 5)
)
X

## compute mean of the censored distribution
mean(X)
## higher moments (variance, skewness, kurtosis) are not implemented yet

## support interval (minimum and maximum)
support(X)

## simulate random variables
random(X, 5)

## histograms of 1,000 simulated observations
x <- random(X, 1000)
hist(x[1, ], main = "uncensored")
hist(x[2, ], main = "left-censored at zero")
hist(x[3, ], main = "interval-censored in [0, 5]")

## probability density function (PDF) and log-density (or log-likelihood)
x <- c(0, 0, 1)
pdf(X, x)
pdf(X, x, log = TRUE)
log_pdf(X, x)

## cumulative distribution function (CDF)
cdf(X, x)

## quantiles
quantile(X, 0.5)

## cdf() and quantile() are inverses (except at censoring points)
cdf(X, quantile(X, 0.5))
quantile(X, cdf(X, 1))

## all methods above can either be applied elementwise or for
## all combinations of X and x, if length(X) = length(x),
## also the result can be assured to be a matrix via drop = FALSE
p <- c(0.05, 0.5, 0.95)
quantile(X, p, elementwise = FALSE)
quantile(X, p, elementwise = TRUE)
quantile(X, p, elementwise = TRUE, drop = FALSE)

## compare theoretical and empirical mean from 1,000 simulated observations
cbind(
  "theoretical" = mean(X),
  "empirical" = rowMeans(random(X, 1000))
)

## evaluate continuous ranked probability score (CRPS) using scoringRules
library("scoringRules")
crps(X, x)
}
