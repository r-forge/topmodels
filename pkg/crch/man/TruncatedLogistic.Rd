\name{TruncatedLogistic}
\alias{TruncatedLogistic}
\alias{mean.TruncatedLogistic}
\alias{variance.TruncatedLogistic}
\alias{skewness.TruncatedLogistic}
\alias{kurtosis.TruncatedLogistic}
\alias{pdf.TruncatedLogistic}
\alias{log_pdf.TruncatedLogistic}
\alias{cdf.TruncatedLogistic}
\alias{quantile.TruncatedLogistic}
\alias{random.TruncatedLogistic}
\alias{crps.TruncatedLogistic}
\alias{support.TruncatedLogistic}

\title{Create a Truncated Logistic Distribution}

\description{
Class and methods for left-, right-, and interval-truncated logistic
distributions using the workflow from the \pkg{distributions3} package.
}

\usage{
TruncatedLogistic(location = 0, scale = 1, left = -Inf, right = Inf)
}

\arguments{
\item{location}{numeric. The location parameter of the underlying \emph{untruncated}
  logistic distribution, typically written \eqn{\mu} in textbooks.
  Can be any real number, defaults to \code{0}.}
\item{scale}{numeric. The scale parameter (standard deviation) of
  the underlying \emph{untruncated} logistic distribution,
  typically written \eqn{\sigma} in textbooks.
  Can be any positive number, defaults to \code{1}.}
\item{left}{numeric. The left truncation point. Can be any real number,
  defaults to \code{-Inf} (untruncated). If set to a finite value, the
  distribution has a point mass at \code{left} whose probability corresponds
  to the cumulative probability function of the \emph{untruncated}
  logistic distribution at this point.}
\item{right}{numeric. The right truncation point. Can be any real number,
  defaults to \code{Inf} (untruncated). If set to a finite value, the
  distribution has a point mass at \code{right} whose probability corresponds
  to 1 minus the cumulative probability function of the \emph{untruncated}
  logistic distribution at this point.}
}

\details{
The constructor function \code{TruncatedLogistic} sets up a distribution
object, representing the truncated logistic probability distribution by the
corresponding parameters: the latent mean \code{location} = \eqn{\mu} and
latent standard deviation \code{scale} = \eqn{\sigma} (i.e., the parameters
of the underlying \emph{untruncated} logistic variable), the \code{left} truncation
point (with \code{-Inf} corresponding to untruncated), and the
\code{right} truncation point (with \code{Inf} corresponding to untruncated). 

The truncated logistic distribution has probability density function (PDF):

\deqn{f(x) = 1/\sigma \lambda((x - \mu)/\sigma) / (\Lambda((right - \mu)/\sigma) - \Lambda((left - \mu)/\sigma))}
for \eqn{left \le x \le right}, and 0 otherwise,
where \eqn{\Lambda} and \eqn{\lambda} are the cumulative distribution function
and probability density function of the standard logistic distribution,
respectively.

All parameters can also be vectors, so that it is possible to define a vector
of truncated logistic distributions with potentially different parameters.
All parameters need to have the same length or must be scalars (i.e.,
of length 1) which are then recycled to the length of the other parameters.

For the \code{TruncatedLogistic} distribution objects there is a wide range
of standard methods available to the generics provided in the \pkg{distributions3}
package: \code{\link[distributions3]{pdf}} and \code{\link[distributions3]{log_pdf}}
for the (log-)density (PDF), \code{\link[distributions3]{cdf}} for the probability
from the cumulative distribution function (CDF), \code{quantile} for quantiles,
\code{\link[distributions3]{random}} for simulating random variables,
\code{\link[scoringRules]{crps}} for the continuous ranked probability score
(CRPS), and \code{\link[distributions3]{support}} for the support interval
(minimum and maximum). Internally, these methods rely on the usual d/p/q/r
functions provided for the truncated logistic distributions in the \pkg{crch}
package, see \code{\link{dtlogis}}, and the \code{\link[scoringRules]{crps_tlogis}}
function from the \pkg{scoringRules} package.

See the examples below for an illustration of the workflow for the class and methods.
}

\value{
A \code{TruncatedLogistic} distribution object.
}

\seealso{\code{\link{dtlogis}}, \code{\link[distributions3]{Logistic}}, \code{\link{CensoredLogistic}},
  \code{\link{TruncatedNormal}}, \code{\link{TruncatedStudentsT}}}

\examples{
\dontshow{ if(!requireNamespace("distributions3")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }
}
## package and random seed
library("distributions3")
set.seed(6020)

## three truncated logistic distributions:
## - untruncated standard logistic
## - left-truncated at zero with latent location = 1 and scale = 1
## - interval-truncated in [0, 5] with latent location = 2 and scale = 1
X <- TruncatedLogistic(
  location = c(   0,   1, 2),
  scale    = c(   1,   1, 1),
  left     = c(-Inf,   0, 0),
  right    = c( Inf, Inf, 5)
)
X

## compute mean of the truncated distribution
mean(X)
## higher moments (variance, skewness, kurtosis) are not implemented yet

## support interval (minimum and maximum)
support(X)

## simulate random variables
random(X, 5)

## histograms of 1,000 simulated observations
x <- random(X, 1000)
hist(x[1, ], main = "untruncated")
hist(x[2, ], main = "left-truncated at zero")
hist(x[3, ], main = "interval-truncated in [0, 5]")

## probability density function (PDF) and log-density (or log-likelihood)
x <- c(0, 0, 1)
pdf(X, x)
pdf(X, x, log = TRUE)
log_pdf(X, x)

## cumulative distribution function (CDF)
cdf(X, x)

## quantiles
quantile(X, 0.5)

## cdf() and quantile() are inverses (except at truncation points)
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
