\name{qresiduals}
\alias{qresiduals}
\alias{qresiduals.default}

\title{(Randomized) Quantile Residuals}

\description{
  Generic function and methods for computing (randomized) quantile residuals.
}

\usage{
qresiduals(object, \dots)

\method{qresiduals}{default}(object, newdata = NULL, trafo = qnorm,
  type = c("random", "quantile"), nsim = 1L, delta = NULL, prob = 0.5, \dots)
}

\arguments{
  \item{object}{an object. For the \code{default} method this needs to
    be either a specification of probabilities (vector or 2-dimensional matrix of
    probabilities) or an object from which the these can be obtained with
    \code{\link{procast}}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{trafo}{function for tranforming residuals from probability scale to
    a different distribution scale (default: Gaussian).}
  \item{type}{character specifying whether - in the case of discrete response distributions - 
    randomized quantile residuals or their corresponding quantiles should be computed.}
  \item{nsim}{numeric. The number of simulated randomized quantile residuals
    per observation (for \code{type = "numeric"}).}
  \item{delta}{numeric. The minimal difference to compute the range of proabilities
    corresponding to each observation according to get (randomized) quantile residuals. 
    For \code{NULL}, the minimal observed difference in the resonse 
    divided by \code{5e-6} is used.}
  \item{prob}{numeric. The probabilities at which quantile residuals should be
    computed (for \code{type = "quantile"}), defaulting to the median.}
  \item{\dots}{further parameters passed to methods.}
}

\details{
(Randomized) quantile residuals have been suggested by Dunn and Smyth (1996).
For regression models with a continuous response distribution this simply computes
theoretical standard normal quantiles corresponding to the probability integral transform of the
fitted distribution. For discrete distributions, a random theoretical normal
quantile is drawn from the range of probabilities corresponding to each observation.
Hence, in \code{qqrplot} the default is to use \code{trafo = qnorm} but other
transformations can also be used, specifically using the uniform probability
scale (via \code{trafo = NULL} or equivalently \code{qunif} or \code{identity}).

The default \code{qresiduals} method can compute randomized quantile residuals
from a vector (which essentially just calls \code{\link[stats]{qnorm}}) or
a 2-column matrix of probabilities. The latter offers to either draw \code{"random"}
samples from the distribution or compute corresponding \code{"quantile"}s such as
the median etc.
}

\note{
Note that there is also a \code{\link[statmod]{qresiduals}} function in the \pkg{statmod}
package that is not generic and always returns a single random quantile residual.
}

\value{
  A vector or matrix of quantile residuals.
}

\references{
Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile Residuals.}
  \emph{Journal of Computational and Graphical Statistics}, \bold{5}, 1--10.
}

\seealso{\code{\link[stats]{qnorm}}, \code{\link{qqrplot}}}

\examples{
## linear regression models (homoscedastic Gaussian response)
m <- lm(dist ~ speed, data = cars)
qresiduals(m)
}

\keyword{regression}
