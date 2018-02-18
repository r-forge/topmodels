\name{newresponse}

\alias{newresponse}
\alias{newresponse.default}

\title{Extract Observed Responses from New Data}

\description{
  Generic function and methods for computing (randomized) quantile residuals.
}

\usage{
newresponse(object, \dots)

\method{newresponse}{default}(object, newdata, na.action = na.pass, \dots)
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
  \item{\dots}{further parameters passed to methods.}
}

\details{
This is a convenience function to accompany \code{\link{p4}}, e.g.,
for computing predicted probabilities at the observed responses
(also known as the probability integral transform, PIT).
}

\value{
  A vector of new responses.
}

\seealso{\code{\link[stats]{terms}}, \code{\link[stats]{model.frame}}}

\examples{
## linear regression models (homoscedastic Gaussian response)
m <- lm(dist ~ speed, data = cars)
newresponse(m)
newresponse(m, newdata = cars[1:3, ])
}

\keyword{regression}
