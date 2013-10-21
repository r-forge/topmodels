\name{predict.crch}
\alias{predict.crch}

\title{Predicted/Fitted Values for CRCH Fits}

\description{
Obtains various types of predictions for \code{crch} models.
}

\usage{
\method{predict}{crch}(object, newdata = NULL, type = c("response", "location",
  "scale", "quantile"), na.action = na.pass, at = 0.5, \dots)
}

\arguments{
  \item{object}{an object of class \code{"crch"}.}
  \item{newdata}{an optional data frame in which to look for variables with
  which to predict.}
  \item{type}{type of prediction: \code{"response"} is equivalent to
  \code{"location"} and returns the location of the predicted distribution.
  \code{"scale"} returns the scale of the predicted distribution and
  \code{"quantile"} returns a matrix of predicted quantiles with quantile
  probabilities \code{at}.}
  \item{na.action}{a function which indicates what should happen when the data
  contain \code{NA}s. Default is na.pass}
  \item{at}{vector of quantile probabilities used for \code{type =
  "quantile"}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
For type \code{"response"}, \code{"location"}, or \code{"scale"} a vector with
either the location or the scale of the predicted distribution.

For type \code{"quantile"} a matrix of predicted quantiles each column
corresponding to an element of \code{at}. 
}

\seealso{\code{\link{crch}}}

\keyword{regression}
