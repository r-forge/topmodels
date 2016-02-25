\name{predict.crch.boost}
\alias{predict.crch.boost}

\title{Predicted/Fitted Values for boosted CRCH Fits}

\description{
Obtains various types of predictions for \code{crch.boost} models.
}

\usage{
\method{predict}{crch.boost}(object, newdata = NULL, mstop = NULL, type = c("response", 
  "location", "scale", "quantile"), na.action = na.pass, at = 0.5, left = NULL, 
  right = NULL, \dots)
}

\arguments{
  \item{object}{an object of class \code{"crch.boost"}.}
  \item{newdata}{an optional data frame in which to look for variables with
  which to predict.}
  \item{mstop}{stopping iteration. Can be either a character (\code{"max"}, 
  \code{"aic"}, \code{"bic"}, \code{"cv"}) or a numeric value. If not \code{NULL},
  newdata has to be supplied.}
  \item{type}{type of prediction: \code{"response"} is equivalent to
  \code{"location"} and returns the location of the predicted distribution.
  \code{"scale"} returns the scale of the predicted distribution and
  \code{"quantile"} returns a matrix of predicted quantiles with quantile
  probabilities \code{at}.}
  \item{na.action}{a function which indicates what should happen when the data
  contain \code{NA}s. Default is na.pass}
  \item{at}{vector of quantile probabilities used for \code{type = "quantile"}.}
  \item{left}{left censoring or truncation point. Only used for \code{type = "quantile"}.
  If NULL, censoring or truncation point is obtained from \code{object}.}
  \item{right}{right censoring or truncation point. Only used for \code{type = "quantile"}.
  If NULL, censoring or truncation point is obtained from \code{object}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
For type \code{"response"}, \code{"location"}, or \code{"scale"} a vector with
either the location or the scale of the predicted distribution.

For type \code{"quantile"} a matrix of predicted quantiles each column
corresponding to an element of \code{at}. 
}

\seealso{\code{\link{crch.boost}},\code{\link{predict.crch}}}

\keyword{regression}
