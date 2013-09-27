\name{predict.hxlr}
\alias{predict.hxlr}

\title{Predictions for Heteroskedastic Extended Logistic Regressions}

\description{
Computes various types of predictions for heteroskedastic extended logistic
regression (HXLR) models.
}

\usage{
\method{predict}{hxlr}(object, newdata = NULL,
  type = c("class", "probability", "cumprob", "location", "scale"),
  thresholds = object$thresholds, na.action = na.pass, \dots)
}

\arguments{
  \item{object}{An object of class \code{"hxlr"}.}
  \item{newdata}{An optional data frame in which to look for variables
    with which to predict.}
  \item{type}{Type of prediction: \code{"probability"} returns data frame with
    category probabilities, \code{"cumprob"} returns cumulative probabilities,
    \code{"location"} and \code{"scale"} return the location and scale of the
    predicted latent distribution respectively, and \code{"class"} returns the
    category with the highest probability. Default is \code{"class"}.}
  \item{thresholds}{Optional thresholds used for defining the thresholds
    (for types \code{"probability"}, \code{"cumprob"}, and \code{"class"}).
    Can be different from thresholds used for fitting. If omitted, the same
    thresholds as for fitting are used.}
  \item{na.action}{A function which indicates what should happen when the data contain ‘NA’s. Default is na.pass}
  \item{\dots}{Currently not used.}
}

\value{
For type \code{"prob"} a matrix with number of intervals (= number of thresholds
+ 1) columns is produced. Each row corresponds to a row in newdata and contains
the predicted probabilities to fall in the corresponding interval.

For type \code{"prob"} a matrix with number of thresholds is produced. Each row
corresponds to a row in \code{newdata} and contains the predicted probabilities to fall
below the corresponding threshold.

For types \code{"class"}, \code{"location"}, and \code{"scale"} a vector is returned. 
}

\seealso{\code{\link{hxlr}}}

\keyword{regression}
