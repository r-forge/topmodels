\name{coef.crch}
\encoding{latin1}
\alias{coef.crch}
\alias{vcov.crch}
\alias{print.crch}
\alias{summary.crch}
\alias{print.summary.crch}
\alias{logLik.crch}
\alias{terms.crch}
\alias{residuals.crch}
\alias{fitted.crch}

\title{Methods for CRCH Objects}

\description{
Methods for extracting information from fitted \code{crch} objects.
}

\usage{
\method{coef}{crch}(object, model = c("full", "location", "scale", "df"), \dots)
\method{vcov}{crch}(object, model = c("full", "location", "scale", "df"), \dots)
\method{terms}{crch}(x, model = c("location", "scale", "full"), 
        \dots)
\method{fitted}{crch}(object, type = c("location", "scale"), 
        \dots)
}

\arguments{
  \item{object, x}{an object of class \code{"crch"}.}
  \item{model}{model for which coefficients shall be returned.}
  \item{type}{type of fitted values.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"crch"}, including methods to the generic functions
  \code{\link[base]{print}} and \code{\link[base]{summary}} which print the estimated
  coefficients along with some further information. The \code{summary} in particular
  supplies partial Wald tests based on the coefficients and the covariance matrix.
  As usual, the \code{summary} method returns an object of class \code{"summary.crch"}
  containing the relevant summary statistics which can subsequently be printed
  using the associated \code{print} method. 

  A \code{\link[stats]{logLik}} method is provided, hence \code{\link[stats]{AIC}}
  can be called to compute information criteria.
}

%\value{
%}

\seealso{\code{\link{crch}}}

\keyword{regression}
