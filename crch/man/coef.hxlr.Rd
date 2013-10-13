\name{coef.hxlr}
\encoding{latin1}
\alias{coef.hxlr}
\alias{vcov.hxlr}
\alias{print.hxlr}
\alias{summary.hxlr}
\alias{print.summary.hxlr}
\alias{logLik.hxlr}
\alias{terms.hxlr}


\title{Methods for HXLR Objects}

\description{
Methods for extracting information from fitted \code{hxlr} 
objects.
}

\usage{
\method{coef}{hxlr}(object, model = c("full", "intercept", 
        "location", "scale"), type = c("CLM", "latent"), \dots)
\method{vcov}{hxlr}(object, model = c("full", "intercept", 
        "location", "scale"), \dots)
\method{terms}{hxlr}(x, model = c("full", "location", "scale"), 
        \dots)
}

\arguments{
  \item{object, x}{an object of class \code{"hxlr"}.}
  \item{model}{model for which coefficients shall be returned.}
  \item{type}{type of coefficients. Default are CLM type coefficients.
    For type \code{"latent"} coefficients are converted in 
    coefficients for location and scale of the latent distribution 
    (analog to \code{crch} models).}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"hxlr"}, including methods to the generic functions
  \code{\link[base]{print}} and \code{\link[base]{summary}} which print the estimated
  coefficients along with some further information. The \code{summary} in particular
  supplies partial Wald tests based on the coefficients and the covariance matrix.
  As usual, the \code{summary} method returns an object of class \code{"summary.hxlr"}
  containing the relevant summary statistics which can subsequently be printed
  using the associated \code{print} method. 

  A \code{\link[stats]{logLik}} method is provided, hence \code{\link[stats]{AIC}}
  can be called to compute information criteria.
}

%\value{
%}

\seealso{\code{\link{hxlr}}}

\keyword{regression}
