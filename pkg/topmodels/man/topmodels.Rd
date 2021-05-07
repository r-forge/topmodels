\name{topmodels}

\alias{topmodels}

\title{Plotting Graphical Evaluation Tools for Probabilistic Models}

\description{
  Plotting methods for probabilistic (regression) model objects supported by \code{procast}.
}

\usage{
topmodels(object, newdata = NULL, na.action = na.pass, which = NULL, 
  ask = dev.interactive(), spar = TRUE, pages = NULL, ...)
}

\arguments{
  \item{object}{An object supported by \code{"procast"}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{na.action}{function determining what should be done with missing
    values in \code{newdata}.  The default is to employ \code{NA}.}
  \item{which}{Character or integer, selects the type of plot: \code{"rootogram"} 
    graphically compares (square roots) of empirical frequencies with fitted
    frequencies from a probability model, \code{"pithist"} compares empirical
    probabilities from fitted models with a uniform distribution,
    \code{"reliagram"} shows reliability diagram for assessing the reliability of a
    fitted probabilistic distributional forecast, \code{"qqrplot"} shows a
    quantile-quantile plot of quantile residuals, \code{"wormplot"} shows a worm
    plot again using quantile resiudals.}
  \item{ask}{For multiple plots, the user is asked to show the next plot.}
  \item{spar}{Should graphical parameters be set?}
  \item{pages}{For numeric one, all plots will be shown on a single page.}
  \item{\dots}{Arguments to be passed to \code{\link{rootogram}}, \code{\link{pithist}},
    \code{\link{reliagram}}, \code{\link{qqrplot}}, and \code{\link{wormplot}}.}
}

\details{
FIXME
}

\value{
FIXME
}

\examples{
m <- lm(dist ~ speed, data = cars)
topmodels(m)
}

\keyword{regression}
