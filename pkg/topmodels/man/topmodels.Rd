\name{topmodels}

\alias{topmodels}

\title{Plotting Graphical Evaluation Tools for Probabilistic Models}

\description{
  Plotting methods for probabilistic (regression) model objects supported by \code{procast}.
}

\usage{
topmodels(object, plot = TRUE, class = NULL, newdata = NULL, na.action = na.pass, 
  which = NULL, ask = dev.interactive(), spar = TRUE, pages = NULL, ...)
}

\arguments{
  \item{object}{An object supported by \code{"procast"}.}
  \item{plot}{Should the \code{plot} or \code{autoplot} method be called to
    draw all chosen plots? Either set \code{plot} expicitly to
    "base" vs. "ggplot2" to choose the type of plot, or for a logical \code{plot} argument it's chosen
    conditional if the package \code{ggplot2} is loaded.}
  \item{class}{Should the invisible return value be either a \code{data.frame}
    or a \code{tibble}. Either set \code{class} expicitly to "data.frame" vs.
    "tibble", or for NULL it's chosen automatically conditional
    if the package \code{tibble} is loaded.}
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
data("CrabSatellites", package = "countreg")
CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]

m1 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
m2 <- glm(satellites ~ width + color, data = CrabSatellites2, family = binomial)

## base graphics
topmodels(m1, pages = 1, plot = "base")
topmodels(m1, pages = 1, nsim = 10, plot = "base")
topmodels(m1, pages = 1, nsim = 30, fill = 2, plot = "base")
topmodels(m1, pages = 1, nsim = 30, ref = 2, plot = "base")
topmodels(m1, pages = 1, nsim = 30, fill =2, add_hist = 2, plot = "base")

## ggplot2 graphics
topmodels(m1, pages = 1, nsim = 30, plot = "ggplot2")
topmodels(m2, pages = 1, nsim = 30, plot = "ggplot2")

}

\keyword{regression}
