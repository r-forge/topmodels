\name{pithist}
\alias{pithist}
\alias{pithist.default}

\alias{plot.pithist}
\alias{lines.pithist}
\alias{c.pithist}
\alias{autoplot.pithist}

\title{PIT Histograms for Assessing Goodness of Fit of Probability Models}

\description{
  PIT histograms graphically compare empirical probabilities from
  fitted models with a uniform distribution.
}

\usage{
pithist(object, \dots)

\method{pithist}{default}(object, newdata = NULL, plot = TRUE, flavor = NULL,
  style = c("histogram", "lines"), type = c("random", "proportional"), nsim = 1L, 
  delta = NULL, freq = FALSE, breaks = NULL, confint = TRUE, 
  confint_level = 0.95, confint_type = c("exact", "approximation"),
  single_graph = FALSE, xlim = c(0, 1), ylim = c(0, NA),
  xlab = "PIT", ylab = if (freq) "Frequency" else "Density", main = NULL, ...)
}

\arguments{
  \item{object}{an object from which probability integral transforms can be
    extracted with \code{\link{procast}}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{plot}{logical. Should the \code{plot} method be called to
    draw the computed PIT histogram?}
  \item{flavor}{Should the rootogram be a \code{base} or \code{ggplot2} style
    graphic, accordingly the invisible return value is either a \code{data.frame}
    or a \code{tibble}. Either set \code{flavor} expicitly to "base" vs.
    "tidyverse", or it's chosen automatically conditional
    if the packages \code{ggplot2} and \code{dplyr} or \code{tibble} are loaded.}
  \item{style}{character specifying the syle of rootogram (see below). FIXME: Description}
  \item{type}{character. In case of discrete distributions should
    the PITs be drawn randomly from the corresponding interval
    or distributed proportionally?}
  \item{nsim}{integer. If \code{type} is \code{"random"} how
    many simulated PITs should be drawn?}
  \item{delta}{numeric. The minimal difference to compute the range of proabilities
    corresponding to each observation according to get (randomized) quantile residuals. 
    For \code{NULL}, the minimal observed difference in the resonse 
    divided by \code{5e-6} is used.}
  \item{freq}{logical. If \code{TRUE}, the PIT histogram is represented by
    frequencies, the \code{counts} component of the result; if
    \code{FALSE}, probability densities, component \code{density}, are
    plotted (so that the histogram has a total area of one).}
  \item{breaks}{numeric. Breaks for the histogram intervals.}
  \item{confint}{logical. Should confident intervals be drawn?}
  \item{confint_level}{numeric. The confidence level required.}
  \item{confint_type}{character. Which type of confidence interval. 
    According to Agresti and Coull (1998) for interval estimation of binomial 
    proportions an approximation can be better than exact.}
  \item{single_graph}{logical. Should all computed extended reliability diagrams be plotted in a 
    single graph?}
  \item{xlim, ylim}{graphical parameters. These may pertain either to the 
    whole plot or just the histogram or just the fitted line.}
  \item{xlab, ylab, main}{graphical parameters.}
  \item{\dots}{further graphical parameters.}
}

\details{
PIT histograms graphically the probability integral transform (PIT),
i.e., observed probabilities from fitted probability models, with
a uniform distribution. It leverages the \code{\link{procast}} generic
and then essentially draws a \code{\link[graphics]{hist}}.

In case of discrete distributions the PIT is either drawn randomly
from the corresponding interval or distributed proportionally in the
histogram (FIXME: not yet implemented).
}

\references{
  Czado C, Gneiting T, Held L (2009). \dQuote{Predictive Model Assessment for Count Data.}
  \emph{Biometrics}, \bold{65}(4), 1254--1261.

  Agresti A, Coull A B (1998). \dQuote{Approximate is Better than ``Exact'' for Interval 
  Estimation of Binomial Proportions.} \emph{The American Statistician}, \bold{52}(2), 119--126.
}

\seealso{\code{\link{procast}}, \code{\link[graphics]{hist}}}

\examples{
require("crch")
m1 <- lm(dist ~ speed, data = cars)
m2 <- crch(dist ~ speed | speed, data = cars)
m3 <- crch(dist ~ speed | speed, left = 30, data = cars)

pit1 <- pithist(m1)
pit2 <- pithist(m2, plot = FALSE)
pit3 <- pithist(m3, plot = FALSE)

plot(pit1, confint = "red", ref = "blue", fill = "lightblue")

plot(c(pit1, pit2, pit3), col = c(1, 2, 3), style = "lines")

plot(c(pit1, pit2), col = c(1, 2), single_graph = TRUE)
lines(pit3, col = 3)
}

\keyword{hplot}
