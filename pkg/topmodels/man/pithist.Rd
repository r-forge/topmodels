\name{pithist}
\alias{pithist}
\alias{pithist.default}

\alias{plot.pithist}
\alias{lines.pithist}
\alias{autoplot.pithist}

\title{PIT Histograms for Assessing Goodness of Fit of Probability Models}

\description{
  PIT histograms graphically compare empirical probabilities from
  fitted models with a uniform distribution.
}

\usage{
pithist(object, \dots)

\method{pithist}{default}(object, newdata = NULL, plot = TRUE, 
  style = c("histogram", "lines"), type = c("random", "proportional"), nsim = 1L,
  freq = FALSE, breaks = NULL, xlim = c(0, 1), ylim = NULL, main = NULL,
  xlab = "PIT", ylab = if(freq) "Frequency" else "Density", 
  confint = TRUE, confint_level = 0.95, confint_type = c("exact", "approximation"), ...)
}

\arguments{
  \item{object}{an object from which probability integral transforms can be
    extracted with \code{\link{procast}}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{plot}{logical. Should the \code{plot} method be called to
    draw the computed PIT histogram?}
  \item{style}{character specifying the syle of rootogram (see below). FIXME: Description}
  \item{type}{character. In case of discrete distributions should
    the PITs be drawn randomly from the corresponding interval
    or distributed proportionally?}
  \item{nsim}{integer. If \code{type} is \code{"random"} how
    many simulated PITs should be drawn?}
  \item{freq}{logical. If \code{TRUE}, the PIT histogram is represented by
    frequencies, the \code{counts} component of the result; if
    \code{FALSE}, probability densities, component \code{density}, are
    plotted (so that the histogram has a total area of one).}
  \item{breaks}{numeric. Breaks for the histogram intervals.}
  \item{xlim, ylim}{graphical parameters. These may pertain either to the 
    whole plot or just the histogram or just the fitted line.}
  \item{xlab, ylab, main}{graphical parameters.}
  \item{confint}{logical. Should confident intervals be drawn?}
  \item{confint_level}{numeric. The confidence level required.}
  \item{confint_type}{character. Which type of confidence interval. 
    According to Agresti and Coull (1998) for interval estimation of binomial 
    proportions an approximation can be better than exact.}
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
## linear regression models (homoscedastic Gaussian response)
m1 <- lm(dist ~ speed, data = cars)
pit1 <- pithist(m1)

## linear regression models (heteroscedastic Gaussian response)
m2 <- crch(dist ~ speed | speed, data = cars)

## censored linear regression models (heteroscedastic Gaussian response)
m3 <- crch(dist ~ speed | speed, left = 30, data = cars)

pit2 <- pithist(m2, plot = FALSE)
pit3 <- pithist(m3, plot = FALSE)

plot(pit1, style = "lines")
lines(c(pit2, pit3), col = c(2, 3))
}

\keyword{hplot}
