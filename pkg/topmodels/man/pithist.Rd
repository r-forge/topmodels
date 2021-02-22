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

\method{pithist}{default}(object, newdata = NULL, type = c("random", "proportional"), nsim = 1L,
  breaks = NULL, plot = TRUE, xlim = c(0, 1), ylim = NULL,
  xlab = "PIT", ylab = if(freq) "Frequency" else "Density", main = NULL,
  border = "black", fill = "lightgray", col = "#B61A51",
  lwd = 2, lty = c(1, 2), freq = FALSE,
  confint = TRUE, confint_level = 0.95, confint_type = c("default", "agresti"), ...)
}

\arguments{
  \item{object}{an object from which probability integral transforms can be
    extracted with \code{\link{procast}}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{type}{character. In case of discrete distributions should
    the PITs be drawn randomly from the corresponding interval
    or distributed proportionally?}
  \item{nsim}{integer. If \code{type} is \code{"random"} how
    many simulated PITs should be drawn?}
  \item{breaks}{numeric. Breaks for the histogram intervals.}
  \item{plot}{logical. Should the \code{plot} method be called to
    draw the computed PIT histogram?}
  \item{xlab, ylab, main}{graphical parameters.}
  \item{xlim, ylim, border, fill, col, lwd, lty}{graphical
    parameters. These may pertain either to the whole plot or just the
    histogram or just the fitted line.}
  \item{confint}{logical. Should confident intervals be drawn?}
  \item{confint_level}{numeric. The confidence level required.}
  \item{confint_type}{character. Which type of confidence interval.}  #FIXME: (ML) Improve documentation.
  \item{freq, \dots}{further arguments passed to
    \code{\link[graphics]{hist}}.}
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
}

\seealso{\code{\link{procast}}, \code{\link[graphics]{hist}}}

\examples{
## linear regression models (homoscedastic Gaussian response)
m <- lm(dist ~ speed, data = cars)
plot(m)
}

\keyword{hplot}
