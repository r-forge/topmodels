\name{reliagram}
\alias{reliagram}
\alias{reliagram.default}

\alias{plot.reliagram}
\alias{lines.reliagram}
\alias{c.reliagram}
\alias{autoplot.reliagram}

\title{Reliagram (Extended Reliability Diagram)}

\description{
  Reliagram (extended reliability diagram) for assessing the reliability
  of a fitted probabilistic distributional forecast.
}

\usage{
reliagram(object, \dots)

\method{reliagram}{default}(object, newdata = NULL, plot = TRUE,
  breaks = seq(0, 1, by = 0.1), probs = 0.5, thresholds = NULL,
  y = NULL, confint = TRUE, confint_level = 0.95, confint_nboot = 250,
  single_graph = FALSE, xlab = "Forecast probability", 
  ylab = "Observed relative frequency", main = NULL, ...) 
}

\arguments{
  \item{object}{an object from which an extended reliability diagram 
    can be extracted with \code{\link{procast}}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{plot}{logical. Should the \code{plot} method be called to
    draw the computed extended reliability diagram?}
  \item{breaks}{numeric vector passed on to \code{\link{cut}} in order to bin 
    the observations and the predicted probabilities.}
  \item{probs}{numeric vector of probabilities with values in [0,1] to calculate single
    or several thresholds. Only used if \code{thresholds} is not specified.}
  \item{thresholds}{numeric vector specifying both where to cut the observations into binary values
      and at which values the predicted probabilities should be calculated (\code{\link{procast}}).}
  \item{y}{numeric vector. Observed or out-of-sample observations. FIXME: (ML) Might be removed.}
  \item{confint}{logical. Should confident intervals be calculated and drawn?}
  \item{confint_level}{numeric. The confidence level required.}
  \item{confint_nboot}{numeric. The number of bootstrap steps.}
  \item{single_graph}{logical. Should all computed extended reliability diagrams be plotted in a single
    graph?}
  \item{xlab, ylab, main}{graphical parameters.}
  \item{\dots}{further graphical parameters.}
}


\details{
TODO: (ML)
}

\references{
  Br{\''o}cker J, Smith L (2007). \dQuote{Increasing the Reliability of Reliability Diagrams}. 
  \emph{Weather and Forecasting}, \bold{22}(3), 651--661. doi:10.1175/WAF993.1.
}

\seealso{\code{\link{procast}}}


\examples{
require("crch")
m1 <- lm(dist ~ speed, data = cars)
m2 <- crch(dist ~ speed | speed, data = cars)
m3 <- crch(dist ~ speed | speed, left = 30, data = cars)

rel1 <- reliagram(m1, probs = c(0.2, 0.6))
rel2 <- reliagram(m2, plot = FALSE)
rel3 <- reliagram(m3, plot = FALSE)

plot(c(rel1, rel2), single_graph = TRUE, col = c(1, 2, 3), lty = c(1, 2, 3), pch = c(1, 2, 3))

plot(rel1, single_graph = TRUE, col = c(1, 2), fill = c(1, 2))
lines(rel3, col = 3, lty = 2, confint = 3)
}

