\name{wormplot}
\alias{wormplot}
\alias{wormplot.default}

\alias{plot.wormplot}
\alias{points.wormplot}
\alias{c.wormplot}
\alias{autoplot.wormplot}

\title{Worm Plots for Quantile Residuals}

\description{
  Visualize goodness of fit of regression models by Worm plots using
  quantile residuals.
}

\usage{
wormplot(object, \dots)

\method{wormplot}{default}(object, newdata = NULL, plot = TRUE, flavor = NULL,
  trafo = qnorm, nsim = 1L, delta = NULL, confint = TRUE, 
  confint_level = 0.95, confint_nsim = 250, confint_seed = 1, single_graph = FALSE, 
  xlab = "Theoretical quantiles", ylab = "Deviation",
  main = NULL, \dots)
}

\arguments{
  \item{object}{an object (for which a \code{\link{qresiduals}} method exists).}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{plot}{logical. Should the \code{plot} method be called to draw the
    computed Q-Q plot?}
  \item{flavor}{Should the rootogram be a \code{base} or \code{ggplot2} style
    graphic, accordingly the invisible return value is either a \code{data.frame}
    or a \code{tibble}. Either set \code{flavor} expicitly to "base" vs.
    "tidyverse", or it's chosen automatically conditional
    if the packages \code{ggplot2} and \code{dplyr} or \code{tibble} are loaded.}
  \item{trafo}{function for tranforming residuals from probability scale to
    a different distribution scale (default: Gaussian).}
  \item{nsim, delta}{arguments passed to \code{qresiduals}.}
  \item{confint}{logical or quantile specification. Should the range of
    quantiles of the randomized quantile residuals be visualized? If \code{TRUE},
    then \code{range = c(0.01, 0.99)} is used.}
  \item{confint_level}{numeric. The confidence level required.}
  \item{confint_nsim}{numeric. The number of simulated quantiles.}
  \item{confint_seed}{numeric. The seed to be set for calculating the confidence interval.}
  \item{single_graph}{logical. Should all computed extended reliability diagrams be plotted in a single
    graph?}
  \item{xlab, ylab, main, \dots}{graphical plotting parameters
    passed to \code{\link[graphics]{plot}} or \code{\link[graphics]{points}},
    respectively.}
}

\details{
FIXME: Q-Q residual draw quantile residuals (by default: transformed to standard
normal scale) against theoretical quantiles from the same distribution.
Alternatively, transformations to other distributions can also be used, specifically
using no transformation at all, i.e., remaining on the uniform scale
(via \code{trafo = NULL} or equivalently \code{qunif} or \code{identity}).

Additional options are offered for models with discrete responses where randomization
of quantiles is needed.
}

\value{
  An list is returned invisibly with:
  \item{normal}{the theoretical normal quantiles,}
  \item{residuals}{the empirical quantile residuals.}
}

\references{
FIXME: Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile Residuals.}
  \emph{Journal of Computational and Graphical Statistics}, \bold{5}, 1--10.
}

\seealso{\code{\link{qresiduals}}, \code{\link[stats]{qqnorm}}}

\examples{
data("CrabSatellites", package = "countreg")
CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]

m1 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
m2 <- glm(satellites ~ width + color, data = CrabSatellites2, family = binomial)
m3 <- lm(dist ~ speed, data = cars)

q1 <- wormplot(m1, nsim = 100, confint = TRUE)
q2 <- wormplot(m2, nsim = 100, confint = TRUE, plot = FALSE)
q3 <- wormplot(m3, nsim = 100, confint = TRUE, plot = FALSE)

plot(c(q1, q2), single_graph = FALSE, fill = c(1, 3), ref = c(2, 2))
points(q3, col = "lightblue")
}

\keyword{hplot}
