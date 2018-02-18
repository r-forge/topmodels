\name{qqrplot}
\alias{qqrplot}

\alias{plot.qqrplot}
\alias{autoplot.qqrplot}

\title{Q-Q Plots for Quantile Residuals}

\description{
  Visualize goodness of fit of regression models by Q-Q plots using
  quantile residuals.
}

\usage{
qqrplot(object, newdata = NULL, trafo = qnorm, type = c("random", "quantile"),
  nsim = 1L, prob = 0.5, range = FALSE, diag = TRUE,
  col = "black", fill = "lightgray", xlim = NULL, ylim = NULL,
  main = "Q-Q residuals plot", xlab = "Theoretical quantiles",
  ylab = "Quantile residuals", \dots)
}

\arguments{
  \item{object}{an object (for which a \code{\link{qresiduals}} method exists).}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{trafo}{function for tranforming residuals from probability scale to
    a different distribution scale (default: Gaussian).}
  \item{type, nsim, prob}{arguments passed to \code{qresiduals}.}
  \item{range}{logical or quantile specification. Should the range of
    quantiles of the randomized quantile residuals be visualized? If \code{TRUE},
    then \code{range = c(0.01, 0.99)} is used.}
  \item{diag}{logical or color specification. Should a diagonal
    reference line be drawn?}
  \item{col, fill, xlim, ylim, main, xlab, ylab, \dots}{graphical plotting parameters
    passed to \code{\link[graphics]{plot}} or \code{\link[graphics]{points}},
    respectively.}
}

\details{
Q-Q residual draw quantile residuals (by default: transformed to standard
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
Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile Residuals.}
  \emph{Journal of Computational and Graphical Statistics}, \bold{5}, 1--10.
}

\seealso{\code{\link{qresiduals}}, \code{\link[stats]{qqnorm}}}

\examples{
## linear regression models (homoscedastic Gaussian response)
m <- lm(dist ~ speed, data = cars)
qqrplot(m)
qqrplot(m, trafo = NULL)
}

\keyword{hplot}
