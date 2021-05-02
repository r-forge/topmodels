\name{rootogram}
\alias{rootogram}

\alias{plot.rootogram}
\alias{autoplot.rootogram}

\alias{+.rootogram}
\alias{c.rootogram}
\alias{rbind.rootogram}

\alias{rootogram.default}
\alias{rootogram_glm}
\alias{rootogram_procast}

\title{Rootograms for Assessing Goodness of Fit of Probability Models}

\description{
  Rootograms graphically compare (square roots) of empirical frequencies
  with fitted frequencies from a probability model.
}

\usage{
rootogram(object, \dots)

\method{rootogram}{default}(object, newdata = NULL, plot = TRUE, 
  style = c("hanging", "standing", "suspended"), 
  scale = c("sqrt", "raw"), breaks = NULL, width = NULL, 
  xlab = NULL, ylab = NULL, main = NULL, \dots)

\method{plot}{rootogram}(x, style = c("hanging", "standing", "suspended"), 
  ref = TRUE, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, 
  main = NULL, border = "black", fill = adjustcolor("black", alpha.f = 0.2), 
  col = "red", lwd = 2, pch = 19, lty = 1, type = NULL, axes = TRUE, \dots)

\method{autoplot}{rootogram}(object, colour = c("black", "#B61A51"),
  fill = "darkgray", size = c(1.2, 4), \dots)
}

\arguments{
  \item{object}{an object from which an rootogram can be extracted 
    with \code{\link{procast}}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{plot}{logical. Should the \code{plot} method be called to
    draw the computed extended reliability diagram?}
  \item{style}{character specifying the syle of rootogram (see below).}
  \item{scale}{character specifying whether raw frequencies or their
    square roots (default) should be drawn.}
  \item{breaks}{numeric. Breaks for the histogram intervals.}
  \item{width}{numeric. Widths of the histogram bars.}
  \item{xlab, ylab, main}{graphical parameters.}
  \item{x}{object of class \code{"rootogram"}.}
  \item{ref}{logical. Should a reference line be plotted?}
  \item{xlim, ylim, border, fill, col, lwd, pch, lty, type, axes}{graphical
    parameters. These may pertain either to the whole plot or just the
    histogram or just the fitted line.}
  \item{colour, size}{graphical parameters passed to \code{geom_line}
    and \code{geom_point}, respectively.}
  \item{\dots}{further graphical parameters passed to the plotting
    function.}
}

\details{
Rootograms graphically compare frequencies of empirical
distributions and fitted probability models. For the observed distribution
the histogram is drawn on a square root scale (hence the name) and superimposed
with a line for the fitted frequencies. The histogram can be \code{"standing"}
on the x-axis (as usual), or \code{"hanging"} from the fitted curve, or
a \code{"suspended"} histogram of deviations can be drawn.

\code{rootogram} is the generic function for generating rootograms from data
or fitted model objects. The workhorse function is the default method (that
computes all necessary coordinates based on observed and fitted frequencies and
the breaks for the histogram intervals) and the associating \code{plot} method
that carries out the actual drawing (using base graphics).

There is a wide range of further \code{rootogram} methods that all take
the following approach: based on a fitted probability model observed and
expected frequencies are computed and then the default method is called.
Currently, there is a method for \code{\link[stats]{glm}}.

Furthermore, there is a \code{numeric} method that uses \code{link[MASS]{fitdistr}}
to obtain a fitted (by maximum likelihood) probability model for a univariate variable. 
For this method, \code{fitted} can either be a character string or a density function
that is passed to \code{fitdistr}. In the latter case, a \code{start} list also
has to be supplied.

In addition to the \code{plot} method for rootogram objects, there are also two
methods that combine two (or more) rootograms: \code{c}/\code{rbind} creates
a set of rootograms that can then be plotted in one go. The \code{+} method adds
up the observed and fitted frequencies from two rootograms (if these use the same
intervals).

The \code{\link[ggplot2]{autoplot}} method creates a \code{\link[ggplot2]{ggplot}}
version of the rootogram.
}

\note{
Note that there is also a \code{\link[vcd]{rootogram}} function in the \pkg{vcd}
package that is similar to the \code{numeric} method provided here. However, it
is much more limited in scope, hence a function has been created here.
}

\value{
  An object of class \code{"rootogram"} inheriting from \code{"data.frame"} with the
  following variables:
  \item{observed}{observed frequencies,}
  \item{expected}{fitted frequencies,}
  \item{x}{histogram interval midpoints on the x-axis,}
  \item{y}{bottom coordinate of the histogram bars,}
  \item{width}{widths of the histogram bars,}
  \item{height}{height of the histogram bars,}
  \item{line}{y-coordinates of the fitted curve.}
  Additionally, \code{style}, \code{scale}, \code{xlab}, \code{ylab}, and
  \code{main} are stored as attributes.
}

\references{
Friendly M (2000), \emph{Visualizing Categorical Data}. SAS Institute, Cary.

Kleiber C, Zeileis A (2016).
  \dQuote{Visualizing Count Data Regressions Using Rootograms.} 
  \emph{The American Statistician}, \bold{70}(3), 296--303.
  \doi{10.1080/00031305.2016.1173590}.

Tukey JW (1977). \emph{Exploratory Data Analysis}. Addison-Wesley, Reading.
}

\seealso{\code{\link[stats]{glm}}}

\examples{
require("crch")
m1 <- glm(dist ~ speed, data = cars)
m2 <- crch(dist ~ speed | speed, data = cars)
plot(c(rootogram(m1, plot = FALSE), rootogram(m2, plot = FALSE)), col = c(1, 2))
}

\keyword{hplot}
