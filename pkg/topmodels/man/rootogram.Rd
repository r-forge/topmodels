\name{rootogram}
\alias{rootogram}

\alias{plot.rootogram}
\alias{autoplot.rootogram}

%\alias{+.rootogram}
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

\method{rootogram}{default}(object, newdata = NULL, 
  plot = TRUE,  flavor = NULL, style = c("hanging", "standing", "suspended"), 
  scale = c("sqrt", "raw"), breaks = NULL, width = NULL, 
  response_type = NULL, xlab = NULL, ylab = NULL, main = NULL, \dots)

\method{plot}{rootogram}(x, ref = TRUE, xlim = NULL, ylim = NULL, xlab = NULL, 
  ylab = NULL, main = NULL, border = "black", fill = adjustcolor("black", alpha.f = 0.2),
  col = 2, lwd = 2, pch = 19, lty = 1, type = NULL, axes = TRUE, box = FALSE, \dots)

\method{autoplot}{rootogram}(object, ref = TRUE,  xlab = NULL, ylab = NULL,
  main = NULL, border = "black", fill = "darkgray", colour = 2, size = 1.2, shape = 19,
  linetype = 1, type = NULL, \dots)
}

\arguments{
  \item{object}{an object from which an rootogram can be extracted 
    with \code{\link{procast}}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  %\item{na.action}{function determining what should be done with missing
  %  values in \code{newdata}.  The default is to employ \code{NA}.}
  \item{plot}{logical. Should the \code{plot} method be called to
    draw the computed extended reliability diagram?}
  \item{flavor}{Should the rootogram be a \code{base} or \code{ggplot2} style
    graphic, accordingly the invisible return value is either a \code{data.frame}
    or a \code{tibble}. Either set \code{flavor} expicitly to "base" vs.
    "tidyverse", or it's chosen automatically conditional
    if the packages \code{ggplot2} and \code{dplyr} or \code{tibble} are loaded.}
  \item{style}{character specifying the syle of rootogram (see below).}
  \item{scale}{character specifying whether raw frequencies or their
    square roots (default) should be drawn.}
  \item{breaks}{numeric. Breaks for the histogram intervals.}
  \item{width}{numeric. Widths of the histogram bars.}
  \item{response_type}{To set the default values for \code{breaks} and \code{widths}.
    Currently different defaults are available for "discrete" and "continous" responses,
    as well as for the special case of a "logseries" response distribution.}
  \item{xlab, ylab, main}{graphical parameters.}
  \item{x}{object of class \code{"rootogram"}.}
  \item{ref}{logical. Should a reference line be plotted?}
  \item{xlim, ylim, border, fill, col, lwd, pch, lty, type, axes, box}{graphical
    parameters. These may pertain either to the whole plot or just the
    histogram or just the fitted line.}
  \item{colour, size, shape, linetype}{graphical parameters passed to \code{geom_line}
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
## plots and output

## number of deaths by horsekicks in Prussian army (Von Bortkiewicz 1898)
deaths <- rep(0:4, c(109, 65, 22, 3, 1))

## fit glm model
m1 <- glm(deaths ~ 1, family = poisson)
rootogram(m1)

## inspect output (without plotting)
r1 <- rootogram(m1, plot = FALSE)
r1

## combine plots
plot(c(r1, r1), col = c(1, 2), ref = 4, lty = c(1, 2))


#-------------------------------------------------------------------------------

## different styles

## artificial data from negative binomial (mu = 3, theta = 2)
## and Poisson (mu = 3) distribution
set.seed(1090)
y <- rnbinom(100, mu = 3, size = 2)
x <- rpois(100, lambda = 3)

## glm method: fitted values via glm()
m2 <- glm(y ~ x, family = poisson)

## correctly specified Poisson model fit
par(mfrow = c(1, 3))
rootogram(m2, style = "standing",  ylim = c(-2.2, 4.8), main = "Standing")
rootogram(m2, style = "hanging",   ylim = c(-2.2, 4.8), main = "Hanging")
rootogram(m2, style = "suspended", ylim = c(-2.2, 4.8), main = "Suspended")
par(mfrow = c(1, 1))

#-------------------------------------------------------------------------------
## linear regression with normal/Gaussian response: anorexia data

data("anorexia", package = "MASS")

m3 <- glm(Postwt ~ Prewt + Treat + offset(Prewt), family = gaussian, data = anorexia)
rootogram(m3, ylim = c(-1, 4))
abline(h = c(-1, 1), col = 4, lty = 2, lwd = 2)
}

\keyword{hplot}
