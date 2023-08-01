# -------------------------------------------------------------------
# Programming outline: Wormplot
# -------------------------------------------------------------------
# - Showing the difference between the empirical quantile and
#   the unit normal quantile following Van Buuren and Fredriks (2001)
#
# Functions:
# - wormplot() generic plus default method
# - Return object of class "wormplot" that is plotted by default
# - But has plot=FALSE so that suitable methods can be added afterwards
# - At least methods: plot(), autoplot()
# -------------------------------------------------------------------


#' Worm Plots for Quantile Residuals
#' 
#' Visualize goodness of fit of regression models by worm plots using quantile
#' residuals. If \code{plot = TRUE}, the resulting object of class
#' \code{"wormplot"} is plotted by \code{\link{plot.qqrplot}} or
#' \code{\link{autoplot.qqrplot}} before it is returned, depending on whether the
#' package \code{ggplot2} is loaded.
#' 
#' Worm plots (de-trended Q-Q plots) draw deviations of quantile residuals (by
#' default: transformed to standard normal scale) and theoretical quantiles from
#' the same distribution against the same theoretical quantiles. For computation,
#' \code{\link{wormplot}} leverages the function \code{\link{qresiduals}}
#' employing the \code{\link{procast}}.
#' 
#' Additional options are offered for models with discrete responses where
#' randomization of quantiles is needed.
#'
#' In addition to the \code{plot} and \code{\link[ggplot2]{autoplot}} method for
#' wormplot objects, it is also possible to combine two (or more) worm plots by
#' \code{c}/\code{rbind}, which creates a set of worm plots that can then be
#' plotted in one go. 
#' 
#' @aliases wormplot wormplot.default
#' @param object an object from which probability integral transforms can be
#' extracted using the generic function \code{\link{procast}}.
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param plot Should the \code{plot} or \code{autoplot} method be called to
#' draw the computed Q-Q plot? Either set \code{plot} expicitly to \code{"base"} vs.
#' \code{"ggplot2"} to choose the type of plot, or for a logical \code{plot} argument
#' it's chosen conditional if the package \code{ggplot2} is loaded.
#' @param class Should the invisible return value be either a \code{data.frame}
#' or a \code{tibble}. Either set \code{class} expicitly to \code{"data.frame"} vs.
#' \code{"tibble"}, or for \code{NULL} it's chosen automatically conditional if the package
#' \code{tibble} is loaded.
#' @param detrend logical. Should the qqrplot be detrended, i.e, plotted as a `wormplot()`?
#' @param scale On which scale should the quantile residuals be shown: on the probability scale 
#' (\code{"uniform"}) or on the normal scale (\code{"normal"}).
#' @param nsim,delta arguments passed to \code{qresiduals}.
#' @param confint logical or character string describing the type for plotting `c("polygon", "line")`.
#' If not set to `FALSE`, the pointwise confidence interval of the (randomized)
#' quantile residuals are visualized.
#' @param simint logical. In case of discrete distributions, should the simulation
#' (confidence) interval due to the randomization be visualized?
#' @param simint_level numeric. The confidence level required for calculating the simulation
#' (confidence) interval due to the randomization.
#' @param simint_nrep numeric. The repetition number of simulated quantiles for calculating the
#' simulation (confidence) interval due to the randomization.
#' @param single_graph logical. Should all computed extended reliability
#' diagrams be plotted in a single graph?
#' @param xlab,ylab,main,\dots graphical parameters passed to
#' \code{\link{plot.qqrplot}} or \code{\link{autoplot.qqrplot}}.
#' @return An object of class \code{"qqrplot"} inheriting from
#' \code{"data.frame"} or \code{"tibble"} conditional on the argument \code{class}
#' with the following variables: \item{x}{theoretical quantiles,}
#' \item{y}{deviations between theoretical and empirical quantiles.} In case of
#' randomized residuals, \code{nsim} different \code{x} and \code{y} values, and
#' lower and upper confidence interval bounds (\code{x_rg_lwr}, \code{y_rg_lwr},
#' \code{x_rg_upr}, \code{y_rg_upr}) can optionally be returned.  Additionally,
#' \code{xlab}, \code{ylab}, \code{main}, and \code{simint_level}, as well as the
#' the (\code{scale}) and wether a \code{detrended} Q-Q residuals plot
#' was computed are stored as attributes.
#' @seealso \code{\link{qqrplot}}, \code{\link{plot.qqrplot}}, \code{\link{qqrplot}}, 
#' \code{\link{qresiduals}}, \code{\link[stats]{qqnorm}}
#' @references van Buuren S and Fredriks M (2001). \dQuote{Worm plot: simple diagnostic
#' device for modelling growth reference curves}. \emph{Statistics in
#' Medicine}, \bold{20}, 1259--1277. \doi{10.1002/sim.746}
#' @keywords hplot
#' @examples
#' 
#' ## speed and stopping distances of cars
#' m1_lm <- lm(dist ~ speed, data = cars)
#' 
#' ## compute and plot wormplot
#' wormplot(m1_lm)
#' 
#' #-------------------------------------------------------------------------------
#' ## determinants for male satellites to nesting horseshoe crabs
#' data("CrabSatellites", package = "countreg")
#' 
#' ## linear poisson model
#' m1_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#' m2_pois <- glm(satellites ~ color, data = CrabSatellites, family = poisson)
#' 
#' ## compute and plot wormplot as base graphic
#' w1 <- wormplot(m1_pois, plot = FALSE)
#' w2 <- wormplot(m2_pois, plot = FALSE)
#' 
#' ## plot combined wormplot as "ggplot2" graphic
#' ggplot2::autoplot(c(w1, w2), single_graph = TRUE, col = c(1, 2), fill = c(1, 2))
#' 
#' @export
wormplot <- function(object, ...) {
  UseMethod("wormplot")
}


#' @rdname wormplot
#' @method wormplot default
#' @export
wormplot.default <- function(object,
                             newdata = NULL,
                             plot = TRUE,
                             class = NULL,
                             detrend = TRUE,
                             scale = c("normal", "uniform"),
                             nsim = 1L,
                             delta = NULL,
                             confint = TRUE,
                             simint = TRUE,
                             simint_level = 0.95,
                             simint_nrep = 250,
                             single_graph = FALSE,
                             xlab = "Theoretical quantiles",
                             ylab = "Deviation",
                             main = NULL,
                             ...) {

  scale <- match.arg(scale)
  if (is.null(main)) {
      # Checking for class name required when called via do.call(wormplot, ...)
      main <- if (inherits(substitute(object), "name")) deparse(substitute(object)) else "object"
  }

  ## TODO: (ML) Like this, only with arg detrend or with call/eval (latter difficulties with tinytest)
  qqrplot(object,
          newdata = newdata,
          plot = plot,
          class = class,
          detrend = detrend,
          scale = scale,
          nsim = nsim,
          delta = delta,
          confint = confint,
          simint = simint,
          simint_level = simint_level,
          simint_nrep = simint_nrep,
          single_graph = single_graph,
          xlab = xlab,
          ylab = ylab,
          main = main,
          ... = ...) 
  
}
