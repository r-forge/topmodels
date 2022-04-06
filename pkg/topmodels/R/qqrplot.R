# -------------------------------------------------------------------
# Programming outline: (Randomized) Q-Q residuals plot
# -------------------------------------------------------------------
#
# - Observed y in-sample or out-of-sample (n x 1)
# - Predicted probabilities F_y(y - eps) and F_y(y) (n x 2)
# - Two columns can be essentially equal -> continuous
#   or different -> (partially) discrete
# - Potentially transform uniform scale to different
#   distribution (default: Gaussian, via qnorm()).
#
# - Plot ordered empirical quantile residuals against
#   theoretical quantiles (from same distribution)
# - To deal with point masses, draw either multiple random
#   draws (enable alpha blending by default) or shade quantiles

# Functions:
# - qqrplot() generic plus default method
# - Return object of class "qqrplot" that is plotted by default
# - But has plot=FALSE so that suitable methods can be added afterwards
# - At least methods: plot(), autoplot()
# -------------------------------------------------------------------


#' Q-Q Plots for Quantile Residuals
#' 
#' Visualize goodness of fit of regression models by Q-Q plots using quantile
#' residuals. If \code{plot = TRUE}, the resulting object of class
#' \code{"qqrplot"} is plotted by \code{\link{plot.qqrplot}} or
#' \code{\link{autoplot.qqrplot}} before it is returned, depending on whether the
#' package \code{ggplot2} is loaded.
#' 
#' Q-Q residuals plots draw quantile residuals (by default: transformed to standard
#' normal scale) against theoretical quantiles from the same distribution.
#' Alternatively, transformations to other distributions can also be used,
#' specifically using no transformation at all, i.e., remaining on the uniform
#' scale (via \code{trafo = NULL} or equivalently \code{qunif} or
#' \code{identity}). For computation, \code{\link{qqrplot}} leverages the function
#' \code{\link{qresiduals}} employing the \code{\link{procast}} generic.
#' 
#' Additional options are offered for models with discrete responses where
#' randomization of quantiles is needed.
#'
#' In addition to the \code{plot} and \code{\link[ggplot2]{autoplot}} method for
#' qqrplot objects, it is also possible to combine two (or more) Q-Q residuals plots by
#' \code{c}/\code{rbind}, which creates a set of Q-Q residuals plots that can then be
#' plotted in one go. 
#' 
#' @aliases qqrplot qqrplot.default c.qqrplot
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
#' @param trafo function for tranforming residuals from probability scale to a
#' different distribution scale (default: Gaussian).
#' @param nsim,delta arguments passed to \code{qresiduals}.
#' @param simint logical. In case of discrete distributions, should the simulation
#' (confidence) interval due to the randomization be visualized?
#' @param simint_level numeric. The confidence level required for calculating the simulation
#' (confidence) interval due to the randomization.
#' @param simint_nrep numeric. The repetition number of simulated quantiles for calculating the
#' simulation (confidence) interval due to the randomization.
#' @param confint logical or character string describing the style for plotting `c("polygon", "line")`.
#' If not set to `FALSE`, the pointwise confidence interval of the (randomized)
#' quantile residuals are visualized.
#' @param ref logical. Should a reference line be plotted?
#' @param xlab,ylab,main,\dots graphical parameters passed to
#' \code{\link{plot.qqrplot}} or \code{\link{autoplot.qqrplot}}.
#' @return An object of class \code{"qqrplot"} inheriting from
#' \code{"data.frame"} or \code{"tibble"} conditional on the argument \code{class}
#' with the following variables: \item{x}{theoretical quantiles,}
#' \item{y}{deviations between theoretical and empirical quantiles.} In case of
#' randomized residuals, \code{nsim} different \code{x} and \code{y} values, and
#' lower and upper confidence interval bounds (\code{simint_expected}, \code{simint_observed_lwr},
#' \code{simint_observed_upr}) can optionally be returned.  Additionally,
#' \code{xlab}, \code{ylab}, \code{main}, and \code{simint_level}, as well as the
#' trafo function (\code{trafo}) and wether a \code{detrended} Q-Q residuals plot
#' was computed are stored as attributes.
#' @seealso \code{\link{plot.qqrplot}}, \code{\link{wormplot}},
#' \code{\link{qresiduals}}, \code{\link[stats]{qqnorm}}
#' @references Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile
#' Residuals.} \emph{Journal of Computational and Graphical Statistics},
#' \bold{5}, 1--10. \doi{10.2307/1390802}
#' @keywords hplot
#' @examples
#' 
#' ## speed and stopping distances of cars
#' m1_lm <- lm(dist ~ speed, data = cars)
#' 
#' ## compute and plot qqrplot
#' qqrplot(m1_lm)
#' 
#' #-------------------------------------------------------------------------------
#' ## determinants for male satellites to nesting horseshoe crabs
#' data("CrabSatellites", package = "countreg")
#' 
#' ## linear poisson model
#' m1_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#' m2_pois <- glm(satellites ~ color, data = CrabSatellites, family = poisson)
#' 
#' ## compute and plot qqrplot as base graphic
#' q1 <- qqrplot(m1_pois, plot = FALSE)
#' q2 <- qqrplot(m2_pois, plot = FALSE)
#' 
#' ## plot combined qqrplot as "ggplot2" graphic
#' ggplot2::autoplot(c(q1, q2), single_graph = TRUE, col = c(1, 2), fill = c(1, 2))
#'
#' ## Use different `trafo`s with confidence intervals
#' qqrplot(m1_pois, trafo = qunif)
#' qqrplot(m1_pois, trafo = qnorm)
#' qqrplot(m1_pois, detrend = TRUE, trafo = qunif, confint = "line")
#' qqrplot(m1_pois, detrend = TRUE, trafo = qnorm, confint = "line")
#' 
#' @export
qqrplot <- function(object, ...) {
  UseMethod("qqrplot")
}


#' @rdname qqrplot
#' @method qqrplot default
#' @export
qqrplot.default <- function(
                            ## computation arguments
                            object,
                            newdata = NULL,
                            plot = TRUE,
                            class = NULL,
                            detrend = FALSE,
                            trafo = qnorm,
                            nsim = 1L,
                            delta = NULL,
                            simint = TRUE,
                            simint_level = 0.95,
                            simint_nrep = 250,

                            ## plotting arguments
                            confint = TRUE,
                            ref = TRUE,
                            xlab = "Theoretical quantiles",
                            ylab = if (!detrend) "Quantile residuals" else "Deviation",
                            main = NULL,
                            ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * `object`, `newdata`, `delta w/i `qresiduals()`
  ## * `simint` w/i `polygon()`
  ## * `delta` w/i `qresiduals()`
  ## * `...` in `plot()` and `autoplot()`
  stopifnot(is.logical(detrend))
  stopifnot(is.null(trafo) | is.function(trafo))
  stopifnot(is.numeric(nsim), length(nsim) == 1)
  stopifnot(
    is.numeric(simint_level),
    length(simint_level) == 1,
    simint_level >= 0,
    simint_level <= 1
  )
  stopifnot(is.numeric(simint_nrep), length(simint_nrep) == 1)
  stopifnot(length(xlab) == 1)
  stopifnot(length(ylab) == 1)
  stopifnot(length(main) == 1 || length(main) == 0)

  ## guess plotting flavor
  if (isFALSE(plot)) {
    plot <- "none"
  } else if (isTRUE(plot)) {
    plot <- if ("ggplot2" %in% .packages()) "ggplot2" else "base"
  } else if (!is.character(plot)) {
    plot <- "base"
  }
  plot <- try(match.arg(plot, c("none", "base", "ggplot2")))
  stopifnot(
    "`plot` must either be logical, or match the arguments 'none', 'base' or 'ggplot2'" =
      !inherits(plot, "try-error")
  )

  ## guess output class
  if (is.null(class)) {
    class <- if ("tibble" %in% .packages()) "tibble" else "data.frame"
  }
  class <- try(match.arg(class, c("tibble", "data.frame")))
  stopifnot(
    "`class` must either be NULL, or match the arguments 'tibble', or 'data.frame'" =
      !inherits(class, "try-error")
  )

  # -------------------------------------------------------------------
  # COMPUTATION OF QUANTILE RESIDUALS
  # -------------------------------------------------------------------
  qres <- qresiduals(object,
    newdata = newdata, trafo = trafo, type = "random", nsim = nsim, delta = delta
  )
  if (is.null(dim(qres))) qres <- matrix(qres, ncol = 1L)

  ## compute corresponding quantiles on the transformed scale (default: normal)
  if (is.null(trafo)) trafo <- identity
  q2q <- function(y) trafo(ppoints(length(y)))[order(order(y))]
  qthe <- apply(qres, 2L, q2q)

  ## compute rg interval
  ## FIXME: (ML) Implement exact method if exists (see "inst/misc/2021_04_16_errorsearch_qqrplot.Rmd")
  ## FIXME: (ML) Return all in the same order w/o x values (same for additional nsim) -> might be an error
  if (!identical(simint, FALSE)) {
    tmp <- qresiduals(object,
      newdata = newdata, trafo = trafo, type = "random", nsim = simint_nrep,
      delta = delta
    )
    simint_prob <- (1 - simint_level) / 2
    simint_prob <- c(simint_prob, 1 - simint_prob)
    qres_rg_lwr <- apply(apply(tmp, 2, sort), 1, quantile, probs = simint_prob[1], na.rm = TRUE)
    qres_rg_upr <- apply(apply(tmp, 2, sort), 1, quantile, probs = simint_prob[2], na.rm = TRUE)
    qthe_rg_lwr <- q2q(qres_rg_lwr)
    qthe_rg_upr <- q2q(qres_rg_upr)

    ## FIXME: (ML) Improve workaround to get simint only for discrete values
    if (isTRUE(all.equal(qres_rg_lwr, qres_rg_upr, tol = .Machine$double.eps^0.4))) {
      qres_rg_lwr <- NULL
      qres_rg_upr <- NULL
      qthe_rg_lwr <- NULL
      qthe_rg_upr <- NULL
      simint <- FALSE
    }
  } else {
    qres_rg_lwr <- NULL
    qres_rg_upr <- NULL
    qthe_rg_lwr <- NULL
    qthe_rg_upr <- NULL
  }

  ## labels
  if (is.null(main)) main <- deparse(substitute(object))

  # -------------------------------------------------------------------
  # OUTPUT AND OPTIONAL PLOTTING
  # -------------------------------------------------------------------
  ## collect everything as data.frame (for detrend TRUE/FALSE)
  if (!detrend) {
    if (any(vapply(
      list(qres_rg_lwr, qres_rg_upr, qthe_rg_lwr, 1),
      FUN = is.null,
      FUN.VALUE = FALSE
    ))) {
      rval <- data.frame(
        observed = qres,
        expected = qthe,
        simint_observed_lwr = NA_real_,
        simint_observed_upr = NA_real_,
        simint_expected = NA_real_
      )
    } else {
      rval <- data.frame(
        observed = qres,
        expected = qthe,
        simint_observed_lwr = qres_rg_lwr,
        simint_observed_upr = qres_rg_upr,
        simint_expected = qthe_rg_lwr
      )
    }
  } else { 
    if (any(vapply(
      list(qres_rg_lwr, qres_rg_upr, qthe_rg_lwr, 1),
      FUN = is.null,
      FUN.VALUE = FALSE
    ))) {
      rval <- data.frame(
        observed = qres - qthe,
        expected = qthe,
        simint_observed_lwr = NA_real_,
        simint_observed_upr = NA_real_,
        simint_expected = NA_real_
      )
    } else {
      rval <- data.frame(
        observed = qres - qthe,
        expected = qthe,
        simint_observed_lwr = qres_rg_lwr - qthe_rg_lwr,
        simint_observed_upr = qres_rg_upr - qthe_rg_upr,
        simint_expected = qthe_rg_lwr
      )
    }
  }
 
  names(rval) <- gsub("(\\.r|\\.q)", "", names(rval))
 
  ## attributes for graphical display
  attr(rval, "detrend") <- detrend
  attr(rval, "trafo") <- trafo

  attr(rval, "simint") <- simint
  attr(rval, "confint") <- confint
  attr(rval, "ref") <- ref
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main

  ## add class
  if (class == "data.frame") {
    class(rval) <- c("qqrplot", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("qqrplot", class(rval))
  }

  ## plot by default
  if (plot == "ggplot2") {
    try(print(ggplot2::autoplot(rval, confint = confint, simint = simint, ...)))
  } else if (plot == "base") {
    try(plot(rval, ...))
  }

  ## return invisibly
  invisible(rval)
}


#' @export
c.qqrplot <- function(...) {
  # -------------------------------------------------------------------
  # GET DATA
  # -------------------------------------------------------------------
  ## list of qqrplots
  rval <- list(...)

  ## set class to tibble if any rval is a tibble
  if (any(do.call("c", lapply(rval, class)) %in% "tbl")) {
    class <- "tibble"
  } else {
    class <- "data.frame"
  }

  ## remove temporary the class (needed below for `c()`)
  ## FIXME: (ML) Rewrite by, e.g., employing `lapply()`
  for (i in 1:length(rval)) class(rval[[i]]) <- class(rval[[i]])[!class(rval[[i]]) %in% "qqrplot"]

  ## convert always to data.frame
  if (class == "tibble") {
    rval <- lapply(rval, as.data.frame)
  }

  ## group sizes
  for (i in seq_along(rval)) {
    if (is.null(rval[[i]]$group)) rval[[i]]$group <- 1L
  }
  n <- lapply(rval, function(r) table(r$group))

  # -------------------------------------------------------------------
  # PREPARE DATA
  # -------------------------------------------------------------------
  ## labels
  xlab <- prepare_arg_for_attributes(rval, "xlab")
  ylab <- prepare_arg_for_attributes(rval, "ylab")
  nam <- names(rval)
  main <- if (is.null(nam)) {
    prepare_arg_for_attributes(rval, "main")
  } else {
    make.unique(rep.int(nam, sapply(n, length)))
  }

  ## parameters
  detrend <- prepare_arg_for_attributes(rval, "detrend", force_single = TRUE)
  trafo <- prepare_arg_for_attributes(rval, "trafo", force_single = FALSE) # check/force below
  simint <- prepare_arg_for_attributes(rval, "simint")
  confint <- prepare_arg_for_attributes(rval, "confint")
  ref <- prepare_arg_for_attributes(rval, "ref")
  n <- unlist(n)

  # -------------------------------------------------------------------
  # CHECK FOR COMPATIBILITY
  # -------------------------------------------------------------------
  if (length(trafo) > 1) {
    if(!all(sapply(2:length(trafo), function(i) identical(trafo[[i-1]], trafo[[i]])))) {
      stop("objects with different `trafo`s are on different scales and hence must not be combined")
    } else {
    trafo <- trafo[[1]]
    }
  }

  # -------------------------------------------------------------------
  # RETURN DATA
  # -------------------------------------------------------------------
  ## combine and return (fill up missing variables with NAs)
  all_names <- unique(unlist(lapply(rval, names)))

  ## get both objects on the same scale
  rval <- lapply(rval, summary.qqrplot, detrend = detrend)
  rval <- lapply(rval, as.data.frame) # remove inner class

  ## combine and return (fill up missing variables with NAs)
  rval <- do.call(
    "rbind.data.frame",
    c(lapply(
      rval,
      function(x) {
        data.frame(c(x, sapply(setdiff(all_names, names(x)), function(y) NA)))
      }
    ),
    make.row.names = FALSE
    )
  )

  rval$group <- if (length(n) < 2L) NULL else rep.int(seq_along(n), n)

  ## add attributes
  attr(rval, "detrend") <- detrend
  attr(rval, "trafo") <- trafo

  attr(rval, "simint") <- simint
  attr(rval, "confint") <- confint
  attr(rval, "ref") <- ref
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main

  ## set class to data.frame or tibble
  if (class == "data.frame") {
    class(rval) <- c("qqrplot", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("qqrplot", class(rval))
  }

  ## return
  return(rval)
}


#' @export
rbind.qqrplot <- c.qqrplot


#' S3 Methods for Plotting Q-Q Residuals Plots
#' 
#' Generic plotting functions for Q-Q residuals plots of the class \code{"qqrplot"}
#' computed by \code{link{qqrplot}}. 
#' 
#' Q-Q residuals plot draw quantile residuals (by default: transformed to standard
#' normal scale) against theoretical quantiles from the same distribution.
#' Alternatively, transformations to other distributions can also be used,
#' specifically using no transformation at all, i.e., remaining on the uniform
#' scale (via \code{trafo = NULL} or equivalently \code{qunif} or
#' \code{identity}).
#'
#' Q-Q residuals plots can be rendered as \code{ggplot2} or base R graphics by using
#' the generics \code{\link[ggplot2]{autoplot}} or \code{\link[graphics]{plot}}. 
#' For a single base R graphically panel, \code{\link{points}} adds an additional Q-Q
#' residuals plot.
#' 
#' @aliases plot.qqrplot points.qqrplot autoplot.qqrplot
#' @param x,object an object of class \code{qqrplot}.
#' @param single_graph logical. Should all computed extended reliability
#' diagrams be plotted in a single graph?
#' @param detrend logical. Should the qqrplot be detrended, i.e, plotted as a `wormplot()`?
#' @param simint logical or quantile specification. Should the simint of
#' quantiles of the randomized quantile residuals be visualized? 
#' @param confint logical or character string describing the style for plotting `c("polygon", "line")`.
#' If not set to `FALSE`, the pointwise confidence interval of the (randomized)
#' quantile residuals are visualized.
#' @param confint_level numeric. The confidence level required.
#' @param ref logical. Should a reference line be plotted?
#' @param ref_identity,ref_probs Should the idenity line be plotted as reference 
#' and otherwise which probabilities should be used for defining the reference line?
#' @param xlab,ylab,main,\dots graphical plotting parameters passed to
#' \code{\link[graphics]{plot}} or \code{\link[graphics]{points}},
#' respectively.
#' @param xlim,ylim,axes,box additional graphical
#' parameters for base plots, whereby \code{x} is a object of class \code{qqrplot}.
#' @param col,pch, graphical parameters for the main part of the base plot.
#' @param colour,fill,alpha,shape,size,stroke graphical parameters passed to \code{ggplot2} 
#' style plots, whereby \code{object} is a object of class \code{qqrplot}.
#' @param legend logical. Should a legend be added in the \code{ggplot2} style
#' graphic?
#' @param theme Which `ggplot2` theme should be used. If not set, \code{\link[ggplot2]{theme_bw}} is employed.
#' @param simint_col,simint_alpha,confint_col,confint_lty,confint_lwd,ref_col,ref_lty,ref_lwd Further graphical parameters for the `confint` and `simint` line/polygon in the base plot.
#' @param simint_fill,confint_colour,confint_fill,confint_size,confint_linetype,confint_alpha,ref_colour,ref_size,ref_linetype Further graphical parameters for the `confint` and `simint` line/polygon using \code{\link[ggplot2]{autoplot}}.
#' @seealso \code{\link{qqrplot}}, \code{\link{wormplot}},
#' \code{\link{qresiduals}}, \code{\link[stats]{qqnorm}}
#' @references Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile
#' Residuals.} \emph{Journal of Computational and Graphical Statistics},
#' \bold{5}, 1--10. \doi{10.2307/1390802}
#' @keywords hplot
#' @examples
#' 
#' ## speed and stopping distances of cars
#' m1_lm <- lm(dist ~ speed, data = cars)
#' 
#' ## compute and plot qqrplot
#' qqrplot(m1_lm)
#' 
#' ## customize colors
#' qqrplot(m1_lm, plot = "base", ref_col = "blue", lty = 2, pch = 20)
#' 
#' ## add separate model
#' if (require("crch", quietly = TRUE)) {
#'   m1_crch <- crch(dist ~ speed | speed, data = cars)
#'   points(qqrplot(m1_crch, plot = FALSE), col = 2, lty = 2, simint = 2)
#' }
#' 
#' #-------------------------------------------------------------------------------
#' if (require("crch")) {
#' 
#'   ## precipitation observations and forecasts for Innsbruck
#'   data("RainIbk", package = "crch")
#'   RainIbk <- sqrt(RainIbk)
#'   RainIbk$ensmean <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, mean)
#'   RainIbk$enssd <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, sd)
#'   RainIbk <- subset(RainIbk, enssd > 0)
#' 
#'   ## linear model w/ constant variance estimation
#'   m2_lm <- lm(rain ~ ensmean, data = RainIbk)
#' 
#'   ## logistic censored model 
#'   m2_crch <- crch(rain ~ ensmean | log(enssd), data = RainIbk, left = 0, dist = "logistic")
#' 
#'   ## compute qqrplots
#'   qq2_lm <- qqrplot(m2_lm, plot = FALSE)
#'   qq2_crch <- qqrplot(m2_crch, plot = FALSE)
#' 
#'   ## plot in single graph
#'   plot(c(qq2_lm, qq2_crch), col = c(1, 2), simint_col = c(1, 2), single_graph = TRUE)
#' }
#' 
#' #-------------------------------------------------------------------------------
#' ## determinants for male satellites to nesting horseshoe crabs
#' data("CrabSatellites", package = "countreg")
#' 
#' ## linear poisson model
#' m3_pois  <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#' 
#' ## compute and plot qqrplot as "ggplot2" graphic
#' qqrplot(m3_pois, plot = "ggplot2")
#'
#' @export
plot.qqrplot <- function(x,
                         single_graph = FALSE,
                         detrend = NULL,
                         simint = NULL,
                         confint = NULL,  # FIXME: (ML) Implement different plotting styles
                         confint_level = 0.95,
                         ref = NULL,
                         ref_identity = TRUE,
                         ref_probs = c(0.25, 0.75),
                         xlim = c(NA, NA),
                         ylim = c(NA, NA),
                         xlab = NULL,
                         ylab = NULL,
                         main = NULL,
                         axes = TRUE,
                         box = TRUE,
                         col = "black",
                         pch = 19,
                         simint_col = "black",
                         simint_alpha = 0.2,
                         confint_col = "black",
                         confint_lty = 2,
                         confint_lwd = 1.25,
                         ref_col = "black",
                         ref_lty = 2,
                         ref_lwd = 1.25,
                         ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## get default arguments
  detrend <- use_arg_from_attributes(x, "detrend", default = FALSE, force_single = TRUE)
  trafo <- use_arg_from_attributes(x, "trafo", default = NULL, force_single = TRUE)
  simint <- use_arg_from_attributes(x, "simint", default = TRUE, force_single = FALSE)
  confint <- use_arg_from_attributes(x, "confint", default = TRUE, force_single = FALSE)
  ref <- use_arg_from_attributes(x, "ref", default = TRUE, force_single = FALSE)

  ## sanity checks
  ## * lengths of all arguments are checked by recycling
  ## * `ref` w/i `abline()`
  ## * `xlab`, `ylab`, `main` and `....` w/i `plot()`
  ## * `col`, `pch` w/i `lines()`
  ## * `simint`, `simint_col` in `polygon()`
  ## * `alpha_min` w/i `set_minimum_transparency()`
  stopifnot(is.logical(single_graph))
  stopifnot(is.logical(simint))
  stopifnot(is.logical(ref))
  stopifnot(is.logical(ref_identity))
  stopifnot(is.numeric(ref_probs), length(ref_probs) == 2)
  stopifnot(length(trafo) <= 1, is.null(trafo) || is.function(trafo))
  stopifnot(length(detrend) <= 1, is.null(detrend) || is.logical(detrend))
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))

  ## get input object on correct scale
  x <- summary(x, detrend = detrend)

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  # -------------------------------------------------------------------
  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING
  # -------------------------------------------------------------------
  ## determine which confint should be plotted
  if (is.logical(confint)) {
    confint <- ifelse(confint,
      "line",
      "none"
    )
  }
  ## FIXME: (ML) Implemnt style polygon
  if (any(confint == "polygon")) {
    confint[confint == "polygon"] <- "line"
    warning("confint style polygon not yet implemented in base plots, set to `confint = 'line'`")
  }
  stopifnot(all(confint %in% c("polygon", "line", "none")))

  ## recycle arguments for plotting to match the number of groups
  if (is.list(xlim)) xlim <- as.data.frame(do.call("rbind", xlim))
  if (is.list(ylim)) ylim <- as.data.frame(do.call("rbind", ylim))
  plot_arg <- data.frame(1:n, simint, ref, confint,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]],
    col, pch, axes, box,
    simint_col, simint_alpha, confint_col, confint_lty, confint_lwd, ref_col, ref_lty, ref_lwd
  )[, -1]

  ## annotation
  ylab_missing <- missing(ylab)
  main_missing <- missing(main)
  if (single_graph) {
    xlab <- use_arg_from_attributes(x, "xlab", default = "Theoretical quantiles", force_single = TRUE)
    ylab <- use_arg_from_attributes(x, "ylab",
      default = if (detrend) "Deviation" else "Quantile residuals",
      force_single = TRUE
    )
    if (is.null(main)) main <- if (detrend) "Worm plot" else "Q-Q residuals plot"

  } else {
    xlab <- use_arg_from_attributes(x, "xlab", default = "Theoretical quantiles", force_single = FALSE)
    ylab <- use_arg_from_attributes(x, "ylab",
      default = if (detrend) "Deviation" else "Quantile residuals",
      force_single = FALSE
    )
    main <- use_arg_from_attributes(x, "main", # FIXME: (ML) Different in pithist()
      default = if (detrend) "Wormplot" else "Q-Q residuals plot", 
      force_single = FALSE
    )
  }

  ## fix `ylabel` according to possible new `detrend`
  if (ylab_missing) {
    ylab[(!detrend & ylab == "Deviation")] <- "Quantile residuals"
    ylab[(detrend & ylab == "Quantile residuals")] <- "Deviation"
  }

  ## fix `main` according to possible new `detrend`
  if (main_missing) {
    main[(!detrend & main == "Wormplot")] <- "Q-Q residuals plot"
    main[(detrend & main == "Q-Q residuals plot")] <- "Wormplot"
  }

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION
  # -------------------------------------------------------------------
  qqrplot_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get xlim and ylim (needs data for all groups) 
    ## TODO: (ML) In case of an object q/ simint but simint = FALSE incorrect ylim
    ylim_idx <- c(is.na(plot_arg$ylim1[j]), is.na(plot_arg$ylim2[j])) 
    xlim_idx <- c(is.na(plot_arg$xlim1[j]), is.na(plot_arg$xlim2[j])) 

    if (any(xlim_idx) && !single_graph) {
      plot_arg[j, c("xlim1", "xlim2")[xlim_idx]] <- 
        range(
          as.matrix(d[grepl("expected|simint_expected", names(d))]), 
          finite = TRUE
        )[xlim_idx]
    } else if (any(xlim_idx) && single_graph) {
      plot_arg[j, c("xlim1", "xlim2")[xlim_idx]] <- 
        range(
          as.matrix(x[grepl("expected|simint_expected", names(x))]), 
          finite = TRUE
        )[xlim_idx]
    }

    if (any(ylim_idx) && !single_graph) {
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <-
        range(
          as.matrix(d[grepl("observed|simint_observed_lwr|simint_observed_upr", names(d))]), 
          finite = TRUE
        )[ylim_idx]
    } else if (any(ylim_idx) && single_graph) {
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <-
        range(
          as.matrix(x[grepl("observed|simint_observed_lwr|simint_observed_upr", names(x))]), 
          finite = TRUE
        )[ylim_idx]
    }

    ## trigger plot
    if (j == 1 || (!single_graph && j > 1)) {
      plot(0, 0,
        type = "n", xlim = c(plot_arg$xlim1[j], plot_arg$xlim2[j]),
        ylim = c(plot_arg$ylim1[j], plot_arg$ylim2[j]),
        xlab = xlab[j], ylab = ylab[j], main = main[j], axes = FALSE, ...
      )
      if (plot_arg$axes[j]) {
        axis(1)
        axis(2)
      }
      if (plot_arg$box[j]) {
        box()
      }
    }

    ## plot simint polygon
    if (isTRUE(plot_arg$simint[j])) {
      idx <- order(d$simint_expected)
      x_pol <- c(d$simint_expected[idx], d$simint_expected[rev(idx)])
      y_pol <- c(d$simint_observed_lwr[idx], d$simint_observed_upr[rev(idx)])
      x_pol[!is.finite(x_pol)] <- 100 * sign(x_pol[!is.finite(x_pol)]) # TODO: (ML) needed?
      y_pol[!is.finite(y_pol)] <- 100 * sign(y_pol[!is.finite(y_pol)]) # TODO: (ML) needed?

      polygon(
        x_pol,
        y_pol,
        col = colorspace::adjust_transparency(plot_arg$simint_col[j], alpha = plot_arg$simint_alpha[j]),
        border = NA
      )
    }

    ## helper function for plotting confint lines
    fun <- function(x, n, level = 0.95, which = c("lower", "upper"), slope = 0, intercept = 1) {
      stopifnot(is.numeric(n), length(n) == 1)
      stopifnot(is.numeric(level), length(level) == 1, level >= 0, level <= 1)
      which <- match.arg(which)

      p <- pnorm(x)
      se <- (1 / dnorm(x)) * (sqrt(p * (1 - p) / n))
      rval <- as.numeric(trafo((1 - level) / 2) * se)

      if (which == "lower") {
        (intercept + slope * x) + rval
      } else {
        (intercept + slope * x) - rval
      }
    }

    ## compute intercept and slope of reference line
    if (j == 1 || (!single_graph && j > 1)) {
      if (!detrend) {
        if (!ref_identity) {
          y_tmp <- quantile(d[grepl("^y$|y_0", names(d))], ref_probs, names = FALSE, na.rm = TRUE)
          x_tmp <- trafo(ref_probs)
          slope <- diff(y_tmp) / diff(x_tmp)
          intercept <- y_tmp[1L] - slope * x_tmp[1L]

        } else { 
          slope = 1
          intercept = 0
        }
      } else {
        slope = 0
        intercept = 0
      }
      
      ## plot reference line
      if (isTRUE(plot_arg$ref[j])) {
        abline(
          a = intercept, 
          b = slope, 
          col = plot_arg$ref_col[j], 
          lty = plot_arg$ref_lty[j], 
          lwd = plot_arg$ref_lwd[j]
        )
      }

      ## plot confidence lines
      if (plot_arg$confint[j] == "line") {
        curve(
          fun(
            x,
            n = NROW(d),
            level = confint_level,
            which = "lower",
            slope = slope,
            intercept = intercept
          ),
          col = plot_arg$confint_col[j],
          lty = plot_arg$confint_lty[j],
          lwd = plot_arg$confint_lwd[j],
          from = plot_arg$xlim1[j],
          to = plot_arg$xlim2[j],
          add = TRUE
        )
        curve(
          fun(
            x,
            n = NROW(d),
            level = confint_level,
            which = "upper",
            slope = slope,
            intercept = intercept
          ),
          col = plot_arg$confint_col[j],
          lty = plot_arg$confint_lty[j],
          lwd = plot_arg$confint_lwd[j],
          from = plot_arg$xlim1[j],
          to = plot_arg$xlim2[j],
          add = TRUE
        )
      }
    }

    ## add qq plot
    for (i in 1L:ncol(d[grepl("^observed$|observed_[0-9]", names(d))])) {
      points.default(
        d[grepl("expected", names(d))][, i],
        d[grepl("observed", names(d))][, i],
        col = plot_arg$col[j], pch = plot_arg$pch[j], ...
      )
    }
  }

  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  ## set up necessary panels
  if (!single_graph && n > 1L) {
    old_pars <- par(mfrow = n2mfrow(n))
    on.exit(par(old_pars), add = TRUE)
  }

  ## draw qqrplots
  for (i in 1L:n) qqrplot_plot(x[x$group == i, ], ...)
}


#' @rdname plot.qqrplot
#' @method points qqrplot
#' @export
points.qqrplot <- function(x,
                           detrend = NULL,
                           simint = FALSE,
                           col = "black",
                           pch = 19,
                           simint_col = "black",
                           simint_alpha = 0.2,
                           ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## get default arguments
  detrend <- use_arg_from_attributes(x, "detrend", default = FALSE, force_single = TRUE)

  ## sanity checks
  ## * lengths of all arguments are checked by recycling
  ## * `col`, `pch` w/i `lines()`
  ## * `simint`, `fill` in `polygon()`
  ## * `alpha_min` w/i `set_minimum_transparency()`

  ## get input object on correct scale
  x <- summary(x, detrend = detrend)

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(1:n, simint, col, pch, simint_col, simint_alpha)[, -1]

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR POINTS
  # -------------------------------------------------------------------
  qqrplot_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## plot simint polygon
    if (isTRUE(plot_arg$simint[j])) {
      idx <- order(d$simint_expected)
      x_pol <- c(d$simint_expected[idx], d$simint_expected[rev(idx)])
      y_pol <- c(d$simint_observed_lwr[idx], d$simint_observed_upr[rev(idx)])
      x_pol[!is.finite(x_pol)] <- 100 * sign(x_pol[!is.finite(x_pol)]) # TODO: (ML) needed?
      y_pol[!is.finite(y_pol)] <- 100 * sign(y_pol[!is.finite(y_pol)]) # TODO: (ML) needed?

      polygon(
        x_pol,
        y_pol,
        col = colorspace::adjust_transparency(plot_arg$simint_col[j], alpha = plot_arg$simint_alpha[j]),
        border = NA
      )
    }

    ## add qq plot
    for (i in 1L:ncol(d[grepl("^observed$|observed_[0-9]", names(d))])) {
      points.default(
        d[grepl("expected", names(d))][, i],
        d[grepl("observed", names(d))][, i],
        col = plot_arg$col[j], pch = plot_arg$pch[j], ...
      )
    }
  }

  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  for (i in 1L:n) {
    qqrplot_plot(x[x$group == i, ], ...)
  }
}


#' @rdname plot.qqrplot
#' @method autoplot qqrplot
#' @exportS3Method ggplot2::autoplot
autoplot.qqrplot <- function(object,
                             single_graph = FALSE,
                             detrend = NULL,
                             simint = NULL,
                             confint = NULL,
                             confint_level = 0.95,
                             ref = NULL,
                             ref_identity = TRUE, 
                             ref_probs = c(0.25, 0.75), 
                             xlim = c(NA, NA),
                             ylim = c(NA, NA),
                             xlab = NULL,
                             ylab = NULL,
                             main = NULL,
                             legend = FALSE,
                             theme = NULL,
                             alpha = NA,
                             colour = "black",
                             fill = NA, 
                             shape = 19,
                             size = 2,
                             stroke = 0.5,
                             simint_fill = "black",
                             simint_alpha = 0.2,
                             confint_colour = NULL,
                             confint_fill = NULL,
                             confint_size = NULL,
                             confint_linetype = NULL,
                             confint_alpha = NULL,
                             ref_colour = "black",
                             ref_size = 0.5,
                             ref_linetype = 2,
                             ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## check if ylab is defined
  ylab_missing <- missing(ylab)

  ## get default arguments
  detrend <- use_arg_from_attributes(object, "detrend", default = FALSE, force_single = TRUE)
  trafo <- use_arg_from_attributes(object, "trafo", default = TRUE, force_single = TRUE)
  simint <- use_arg_from_attributes(object, "simint", default = TRUE, force_single = TRUE)
  confint <- use_arg_from_attributes(object, "confint", default = TRUE, force_single = TRUE)
  ref <- use_arg_from_attributes(object, "ref", default = TRUE, force_single = TRUE)
  xlab <- use_arg_from_attributes(object, "xlab", default = "Theoretical quantiles", force_single = TRUE)
  ylab <- use_arg_from_attributes(object, "ylab",
    default = if (detrend) "Deviation" else "Quantile residuals",
    force_single = TRUE
  )

  ## get base style arguments
  add_arg <- list(...)
  if (!is.null(add_arg$pch)) shape <- add_arg$pch
  if (!is.null(add_arg$lwd)) size <- add_arg$lwd

  ## fix `ylab` according to possible new `detrend`
  if (ylab_missing) {
    if (detrend && ylab == "Quantile residuals") ylab <- "Deviation"
    if (!detrend && ylab == "Deviation") ylab <- "Quantile residuals"
  }

  ## sanity checks
  stopifnot(is.logical(single_graph))
  stopifnot(is.logical(simint))
  stopifnot(is.logical(ref))
  stopifnot(is.logical(ref_identity))
  stopifnot(is.numeric(ref_probs), length(ref_probs) == 2)
  stopifnot(length(trafo) <= 1, is.null(trafo) || is.function(trafo))
  stopifnot(length(detrend) <= 1, is.null(detrend) || is.logical(detrend))
  stopifnot(is.logical(legend))
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))

  ## get input object on correct scale
  object <- summary(object, detrend = detrend)

  ## convert data always to data.frame
  object <- as.data.frame(object)

  ## determine grouping
  if (is.null(object$group)) object$group <- 1L
  n <- max(object$group)

  ## get title (must be done before handling of `main`)
  if (!is.null(main)) {
    title <- main[1]
    object$title <- factor(title)
  }

  ## get main and into the right length (must be done after handling of `title`)
  main <- use_arg_from_attributes(object, "main", default = "model", force_single = FALSE)
  stopifnot(is.character(main))
  main <- make.names(rep_len(main, n), unique = TRUE)

  ## prepare grouping
  object$group <- factor(object$group, levels = 1L:n, labels = main)

  # ------------------------------------------------------------------- 
  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING 
  # ------------------------------------------------------------------- 
  ## get confint
  if (isFALSE(confint)) {
    confint <- "none"
  } else if (isTRUE(confint)) {
    confint <- "line"
  }
  confint <- match.arg(confint, c("polygon", "line", "none"))

  ## Get a long data.frame with all x and y simulations
  ## FIXME: (ML) This must be done in base and somehow nicer
  object_long <- tidyr::pivot_longer(object,
    cols = names(object)[grepl("^expected$|expected_[0-9]", names(object))],
    names_to = "expected_sim", values_to = "expected"
  )
  object_long <- tidyr::pivot_longer(object_long,
    cols = names(object_long)[grepl("^observed$|observed_[0-9]", names(object_long))],
    names_to = "observed_sim", values_to = "observed"
  )
  object_long <- object_long[which(gsub("expected", "", object_long$expected_sim) == gsub("observed", "", object_long$observed_sim)), ]
  object_long$observed_sim <- NULL
  object_long <- as.data.frame(object_long)

  ## FIXME: (ML) Improve somehow (see fixme above)
  names(object)[grepl("^expected_1$|^observed_1$", names(object))] <- c("expected", "observed")

  ## set plotting aes
  aes_ref <- set_aes_helper_geoms(
    GeomQqrplotRef$default_aes,
    list(
      colour = ref_colour,
      size = ref_size,
      linetype = ref_linetype
    )
  )

  aes_confint <- set_aes_helper_geoms(
    set_default_aes_qqrplot_confint(confint),
    list(
      colour = confint_colour,
      fill = confint_fill,
      size = confint_size,
      linetype = confint_linetype,
      alpha = confint_alpha
    )
  )

  aes_simint <- set_aes_helper_geoms(
    GeomQqrplotSimint$default_aes,
    list(
      fill = simint_fill,
      alpha = simint_alpha
    )
  )


  ## recycle arguments for plotting to match the number of groups (for `scale_<...>_manual()`)
  plot_arg <- data.frame(
    1:n,
    alpha, colour, fill, shape, size
  )[, -1]

  # -------------------------------------------------------------------
  # MAIN PLOTTING
  # -------------------------------------------------------------------
  ## actual plotting
  rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "expected", y = "observed")) 

  ## add ref
  if (ref) {
    rval <- rval +
      geom_qqrplot_ref(
        detrend = detrend,
        identity = ref_identity, 
        probs = ref_probs, 
        trafo = trafo,
        colour = aes_ref$colour,
        size = aes_ref$size,
        linetype = aes_ref$linetype
      )
  }

  ## add conf
  if (confint != "none") {
    rval <- rval +
      geom_qqrplot_confint(
        detrend = detrend,
        level = confint_level,
        identity = ref_identity, 
        probs = ref_probs, 
        trafo = trafo,
        style = confint,
        xlim = xlim,
        colour = aes_confint$colour,
        fill = aes_confint$fill,
        size = aes_confint$size,
        linetype = aes_confint$linetype,
        alpha = aes_confint$alpha
      )
  }

  ## add simint
  if (simint) {
    rval <- rval +
      geom_qqrplot_simint(
        ggplot2::aes_string(
          x = "simint_expected", 
          ymin = "simint_observed_lwr", 
          ymax = "simint_observed_upr",
          group = "group"
        ),
        na.rm = TRUE,
        fill = aes_simint$fill,
        alpha = aes_simint$alpha
      )
  }

  ## add points
  rval <- rval +
    geom_qqrplot(ggplot2::aes_string(x = "expected", y = "observed", 
      alpha = "group", colour = "group", fill = "group", 
      shape = "group", size = "group"), data = object_long, stroke = stroke
    )
  ## FIXME: (ML) alpha is not correctly represented in the legend 
  ##  (compare: https://stackoverflow.com/q/69634268/6583972?sem=2)

  ## set the colors, shapes, etc.
  rval <- rval +
    ggplot2::scale_alpha_manual(values = plot_arg$alpha, na.value = NA) +
    ggplot2::scale_colour_manual(values = plot_arg$colour, na.value = NA) +
    ggplot2::scale_fill_manual(values = plot_arg$fill, na.value = NA) + 
    ggplot2::scale_shape_manual(values = plot_arg$shape, na.value = NA) +
    ggplot2::scale_size_manual(values = plot_arg$size, na.value = NA) 

  ## annotation
  rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

  ## add legend
  if (legend) {
    rval <- rval + 
      ggplot2::labs(alpha = "Model", colour = "Model", fill = "Model", shape = "Model", 
        size = "Model") +
      ggplot2::guides(alpha = "legend", colour = "legend", fill = "legend", shape = "legend", 
        size = "legend")
  } else {
    rval <- rval + 
      ggplot2::guides(alpha = "none", colour = "none", fill = "none", shape = "none", size = "none")
  }

  ## set x and y limits
  rval <- rval + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE)

  ## set ggplot2 theme
  if (is.character(theme)) {
    theme_tmp <- try(eval(parse(text = theme)), silent = TRUE)
    if (inherits(theme_tmp, "try-error") && !grepl("^ggplot2::", theme)) {
      theme_tmp <- try(eval(parse(text = paste0("ggplot2::", theme))), silent = TRUE)
    }
    theme <- theme_tmp
    if (!is.function(theme)) {
      warning("`theme` must be a ggplot2 theme, theme-generating function or valid 'character string'")
      theme <- ggplot2::theme_bw()
    }
  }

  if (is.function(theme)) {
    theme <- try(theme(), silent = TRUE)
    if (inherits(theme, "try-error") || !inherits(theme, "theme")) {
      warning("`theme` must be a ggplot2 theme, theme-generating function or valid 'character string'")
      theme <- ggplot2::theme_bw()
    }
  }

  if (inherits(theme, "theme")) {
    rval <- rval + theme
  } else if (isTRUE(all.equal(ggplot2::theme_get(), ggplot2::theme_gray()))) {
    rval <- rval + ggplot2::theme_bw()
  }

  # -------------------------------------------------------------------
  # GROUPING (IF ANY) AND RETURN PLOT
  # -------------------------------------------------------------------
  ## grouping
  if (!single_graph && n > 1L) {
    rval <- rval + ggplot2::facet_grid(group ~ .)
  } else if (!is.null(object$title)) {
    rval <- rval + ggplot2::facet_wrap(title ~ .)
  }

  ## return ggplot object
  rval
}

# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_qqrplot()`
# -------------------------------------------------------------------

#' \code{geom_*} and \code{stat_*} for Producing Quantile Residual Q-Q Plots with `ggplot2`
#' 
#' Various \code{geom_*} and \code{stat_*} used within
#' \code{\link[ggplot2]{autoplot}} for producing quantile residual Q-Q plots.
#' 
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#' @param identity logical, should the identity line be plotted or a theoretical line
#' which passes through \code{probs} quantiles computed by \code{trafo}. 
#' @param trafo function for calculating reference line through first and third
#' quartile of theoretical distribution (default: Gaussian \code{qnorm}).
#' @param probs numeric vector of length two, representing probabilities of reference
#' line used in \code{trafo}.
#' @param detrend logical. Should the qqrplot be detrended, i.e, plotted as a `wormplot()`?
#' @param level numeric. The confidence level required.
#' @param xlim The x limits for computing the confidence intervals.
#' @param n The number of points used to compute the confidence intervals, the more the smoother.
#' @param style The confidence style.
#' @examples
#' if (require("ggplot2")) {
#'   ## Fit model
#'   data("CrabSatellites", package = "countreg")
#'   m1_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#'   m2_pois <- glm(satellites ~ color, data = CrabSatellites, family = poisson)
#'   
#'   ## Compute qqrplot
#'   q1 <- qqrplot(m1_pois, plot = FALSE)
#'   q2 <- qqrplot(m2_pois, plot = FALSE)
#'   
#'   d <- c(q1, q2) 
#'   
#'   ## Get label names
#'   xlab <- unique(attr(d, "xlab"))
#'   ylab <- unique(attr(d, "ylab"))
#'   main <- attr(d, "main")
#'   main <- make.names(main, unique = TRUE)
#'   d$group <- factor(d$group, labels = main)
#'   
#'   ## Polygon CI around identity line used as reference 
#'   gg1 <- ggplot(data = d, aes(x = expected, y = observed, na.rm = TRUE)) + 
#'     geom_qqrplot_ref() + 
#'     geom_qqrplot_confint(fill = "red") + 
#'     geom_qqrplot() + 
#'     geom_qqrplot_simint(
#'       aes(
#'         x = simint_expected, 
#'         ymin = simint_observed_lwr, 
#'         ymax = simint_observed_upr,
#'         group = group
#'       )
#'     ) + 
#'     xlab(xlab) + ylab(ylab)
#'
#'   gg1
#'   gg1 + facet_wrap(~group)
#'   
#'   ## Polygon CI around robust reference line
#'   gg2 <- ggplot(data = d, aes(x = expected, y = observed, na.rm = TRUE)) + 
#'     geom_qqrplot_ref(identity = FALSE, trafo = attr(d, "trafo")) + 
#'     geom_qqrplot_confint(identity = FALSE, trafo = attr(d, "trafo"), style = "line") + 
#'     geom_qqrplot() + 
#'     geom_qqrplot_simint(
#'       aes(
#'         x = simint_expected, 
#'         ymin = simint_observed_lwr, 
#'         ymax = simint_observed_upr,
#'         group = group
#'       )
#'     ) + 
#'     xlab(xlab) + ylab(ylab)
#'
#'   gg2
#'   gg2 + facet_wrap(~group)
#' 
#'   ## Use different `trafo`s with confidence intervals
#'   q1 <- qqrplot(m1_pois, trafo = qunif, plot = FALSE)
#'   q2 <- qqrplot(m2_pois, plot = FALSE)
#'   
#'   gg3 <- ggplot(data = q1, aes(x = expected, y = observed, na.rm = TRUE)) +
#'     geom_qqrplot_ref() +
#'     geom_qqrplot_confint(fill = "red", trafo = qunif) +
#'     geom_qqrplot()
#'   gg3
#'   
#'   gg4 <- ggplot(data = q2, aes(x = expected, y = observed, na.rm = TRUE)) +
#'     geom_qqrplot_ref() +
#'     geom_qqrplot_confint(fill = "red", trafo = qnorm) +
#'     geom_qqrplot()
#'   gg4
#' } 
#' @export
geom_qqrplot <- function(mapping = NULL, data = NULL, stat = "identity",
                            position = "identity", na.rm = FALSE, 
                            show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    geom = GeomQqrplot, mapping = mapping,  
    data = data, stat = stat, position = position, 
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' @rdname geom_qqrplot
#' @format NULL
#' @usage NULL
#' @export
GeomQqrplot <- ggplot2::ggproto("GeomQqrplot", ggplot2::Geom,
  required_aes = c("x", "y"),
  non_missing_aes = c("size", "shape", "colour"), # TODO: (ML) what is that for?
  default_aes = ggplot2::aes(
    shape = 19, colour = "black", size = 2,
    fill = NA, alpha = NA, stroke = 0.5
  ),

  setup_params = function(data, params) {
    n <- nrow(data)
    if (n > 100 && n <= 200) {
      params$alpha <- 0.3
    } else if (n > 200) {
      params$alpha <- 0.15
    } else {
      params$alpha <- 1
    }
    params
  },

  draw_panel = function(data, panel_scales, coord, alpha) {
    if (is.character(data$shape)) {
      data$shape <- translate_shape_string(data$shape)
    }

    ## Transform the data first
    coords <- coord$transform(data, panel_scales)

    ## Get alpha conditional on number of data points
    n <- nrow(data)
    if (any(is.na(coords$alpha))) {
      coords$alpha <- alpha
    }

    ## Construct a grid grob
    grid::pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = coords$shape,
      gp = grid::gpar(
        col = colorspace::adjust_transparency(coords$colour, coords$alpha),
        fill = colorspace::adjust_transparency(coords$fill, coords$alpha),
        # Stroke is added around the outside of the point
        fontsize = coords$size * ggplot2::.pt + coords$stroke * ggplot2::.stroke / 2,
        lwd = coords$stroke * ggplot2::.stroke / 2
      )
    )
  },

  draw_key = function(data, params, size) {
    if (is.na(data$alpha)) { 
      data$alpha <- params$alpha
    } 
    ggplot2::draw_key_point(data, params, size)
  }
)

# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_qqrplot_simint()`
# -------------------------------------------------------------------

#' @rdname geom_qqrplot
#' @export
stat_qqrplot_simint <- function(mapping = NULL, data = NULL, geom = "qqrplot_simint",
                             position = "identity", na.rm = FALSE, 
                             show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatQqrplotSimint, 
    data = data, 
    mapping = mapping, 
    geom = geom, 
    position = position, 
    show.legend = show.legend, 
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' @rdname geom_qqrplot
#' @format NULL
#' @usage NULL
#' @export
StatQqrplotSimint <- ggplot2::ggproto("StatQqrplotSimint", ggplot2::Stat,
## TODO: (ML) Alternative to use `stat = "identity"` in `geom_qqrplot_simint()` and write `setup_data()`
##            fails as here aes `x_lwr`, ... are unknown and ignored
  compute_group = function(data, scales) {
    
    ## Manipulate object
    nd <- data.frame(
      x = c(data$x, rev(data$x)),
      y = c(data$ymin, rev(data$ymax))
    )
    nd
  },

  # Tells us what we need
  required_aes = c("x", "ymin", "ymax")
)


#' @rdname geom_qqrplot
#' @export
geom_qqrplot_simint <- function(mapping = NULL, data = NULL, stat = "qqrplot_simint",
                             position = "identity", na.rm = FALSE, 
                             show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    geom = GeomQqrplotSimint, mapping = mapping,  
    data = data, stat = stat, position = position, 
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' @rdname geom_qqrplot
#' @format NULL
#' @usage NULL
#' @export
GeomQqrplotSimint <- ggplot2::ggproto("GeomQqrplotSimint", ggplot2::GeomPolygon,
  default_aes = ggplot2::aes(colour = "NA", fill = "black", size = 0.5, linetype = 1,
  alpha = 0.2, subgroup = NULL)
)

# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_qqrplot_ref()`
# -------------------------------------------------------------------

#' @rdname geom_qqrplot
#' @export
stat_qqrplot_ref <- function(mapping = NULL, data = NULL, geom = "qqrplot_ref",
                         position = "identity", na.rm = FALSE, 
                         show.legend = NA, inherit.aes = TRUE, 
                         detrend = FALSE, identity = TRUE, probs = c(0.25, 0.75), trafo = qnorm, ...) {
  ggplot2::layer(
    stat = StatQqrplotRef, 
    data = data, 
    mapping = mapping, 
    geom = geom, 
    position = position, 
    show.legend = show.legend, 
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm, 
      detrend = detrend,
      identity = identity,
      probs = probs,
      trafo = trafo,
      ...
    )
  )
}


#' @rdname geom_qqrplot
#' @format NULL
#' @usage NULL
#' @export
StatQqrplotRef <- ggplot2::ggproto("StatQqrplotRef", ggplot2::Stat,

  compute_group = function(data, scales, detrend, identity, probs, trafo) {
    ## Manipulate object depending on arguments `detrend` and `identity`
    if (!detrend) {
      if (!identity) {
        stopifnot(is.numeric(probs), length(probs) == 2)
        stopifnot(is.function(trafo))

        y_tmp <- quantile(data$y, probs, names = FALSE, na.rm = TRUE)
        x_tmp <- trafo(probs)
        slope <- diff(y_tmp) / diff(x_tmp)
        intercept <- y_tmp[1L] - slope * x_tmp[1L]
        nd <- data.frame(
          slope = slope,
          intercept = intercept
        )

      } else { 
        nd <- data.frame(
          slope = 1,
          intercept = 0
        )
      }
      nd
    } else {
      nd <- data.frame(
        slope = 0,
        intercept = 0
      )
      nd
    }
  },

  # Tells us what we need
  required_aes = c("x", "y")
)


#' @rdname geom_qqrplot
#' @export
geom_qqrplot_ref <- function(mapping = NULL, data = NULL, stat = "qqrplot_ref",
                         position = "identity", na.rm = FALSE, 
                         show.legend = NA, inherit.aes = TRUE, detrend = FALSE, identity = TRUE,
                         probs = c(0.25, 0.75), trafo = qnorm, ...) {
  ggplot2::layer(
    geom = GeomQqrplotRef, 
    mapping = mapping, 
    data = data, 
    stat = stat, 
    position = position, 
    show.legend = show.legend, 
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm, 
      detrend = detrend,
      identity = identity,
      probs = probs,
      trafo = trafo,
      ...
    )
  )
}


#' @rdname geom_qqrplot
#' @format NULL
#' @usage NULL
#' @export
GeomQqrplotRef <- ggplot2::ggproto("GeomQqrplotRef", ggplot2::GeomAbline, 
  # FIXME: (ML) Maybe change it to a GeomPath to be plotted equivalent to `geom_qqrplot_confint()`
  default_aes = ggplot2::aes(colour = "black", size = 0.5, linetype = 2,
  alpha = NA)
)

# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_qqrplot_confint()`
# -------------------------------------------------------------------


#' @rdname geom_qqrplot
#' @format NULL
#' @usage NULL
#' @export
stat_qqrplot_confint <- function(mapping = NULL, data = NULL, geom = "qqrplot_confint", 
                             position = "identity", na.rm = FALSE,
                             show.legend = NA, inherit.aes = TRUE,
                             xlim = NULL, n = 101, 
                             detrend = FALSE, level = 0.95, 
                             identity = TRUE, probs = c(0.25, 0.75), trafo = qnorm, 
                             style = c("polygon", "line"), ...) {

  style <- match.arg(style)

  ggplot2::layer(
    geom = geom, 
    stat = StatQqrplotConfint,
    data = data,
    mapping = mapping,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      xlim = xlim,
      n = n,
      detrend = detrend,
      level = level,
      identity = identity,
      probs = probs,
      trafo = trafo,
      style = style,
      ...
    )
  )
}


#' @rdname geom_qqrplot
#' @format NULL
#' @usage NULL
#' @export
StatQqrplotConfint <- ggplot2::ggproto("StatQqrplotConfint", ggplot2::Stat,

  required_aes = c("x", "y"),

  compute_group = function(data, 
                           scales, 
                           xlim = NULL, 
                           n = 101, 
                           detrend = FALSE,
                           level = 0.95,
                           identity = TRUE, 
                           probs = c(0.25, 0.75), 
                           trafo = qnorm,
                           style = "polygon") {

    fun <- function(x, n, level = 0.95, which = c("lower", "upper")) {
      stopifnot(is.numeric(n), length(n) == 1)
      stopifnot(is.numeric(level), length(level) == 1, level >= 0, level <= 1)
      which <- match.arg(which)

      ## FIXME: (ML) Update once `distributions3` or alternative is working
      if (identical(trafo, stats::qunif) | identical(trafo, base::identity)) {
        dFun <- dunif
        pFun <- punif
      } else if (identical(trafo, stats::qnorm)) {
        dFun <- dnorm
        pFun <- pnorm
      } else {
        stop("Appropriate `dFun` and `pFun` are not yet supported.")
      }

      p <- pFun(x) # alternative: p <- ppoints(x)
      se <- (1 / dFun(x)) * (sqrt(p * (1 - p) / n))
      rval <- as.numeric(qnorm((1 - level) / 2) * se)

      if (which == "lower") {
        -rval
      } else {
        rval
      }
    }

    ## Copied and modified from `StatFunction$compute_group()`
    if (is.null(scales$x)) {
      simint <- if(is.null(xlim)) c(0, 1) else xlim
      xseq <- seq(simint[1], simint[2], length.out = n)
      x_trans <- xseq
    } else {
      simint <- if(is.null(xlim)) scales$x$dimension() else xlim

      ## Make sure simint is not NA and add default ggplot2 expansion
      simint[is.na(simint)] <- scales$x$dimension()[is.na(simint)]
      simint <- simint + c(-1, 1) * diff(simint) * 0.05 
      ## FIXME: (ML) Better idea how to get the scales of the plot?
      xseq <- seq(simint[1], simint[2], length.out = n) # alternative: xseq <- trafo(ppoints(xseq))

      if (scales$x$is_discrete()) {
        x_trans <- xseq
      } else {
        # For continuous scales, need to back transform from transformed simint
        # to original values
        x_trans <- scales$x$trans$inverse(xseq)
      }
    }

    y_out1 <- do.call(fun, c(list(quote(x_trans)), list(n = length(data$x), level = level, which = "upper")))
    if (!is.null(scales$y) && !scales$y$is_discrete()) {
      # For continuous scales, need to apply transform
      y_out1 <- scales$y$trans$transform(y_out1)
    }
    y_out2 <- do.call(fun, c(list(quote(x_trans)), list(n = length(data$x), level = level, which = "lower")))
    if (!is.null(scales$y) && !scales$y$is_discrete()) {
      # For continuous scales, need to apply transform
      y_out2 <- scales$y$trans$transform(y_out2)
    }

    # Must make sure that is not NA for specific trafo (due to extension of plot simint)
    idx_na <- is.na(y_out1) | is.na(y_out2)

    ## Employing StatQqrplotRef Method
    intercept <- StatQqrplotRef$compute_group(data = data,
                                          scales = scales,
                                          detrend = detrend,
                                          identity = identity,
                                          probs = probs,
                                          trafo = trafo)$intercept
    slope <- StatQqrplotRef$compute_group(data = data,
                                      scales = scales,
                                      detrend = detrend,
                                      identity = identity,
                                      probs = probs,
                                      trafo = trafo)$slope

    if (style == "line") {
      ## prepare long format with group variable
      d <- as.data.frame(tidyr::pivot_longer(
        data.frame(
          x_noaes = x_trans,  
          y1 = (intercept + slope * xseq) + y_out1,
          y2 = (intercept + slope * xseq) + y_out2
        )[!idx_na, ],
        cols = c(y1, y2),
        names_to = "topbottom",
        values_to = "y_noaes",
        names_prefix = "y"
      ))
      rbind(subset(d, subset = topbottom == 1), c(NA, NA, NA), subset(d, subset = topbottom == 2))

    } else {
      ## prepare short format
      data.frame(
        x_noaes = c(x_trans, rev(x_trans)),
        y_noaes = c(
          (intercept + slope * xseq) + y_out2,
          rev((intercept + slope * xseq) + y_out1)
        )
      )[!idx_na, ]
    }
     
  }
)


#' @rdname geom_qqrplot
#' @export
geom_qqrplot_confint <- function(mapping = NULL, data = NULL, stat = "qqrplot_confint",
                            position = "identity", na.rm = FALSE,
                            show.legend = NA, inherit.aes = TRUE,
                            xlim = NULL, n = 101, 
                            detrend = FALSE, level = 0.95, 
                            identity = TRUE, probs = c(0.25, 0.75), trafo = qnorm, 
                            style = c("polygon", "line"), ...) {
  style <- match.arg(style)

  ggplot2::layer(
    geom = GeomQqrplotConfint,
    mapping = mapping,
    data = data,
    stat = stat,
    position = ggplot2::PositionIdentity,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      xlim = xlim,
      n = n,
      detrend = detrend,
      level = level,
      identity = identity,
      probs = probs,
      trafo = trafo,
      style = style,
      ...
    )
  )
}


#' @rdname geom_qqrplot
#' @export
GeomQqrplotConfint <- ggplot2::ggproto("GeomQqrplotConfint", ggplot2::Geom,

  required_aes = c("x_noaes", "y_noaes"),  # TODO: (ML) For ref = identity, would theoret. work w/o `x` and `y`

  # FIXME: (ML) Does not vary for style; this is a copy of `GeomPolygon$handle_na()`
  handle_na = function(data, params) {
    data
  },

  ## Setting up all defaults needed for `GeomPolygon` and `GeomPath`
  default_aes = ggplot2::aes(
    colour = NA,
    fill = NA,
    size = 0.5,
    linetype = NA,
    alpha = NA,
    subgroup = NULL
  ),

  draw_panel = function(data, panel_params, coord,
                        rule = "evenodd", # polygon arguments
                        lineend = "butt", linejoin = "round", # line arguments
                        linemitre = 10, na.rm = FALSE, arrow = NULL, # line arguments
                        style = c("polygon", "line")) {
    style <- match.arg(style)

    ## Swap NAs in `default_aes` with own defaults 
    data <- my_modify_list(data, set_default_aes_qqrplot_confint(style), force = FALSE)
    data$x <- data$x_noaes
    data$y <- data$y_noaes

    if (style == "polygon") {
      ggplot2::GeomPolygon$draw_panel(data, panel_params, coord, rule)
    } else {
      ggplot2::GeomPath$draw_panel(data, panel_params, coord,
                          arrow, lineend, linejoin, linemitre, na.rm)
    }

  },

  draw_key = function(data, params, size) {
    ## Swap NAs in `default_aes` with own defaults 
    data <- my_modify_list(data, set_default_aes_qqrplot_confint(params$style), force = FALSE)
    if (params$style == "polygon") {
      draw_key_polygon(data, params, size)
    } else {
      draw_key_path(data, params, size)
    }
  }
)


## Helper function inspired by internal from `ggplot2` defined in `geom-sf.R`
set_default_aes_qqrplot_confint <- function(style) {
  if (style == "line") {
    my_modify_list(ggplot2::GeomPath$default_aes, list(colour = "black", fill = NA, size = 0.5, linetype = 2, alpha = NA), 
      force = TRUE)
  } else {
    my_modify_list(ggplot2::GeomPolygon$default_aes, list(colour = "NA", fill = "black", size = 0.5, 
      linetype = 1, alpha = 0.2, subgroup = NULL), force = TRUE)
  }
}


# -------------------------------------------------------------------
# HELPER FUNCTIONS FOR `qqrplot`
# -------------------------------------------------------------------

#' @export
summary.qqrplot <- function(object,
                            detrend = NULL,
                            ...) {

  ## get arg `freq`
  detrend_object <- attr(object, "detrend")
  detrend <- use_arg_from_attributes(object, "detrend", default = FALSE, force_single = TRUE)

  stopifnot(is.logical(detrend), is.logical(detrend_object))

  if (detrend != detrend_object && detrend_object) {
    rval <- transform(object,
      observed = object$observed + object$expected,
      simint_observed_lwr = object$simint_observed_lwr + object$simint_expected,
      simint_observed_upr = object$simint_observed_upr + object$simint_expected
    )
  } else if (detrend != detrend_object && !detrend_object) {
    rval <- transform(object,
      observed = object$observed - object$expected,
      simint_observed_lwr = object$simint_observed_lwr - object$simint_expected,
      simint_observed_upr = object$simint_observed_upr - object$simint_expected
    )
  } else {
    rval <- object
  }

  ## set attributes
  attr(rval, "simint") <- attr(object, "simint")
  attr(rval, "confint") <- attr(object, "confint")
  attr(rval, "ref") <- attr(object, "ref")
  attr(rval, "xlab") <- attr(object, "xlab")
  attr(rval, "ylab") <- attr(object, "ylab")
  attr(rval, "main") <- attr(object, "main")

  attr(rval, "trafo") <- attr(object, "trafo")
  attr(rval, "detrend") <- detrend


  ## return as `data.frame` or `tibble`
  if ("data.frame" %in% class(object)) {
    class(rval) <- c("qqrplot", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("qqrplot", class(rval))
  }

  rval
}

#' @export
print.qqrplot <- function(x, ...) {

  ## get arg `type`, `style` and `freq`
  detrend <- use_arg_from_attributes(x, "detrend", default = NULL, force_single = TRUE)

  ## return custom print statement
  if (is.null(detrend)) {
    cat("A `qqrplot` object without mandatory attribute `detrend`\n\n")
  } else {
    cat(
      sprintf(
        "A `qqrplot` object with `detrend = \"%s\"`\n\n",
        detrend
      )
    )
  }

  ## call next print method
  NextMethod()
}
