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
#' @param confint logical or character string describing the style for plotting `c("polygon", "line")`.
#' If not set to `FALSE`, the pointwise confidence interval of the (randomized)
#' quantile residuals are visualized.
#' @param range logical or quantile specification. In case of discrete distributions, should the range 
#' (confidence interval) of values due to the randomization be visualized? If \code{TRUE}, 
#' then \code{range = c(0.01, 0.99)} is used.
#' @param range_level numeric. The confidence level required for calculating the range of values 
#' due to the randomization. 
#' @param range_nsim numeric. The number of simulated quantiles for calculating the range of values
#' due to the randomization. 
#' @param range_seed numeric. The seed to be set for calculating the range of values
#' due to the randomization. 
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
#' \code{xlab}, \code{ylab}, \code{main}, and \code{range_level}, as well as the
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
#' @export
qqrplot <- function(object, ...) {
  UseMethod("qqrplot")
}


#' @rdname qqrplot
#' @method qqrplot default
#' @export
qqrplot.default <- function(object,
                            newdata = NULL,
                            plot = TRUE,
                            class = NULL,
                            detrend = FALSE,
                            trafo = qnorm,
                            nsim = 1L,
                            delta = NULL,
                            confint = TRUE,
                            range = TRUE,
                            range_level = 0.95,
                            range_nsim = 250,
                            range_seed = 1,
                            single_graph = FALSE,
                            xlab = "Theoretical quantiles",
                            ylab = if (!detrend) "Quantile residuals" else "Deviation",
                            main = NULL,
                            ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * `object`, `newdata`, `delta w/i `qresiduals()`
  ## * `range` w/i `polygon()`
  ## * `delta` w/i `qresiduals()`
  ## * `...` in `plot()` and `autoplot()`
  stopifnot(is.logical(detrend))
  stopifnot(is.null(trafo) | is.function(trafo))
  stopifnot(is.numeric(nsim), length(nsim) == 1)
  stopifnot(
    is.numeric(range_level),
    length(range_level) == 1,
    range_level >= 0,
    range_level <= 1
  )
  stopifnot(is.numeric(range_nsim), length(range_nsim) == 1)
  stopifnot(is.numeric(range_seed), length(range_seed) == 1)
  stopifnot(is.logical(single_graph))
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
  if (!identical(range, FALSE)) {
    set.seed(range_seed)
    tmp <- qresiduals(object,
      newdata = newdata, trafo = trafo, type = "random", nsim = range_nsim,
      delta = delta
    )
    range_prob <- (1 - range_level) / 2
    range_prob <- c(range_prob, 1 - range_prob)
    qres_rg_lwr <- apply(apply(tmp, 2, sort), 1, quantile, probs = range_prob[1], na.rm = TRUE)
    qres_rg_upr <- apply(apply(tmp, 2, sort), 1, quantile, probs = range_prob[2], na.rm = TRUE)
    qthe_rg_lwr <- q2q(qres_rg_lwr)
    qthe_rg_upr <- q2q(qres_rg_upr)

    ## FIXME: (ML) Improve workaround to get range only for discrete values
    if (isTRUE(all.equal(qres_rg_lwr, qres_rg_upr, tol = .Machine$double.eps^0.4))) {
      qres_rg_lwr <- NULL
      qres_rg_upr <- NULL
      qthe_rg_lwr <- NULL
      qthe_rg_upr <- NULL
      range <- FALSE
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
        x = qthe,
        y = qres
      )
    } else {
      rval <- data.frame(
        x = qthe,
        y = qres,
        y_rg_lwr = qres_rg_lwr,
        y_rg_upr = qres_rg_upr,
        x_rg_lwr = qthe_rg_lwr,
        x_rg_upr = qthe_rg_upr
      )
    }
  } else { 
    if (any(vapply(
      list(qres_rg_lwr, qres_rg_upr, qthe_rg_lwr, 1),
      FUN = is.null,
      FUN.VALUE = FALSE
    ))) {
      rval <- data.frame(
        x = qthe,
        y = qres - qthe
      )
    } else {
      rval <- data.frame(
        x = qthe,
        y = qres - qthe,
        y_rg_lwr = qres_rg_lwr - qthe_rg_lwr,
        y_rg_upr = qres_rg_upr - qthe_rg_upr,
        x_rg_lwr = qthe_rg_lwr,
        x_rg_upr = qthe_rg_lwr
      )
    }
  }

  names(rval) <- gsub("(\\.r|\\.q)", "", names(rval))

  ## attributes for graphical display
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "range_level") <- ifelse(range, range_level, NA)
  attr(rval, "trafo") <- trafo
  attr(rval, "detrend") <- detrend

  ## add class
  if (class == "data.frame") {
    class(rval) <- c("qqrplot", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("qqrplot", class(rval))
  }

  ## plot by default
  if (plot == "ggplot2") {
    try(print(ggplot2::autoplot(rval, confint = confint, range = range, ...)))
  } else if (plot == "base") {
    try(plot(rval, confint = confint, range = range, ...))
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
  xlab <- unlist(lapply(rval, function(r) attr(r, "xlab")))
  ylab <- unlist(lapply(rval, function(r) attr(r, "ylab")))
  nam <- names(rval)
  main <- if (is.null(nam)) {
    as.vector(sapply(rval, function(r) attr(r, "main")))
  } else {
    make.unique(rep.int(nam, sapply(n, length)))
  }

  ## parameters
  detrend <- unlist(lapply(rval, function(r) attr(r, "detrend")))
  range_level <- unlist(lapply(rval, function(r) attr(r, "range_level")))
  trafo <- unlist(lapply(rval, function(r) attr(r, "trafo")))
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

  if (length(detrend) > 1) {
    if(!all(sapply(2:length(detrend), function(i) identical(detrend[[i-1]], detrend[[i]])))) {
      stop("objects with different `detrend`s are on different scales and hence must not be combined")
    } else {
    detrend <- detrend[[1]]
    }
  }

  # -------------------------------------------------------------------
  # RETURN DATA
  # -------------------------------------------------------------------
  ## combine and return (fill up missing variables with NAs)
  all_names <- unique(unlist(lapply(rval, names)))
  if (any(grepl("x_1", all_names)) & any(grepl("^x$", all_names))) {
    for (i in 1:length(rval)) {
      names(rval[[i]])[grepl("^x$", names(rval[[i]]))] <- "x_1"
      names(rval[[i]])[grepl("^y$", names(rval[[i]]))] <- "y_1"
    }
    all_names <- unique(unlist(lapply(rval, names)))
  }

  rval <- do.call(
    "rbind.data.frame",
    c(lapply(
      rval,
      function(x) data.frame(c(x, sapply(setdiff(all_names, names(x)), function(y) NA)))
    ),
    make.row.names = FALSE
    )
  )

  rval$group <- if (length(n) < 2L) NULL else rep.int(seq_along(n), n)
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "detrend") <- detrend
  attr(rval, "range_level") <- range_level
  attr(rval, "trafo") <- trafo

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
#' @param confint logical or character string describing the style for plotting `c("polygon", "line")`.
#' If not set to `FALSE`, the pointwise confidence interval of the (randomized)
#' quantile residuals are visualized.
#' @param range logical or quantile specification. Should the range of
#' quantiles of the randomized quantile residuals be visualized? If
#' \code{TRUE}, then \code{range = c(0.01, 0.99)} is used.
#' @param xlab,ylab,main,\dots graphical plotting parameters passed to
#' \code{\link[graphics]{plot}} or \code{\link[graphics]{points}},
#' respectively.
#' @param ref,xlim,ylim,col,fill,alpha_min,pch,axes,box additional graphical
#' parameters for base plots, whereby \code{x} is a object of class \code{qqrplot}.
#' @param alpha,colour,shape,size,stroke,legend,detrend,identity,trafo,probs 
#' graphical parameters passed to \code{ggplot2} style plots, whereby
#' \code{object} is a object of class \code{qqrplot}.
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
#' qqrplot(m1_lm, plot = "base", ref = "blue", lty = 2, pch = 20)
#' 
#' ## add separate model
#' if (require("crch", quietly = TRUE)) {
#'   m1_crch <- crch(dist ~ speed | speed, data = cars)
#'   points(qqrplot(m1_crch, plot = FALSE), col = 2, lty = 2, range = 2)
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
#'   plot(c(qq2_lm, qq2_crch), col = c(1, 2), range = c(1, 2), ref = 3, single_graph = TRUE)
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
                         confint = TRUE,  # FIXME: (ML) Implement different plotting styles
                         range = TRUE,
                         ref = TRUE,
                         identity = TRUE,
                         probs = c(0.25, 0.75),
                         trafo = NULL,
                         xlim = c(NA, NA),
                         ylim = c(NA, NA),
                         xlab = NULL,
                         ylab = NULL,
                         main = NULL,
                         col = adjustcolor("black", alpha.f = 0.4),
                         fill = adjustcolor("black", alpha.f = 0.2),
                         alpha_min = 0.2,
                         pch = 19,
                         axes = TRUE,
                         box = TRUE,
                         ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * lengths of all arguments are checked by recycling
  ## * `ref` w/i `abline()`
  ## * `xlab`, `ylab`, `main` and `....` w/i `plot()`
  ## * `col`, `pch` w/i `lines()`
  ## * `range`, `fill` in `polygon()`
  ## * `alpha_min` w/i `set_minimum_transparency()`
  stopifnot(is.logical(single_graph))
  stopifnot(is.logical(identity))
  stopifnot(is.numeric(probs), length(probs) == 2)
  stopifnot(length(trafo) <= 1, is.null(trafo) || is.function(trafo))
  stopifnot(length(detrend) <= 1, is.null(detrend) || is.logical(detrend))
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))

  ## get confint
  if (isFALSE(confint)) {
    confint <- "none"
  } else if (isTRUE(confint)) {
    #confint <- "polygon"
    confint <- "line"
  }
  ## FIXME: (ML) Implemnt style polygon
  if (confint == "polygon") {
    confint <- "line"
    warning("confint style polygon not yet implemented in base plots, set to `confint = 'line'`")
  }
  confint <- match.arg(confint, c("polygon", "line", "none"))

  ## get detrend
  if (!is.null(attr(x, "detrend")) && !is.null(detrend)) {
    detrend <- attr(x, "detrend")
    warning(paste0(
      "argument `detrend` is overwritten by x's attribute `detrend = ",
      detrend,
      "`"
    ))
  } else if (!is.null(attr(x, "detrend")) && is.null(detrend)) {
    detrend <- attr(x, "detrend")
  }

  ## get trafo
  if (!is.null(attr(x, "trafo")) && !is.null(trafo)) {
    trafo <- attr(x, "trafo")
    warning("argument `trafo` is overwritten by x's attribute `trafo`")
  } else if (!is.null(attr(x, "trafo")) && is.null(trafo)) {
    trafo <- attr(x, "trafo")
  }

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  if (is.list(xlim)) xlim <- as.data.frame(do.call("rbind", xlim))
  if (is.list(ylim)) ylim <- as.data.frame(do.call("rbind", ylim))
  plot_arg <- data.frame(1:n, range, ref,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]],
    col, fill, alpha_min, pch, axes, box
  )[, -1]

  ## annotation
  if (single_graph) {
    if (is.null(xlab)) xlab <- "Theoretical quantiles"
    if (is.null(ylab)) ylab <- if (!detrend) "Quantile residuals" else "Deviation"
    if (is.null(main)) main <- if (!detrend) "Q-Q residuals plot" else "Worm plot"  # FIXME: (ML) Achim prefers other title
  } else {
    if (is.null(xlab)) xlab <- TRUE
    if (is.null(ylab)) ylab <- TRUE
    if (is.null(main)) main <- TRUE
    xlab <- rep(xlab, length.out = n)
    ylab <- rep(ylab, length.out = n)
    main <- rep(main, length.out = n)
    if (is.logical(xlab)) xlab <- ifelse(xlab, attr(x, "xlab"), "")
    if (is.logical(ylab)) ylab <- ifelse(ylab, attr(x, "ylab"), "")
    if (is.logical(main)) main <- ifelse(main, attr(x, "main"), "")
  }

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION
  # -------------------------------------------------------------------
  qqrplot_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get xlim and ylim conditional on range and on single_graph
    if (single_graph) {
      if (
        !identical(plot_arg$range[j], FALSE) &&
          any(!is.na(attr(d, "range_level")))
      ) {
        if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) {
          tmp <- range(as.matrix(x[grepl("x", names(x))]), finite = TRUE)
          plot_arg[j, c("xlim1", "xlim2")] <- tmp
        }
        if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) {
          tmp <- range(as.matrix(x[grepl("y", names(x))]), finite = TRUE)
          plot_arg[j, c("ylim1", "ylim2")] <- tmp
        }
      } else {
        if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) {
          tmp <- range(as.matrix(x[grepl("^x$|x_[0-9]", names(x))]), finite = TRUE)
          plot_arg[j, c("xlim1", "xlim2")] <- tmp
        }
        if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) {
          tmp <- range(as.matrix(x[grepl("^y$|y_[0-9]", names(x))]), finite = TRUE)
          plot_arg[j, c("ylim1", "ylim2")] <- tmp
        }
      }
    } else { 
      if (
        !identical(plot_arg$range[j], FALSE) &&
          !is.na(attr(d, "range_level")[j])
      ) {
        if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) {
          tmp <- range(as.matrix(d[grepl("x", names(d))]), finite = TRUE)
          plot_arg[j, c("xlim1", "xlim2")] <- tmp
        }
        if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) {
          tmp <- range(as.matrix(d[grepl("y", names(d))]), finite = TRUE)
          plot_arg[j, c("ylim1", "ylim2")] <- tmp
        }
      } else {
        if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) {
          tmp <- range(as.matrix(d[grepl("^x$|x_[0-9]", names(d))]), finite = TRUE)
          plot_arg[j, c("xlim1", "xlim2")] <- tmp
        }
        if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) {
          tmp <- range(as.matrix(d[grepl("^y$|y_[0-9]", names(d))]), finite = TRUE)
          plot_arg[j, c("ylim1", "ylim2")] <- tmp
        }
      }
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

    ## plot range polygon
    if (!identical(plot_arg$range[j], FALSE) && !is.na(attr(d, "range_level")[j])) {
      if (isTRUE(plot_arg$range[j])) plot_arg$range[j] <- plot_arg$fill[j]

      idx_upr <- order(d$x_rg_upr)
      idx_lwr <- order(d$x_rg_lwr)
      x_pol <- c(d$x_rg_lwr[idx_lwr], d$x_rg_upr[rev(idx_upr)])
      y_pol <- c(d$y_rg_lwr[idx_lwr], d$y_rg_upr[rev(idx_upr)])
      x_pol[!is.finite(x_pol)] <- 100 * sign(x_pol[!is.finite(x_pol)]) # TODO: (ML) needed?
      y_pol[!is.finite(y_pol)] <- 100 * sign(y_pol[!is.finite(y_pol)]) # TODO: (ML) needed?

      polygon(
        x_pol,
        y_pol,
        col = set_minimum_transparency(plot_arg$range[j], alpha_min = plot_arg$alpha_min[j]),
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
        if (!identity) {
          y_tmp <- quantile(d[grepl("^y$|y_0", names(d))], probs, names = FALSE, na.rm = TRUE)
          x_tmp <- trafo(probs)
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
      if (!identical(plot_arg$ref[j], FALSE)) {
        if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
        abline(a = intercept, b = slope, col = plot_arg$ref[j], lty = 2, lwd = 1.25)
      }

      ## plot confidence lines
      if (!identical(plot_arg$confint[j], FALSE)) {
        if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- "black"
        curve(
          fun(
            x,
            n = NROW(d),
            level = 0.95,
            which = "lower",
            slope = slope,
            intercept = intercept
          ),
          lty = 2,
          lwd = 1.25,
          col = plot_arg$ref[j],
          from = plot_arg$xlim1[j],
          to = plot_arg$xlim2[j],
          add = TRUE
        )
        curve(
          fun(
            x,
            n = NROW(d),
            level = 0.95,
            which = "upper",
            slope = slope,
            intercept = intercept
          ),
          lty = 2,
          lwd = 1.25,
          col = plot_arg$ref[j],
          from = plot_arg$xlim1[j],
          to = plot_arg$xlim2[j],
          add = TRUE
        )
      }
    }

    ## add qq plot
    for (i in 1L:ncol(d[grepl("^y$|y_[0-9]", names(d))])) {
      points.default(
        d[grepl("x", names(d))][, i],
        d[grepl("y", names(d))][, i],
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
                           range = FALSE,
                           col = adjustcolor("black", alpha.f = 0.4),
                           fill = adjustcolor("black", alpha.f = 0.2),
                           alpha_min = 0.2,
                           pch = 19,
                           ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * lengths of all arguments are checked by recycling
  ## * `col`, `pch` w/i `lines()`
  ## * `range`, `fill` in `polygon()`
  ## * `alpha_min` w/i `set_minimum_transparency()`

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(1:n, range, col, fill, alpha_min, pch)[, -1]

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR POINTS
  # -------------------------------------------------------------------
  qqrplot_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## plot range polygon
    if (!identical(plot_arg$range[j], FALSE) && !is.na(attr(d, "range_level")[j])) {
      if (isTRUE(plot_arg$range[j])) plot_arg$range[j] <- plot_arg$fill[j]

      idx_upr <- order(d$x_rg_upr)
      idx_lwr <- order(d$x_rg_lwr)
      x_pol <- c(d$x_rg_lwr[idx_lwr], d$x_rg_upr[rev(idx_upr)])
      y_pol <- c(d$y_rg_lwr[idx_lwr], d$y_rg_upr[rev(idx_upr)])
      x_pol[!is.finite(x_pol)] <- 100 * sign(x_pol[!is.finite(x_pol)]) # TODO: (ML) needed?
      y_pol[!is.finite(y_pol)] <- 100 * sign(y_pol[!is.finite(y_pol)]) # TODO: (ML) needed?

      polygon(
        x_pol,
        y_pol,
        col = set_minimum_transparency(plot_arg$range[j], alpha_min = plot_arg$alpha_min[j]),
        border = NA
      )
    }

    ## add qq plot
    for (i in 1L:ncol(d[grepl("^y$|y_[0-9]", names(d))])) {
      points(
        d[grepl("x", names(d))][, i],
        d[grepl("y", names(d))][, i],
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
                             confint = TRUE,
                             range = TRUE,
                             ref = TRUE,
                             identity = TRUE, 
                             probs = c(0.25, 0.75), 
                             trafo = NULL,
                             xlim = c(NA, NA),
                             ylim = c(NA, NA),
                             xlab = NULL,
                             ylab = NULL,
                             main = NULL,
                             alpha = NA,
                             colour = "black",
                             fill = NA, 
                             shape = 19,
                             size = 2,
                             stroke = 0.5,
                             legend = FALSE,
                             ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## get base style arguments
  add_arg <- list(...)
  if (!is.null(add_arg$pch)) shape <- add_arg$pch
  if (!is.null(add_arg$lwd)) size <- add_arg$lwd
  #if (!is.null(add_arg$lty)) linetype <- add_arg$lty # TODO: (ML) Currently not needed as no linetype is used

  ## sanity checks
  stopifnot(is.logical(single_graph))
  stopifnot(is.logical(range))
  stopifnot(is.logical(ref))
  stopifnot(is.logical(identity))
  stopifnot(is.numeric(probs), length(probs) == 2)
  stopifnot(length(trafo) <= 1, is.null(trafo) || is.function(trafo))
  stopifnot(length(detrend) <= 1, is.null(detrend) || is.logical(detrend))
  stopifnot(is.logical(legend))
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))

  ## get confint
  if (isFALSE(confint)) {
    confint <- "none"
  } else if (isTRUE(confint)) {
    confint <- "polygon"
  }
  confint <- match.arg(confint, c("polygon", "line", "none"))

  ## get detrend
  if (!is.null(attr(object, "detrend")) && !is.null(detrend)) {
    detrend <- attr(object, "detrend")
    warning(paste0(
      "argument `detrend` is overwritten by object's attribute `detrend = ",
      detrend,
      "`"
    ))
  } else if (!is.null(attr(object, "detrend")) && is.null(detrend)) {
    detrend <- attr(object, "detrend")
  }

  ## get trafo
  if (!is.null(attr(object, "trafo")) && !is.null(trafo)) {
    trafo <- attr(object, "trafo")
    warning("argument `trafo` is overwritten by object's attribute `trafo`")
  } else if (!is.null(attr(object, "trafo")) && is.null(trafo)) {
    trafo <- attr(object, "trafo")
  }

  ## convert data always to data.frame
  object <- as.data.frame(object)

  ## determine grouping
  if (is.null(object$group)) object$group <- 1L
  n <- max(object$group)

  ## get title
  if (!is.null(main)) {
    title <- main[1]
    object$title <- factor(title)
  }

  ## get annotations in the right lengths
  if (is.null(xlab)) xlab <- attr(object, "xlab")
  xlab <- paste(unique(xlab), collapse = "/")
  if (is.null(ylab)) ylab <- attr(object, "ylab")
  ylab <- paste(unique(ylab), collapse = "/")
  if (is.null(main)) main <- attr(object, "main")
  main <- make.names(rep_len(main, n), unique = TRUE)

  ## prepare grouping
  object$group <- factor(object$group, levels = 1L:n, labels = main)

  # -------------------------------------------------------------------
  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING
  # -------------------------------------------------------------------
  ## Get a long data.frame with all x and y simulations
  ## FIXME: (ML) This must be done in base and somehow nicer
  object <- tidyr::pivot_longer(object,
    cols = names(object)[grepl("^x$|x_[0-9]", names(object))],
    names_to = "x_sim", values_to = "x"
  )
  object <- tidyr::pivot_longer(object,
    cols = names(object)[grepl("^y$|y_[0-9]", names(object))],
    names_to = "y_sim", values_to = "y"
  )
  object <- object[which(gsub("x", "", object$x_sim) == gsub("y", "", object$y_sim)), ]
  object$y_sim <- NULL
  object <- as.data.frame(object)

  ## recycle arguments for plotting to match the number of groups (for `scale_<...>_manual()`)
  plot_arg <- data.frame(
    1:n,
    alpha, colour, fill, shape, size
  )[, -1]

  # -------------------------------------------------------------------
  # MAIN PLOTTING
  # -------------------------------------------------------------------
  ## actual plotting
  rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y")) 

  ## add ref
  if (ref) {
    rval <- rval +
      geom_qqr_ref(
        detrend = detrend,
        identity = identity, 
        probs = probs, 
        trafo = trafo
      )
  }

  ## add conf
  if (confint != "none") {
    rval <- rval +
      geom_qqr_confint(
        detrend = detrend,
        identity = identity, 
        probs = probs, 
        trafo = trafo,
        style = confint,
        xlim = xlim
      )
  }

  ## add range
  if (range && all(c("x_rg_lwr", "x_rg_upr", "y_rg_lwr", "y_rg_upr") %in% names(object))) {
    rval <- rval +
      geom_qqr_range(
        ggplot2::aes_string(
          x_lwr = "x_rg_lwr", 
          x_upr = "x_rg_upr", 
          y_lwr = "y_rg_lwr", 
          y_upr = "y_rg_upr",
          group = "group"
        )
      )
  }

  ## add points
  rval <- rval +
    geom_qqr_point(ggplot2::aes_string(alpha = "group", colour = "group", fill = "group", 
      shape = "group", size = "group"), stroke = stroke)
  ## FIXME: (ML) alpha is not correctly represented in the legend 
  ##  (compare: https://stackoverflow.com/q/69634268/6583972?sem=2)

  ## set the colors, shapes, etc.
  rval <- rval +
    ggplot2::scale_alpha_manual(values = plot_arg$alpha) +
    ggplot2::scale_colour_manual(values = plot_arg$colour) +
    ggplot2::scale_fill_manual(values = plot_arg$fill) + 
    ggplot2::scale_shape_manual(values = plot_arg$shape) +
    ggplot2::scale_size_manual(values = plot_arg$size) 

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
#' @param detrend Fix description.
#' @param xlim Fix description.
#' @param n Fix description.
#' @param style Fix description.
#' @examples
#' require("ggplot2")
#' ## Fit model
#' data("CrabSatellites", package = "countreg")
#' m1_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#' m2_pois <- glm(satellites ~ color, data = CrabSatellites, family = poisson)
#' 
#' ## Compute qqrplot
#' q1 <- qqrplot(m1_pois, plot = FALSE)
#' q2 <- qqrplot(m2_pois, plot = FALSE)
#' 
#' d <- c(q1, q2) 
#' 
#' ## Get label names
#' xlab <- unique(attr(d, "xlab"))
#' ylab <- unique(attr(d, "ylab"))
#' main <- attr(d, "main")
#' main <- make.names(main, unique = TRUE)
#' d$group <- factor(d$group, labels = main)
#' 
#' ## Polygon CI around identity line used as reference 
#' gg1 <- ggplot(data = d, aes(x, y, na.rm = TRUE)) + 
#'   geom_qqr_ref() + 
#'   geom_qqr_confint(fill = "red") + 
#'   geom_qqr_point() + 
#'   geom_qqr_range(
#'     aes(
#'       x_lwr = x_rg_lwr, 
#'       x_upr = x_rg_upr, 
#'       y_lwr = y_rg_lwr, 
#'       y_upr = y_rg_upr,
#'       group = group
#'     )
#'   ) + 
#'   xlab(xlab) + ylab(ylab)
#'
#' gg1
#' gg1 + facet_wrap(~group)
#' 
#' ## Polygon CI around robust reference line
#' gg2 <- ggplot(data = d, aes(x, y, na.rm = TRUE)) + 
#'   geom_qqr_ref(identity = FALSE, trafo = attr(d, "trafo")) + 
#'   geom_qqr_confint(identity = FALSE, trafo = attr(d, "trafo"), style = "line") + 
#'   geom_qqr_point() + 
#'   geom_qqr_range(
#'     aes(
#'       x_lwr = x_rg_lwr, 
#'       x_upr = x_rg_upr, 
#'       y_lwr = y_rg_lwr, 
#'       y_upr = y_rg_upr,
#'       group = group
#'     )
#'   ) + 
#'   xlab(xlab) + ylab(ylab)
#'
#' gg2
#' gg2 + facet_wrap(~group)
#' 
#' @export
geom_qqr_point <- function(mapping = NULL, data = NULL, stat = "identity",
                            position = "identity", na.rm = FALSE, 
                            show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    geom = GeomQqrPoints, mapping = mapping,  
    data = data, stat = stat, position = position, 
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' @rdname geom_qqr_point
#' @format NULL
#' @usage NULL
#' @export
GeomQqrPoints <- ggplot2::ggproto("GeomQqrPoints", ggplot2::Geom,
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


#' @rdname geom_qqr_point
#' @export
stat_qqr_range <- function(mapping = NULL, data = NULL, geom = "qqr_range",
                             position = "identity", na.rm = FALSE, 
                             show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatQqrRange, 
    data = data, 
    mapping = mapping, 
    geom = geom, 
    position = position, 
    show.legend = show.legend, 
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' @rdname geom_qqr_point
#' @format NULL
#' @usage NULL
#' @export
StatQqrRange <- ggplot2::ggproto("StatQqrRange", ggplot2::Stat,
## TODO: (ML) Alternative to use `stat = "identity"` in `geom_qqr_range()` and write `setup_data()`
##            fails as here aes `x_lwr`, ... are unknown and ignored
  compute_group = function(data, scales) {
    ## Manipulate object
    nd <- data.frame(
      x = c(data$x_lwr, rev(data$x_upr)),
      y = c(data$y_lwr, rev(data$y_upr))
    )
    nd
  },

  # Tells us what we need
  required_aes = c("x_lwr", "x_upr", "y_lwr", "y_upr")
)


#' @rdname geom_qqr_point
#' @export
geom_qqr_range <- function(mapping = NULL, data = NULL, stat = "qqr_range",
                             position = "identity", na.rm = FALSE, 
                             show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    geom = GeomQqrRange, mapping = mapping,  
    data = data, stat = stat, position = position, 
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' @rdname geom_qqr_point
#' @format NULL
#' @usage NULL
#' @export
GeomQqrRange <- ggplot2::ggproto("GeomQqrRange", ggplot2::GeomPolygon,
  default_aes = ggplot2::aes(colour = "NA", fill = "black", size = 0.5, linetype = 1,
  alpha = 0.2, subgroup = NULL)
)


#' @rdname geom_qqr_point
#' @export
stat_qqr_ref <- function(mapping = NULL, data = NULL, geom = "qqr_ref",
                         position = "identity", na.rm = FALSE, 
                         show.legend = NA, inherit.aes = TRUE, 
                         detrend = FALSE, identity = TRUE, probs = c(0.25, 0.75), trafo = qnorm, ...) {
  ggplot2::layer(
    stat = StatQqrRef, 
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


#' @rdname geom_qqr_point
#' @format NULL
#' @usage NULL
#' @export
StatQqrRef <- ggplot2::ggproto("StatQqrRef", ggplot2::Stat,

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


#' @rdname geom_qqr_point
#' @export
geom_qqr_ref <- function(mapping = NULL, data = NULL, stat = "qqr_ref",
                         position = "identity", na.rm = FALSE, 
                         show.legend = NA, inherit.aes = TRUE, detrend = FALSE, identity = TRUE,
                         probs = c(0.25, 0.75), trafo = qnorm, ...) {
  ggplot2::layer(
    geom = GeomQqrRef, 
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


#' @rdname geom_qqr_point
#' @format NULL
#' @usage NULL
#' @export
GeomQqrRef <- ggplot2::ggproto("GeomQqrRef", ggplot2::GeomAbline, 
  # FIXME: (ML) Maybe change it to a GeomPath to be plotted equivalent to `geom_qqr_confint()`
  default_aes = ggplot2::aes(colour = "black", size = 0.5, linetype = 2,
  alpha = NA)
)


#' @rdname geom_qqr_point
#' @format NULL
#' @usage NULL
#' @export
stat_qqr_confint <- function(mapping = NULL, data = NULL, geom = "qqr_confint", 
                             position = "identity", na.rm = FALSE,
                             show.legend = NA, inherit.aes = TRUE,
                             xlim = NULL, n = 101, 
                             detrend = FALSE, identity = TRUE, probs = c(0.25, 0.75), trafo = qnorm, 
                             style = c("polygon", "line"), ...) {

  style <- match.arg(style)

  ggplot2::layer(
    geom = geom, 
    stat = StatQqrConfint,
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
      identity = identity,
      probs = probs,
      trafo = trafo,
      style = style,
      ...
    )
  )
}


#' @rdname geom_qqr_point
#' @format NULL
#' @usage NULL
#' @export
StatQqrConfint <- ggplot2::ggproto("StatQqrConfint", ggplot2::Stat,

  compute_group = function(data, 
                           scales, 
                           xlim = NULL, 
                           n = 101, 
                           detrend = FALSE,
                           identity = TRUE, 
                           probs = c(0.25, 0.75), 
                           trafo = qnorm,
                           style = "polygon") {

    fun <- function(x, n, level = 0.95, which = c("lower", "upper")) {
      stopifnot(is.numeric(n), length(n) == 1)
      stopifnot(is.numeric(level), length(level) == 1, level >= 0, level <= 1)
      which <- match.arg(which)

      p <- pnorm(x)
      se <- (1 / dnorm(x)) * (sqrt(p * (1 - p) / n))
      rval <- as.numeric(trafo((1 - level) / 2) * se)

      if (which == "lower") {
        rval
      } else {
        -rval
      }
    }

    ## Copied and modified from `StatFunction$compute_group()`
    if (is.null(scales$x)) {
      range <- if(is.null(xlim)) c(0, 1) else xlim
      xseq <- seq(range[1], range[2], length.out = n)
      x_trans <- xseq
    } else {
      range <- if(is.null(xlim)) scales$x$dimension() else xlim

      ## Make sure range is not NA and add default ggplot2 expansion
      range[is.na(range)] <- scales$x$dimension()[is.na(range)]
      range <- range + c(-1, 1) * diff(range) * 0.05 
      ## FIXME: (ML) Better idea how to get the scales of the plot?
      xseq <- seq(range[1], range[2], length.out = n)

      if (scales$x$is_discrete()) {
        x_trans <- xseq
      } else {
        # For continuous scales, need to back transform from transformed range
        # to original values
        x_trans <- scales$x$trans$inverse(xseq)
      }
    }

    y_out1 <- do.call(fun, c(list(quote(x_trans)), list(n = length(data$x), level = 0.95, which = "upper")))
    if (!is.null(scales$y) && !scales$y$is_discrete()) {
      # For continuous scales, need to apply transform
      y_out1 <- scales$y$trans$transform(y_out1)
    }
    y_out2 <- do.call(fun, c(list(quote(x_trans)), list(n = length(data$x), level = 0.95, which = "lower")))
    if (!is.null(scales$y) && !scales$y$is_discrete()) {
      # For continuous scales, need to apply transform
      y_out2 <- scales$y$trans$transform(y_out2)
    }

    ## Employgin StatQqrRef Method
    intercept <- StatQqrRef$compute_group(data = data,
                                          scales = scales,
                                          detrend = detrend,
                                          identity = identity,
                                          probs = probs,
                                          trafo = trafo)$intercept
    slope <- StatQqrRef$compute_group(data = data,
                                      scales = scales,
                                      detrend = detrend,
                                      identity = identity,
                                      probs = probs,
                                      trafo = trafo)$slope

    if (style == "line") {
      ## prepare long format with group variable
      d <- as.data.frame(tidyr::pivot_longer(
        data.frame(
          x_noaes = xseq,
          y1 = (intercept + slope * xseq) + y_out1,
          y2 = (intercept + slope * xseq) + y_out2
        ),
        cols = c(y1, y2),
        names_to = "topbottom",
        values_to = "y_noaes",
        names_prefix = "y"
      ))
      rbind(subset(d, subset = topbottom == 1), c(NA, NA, NA), subset(d, subset = topbottom == 2))

    } else {
      ## prepare short format
      data.frame(
        x_noaes = c(xseq, rev(xseq)),
        y_noaes = c(
          (intercept + slope * xseq) + y_out2,
          rev((intercept + slope * xseq) + y_out1)
        )
      )
    }
     
  },

  # Tells us what we need
  required_aes = c("x", "y")
)


#' @rdname geom_qqr_point
#' @export
geom_qqr_confint <- function(mapping = NULL, data = NULL, stat = "qqr_confint",
                            position = "identity", na.rm = FALSE,
                            show.legend = NA, inherit.aes = TRUE,
                            xlim = NULL, n = 101, 
                            detrend = FALSE, identity = TRUE, probs = c(0.25, 0.75), trafo = qnorm, 
                            style = c("polygon", "line"), ...) {
  style <- match.arg(style)

  ggplot2::layer(
    geom = GeomQqrConfint,
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
      identity = identity,
      probs = probs,
      trafo = trafo,
      style = style,
      ...
    )
  )
}


#' @rdname geom_qqr_point
#' @export
GeomQqrConfint <- ggplot2::ggproto("GeomQqrConfint", ggplot2::Geom,

  required_aes = c("x_noaes", "y_noaes"),

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
    data <- my_modify_list(data, qqr_default_aesthetics(style), force = FALSE)
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
    data <- my_modify_list(data, qqr_default_aesthetics(params$style), force = FALSE)
    if (params$style == "polygon") {
      draw_key_polygon(data, params, size)
    } else {
      draw_key_path(data, params, size)
    }
  }
)


# Helper function inspired by internal from `ggplot2` defined in `performance.R`
my_modify_list <- function(old, new, force = FALSE) {

  if (force) {
    for (i in names(new)) old[[i]] <- new[[i]]
  } else {
    for (i in names(new)) old[[i]] <- if (all(is.na(old[[i]]))) new[[i]] else old[[i]]
  }

  old
}


## Helper function inspired by internal from `ggplot2` defined in `geom-sf.R`
qqr_default_aesthetics <- function(style) {
  if (style == "line") {
    my_modify_list(ggplot2::GeomPath$default_aes, list(colour = "black", size = 0.5, linetype = 2, alpha = NA), 
      force = TRUE)
  } else {
    my_modify_list(ggplot2::GeomPolygon$default_aes, list(colour = "NA", fill = "black", size = 0.5, 
      linetype = 1, alpha = 0.2, subgroup = NULL), force = TRUE)
  }
}

