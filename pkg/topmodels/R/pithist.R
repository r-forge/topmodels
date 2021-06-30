# -------------------------------------------------------------------
# Programming outline: PIT histogram
# -------------------------------------------------------------------
# - Observed y in-sample or out-of-sample (n x 1)
# - Predicted probabilities F_y(y - eps) and F_y(y) (n x 2)
# - Two columns can be essentially equal -> continuous
#   or different -> (partially) discrete
# - Breaks for predicted probabilities in [0, 1] (m x 1)
#
# - Cut probabilities at breaks -> (m-1) groups and draw histogram
# - In case of point masses either use a random draw
#   or distribute evenly across relevant intervals
# - Random draws could be drawn by hist() (current solution)
#   but proportional distribution requires drawing rectangles by hand
# - Add confidence interval as well.
# - Instead of shaded rectangles plus reference line and CI lines
#   support shaded CI plus step lines

# Functions:
# - pithist() generic plus default method
# - Return object of class "pithist" that is plotted by default
# - But has plot=FALSE so that suitable methods can be added afterwards
# - At least methods: plot(), autoplot(), lines(), possibly c(), +
# -------------------------------------------------------------------


#' PIT Histograms for Assessing Goodness of Fit of Probability Models
#' 
#' PIT histograms graphically compare empirical probabilities from fitted
#' models with a uniform distribution.
#' 
#' PIT histograms graphically the probability integral transform (PIT), i.e.,
#' observed probabilities from fitted probability models, with a uniform
#' distribution. It leverages the \code{\link{procast}} generic and then
#' essentially draws a \code{\link[graphics]{hist}}.
#' 
#' In case of discrete distributions the PIT is either drawn randomly from the
#' corresponding interval or distributed proportionally in the histogram
#' (FIXME: not yet implemented).
#' 
#' @aliases pithist pithist.default c.pithist
#' @param object an object from which probability integral transforms can be
#' extracted with \code{\link{procast}}.
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param plot Should the \code{plot} or \code{autoplot} method be called to
#' draw the computed extended reliability diagram? Either set \code{plot}
#' expicitly to "base" vs. "ggplot2" to choose the type of plot, or for a
#' logical \code{plot} argument it's chosen conditional if the package
#' \code{ggplot2} is loaded.
#' @param class Should the invisible return value be either a \code{data.frame}
#' or a \code{tibble}. Either set \code{class} expicitly to "data.frame" vs.
#' "tibble", or for NULL it's chosen automatically conditional if the package
#' \code{tibble} is loaded.
#' @param style character specifying the syle of rootogram (see below). FIXME:
#' Description
#' @param type character. In case of discrete distributions should the PITs be
#' drawn randomly from the corresponding interval or distributed
#' proportionally?
#' @param nsim integer. If \code{type} is \code{"random"} how many simulated
#' PITs should be drawn?
#' @param delta numeric. The minimal difference to compute the range of
#' proabilities corresponding to each observation according to get (randomized)
#' quantile residuals.  For \code{NULL}, the minimal observed difference in the
#' resonse divided by \code{5e-6} is used.
#' @param freq logical. If \code{TRUE}, the PIT histogram is represented by
#' frequencies, the \code{counts} component of the result; if \code{FALSE},
#' probability densities, component \code{density}, are plotted (so that the
#' histogram has a total area of one).
#' @param breaks numeric. Breaks for the histogram intervals.
#' @param confint logical. Should confident intervals be drawn?
#' @param confint_level numeric. The confidence level required.
#' @param confint_type character. Which type of confidence interval.  According
#' to Agresti and Coull (1998) for interval estimation of binomial proportions
#' an approximation can be better than exact.
#' @param single_graph logical. Should all computed extended reliability
#' diagrams be plotted in a single graph?
#' @param xlim,ylim graphical parameters. These may pertain either to the whole
#' plot or just the histogram or just the fitted line.
#' @param xlab,ylab,main graphical parameters.
#' @param \dots further graphical parameters.
#' @seealso \code{\link{procast}}, \code{\link[graphics]{hist}}
#' @references Czado C, Gneiting T, Held L (2009). \dQuote{Predictive Model
#' Assessment for Count Data.} \emph{Biometrics}, \bold{65}(4), 1254--1261.
#' 
#' Agresti A, Coull A B (1998). \dQuote{Approximate is Better than ``Exact''
#' for Interval Estimation of Binomial Proportions.} \emph{The American
#' Statistician}, \bold{52}(2), 119--126.
#' @keywords hplot
#' @examples
#' 
#' require("crch")
#' m1 <- lm(dist ~ speed, data = cars)
#' m2 <- crch(dist ~ speed | speed, data = cars)
#' m3 <- crch(dist ~ speed | speed, left = 30, data = cars)
#' 
#' pit1 <- pithist(m1)
#' pit2 <- pithist(m2, plot = FALSE)
#' pit3 <- pithist(m3, plot = FALSE)
#' 
#' plot(pit1, confint = "red", ref = "blue", fill = "lightblue")
#' 
#' plot(c(pit1, pit2), col = c(1, 2), single_graph = TRUE, style = "lines")
#' lines(pit3, col = 3)
#' 
#' @export
pithist <- function(object, ...) {
  UseMethod("pithist")
}


#' @rdname pithist
#' @method pithist default
#' @export
pithist.default <- function(object,
                            newdata = NULL,
                            plot = TRUE,
                            class = NULL,
                            style = c("histogram", "lines"),
                            type = c("random", "proportional"), # FIXME: (ML) not yet implemented
                            nsim = 1L,
                            delta = NULL,
                            freq = FALSE,
                            breaks = NULL,
                            confint = TRUE,
                            confint_level = 0.95,
                            confint_type = c("exact", "approximation"),
                            single_graph = FALSE,
                            xlim = c(0, 1),
                            ylim = c(0, NA),
                            xlab = "PIT",
                            ylab = if (freq) "Frequency" else "Density",
                            main = NULL,
                            ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * `object` and `newdata` w/i `newrepsone()`
  ## * `delta w/i `qresiduals()`
  ## * `confint` w/i `abline()`
  ## * `...` w/i `plot()` and `autoplot()`
  stopifnot(is.numeric(nsim), length(nsim) == 1)
  stopifnot(is.logical(freq))
  stopifnot(is.null(breaks) || (is.numeric(breaks) && is.null(dim(breaks))))
  stopifnot(
    is.numeric(confint_level),
    length(confint_level) == 1,
    confint_level >= 0,
    confint_level <= 1
  )
  stopifnot(is.logical(single_graph))
  stopifnot(length(xlim) == 2 && (all(is.na(xlim)) || is.numeric(xlim)))
  stopifnot(length(ylim) == 2 && (all(is.na(ylim)) || is.numeric(ylim)))
  stopifnot(length(xlab) == 1)
  stopifnot(length(ylab) == 1)
  stopifnot(length(main) == 1 | length(main) == 0)

  ## match arguments
  style <- match.arg(style)
  confint_type <- match.arg(confint_type)
  type <- match.arg(type)

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
  # COMPUTATION OF PIT
  # -------------------------------------------------------------------
  ## either compute proportion exactly (to do...) or approximate by simulation
  if (type == "proportional") {
    ## FIXME: (ML)
    ## * Implement proportional over the inteverals (e.g., below censoring point)
    ## * confusing naming, as `type` in `qresiduals()` must be `random` or `quantile`
    stop("not yet implemented")
  } else {
    p <- qresiduals.default(object,
      newdata = newdata, trafo = NULL, type = "random",
      nsim = nsim, delta = delta
    )
  }

  ## get breaks
  if (is.null(breaks)) breaks <- c(4, 10, 20, 25)[cut(NROW(p), c(0, 50, 5000, 1000000, Inf))]
  if (length(breaks) == 1L) breaks <- seq(0, 1, length.out = breaks + 1L)
  ## FIXME: (ML) Maybe use xlim instead or `0` and `1`

  ## compute ci interval
  if (confint_type == "exact") {
    ci <- get_confint(NROW(p), length(breaks) - 1, confint_level, freq)
  } else {
    ci <- get_confint_agresti(
      NROW(p) / (length(breaks) - 1),
      NROW(p),
      confint_level,
      length(breaks) - 1,
      freq
    )
  }

  ## perfect prediction
  pp <- ifelse(freq, NROW(p) / (length(breaks) - 1), 1)

  ## labels
  if (is.null(main)) main <- deparse(substitute(object))

  # -------------------------------------------------------------------
  # OUTPUT AND OPTIONAL PLOTTING
  # -------------------------------------------------------------------
  ## collect everything as data.frame
  tmp_hist <- hist(p, breaks = breaks, plot = FALSE)
  ## TODO: (ML) Maybe get rid of `hist()`
  if (freq) {
    rval <- data.frame(
      x = tmp_hist$mids,
      y = tmp_hist$counts,
      width = diff(tmp_hist$breaks),
      ci_lwr = ci[1],
      ci_upr = ci[2],
      ref = pp
    )
  } else {
    rval <- data.frame(
      x = tmp_hist$mids,
      y = tmp_hist$density,
      width = diff(tmp_hist$breaks),
      ci_lwr = ci[1],
      ci_upr = ci[2],
      ref = pp
    )
  }

  ## attributes for graphical display
  attr(rval, "freq") <- freq
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "confint_level") <- ifelse(confint, confint_level, NA)

  ## add class
  if (class == "data.frame") {
    class(rval) <- c("pithist", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("pithist", class(rval))
  }

  ## plot by default
  if (plot == "ggplot2") {
    try(print(ggplot2::autoplot(rval,
      style = style, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, confint = confint, ...
    )))
  } else if (plot == "base") {
    try(plot(rval,
      style = style, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, confint = confint, ...
    ))
  }

  ## return invisibly
  invisible(rval)
}


#' @export
c.pithist <- function(...) {
  # -------------------------------------------------------------------
  # GET DATA
  # -------------------------------------------------------------------
  ## list of pithists
  rval <- list(...)

  ## set class to tibble if any rval is a tibble
  if (any(do.call("c", lapply(rval, class)) %in% "tbl")) {
    class <- "tibble"
  } else {
    class <- "data.frame"
  }

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
  ## check if all of same `freq`
  freq <- unlist(lapply(rval, function(r) attr(r, "freq")))
  stopifnot(length(unique(freq)) == 1)

  ## labels
  xlab <- unlist(lapply(rval, function(r) attr(r, "xlab")))
  ylab <- unlist(lapply(rval, function(r) attr(r, "ylab")))
  confint_level <- unlist(lapply(rval, function(r) attr(r, "confint_level")))
  nam <- names(rval)
  main <- if (is.null(nam)) {
    as.vector(sapply(rval, function(r) attr(r, "main")))
  } else {
    make.unique(rep.int(nam, sapply(n, length)))
  }
  n <- unlist(n)

  # -------------------------------------------------------------------
  # RETURN DATA
  # -------------------------------------------------------------------
  ## combine and return
  rval <- do.call("rbind.data.frame", rval)
  rval$group <- if (length(n) < 2L) NULL else rep.int(seq_along(n), n)
  attr(rval, "freq") <- freq
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "confint_level") <- confint_level

  ## set class to data.frame or tibble
  if (class == "data.frame") {
    class(rval) <- c("pithist", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("pithist", class(rval))
  }

  ## return
  return(rval)
}


#' @export
rbind.pithist <- c.pithist


#' Plotting a PIT Histogram
#' 
#' PIT histograms graphically compare empirical probabilities from fitted
#' models with a uniform distribution.
#' 
#' PIT histograms graphically the probability integral transform (PIT), i.e.,
#' observed probabilities from fitted probability models, with a uniform
#' distribution. It leverages the \code{\link{procast}} generic and then
#' essentially draws a \code{\link[graphics]{hist}}.
#' 
#' In case of discrete distributions the PIT is either drawn randomly from the
#' corresponding interval or distributed proportionally in the histogram
#' (FIXME: not yet implemented).
#' 
#' @aliases plot.pithist lines.pithist autoplot.pithist
#' @param object,x an object of class \code{pithist}.
#' @param single_graph logical. Should all computed extended reliability
#' diagrams be plotted in a single graph?
#' @param style character specifying the syle of rootogram (see below). FIXME:
#' Description
#' @param confint logical. Should confident intervals be drawn?
#' @param xlim,ylim graphical parameters. These may pertain either to the whole
#' plot or just the histogram or just the fitted line.
#' @param xlab,ylab,main graphical parameters.
#' @param \dots further graphical parameters.
#' @param ref,col,fill,border,alpha_min,lwd,lty,axes,box additional graphical
#' parameters for base plots, whereby \code{x} is a object of class \code{pithist}.
#' @param colour,size,linetype,legend graphical parameters passed for 
#' \code{ggplot2} style plots, whereby \code{object} is a object of class \code{pithist}.
#' @seealso \code{\link{procast}}, \code{\link[graphics]{hist}}
#' @references Czado C, Gneiting T, Held L (2009). \dQuote{Predictive Model
#' Assessment for Count Data.} \emph{Biometrics}, \bold{65}(4), 1254--1261.
#' 
#' Agresti A, Coull A B (1998). \dQuote{Approximate is Better than ``Exact''
#' for Interval Estimation of Binomial Proportions.} \emph{The American
#' Statistician}, \bold{52}(2), 119--126.
#' @keywords hplot
#' @examples
#' 
#' require("crch")
#' m1 <- lm(dist ~ speed, data = cars)
#' m2 <- crch(dist ~ speed | speed, data = cars)
#' m3 <- crch(dist ~ speed | speed, left = 30, data = cars)
#' 
#' pit1 <- pithist(m1)
#' pit2 <- pithist(m2, plot = FALSE)
#' pit3 <- pithist(m3, plot = FALSE)
#' 
#' plot(pit1, confint = "red", ref = "blue", fill = "lightblue")
#' 
#' plot(c(pit1, pit2, pit3), col = c(1, 2, 3), style = "lines")
#' 
#' plot(c(pit1, pit2), col = c(1, 2), single_graph = TRUE)
#' lines(pit3, col = 3)
#' 
#' @export
plot.pithist <- function(x,
                         single_graph = FALSE,
                         style = c("histogram", "lines"),
                         confint = TRUE,
                         ref = TRUE,
                         xlim = c(0, 1),
                         ylim = c(0, NA),
                         xlab = NULL,
                         ylab = NULL,
                         main = NULL,
                         col = "black",
                         fill = adjustcolor("black", alpha.f = 0.2),
                         border = "black",
                         alpha_min = 0.2,
                         lwd = NULL,
                         lty = 1,
                         axes = TRUE,
                         box = TRUE,
                         ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * lengths of all arguments are checked by recycling
  ## * `ref` and `confint` w/i `abline()`
  ## * `xlab`, `ylab`, `main`, `col`, `fill`, `lwd`, `lty` and `...` w/i `plot()`
  ## * `alpha_min` w/i `set_minimum_transparency()`
  stopifnot(is.logical(single_graph))
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))
  if (single_graph) {
    stopifnot(
      "for `single_graph` all `freq` in attr of object `x` must be of the same type" =
        length(unique(attr(x, "freq"))) == 1
    )
  }
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## set style
  style <- match.arg(style)
  if (n > 1 && single_graph && style == "histogram") {
    message(" * For several histograms in a single graph solely line style histograms can be plotted. \n * For proper usage, set `style` = 'lines' when numbers of histograms greater one and `single_graph` = TRUE.")
    style <- "lines"
  }

  ## set lwd
  if (is.null(lwd)) lwd <- if (style == "histogram") 1.5 else 2

  ## recycle arguments for plotting to match the number of groups
  if (is.list(xlim)) xlim <- as.data.frame(do.call("rbind", xlim))
  if (is.list(ylim)) ylim <- as.data.frame(do.call("rbind", ylim))
  plot_arg <- data.frame(1:n, confint, ref,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]],
    border, col, fill, alpha_min, lwd, lty, axes, box
  )[, -1]

  ## annotation
  if (single_graph) {
    if (is.null(xlab)) xlab <- "PIT"
    if (is.null(ylab)) ylab <- if (all(attr(x, "freq"))) "Frequency" else "Density"
    if (is.null(main)) main <- "PIT histogram"
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
  # MAIN PLOTTING FUNCTION FOR 'HISTOGRAM-STYLE PITHIST'
  # -------------------------------------------------------------------
  pithist_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## rect elements
    xleft <- d$x - d$width / 2
    xright <- d$x + d$width / 2
    y <- d$y

    ## prepare confint lines and check if they consist of unique values
    if (!identical(plot_arg$confint[j], FALSE)) {
      ci_lwr <- unique(d$ci_lwr)
      ci_upr <- unique(d$ci_upr)
      stopifnot(
        "`ci_lwr` and `ci_upr` in attr of object `x` must consist of unique values per group index" =
          length(ci_lwr) == 1 & length(ci_upr) == 1
      )
    } else {
      ci_lwr <- NULL
      ci_upr <- NULL
    }

    ## prepare ref line and check if it consists of unique values
    pp <- unique(d$ref)
    stopifnot(
      "`ref` in attr of object `x` must consist of unique values per group index" =
        length(pp) == 1
    )

    ## get xlim and ylim
    ylim_idx <- c(is.na(plot_arg$ylim1[j]), is.na(plot_arg$ylim2[j]))
    xlim_idx <- c(is.na(plot_arg$xlim1[j]), is.na(plot_arg$xlim2[j]))
    if (any(xlim_idx)) {
      plot_arg[j, c("xlim1", "xlim2")[xlim_idx]] <- range(c(xleft, xright))[xlim_idx]
    }
    if (any(ylim_idx)) {
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <- range(c(0, y, ci_lwr, ci_upr))[ylim_idx]
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

    if (plot_arg$fill[j] == adjustcolor("black", alpha.f = 0.2) &&
      plot_arg$col[j] != "black") {
      message(" * As the argument `col` is set but no argument `fill` is specified, \n   the former is used for colorizing the PIT histogram. \n * For proper usage, solely provide `fill` for histogram style plots.")
      plot_arg$fill[j] <- plot_arg$col[j]
    }

    ## plot pithist
    rect(xleft, 0, xright, y,
      border = plot_arg$border[j], col = plot_arg$fill[j],
      lty = plot_arg$lty[j]
    )

    ## plot ref line
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- 2  # red
      abline(h = pp, col = plot_arg$ref[j], lty = 1, lwd = plot_arg$lwd[j])
    }

    ## plot confint lines
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- 2  # red
      abline(h = ci_lwr, col = plot_arg$confint[j], lty = 2, lwd = plot_arg$lwd[j])
      abline(h = ci_upr, col = plot_arg$confint[j], lty = 2, lwd = plot_arg$lwd[j])
    }
  }


  # -------------------------------------------------------------------
  # FUNCTION TO TRIGGER FIGURE AND PLOT CONFINT FOR 'LINE-STYLE PITHIST'
  # -------------------------------------------------------------------
  pitlines_trigger <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## step elements
    z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
    y <- c(d$y, d$y[NROW(d)])

    ## prepare confint lines and check if they consist of unique values
    if (!identical(plot_arg$confint[j], FALSE)) {
      ci_lwr <- unique(d$ci_lwr)
      ci_upr <- unique(d$ci_upr)
      stopifnot(
        "`ci_lwr` and `ci_upr` in attr of object `x` must consist of unique values per group index" =
          length(ci_lwr) == 1 & length(ci_upr) == 1
      )
    } else {
      ci_lwr <- NULL
      ci_upr <- NULL
    }

    ## get xlim and ylim
    ylim_idx <- c(is.na(plot_arg$ylim1[j]), is.na(plot_arg$ylim2[j]))
    xlim_idx <- c(is.na(plot_arg$xlim1[j]), is.na(plot_arg$xlim2[j]))
    if (any(xlim_idx)) {
      plot_arg[j, c("xlim1", "xlim2")[xlim_idx]] <- range(z)[xlim_idx]
    }
    if (any(ylim_idx) && !single_graph) {
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <- range(c(y, ci_lwr, ci_upr))[ylim_idx]
    }
    if (any(ylim_idx) && single_graph) {
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <- range(c(x$y, x$ci_lwr, x$ci_upr))[ylim_idx]
    }

    ## trigger plot
    if (j == 1 || (!single_graph && j > 1)) {
      plot(0, 0,
        type = "n", xlim = c(plot_arg$xlim1[j], plot_arg$xlim2[j]),
        ylim = c(plot_arg$ylim1[j], plot_arg$ylim2[j]),
        xlab = xlab[j], ylab = ylab[j], xaxs = "i", main = main[j], axes = FALSE, ...
      )
      if (plot_arg$axes[j]) {
        axis(1)
        axis(2)
      }
      if (plot_arg$box[j]) {
        box()
      }
    }

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]
      polygon(
        c(0, 1, 1, 0),
        c(ci_lwr, ci_lwr, ci_upr, ci_upr),
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min = plot_arg$alpha_min[j]),
        border = NA
      )
    }
  }

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR 'LINE_STYLE PITHIST'
  # -------------------------------------------------------------------
  pitlines_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## step elements
    z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
    y <- c(d$y, d$y[NROW(d)])

    ## prepare ref line and check if it consists of unique values
    pp <- unique(d$ref)
    stopifnot(
      "`ref` in attr of object `x` must consist of unique values per group index" =
        length(pp) == 1
    )

    ## plot ref line
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
      segments(x0 = 0, y0 = pp, x1 = 1, y1 = pp, col = plot_arg$ref[j], lty = 2, lwd = 1)
    }

    ## plot stepfun
    lines(y ~ z, type = "s", lwd = plot_arg$lwd[j], lty = plot_arg$lty[j], col = plot_arg$col[j])
  }

  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  ## set up necessary panels
  if (!single_graph && n > 1L) {
    old_pars <- par(mfrow = n2mfrow(n))
    on.exit(par(old_pars), add = TRUE)
  }

  ## draw polygons first
  if (single_graph || n == 1) {
    if (style == "histogram") {
      for (i in 1L:n) pithist_plot(x[x$group == i, ], ...)
    } else {
      for (i in 1L:n) pitlines_trigger(x[x$group == i, ], ...)
      for (i in 1L:n) pitlines_plot(x[x$group == i, ], ...)
    }
  } else {
    if (style == "histogram") {
      for (i in 1L:n) pithist_plot(x[x$group == i, ], ...)
    } else {
      for (i in 1L:n) {
        pitlines_trigger(x[x$group == i, ], ...)
        pitlines_plot(x[x$group == i, ], ...)
      }
    }
  }
}


#' @rdname plot.pithist
#' @method lines pithist
#' @export
lines.pithist <- function(x,
                          confint = FALSE,
                          ref = FALSE,
                          col = "black",
                          fill = adjustcolor("black", alpha.f = 0.2),
                          alpha_min = 0.2,
                          lwd = 2,
                          lty = 1,
                          ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * lengths of all arguments are checked by recycling
  ## * `ref` and `confint` w/i `abline()`
  ## * `col`, `fill`, `lwd`, `lty` and `...` w/i `lines()`
  ## * `alpha_min` w/i `set_minimum_transparency()`
  stopifnot(
    "all `freq` in attr of object `x` must be of the same type" =
      length(unique(attr(x, "freq"))) == 1
  )

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(
    1:n, confint, ref, col, fill, alpha_min, lwd, lty
  )[, -1]

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR LINES
  # -------------------------------------------------------------------
  pitlines_lines <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## step elements
    x <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
    y <- c(d$y, d$y[NROW(d)])

    ## prepare confint lines and check if they consist of unique values
    if (!identical(plot_arg$confint[j], FALSE)) {
      ci_lwr <- unique(d$ci_lwr)
      ci_upr <- unique(d$ci_upr)
      stopifnot(
        "`ci_lwr` and `ci_upr` in attr of object `x` must consist of unique values per group index" =
          length(ci_lwr) == 1 & length(ci_upr) == 1
      )
    } else {
      ci_lwr <- NULL
      ci_upr <- NULL
    }

    ## prepare ref line and check if it consists of unique values
    pp <- unique(d$ref)
    stopifnot(
      "`ref` in attr of object `x` must consist of unique values per group index" =
        length(pp) == 1
    )

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]
      polygon(
        c(0, 1, 1, 0),
        c(ci_lwr, ci_lwr, ci_upr, ci_upr),
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min = plot_arg$alpha_min[j]),
        border = NA
      )
    }

    ## plot ref line
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- 2  # red
      segments(x0 = 0, y0 = pp, x1 = 1, y1 = pp, col = plot_arg$ref[j], lty = 2, lwd = 1.5)
    }

    ## plot stepfun
    lines.default(y ~ x, type = "s", lwd = plot_arg$lwd[j], lty = plot_arg$lty[j], col = plot_arg$col[j])
  }

  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  for (i in 1L:n) {
    pitlines_lines(x[x$group == i, ], ...)
  }
}


#' @rdname plot.pithist
#' @method autoplot pithist
#' @exportS3Method ggplot2::autoplot
autoplot.pithist <- function(object,
                             single_graph = FALSE,
                             style = c("histogram", "lines"),
                             confint = TRUE,
                             ref = TRUE,
                             xlim = c(0, 1),
                             ylim = c(0, NA),
                             xlab = NULL,
                             ylab = NULL,
                             main = NULL,
                             colour = "black",
                             fill = "darkgray",
                             border = "black",
                             alpha_min = 0.2,
                             size = NULL,
                             linetype = 1,
                             legend = FALSE,
                             ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## get base style arguments
  add_arg <- list(...)
  if (!is.null(add_arg$lwd)) size <- add_arg$lwd
  if (!is.null(add_arg$lty)) linetype <- add_arg$lty

  ## sanity checks
  stopifnot(is.logical(single_graph))
  if (single_graph) {
    stopifnot(
      "for `single_graph` all `freq` in attr of `object` must be of the same type" =
        length(unique(attr(object, "freq"))) == 1
    )
  }
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))

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
  ## determine which style should be plotted
  style <- match.arg(style)
  if (n > 1 && single_graph && style == "histogram") {
    message(" * For several histograms in a single graph solely line style histograms can be plotted. \n * For proper usage, set `style` = 'lines' when numbers of histograms greater one and `single_graph` = TRUE.")
    style <- "lines"
  }

  ## set size
  if (is.null(size)) size <- if (style == "histogram") 0.7 else 1

  ## set color to 2 (red) or NA for not plotting
  if (is.logical(ref)) ref <- ifelse(ref, 2, NA)

  ## only needed for `style == "lines"`
  ## stat helper function to get left/right points from respective mid points
  calc_pit_points <- ggplot2::ggproto("calc_pit_points", ggplot2::Stat,

    # required as we operate on groups (facetting)
    compute_group = function(data, scales) {
      ## manipulate object
      nd <- data.frame(
        x = c(data$x - data$width / 2, data$x[NROW(data)] + data$width[NROW(data)] / 2),
        y = c(data$y, data$y[NROW(data)])
      )
      nd
    },

    # tells us what we need
    required_aes = c("x", "y")
  )

  ## recycle arguments for plotting to match the number of groups (for `scale_<...>_manual()`)
  plot_arg <- data.frame(
    1:n,
    fill, colour, size, linetype, confint, alpha_min
  )[, -1]

  ## prepare fill and confint color depending on style
  if (style == "histogram") {

    ## set color to 2 (red) or NA for not plotting
    if (is.logical(confint)) confint <- ifelse(confint, 2, NA)

    ## check if colour and no fill is set
    if (all(fill == "darkgray") && any(colour != "black")) {
      message(" * As the argument `colour` is set but no argument `fill` is specified, \n   the former is used for colorizing the PIT histogram. \n * For proper usage, solely provide `fill` for histogram style plots.")
      plot_arg$fill <- plot_arg$colour
    }
  } else {
    if (is.logical(plot_arg$confint)) {

      ## use fill and set alpha
      plot_arg$fill <- sapply(seq_along(plot_arg$fill), function(idx) {
        set_minimum_transparency(plot_arg$fill[idx], alpha_min = plot_arg$alpha_min[idx])
      })

      ## set color to NA for not plotting
      plot_arg$fill[!plot_arg$confint] <- NA
    } else {
      ## use confint and set alpha
      plot_arg$fill <- sapply(seq_along(plot_arg$confint), function(idx) {
        set_minimum_transparency(plot_arg$confint[idx], alpha_min = plot_arg$alpha_min[idx])
      })
    }
  }

  ## recycle arguments for plotting to match the length (rows) of the object (for geom w/ aes)
  plot_arg2 <- data.frame(1:n, border, colour, ref, confint)[, -1]
  plot_arg2 <- as.data.frame(lapply(plot_arg2, rep, table(object$group)))

  # -------------------------------------------------------------------
  # MAIN PLOTTING
  # -------------------------------------------------------------------
  if (style == "histogram") {

    ## actual plotting
    rval <- ggplot2::ggplot(
      object,
      ggplot2::aes_string(x = "x", y = "y / 2", width = "width", height = "y")
    ) +
      ggplot2::geom_tile(ggplot2::aes_string(fill = "group"),
        colour = plot_arg2$border
      ) +
      ggplot2::geom_hline(ggplot2::aes_string(yintercept = "ref", size = "group"),
        colour = plot_arg2$ref
      ) +
      ggplot2::geom_hline(ggplot2::aes_string(yintercept = "ci_lwr", size = "group"),
        colour = plot_arg2$confint, linetype = 2
      ) +
      ggplot2::geom_hline(ggplot2::aes_string(yintercept = "ci_upr", size = "group"),
        colour = plot_arg2$confint, linetype = 2
      )

    ## set the colors, shapes, etc. for the groups
    rval <- rval +
      ggplot2::scale_fill_manual(values = plot_arg$fill) +
      ggplot2::scale_size_manual(values = plot_arg$size)

    ## annotation
    rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

    ## add legend
    if (legend) {
      rval <- rval + ggplot2::labs(fill = "Model") +
        ggplot2::guides(fill = "legend", size = "none")
    } else {
      rval <- rval + ggplot2::guides(fill = "none", size = "none")
    }

    ## set x and y limits
    rval <- rval + ggplot2::scale_x_continuous(limits = xlim, expand = c(0.01, 0.01))
    rval <- rval + ggplot2::scale_y_continuous(limits = ylim, expand = c(0.01, 0.01))
  } else {

    ## actual plotting
    rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y", width = "width")) +
      ggplot2::geom_rect(ggplot2::aes_string(ymin = "ci_lwr", ymax = "ci_upr", fill = "group"),
        xmin = 0, xmax = 1, colour = NA, show.legend = FALSE
      ) +
      ggplot2::geom_step(ggplot2::aes_string(colour = "group", size = "group", linetype = "group"),
        stat = calc_pit_points
      )

    ## set the colors, shapes, etc.
    rval <- rval +
      ggplot2::scale_colour_manual(values = plot_arg$colour) +
      ggplot2::scale_fill_manual(values = plot_arg$fill) +
      ggplot2::scale_size_manual(values = plot_arg$size) +
      ggplot2::scale_linetype_manual(values = plot_arg$linetype)

    ## add annotation
    rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

    ## add legend
    if (legend) {
      rval <- rval + ggplot2::labs(colour = "Model") +
        ggplot2::guides(colour = "legend", size = "none", linetype = "none")
    } else {
      rval <- rval + ggplot2::guides(colour = "none", size = "none", linetype = "none")
    }

    ## set x and y limits
    rval <- rval + ggplot2::scale_x_continuous(limits = xlim, expand = c(0.01, 0.01))
    rval <- rval + ggplot2::scale_y_continuous(limits = ylim, expand = c(0.01, 0.01))
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


get_confint <- function(n, bins, level, freq) {
  ## helper function to calculate CI employing `qbinom()`
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  rval <- qbinom(a, size = n, prob = 1 / bins)
  if (!freq) rval <- rval / (n / bins)
  rval
}


get_confint_agresti <- function(x, n, level, bins, freq) {
  ## helper function to calculate an approximated CI according to Agresti & Coull (1998)
  ## doi=10.1080/00031305.1998.10480550
  rval <- add4ci(x, n, level)$conf.int * n
  if (!freq) rval <- rval / (n / bins)
  rval
}


add4ci <- function(x, n, conf.level) {
  ## copy of `add4ci` package from package `PropCIs` by Ralph Scherer (licensed under GPL-2/GPL-3)
  ptilde <- (x + 2) / (n + 4)
  z <- abs(qnorm((1 - conf.level) / 2))
  stderr <- sqrt(ptilde * (1 - ptilde) / (n + 4))
  ul <- ptilde + z * stderr
  ll <- ptilde - z * stderr
  if (ll < 0) {
    ll <- 0
  }
  if (ul > 1) {
    ul <- 1
  }
  cint <- c(ll, ul)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint, estimate = ptilde)
  class(rval) <- "htest"
  return(rval)
}
