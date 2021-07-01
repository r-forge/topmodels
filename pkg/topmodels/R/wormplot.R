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
#' Visualize goodness of fit of regression models by Worm plots using quantile
#' residuals.
#' 
#' FIXME: Q-Q residual draw quantile residuals (by default: transformed to
#' standard normal scale) against theoretical quantiles from the same
#' distribution. Alternatively, transformations to other distributions can also
#' be used, specifically using no transformation at all, i.e., remaining on the
#' uniform scale (via \code{trafo = NULL} or equivalently \code{qunif} or
#' \code{identity}).
#' 
#' Additional options are offered for models with discrete responses where
#' randomization of quantiles is needed.
#' 
#' @aliases wormplot wormplot.default c.wormplot
#' @param object an object (for which a \code{\link{qresiduals}} method
#' exists).
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param plot Should the \code{plot} or \code{autoplot} method be called to
#' draw the computed worm plot? Either set \code{plot} expicitly to "base" vs.
#' "ggplot2" to choose the type of plot, or for a logical \code{plot} argument
#' it's chosen conditional if the package \code{ggplot2} is loaded.
#' @param class Should the invisible return value be either a \code{data.frame}
#' or a \code{tibble}. Either set \code{class} expicitly to "data.frame" vs.
#' "tibble", or for NULL it's chosen automatically conditional if the package
#' \code{tibble} is loaded.
#' @param trafo function for tranforming residuals from probability scale to a
#' different distribution scale (default: Gaussian).
#' @param nsim,delta arguments passed to \code{qresiduals}.
#' @param confint logical or quantile specification. Should the range of
#' quantiles of the randomized quantile residuals be visualized? If
#' \code{TRUE}, then \code{range = c(0.01, 0.99)} is used.
#' @param confint_level numeric. The confidence level required.
#' @param confint_nsim numeric. The number of simulated quantiles.
#' @param confint_seed numeric. The seed to be set for calculating the
#' confidence interval.
#' @param single_graph logical. Should all computed extended reliability
#' diagrams be plotted in a single graph?
#' @param xlab,ylab,main,\dots graphical plotting parameters passed to
#' \code{\link[graphics]{plot}} or \code{\link[graphics]{points}},
#' respectively.
#' @return An list is returned invisibly with: \item{normal}{the theoretical
#' normal quantiles,} \item{residuals}{the empirical quantile residuals.}
#' @seealso \code{\link{qresiduals}}, \code{\link[stats]{qqnorm}}
#' @references FIXME: Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile
#' Residuals.} \emph{Journal of Computational and Graphical Statistics},
#' \bold{5}, 1--10.
#' @keywords hplot
#' @examples
#' 
#' data("CrabSatellites", package = "countreg")
#' CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]
#' 
#' m1 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#' m2 <- glm(satellites ~ width + color, data = CrabSatellites2, family = binomial)
#' m3 <- lm(dist ~ speed, data = cars)
#' 
#' w1 <- wormplot(m1, nsim = 100, confint = TRUE)
#' w2 <- wormplot(m2, nsim = 100, confint = TRUE, plot = FALSE)
#' w3 <- wormplot(m3, nsim = 100, confint = TRUE, plot = FALSE)
#' 
#' plot(c(w1, w2), single_graph = FALSE, fill = c(1, 3), ref = c(2, 2))
#' points(w3, col = "lightblue")
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
                             trafo = qnorm,
                             nsim = 1L,
                             delta = NULL,
                             confint = TRUE,
                             confint_level = 0.95,
                             confint_nsim = 250,
                             confint_seed = 1,
                             single_graph = FALSE,
                             xlab = "Theoretical quantiles",
                             ylab = "Deviation",
                             main = NULL,
                             ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * `object`, `newdata`, `delta w/i `qresiduals()`
  ## * `confint` w/i `polygon()`
  ## * `delta` w/i `qresiduals()`
  ## * `...` in `plot()` and `autoplot()`
  stopifnot(is.null(trafo) | is.function(trafo))
  stopifnot(is.numeric(nsim), length(nsim) == 1)
  stopifnot(
    is.numeric(confint_level),
    length(confint_level) == 1,
    confint_level >= 0,
    confint_level <= 1
  )
  stopifnot(is.logical(single_graph))
  stopifnot(length(xlab) == 1)
  stopifnot(length(ylab) == 1)
  stopifnot(length(main) == 1 | length(main) == 0)

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

  ## compute ci interval
  ## FIXME: (ML) Implement exact method if exists (see "inst/misc/2021_04_16_errorsearch_qqrplot.Rmd")
  if (!identical(confint, FALSE)) {
    set.seed(confint_seed)
    tmp <- qresiduals(object,
      newdata = newdata, trafo = trafo, type = "random", nsim = confint_nsim,
      delta = delta
    )
    confint_prob <- (1 - confint_level) / 2
    confint_prob <- c(confint_prob, 1 - confint_prob)
    qres_ci_lwr <- apply(apply(tmp, 2, sort), 1, quantile, prob = confint_prob[1], na.rm = TRUE)
    qres_ci_upr <- apply(apply(tmp, 2, sort), 1, quantile, prob = confint_prob[2], na.rm = TRUE)
    qthe_ci_lwr <- q2q(qres_ci_lwr)
    qthe_ci_upr <- q2q(qres_ci_upr)

    ## FIXME: (ML) Improve workaround to get CI only for discrete values
    if (isTRUE(all.equal(qres_ci_lwr, qres_ci_upr, tol = .Machine$double.eps^0.4))) {
      qres_ci_lwr <- NULL
      qres_ci_upr <- NULL
      qthe_ci_lwr <- NULL
      qthe_ci_upr <- NULL
      confint <- FALSE
    }
  } else {
    qres_ci_lwr <- NULL
    qres_ci_upr <- NULL
    qthe_ci_lwr <- NULL
    qthe_ci_upr <- NULL
  }

  ## setup function to calculate ref lines:
  ## * ci interval according to Van Buuren and Fredriks (2001) p. 1276
  ## * all data should be prepared in `wormplot.default()`, but length does not fit to rval
  ## FIXME: (ML) Adapt for other trafos (density and quantile functions): Order statistics.
  ref_fun <- function(x, n, level = 0.95, which = c("lower", "upper")) {
    stopifnot(is.numeric(n), length(n) == 1)
    stopifnot(is.numeric(level), length(level) == 1, level >= 0, level <= 1)
    which <- match.arg(which)

    p <- pnorm(x)
    se <- (1 / dnorm(x)) * (sqrt(p * (1 - p) / n))
    rval <- as.numeric(qnorm((1 - level) / 2) * se)

    if (which == "lower") {
      rval
    } else {
      -rval
    }
  }

  ## labels
  if (is.null(main)) main <- deparse(substitute(object))

  # -------------------------------------------------------------------
  # OUTPUT AND OPTIONAL PLOTTING
  # -------------------------------------------------------------------
  ## collect everything as data.frame
  if (any(vapply(
    list(qres_ci_lwr, qres_ci_upr, qthe_ci_lwr, 1),
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
      y_ci_lwr = qres_ci_lwr - qthe_ci_lwr,
      y_ci_upr = qres_ci_upr - qthe_ci_upr,
      x_ci_lwr = qthe_ci_lwr,
      x_ci_upr = qthe_ci_lwr
    )
  }
  names(rval) <- gsub("(\\.r|\\.q)", "", names(rval))

  ## attributes for graphical display
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "ref_fun") <- ref_fun
  attr(rval, "confint_level") <- ifelse(confint, confint_level, NA)

  if (class == "data.frame") {
    class(rval) <- c("wormplot", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("wormplot", class(rval))
  }

  ## plot by default
  if (plot == "ggplot2") {
    try(print(ggplot2::autoplot(rval, confint = confint, ...)))
  } else if (plot == "base") {
    try(plot(rval, confint = confint, ...))
  }

  ## return coordinates invisibly
  invisible(rval)
}


#' @export
c.wormplot <- function(...) {
  # -------------------------------------------------------------------
  # GET DATA
  # -------------------------------------------------------------------
  ## list of wormplots
  rval <- list(...)

  ## set class to tibble if any rval is a tibble
  if (any(do.call("c", lapply(rval, class)) %in% "tbl")) {
    class <- "tibble"
  } else {
    class <- "data.frame"
  }

  ## remove temporary the class (needed below for `c()`)
  ## FIXME: (ML) Rewrite by, e.g., employing `lapply()`
  for (i in 1:length(rval)) class(rval[[i]]) <- class(rval[[i]])[!class(rval[[i]]) %in% "wormplot"]

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
  confint_level <- unlist(lapply(rval, function(r) attr(r, "confint_level")))
  ref_fun <- unlist(lapply(rval, function(r) attr(r, "ref_fun")))
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
  attr(rval, "confint_level") <- confint_level
  attr(rval, "ref_fun") <- ref_fun

  ## set class to data.frame or tibble
  if (class == "data.frame") {
    class(rval) <- c("wormplot", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("wormplot", class(rval))
  }

  ## return
  return(rval)
}


#' @export
rbind.wormplot <- c.wormplot


#' Plotting Worm Plots for Quantile Residuals
#' 
#' Visualize goodness of fit of regression models by Worm plots using quantile
#' residuals.
#' 
#' FIXME: Q-Q residual draw quantile residuals (by default: transformed to
#' standard normal scale) against theoretical quantiles from the same
#' distribution. Alternatively, transformations to other distributions can also
#' be used, specifically using no transformation at all, i.e., remaining on the
#' uniform scale (via \code{trafo = NULL} or equivalently \code{qunif} or
#' \code{identity}).
#' 
#' Additional options are offered for models with discrete responses where
#' randomization of quantiles is needed.
#' 
#' @aliases plot.wormplot points.wormplot autoplot.wormplot
#' @param x,object an object of class \code{reliagram}.
#' @param single_graph logical. Should all computed extended reliability
#' diagrams be plotted in a single graph?
#' @param confint logical or quantile specification. Should the range of
#' quantiles of the randomized quantile residuals be visualized? If
#' \code{TRUE}, then \code{range = c(0.01, 0.99)} is used.
#' @param xlab,ylab,main,\dots graphical plotting parameters passed to
#' \code{\link[graphics]{plot}} or \code{\link[graphics]{points}},
#' respectively.
#' @param ref,xlim,ylim,col,fill,alpha_min,pch,axes,box additional graphical
#' parameters for base plots, whereby \code{x} is a object of class \code{pithist}.
#' @param colour,size,shape,linetype,legend graphical parameters passed for 
#' \code{ggplot2} style plots, whereby \code{object} is a object of class \code{pithist}.
#' @return An list is returned invisibly with: \item{normal}{the theoretical
#' normal quantiles,} \item{residuals}{the empirical quantile residuals.}
#' @seealso \code{\link{qresiduals}}, \code{\link[stats]{qqnorm}}
#' @references FIXME: Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile
#' Residuals.} \emph{Journal of Computational and Graphical Statistics},
#' \bold{5}, 1--10.
#' @keywords hplot
#' @examples
#' 
#' ## speed and stopping distances of cars
#' m1_lm <- lm(dist ~ speed, data = cars)
#' 
#' ## compute and plot wormplot
#' wormplot(m1_lm)
#' 
#' ## customize colors
#' wormplot(m1_lm, ref = "blue", lty = 2, pch = 20)
#' 
#' ## add separate model
#' if (require("crch", quietly = TRUE)) {
#'   m1_crch <- crch(dist ~ speed | speed, data = cars)
#'   points(wormplot(m1_crch, plot = FALSE), col = 2, lty = 2, confint = 2)
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
#'   ## compute wormplots
#'   wp2_lm <- wormplot(m2_lm, plot = FALSE)
#'   wp2_crch <- wormplot(m2_crch, plot = FALSE)
#' 
#'   ## plot in single graph
#'   plot(c(wp2_lm, wp2_crch), col = c(1, 2), confint = c(1, 2), ref = 3, single_graph = TRUE)
#' }
#'
#' #-------------------------------------------------------------------------------
#' ## determinants for male satellites to nesting horseshoe crabs
#' data("CrabSatellites", package = "countreg")
#' 
#' ## linear poisson model
#' m3_pois  <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#' 
#' ## compute and plot rootogram as ggplot2 graphic
#' wormplot(m3_pois, plot = "ggplot2")
#'
#' @export
plot.wormplot <- function(x,
                          single_graph = FALSE,
                          confint = TRUE,
                          ref = TRUE,
                          xlim = c(NA, NA),
                          ylim = c(NA, NA),
                          xlab = NULL,
                          ylab = NULL,
                          main = NULL,
                          col = adjustcolor("black", alpha.f = 0.4),
                          fill = adjustcolor("black", alpha.f = 0.2),
                          alpha_min = 0.2, # single or n values
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
  ## * `confint`, `fill` in `polygon()`
  ## * `alpha_min` w/i `set_minimum_transparency()`
  stopifnot(is.logical(single_graph))
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  if (is.list(xlim)) xlim <- as.data.frame(do.call("rbind", xlim))
  if (is.list(ylim)) ylim <- as.data.frame(do.call("rbind", ylim))
  plot_arg <- data.frame(1:n, confint, ref,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]],
    col, fill, alpha_min, pch, axes, box
  )[, -1]

  ## annotation
  if (single_graph) {
    if (is.null(xlab)) xlab <- "Theoretical quantiles"
    if (is.null(ylab)) ylab <- "Deviation"
    if (is.null(main)) main <- "Worm plot"
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
  wormplot_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get xlim and ylim conditional on confint
    if (
      !identical(plot_arg$confint[j], FALSE) &&
        !is.na(attr(d, "confint_level")[j])
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

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]

      x_pol <- c(sort(d$x_ci_lwr), sort(d$x_ci_upr, decreasing = TRUE))
      y_pol <- c(sort(d$y_ci_lwr), sort(d$y_ci_upr, decreasing = TRUE))
      x_pol[!is.finite(x_pol)] <- 100 * sign(x_pol[!is.finite(x_pol)]) # TODO: (ML) needed?
      y_pol[!is.finite(y_pol)] <- 100 * sign(y_pol[!is.finite(y_pol)]) # TODO: (ML) needed?

      polygon(
        x_pol,
        y_pol,
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min = plot_arg$alpha_min[j]),
        border = NA
      )
    }

    ## add qq plot
    for (i in 1L:ncol(d[grepl("^y$|y_[0-9]", names(d))])) {
      points.default(
        d[grepl("x", names(d))][, i],
        d[grepl("y", names(d))][, i],
        col = plot_arg$col[j], pch = plot_arg$pch[j], ...
      )
    }

    ## plot ref lines
    if (j == 1 || (!single_graph && j > 1)) {
      if (!identical(plot_arg$ref[j], FALSE)) {
        if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
        ref_fun <- attr(d, "ref_fun")
        ref_fun <- if (is.list(ref_fun)) ref_fun else list(ref_fun)
        curve(
          ref_fun[[j]](
            x,
            n = NROW(d),
            level = 0.95,
            which = "lower"
          ),
          lty = 2,
          lwd = 1.25,
          col = plot_arg$ref[j],
          from = plot_arg$xlim1[j],
          to = plot_arg$xlim2[j],
          add = TRUE
        )
        curve(
          ref_fun[[j]](
            x,
            n = NROW(d),
            level = 0.95,
            which = "upper"
          ),
          lty = 2,
          lwd = 1.25,
          col = plot_arg$ref[j],
          from = plot_arg$xlim1[j],
          to = plot_arg$xlim2[j],
          add = TRUE
        )
        abline(h = 0, col = plot_arg$ref[j], lty = 2, lwd = 1.25)
      }
    }
  }

  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  ## set up necessary panels
  if (n > 1L) {
    old_pars <- par(mfrow = n2mfrow(n))
    on.exit(par(old_pars), add = TRUE)
  }

  ## draw wormplots
  for (i in 1L:n) wormplot_plot(x[x$group == i, ], ...)
}


#' @rdname plot.wormplot
#' @method points wormplot
#' @export
points.wormplot <- function(x,
                            confint = FALSE,
                            ref = FALSE,
                            col = "black",
                            fill = adjustcolor("black", alpha.f = 0.4),
                            alpha_min = 0.2,
                            pch = 19,
                            ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * lengths of all arguments are checked by recycling
  ## * `ref` w/i `abline()`
  ## * `col`, `pch` w/i `lines()`
  ## * `confint`, `fill` in `polygon()`
  ## * `alpha_min` w/i `set_minimum_transparency()`

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(1:n, confint, ref, col, fill, alpha_min, pch)[, -1]

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR POINTS
  # -------------------------------------------------------------------
  wormplot_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get xlim and ylim conditional on confint
    if (
      !identical(plot_arg$confint[j], FALSE) &&
        !is.na(attr(d, "confint_level")[j])
    ) {
      if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) {
        xlim <- range(as.matrix(d[grepl("x", names(d))]), finite = TRUE)
      }
      if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) {
        ylim <- range(as.matrix(d[grepl("y", names(d))]), finite = TRUE)
      }
    } else {
      if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) {
        xlim <- range(as.matrix(d[grepl("^x$|x_[0-9]", names(d))]), finite = TRUE)
      }
      if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) {
        ylim <- range(as.matrix(d[grepl("^y$|y_[0-9]", names(d))]), finite = TRUE)
      }
    }

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]

      x_pol <- c(sort(d$x_ci_lwr), sort(d$x_ci_upr, decreasing = TRUE))
      y_pol <- c(sort(d$y_ci_lwr), sort(d$y_ci_upr, decreasing = TRUE))
      x_pol[!is.finite(x_pol)] <- 100 * sign(x_pol[!is.finite(x_pol)]) # TODO: (ML) needed?
      y_pol[!is.finite(y_pol)] <- 100 * sign(y_pol[!is.finite(y_pol)]) # TODO: (ML) needed?

      polygon(
        x_pol,
        y_pol,
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min = plot_arg$alpha_min[j]),
        border = NA
      )
    }

    ## plot reference diagonal
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
      abline(h = 0, col = plot_arg$ref[j], lty = 2, lwd = 1.25)
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
    wormplot_plot(x[x$group == i, ], ...)
  }
}


#' @rdname plot.wormplot
#' @method autoplot wormplot
#' @exportS3Method ggplot2::autoplot
autoplot.wormplot <- function(object,
                              single_graph = FALSE,
                              confint = TRUE,
                              ref = TRUE,
                              xlim = c(NA, NA),
                              ylim = c(NA, NA),
                              xlab = NULL,
                              ylab = NULL,
                              main = NULL,
                              colour = adjustcolor("black", alpha.f = 0.4),
                              fill = adjustcolor("black", alpha.f = 0.2),
                              alpha_min = 0.2,
                              size = 2,
                              shape = 19,
                              linetype = 1,
                              legend = FALSE,
                              ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## get base style arguments
  add_arg <- list(...)
  if (!is.null(add_arg$pch)) shape <- add_arg$pch
  if (!is.null(add_arg$lwd)) size <- add_arg$lwd
  if (!is.null(add_arg$lty)) linetype <- add_arg$lty

  ## sanity checks
  stopifnot(is.logical(single_graph))
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
  ## Get ref_fun
  ## FIXME: (ML) Uses always (only) the first `ref_fun()`
  ref_fun <- attr(object, "ref_fun")
  ref_fun <- if (is.list(ref_fun)) ref_fun[[1]] else ref_fun

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

  ## stat helper function to get ci polygon
  calc_confint_polygon <- ggplot2::ggproto("calc_confint_polygon", ggplot2::Stat,

    # Required as we operate on groups (facetting)
    compute_group = function(data, scales) {
      ## Manipulate object
      nd <- data.frame(
        x = c(data$x_ci_lwr, rev(data$x_ci_upr)),
        y = c(data$y_ci_lwr, rev(data$y_ci_upr))
      )
      nd
    },

    # Tells us what we need
    required_aes = c("x_ci_lwr", "x_ci_upr", "y_ci_lwr", "y_ci_upr")
  )

  ## recycle arguments for plotting to match the number of groups (for `scale_<...>_manual()`)
  plot_arg <- data.frame(
    1:n,
    fill, colour, size, linetype, confint, alpha_min, shape
  )[, -1]

  ## prepare fill color for confint (must be done on vector to match args)
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

  ## recycle arguments for plotting to match the length (rows) of the object (for geom w/ aes)
  plot_arg2 <- data.frame(1:n, ref)[, -1, drop = FALSE]
  plot_arg2 <- as.data.frame(lapply(plot_arg2, rep, table(object$group)))

  # -------------------------------------------------------------------
  # MAIN PLOTTING
  # -------------------------------------------------------------------
  ## actual plotting
  rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y"))

  ## add conf
  if (all(c("x_ci_lwr", "x_ci_upr", "y_ci_lwr", "y_ci_upr") %in% names(object))) {
    rval <- rval +
      ggplot2::geom_polygon(ggplot2::aes_string(
        x_ci_lwr = "x_ci_lwr", x_ci_upr = "x_ci_upr",
        y_ci_lwr = "y_ci_lwr", y_ci_upr = "y_ci_upr", fill = "group"
      ),
      stat = calc_confint_polygon, show.legend = FALSE, na.rm = TRUE
      )
  }

  ## add ref lines
  if (!identical(ref, FALSE)) {
    if (isTRUE(ref)) plot_arg2$ref <- "black"
    ## FIXME: (ML) no group specific colors work for `geom_function()`
    rval <- rval +
      ggplot2::geom_function(
        fun = ref_fun,
        args = list(n = NROW(object), level = 0.95, which = "lower"),
        linetype = 2, col = unique(plot_arg2$ref)[1]
      ) +
      ggplot2::geom_function(
        fun = ref_fun,
        args = list(n = NROW(object), level = 0.95, which = "upper"),
        linetype = 2, col = unique(plot_arg2$ref)[1]
      ) +
      ggplot2::geom_hline(ggplot2::aes_string(yintercept = 0), col = plot_arg2$ref, linetype = 2)

    ## do not include ref lines in ylim
    ylim[1] <- ifelse(
      is.na(ylim)[1],
      min(c(object$y[is.finite(object$y)], object$y_ci_lwr[is.finite(object$y_ci_lwr)])),
      ylim[1]
    )
    ylim[2] <- ifelse(
      is.na(ylim)[2],
      max(c(object$y[is.finite(object$y)], object$y_ci_up[is.finite(object$y_ci_upr)])),
      ylim[2]
    )
  }

  ## add points
  rval <- rval + 
    ggplot2::geom_point(ggplot2::aes_string(colour = "group", shape = "group", size = "group"))

  ## set the colors, shapes, etc.
  rval <- rval +
    ggplot2::scale_colour_manual(values = plot_arg$colour) +
    ggplot2::scale_fill_manual(values = plot_arg$fill) +
    ggplot2::scale_shape_manual(values = plot_arg$shape) +
    ggplot2::scale_size_manual(values = plot_arg$size) +
    ggplot2::scale_linetype_manual(values = plot_arg$linetype)

  ## annotation
  rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

  ## add legend
  if (legend) {
    rval <- rval + ggplot2::labs(colour = "Model") +
      ggplot2::guides(colour = "legend", shape = "none", size = "none")
  } else {
    rval <- rval + ggplot2::guides(colour = "none", shape = "none", size = "none")
  }

  ## set x and y limits (zoom only due to ref_fun)
  rval <- rval + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = 0.01)

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
