## Programming outline: (Randomized) Q-Q residuals plot
##
## - Observed y in-sample or out-of-sample (n x 1)
## - Predicted probabilities F_y(y - eps) and F_y(y) (n x 2)
## - Two columns can be essentially equal -> continuous
##   or different -> (partially) discrete
## - Potentially transform uniform scale to different
##   distribution (default: Gaussian, via qnorm()).
##
## - Plot ordered empirical quantile residuals against
##   theoretical quantiles (from same distribution)
## - To deal with point masses, draw either multiple random
##   draws (enable alpha blending by default) or shade quantiles

## Functions:
## - qqrplot() generic plus default method
## - Return object of class "qqrplot" that is plotted by default
## - But has plot=FALSE so that suitable methods can be added afterwards
## - At least methods: plot(), autoplot()

qqrplot <- function(object, ...) {
  UseMethod("qqrplot")
}


qqrplot.default <- function(object,
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
                            ylab = "Quantile residuals",
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
  stopifnot(is.numeric(confint_nsim), length(confint_nsim) == 1)
  stopifnot(is.numeric(confint_seed), length(confint_seed) == 1)
  stopifnot(is.logical(single_graph))
  stopifnot(length(xlab) == 1)
  stopifnot(length(ylab) == 1)
  stopifnot(length(main) == 1 || length(main) == 0)

  ## guess plotting flavor
  if (isFALSE(plot)) {
    plot <- "none"
  } else if (isTRUE(plot)) {
    plot <- if ("package:ggplot2" %in% search()) "ggplot2" else "base"
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
    class <- if ("package:tibble" %in% search()) "tibble" else "data.frame"
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
      y = qres
    )
  } else {
    rval <- data.frame(
      x = qthe,
      y = qres,
      y_ci_lwr = qres_ci_lwr,
      y_ci_upr = qres_ci_upr,
      x_ci_lwr = qthe_ci_lwr,
      x_ci_upr = qthe_ci_upr
    )
  }
  names(rval) <- gsub("(\\.r|\\.q)", "", names(rval))

  ## attributes for graphical display
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "confint_level") <- ifelse(confint, confint_level, NA)

  ## add class
  if (class == "data.frame") {
    class(rval) <- c("qqrplot", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("qqrplot", class(rval))
  }

  ## plot by default
  if (plot == "ggplot2") {
    try(print(ggplot2::autoplot(rval, confint = confint, ...)))
  } else if (plot == "base") {
    try(plot(rval, confint = confint, ...))
  }

  ## return invisibly
  invisible(rval)
}


c.qqrplot <- rbind.qqrplot <- function(...) {
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


plot.qqrplot <- function(x,
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
    if (is.null(ylab)) ylab <- "Quantile residuals"
    if (is.null(main)) main <- "Q-Q residuals plot" # FIXME: (ML) Achim prefers other title
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

    ## plot reference diagonal
    if (j == 1 || (!single_graph && j > 1)) {
      if (!identical(plot_arg$ref[j], FALSE)) {
        if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
        abline(0, 1, col = plot_arg$ref[j], lty = 2, lwd = 1.25)
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
  if (n > 1L) par(mfrow = n2mfrow(n))

  ## draw qqrplots
  for (i in 1L:n) qqrplot_plot(x[x$group == i, ], ...)
}


points.qqrplot <- function(x,
                           confint = FALSE,
                           ref = FALSE,
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
  qqrplot_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

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
      abline(0, 1, col = plot_arg$ref[j], lty = 2, lwd = 1.25)
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


autoplot.qqrplot <- function(object,
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

  ## set color to NA for not plotting
  if (is.logical(ref)) ref <- ifelse(ref, 1, NA)

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
    fill, colour, size, ref, linetype, confint, alpha_min, shape
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
  rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y")) +
    ggplot2::geom_abline(ggplot2::aes_string(intercept = 0, slope = 1),
      linetype = 2, colour = plot_arg2$ref
  )

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

  ## add points
  rval <- rval +
    ggplot2::geom_point(ggplot2::aes_string(colour = "group", shape = "group", size = "group"),
      show.legend = FALSE
    )

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
      ggplot2::guides(colour = "legend", size = "none", linetype = "none")
  } else {
    rval <- rval + ggplot2::guides(colour = "none", size = "none", linetype = "none")
  }

  ## set x and y limits
  rval <- rval + ggplot2::scale_x_continuous(limits = xlim, expand = c(0.01, 0.01))
  rval <- rval + ggplot2::scale_y_continuous(limits = ylim, expand = c(0.01, 0.01))

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
