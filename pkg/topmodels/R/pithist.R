## PIT histogram
##
## - Observed y in-sample or out-of-sample (n x 1)
## - Predicted probabilities F_y(y - eps) and F_y(y) (n x 2)
## - Two columns can be essentially equal -> continuous
##   or different -> (partially) discrete
## - Breaks for predicted probabilities in [0, 1] (m x 1)
##
## - Cut probabilities at breaks -> (m-1) groups and draw histogram
## - TODO: In case of point masses either use a random draw
##   or distribute evenly across relevant intervals 
## - TODO: Random draws could be drawn by hist() (current solution)
##   but proportional distribution requires drawing rectangles by hand
## - TODO: add confidence interval as well.
## - TODO: Instead of shaded rectangles plus reference line and CI lines
##   support shaded CI plus step lines

## Functions:
## - pithist() generic plus default method
## - Return object of class "pithist" that is plotted by default
## - But has plot=FALSE so that suitable methods can be added afterwards
## - At least methods: plot(), autoplot(), lines(), possibly c(), +

pithist <- function(object, ...) {
  UseMethod("pithist")
}


pithist.default <- function(object,
                            newdata = NULL,
                            plot = TRUE,
                            style = c("histogram", "lines"),
                            type = c("random", "proportional"),  # FIXME: (ML) see below
                            nsim = 1L,
                            delta = NULL,
                            freq = FALSE,
                            breaks = NULL,
                            confint = TRUE,
                            confint_level = 0.95,
                            confint_type = c("exact", "approximation"),
                            single_graph = FALSE,
                            xlim = c(0, 1),
                            ylim = NULL,
                            xlab = "PIT",
                            ylab = if (freq) "Frequency" else "Density",
                            main = NULL,
                            ...) {

  ## sanity checks
  ## `object` and `newdata` w/i `newrepsone()`; `delta w/i `qresiduals()`, `confint` in `abline()`
  stopifnot(is.logical(plot))
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
  stopifnot(is.numeric(xlim), length(xlim) == 2)
  stopifnot(is.null(ylim) || (is.numeric(ylim) && length(ylim) == 2))
  stopifnot(length(xlab) == 1)
  stopifnot(length(ylab) == 1)
  stopifnot(length(main) == 1 | length(main) == 0)

  ## match arguments
  style <- match.arg(style)
  confint_type <- match.arg(confint_type)

  ## either compute proportion exactly (to do...) or approximate by simulation
  type <- match.arg(type)
  if (type == "proportional") {  
    # FIXME: (ML) implement proportional over the inteverals (e.g., below censoring point)
    # FIXME: (ML) confusing naming, as `type` in `qresiduals()` must be `random` or `quantile`
    stop("not yet implemented")
  } else {
    p <- qresiduals.default(object, newdata = newdata, trafo = NULL, type = "random", 
      nsim = nsim, delta = delta)
  }

  ## get breaks
  if (is.null(breaks)) breaks <- c(4, 10, 20, 25)[cut(NROW(p), c(0, 50, 5000, 1000000, Inf))]
  if (length(breaks) == 1L) breaks <- seq(xlim[1L], xlim[2L], length.out = breaks + 1L)

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

  ## collect everything as data.frame
  tmp_hist <- hist(p, breaks = breaks, plot = FALSE) # TODO: (ML) Maybe get rid of `hist()`
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
  class(rval) <- c("pithist", "data.frame")

  ## plot by default
  if (plot) {
    try(
      plot(rval,
        style = style, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
        confint = confint, ...
      )
    )
  }

  ## return invisibly
  invisible(rval)
}


c.pithist <- rbind.pithist <- function(...) {

  ## list of pithists
  rval <- list(...)

  ## group sizes
  for (i in seq_along(rval)) {
    if (is.null(rval[[i]]$group)) rval[[i]]$group <- 1L
  }
  n <- lapply(rval, function(r) table(r$group))

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

  ## combine and return
  rval <- do.call("rbind.data.frame", rval)
  rval$group <- if (length(n) < 2L) NULL else rep.int(seq_along(n), n)
  attr(rval, "freq") <- freq
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "confint_level") <- confint_level
  class(rval) <- c("pithist", "data.frame")
  return(rval)
}


plot.pithist <- function(x,
                         single_graph = FALSE,
                         style = c("histogram", "lines"),
                         confint = TRUE,
                         ref = TRUE,
                         xlim = c(0, 1),
                         ylim = NULL,
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

  ## sanity checks
  ## lengths of all arguments are checked by recycling; `ref` and `confint` w/i `abline()`;
  ## `xlim`, `ylim`, xlab`, `ylab`, `main`, `col`, `fill`, `lwd`, `lty` and `...` w/i `plot()`
  ## `alpha_min` w/i colorspace fun 
  stopifnot(is.logical(single_graph))
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))
  if (single_graph) stopifnot(
    "for `single_graph` all `freq` in attr of object `x` must be of the same type" = 
    length(unique(attr(x, "freq"))) == 1
  )

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## set style
  style <- match.arg(style)
  if (n > 1 && single_graph && style == "histogram"){
    style <- "lines"
  }

  ## set lwd
  if (is.null(lwd)) lwd <- if (style == "histogram") 1 else 2

  ## recycle arguments for plotting to match the number of groups
  if (is.null(ylim)) ylim <- c(NA, NA)
  plot_arg <- data.frame(1:n, confint, ref,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]], 
    border, col, fill, alpha_min, lwd, lty, axes, box
  )[, -1]

  ## annotation
  if (single_graph) {
    if (is.null(xlab)) xlab <- "PIT"
    if (is.null(ylab)) ylab <- if (all(attr(x, "freq"))) "Frequency" else "Density"
    if (is.null(main)) main <- "PIT Histogram"
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

  ## function to plot pithist (style = "histogram")
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
    if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) xlim <- range(c(xleft, xright))
    if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) ylim <- range(c(0, y, ci_lwr, ci_upr))

    ## trigger plot
    if (j == 1 || (!single_graph && j > 1)) {
      plot(0, 0,
        type = "n", xlim = xlim, ylim = ylim,
        xlab = xlab[j], ylab = ylab[j], main = main[j], axes = FALSE, ...
      )
      if(plot_arg$axes[j]) {
        axis(1)
        axis(2)
      }
      if(plot_arg$box[j]) {
        box()
      }
    }

    if (plot_arg$fill[j] == adjustcolor("black", alpha.f = 0.2) && 
      plot_arg$col[j] != "black") {
      message(" * As the argument `col` is set but no argument `fill` is specified, \n   the former is used for colorizing the PIT histogram. \n * For proper usage, solely provide `fill` for histogram style plots.")
      plot_arg$fill[j] <- plot_arg$col[j]
    }

    ## plot pithist
    rect(xleft, 0, xright, y, border = plot_arg$border[j], col = plot_arg$fill[j], 
      lty = plot_arg$lty[j], lwd = plot_arg$lwd[j])

    ## plot ref line 
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
      abline(h = pp, col = plot_arg$ref[j], lty = 1, lwd = 1.5)
    }

    ## plot confint lines
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- "black"
      abline(h = ci_lwr, col = plot_arg$confint[j], lty = 2, lwd = 1.5)
      abline(h = ci_upr, col = plot_arg$confint[j], lty = 2, lwd = 1.5)
    }

  }

  ## function to trigger figure and plot confint (style = "lines")
  pitlines_trigger <- function(d, ...) {

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

    ## get xlim and ylim
    if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) xlim <- range(x)
    if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) ylim <- range(c(0, y, ci_lwr, ci_upr))

    ## trigger plot
    if (j == 1 || (!single_graph && j > 1)) {
      plot(0, 0,
        type = "n", xlim = xlim, ylim = ylim,
        xlab = xlab[j], ylab = ylab[j], xaxs = "i", main = main[j], axes = FALSE, ...
      )
      if(plot_arg$axes[j]) {
        axis(1)
        axis(2)
      }
      if(plot_arg$box[j]) {
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

  ## function to trigger figure and plot confint (style = "lines")
  pitlines_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## step elements
    x <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
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
    lines(y ~ x, type = "s", lwd = plot_arg$lwd[j], lty = plot_arg$lty[j], col = plot_arg$col[j])
  }

  ## draw plots
  if (!single_graph && n > 1L) par(mfrow = n2mfrow(n))

  ## draw polygons first
  if (single_graph || n == 1) {
    if (style == "histogram") {
      for (i in 1L:n) pithist_plot(x[x$group == i, ], ...)
    } else {
      for (i in 1L:n) pitlines_trigger(x[x$group == i, ], ...)
      for (i in 1L:n) pitlines_plot(x[x$group == i, ], ...)
    }
  } else {
    for (i in 1L:n) {
      pitlines_trigger(x[x$group == i, ], ...)
      pitlines_plot(x[x$group == i, ], ...)
    }
  }
}


lines.pithist <- function(x,
                          confint = FALSE,
                          ref = FALSE,
                          col = "black",
                          fill = adjustcolor("black", alpha.f = 0.2),
                          alpha_min = 0.2, 
                          lwd = 2,
                          lty = 1,
                          ...) {

  ## sanity checks
  ## lengths of all arguments are checked by recycling; `ref` and `confint` w/i `abline()`;
  ## `col`, `fill`, `lwd`, `lty` and `...` w/i `lines()`
  ## `alpha_min` w/i colorspace fun 
  stopifnot(
    "all `freq` in attr of object `x` must be of the same type" =
    length(unique(attr(x, "freq"))) == 1
  )

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(
    1:n, confint, ref, col, fill, alpha_min, lwd, lty
  )[, -1]

  ## plotting function: stepwise pithist
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
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
      segments(x0 = 0, y0 = pp, x1 = 1, y1 = pp, col = plot_arg$ref[j], lty = 2, lwd = 1.5)
    }

    ## plot stepfun
    lines.default(y ~ x, type = "s", lwd = plot_arg$lwd[j], lty = plot_arg$lty[j], col = plot_arg$col[j])
  }

  ## draw plots
  for (i in 1L:n) {
    pitlines_lines(x[x$group == i, ], ...)
  }
}


## ggplot2 interface
autoplot.pithist <- function(object,
                             style = c("histogram", "lines"),
                             grid = TRUE,
                             confint = TRUE,
                             colour = "black",
                             fill = "darkgray",
                             border = "black",
                             linetype = 1,
                             size = 1.0,
                             ylim = c(0, NA),
                             xlim = c(0, 1),
                             ...) {

  ## match arguments
  style <- match.arg(style)

  ## determine grouping
  class(object) <- "data.frame"
  if (is.null(object$group)) object$group <- 1L
  n <- max(object$group)
  object$group <- factor(object$group,
    levels = 1L:n,
    labels = make.names(attr(object, "main"), unique = TRUE)
  )

  ## get CI and perfect prediction
  ci_lwr <- if (confint) unique(object$ci_lwr) else NULL
  ci_upr <- if (confint) unique(object$ci_upr) else NULL
  pp <- unique(object$ref)
  stopifnot(
    "`pp` in attr of object `x` must consist of unique values per group index" =
    length(pp) == 1
  )

  if (style == "histogram") {

    if (length(colour) == 1L) {
      colour <- rep(colour[1], 3L)
    } else if (length(colour) == 2L) {
      colour <- c(colour[1], rep(colour[2], 2L))
    } else {
      colour <- colour[1:3]
    }

    if (length(linetype) == 1L) {
      linetype <- rep(linetype[1], 3L)
    } else if (length(linetype) == 2L) {
      linetype <- c(linetype[1], rep(linetype[2], 2L))
    } else {
      linetype <- linetype[1:3]
    }

    rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y / 2", width = "width", height = "y")) +
      ggplot2::geom_tile(colour = border, fill = fill) +
      ggplot2::geom_hline(yintercept = pp, colour = colour[1], linetype = linetype[1], size = size) +
      ggplot2::geom_hline(yintercept = ci_lwr, colour = colour[2], linetype = linetype[2], size = size) +
      ggplot2::geom_hline(yintercept = ci_upr, colour = colour[3], linetype = linetype[3], size = size)

    ## grouping (if any)
    if (n > 1L) rval <- rval + ggplot2::facet_grid(group ~ .) + ggplot2::labs(colour = "Model")

  } else {

    if (length(colour) < n) {
      group_colours <- rep(colour[1], n)
    } else {
      group_colours <- colour
    }
    names(group_colours) <- levels(object$group)

    if (length(linetype) < n) {
      group_linetypes <- rep(linetype[1], n)
    } else {
      group_linetypes <- linetype
    }
    names(group_linetypes) <- levels(object$group)

    ## stat helper function to get left/right points from respective mid points
    calc_pit_points <- ggplot2::ggproto("calc_pit_points", ggplot2::Stat,

      # Required as we operate on groups (facetting)
      compute_group = function(data, scales) {
        ## Manipulate object  #TODO: Proably not "very nice" code
        nd <- data.frame(
          x = c(data$x - data$width / 2, data$x[NROW(data)] + data$width[NROW(data)] / 2),
          y = c(data$y, data$y[NROW(data)])
        )
        nd
      },

      # Tells us what we need
      required_aes = c("x", "y")
    )

    if (grid == TRUE & length(unique(colour)) == 1L & length(unique(linetype)) == 1L) {
      rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y", width = "width")) +
        ggplot2::geom_rect(
          xmin = 0, xmax = 1, ymin = ci_lwr, ymax = ci_upr,
          fill = fill, alpha = 0.5, colour = NA
        ) +
        ggplot2::geom_step(stat = calc_pit_points, linetype = linetype, size = size, colour = colour)

      if (!all(is.na(ylim))) {
        rval <- rval + ggplot2::ylim(ylim)
      }

      if (!all(is.na(xlim))) {
        rval <- rval + ggplot2::xlim(xlim)
      }

      ## grouping (if any)
      if (n > 1L) {
        rval <- rval + ggplot2::facet_grid(group ~ .) + ggplot2::labs(colour = "Model") + ggplot2::labs(linetype = "Model")
      }
  
    } else if (grid == TRUE & (length(unique(colour)) > 1L | length(unique(linetype)) > 1L)) {
      rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y", width = "width", linetype = "group", colour = "group")) +
        ggplot2::geom_rect(
          xmin = 0, xmax = 1, ymin = ci_lwr, ymax = ci_upr,
          fill = fill, alpha = 0.5, colour = NA
        ) +
        ggplot2::geom_step(stat = calc_pit_points, size = size) +
        ggplot2::scale_colour_manual(values = group_colours) +
        ggplot2::scale_linetype_manual(values = group_linetypes)

      if (!all(is.na(ylim))) {
        rval <- rval + ggplot2::ylim(ylim)
      }

      if (!all(is.na(xlim))) {
        rval <- rval + ggplot2::xlim(xlim)
      }

      ## grouping (if any)
      if (n > 1L) {
        rval <- rval + ggplot2::facet_grid(group ~ .) + ggplot2::labs(colour = "Model") + ggplot2::labs(linetype = "Model")
      }
  
    } else if (grid == FALSE & length(unique(colour)) == 1L & length(unique(linetype)) == 1L) {
      rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y", width = "width", colour = "group")) +
        ggplot2::geom_rect(
          xmin = 0, xmax = 1, ymin = ci_lwr, ymax = ci_upr,
          fill = fill, alpha = 0.5, colour = NA
        ) +
        ggplot2::geom_step(stat = calc_pit_points, linetype = linetype[1], size = size) + 
        ggplot2::scale_colour_manual(values = group_colours) + 
        ggplot2::scale_linetype_manual(values = group_linetypes) + 
        ggplot2::labs(colour = "Model") + 
        ggplot2::theme(legend.position = "none")

      if (!all(is.na(ylim))) {
        rval <- rval + ggplot2::ylim(ylim)
      }

      if (!all(is.na(xlim))) {
        rval <- rval + ggplot2::xlim(xlim)
      }
    } else {  # grid == FALSE  & (length(unique(colour)) > 1L | length(unique(linetype)) > 1L)
      rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y", width = "width", linetype = "group", colour = "group")) +
        ggplot2::geom_rect(
          xmin = 0, xmax = 1, ymin = ci_lwr, ymax = ci_upr,
          fill = fill, alpha = 0.5, colour = NA
        ) +
        ggplot2::geom_step(stat = calc_pit_points, size = size) + 
        ggplot2::scale_colour_manual(values = group_colours) + 
        ggplot2::scale_linetype_manual(values = group_linetypes) + 
        ggplot2::labs(colour = "Model") + 
        ggplot2::labs(linetype = "Model")

      if (!all(is.na(ylim))) {
        rval <- rval + ggplot2::ylim(ylim)
      }

      if (!all(is.na(xlim))) {
        rval <- rval + ggplot2::xlim(xlim)
      }
    }
  }

  ## annotation
  rval <- rval + ggplot2::xlab(paste(unique(attr(object, "xlab")), collapse = "/")) +
    ggplot2::ylab(paste(unique(attr(object, "ylab")), collapse = "/"))

  ## return with annotation
  rval
}


## helper function to calculate CI employing `qbinom()`
get_confint <- function(n, bins, level, freq) {
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  rval <- qbinom(a, size = n, prob = 1 / bins)
  if (!freq) rval <- rval / (n / bins)
  rval
}


## helper function to calculate an approximated CI according to Agresti & Coull (1998)
## doi=10.1080/00031305.1998.10480550
get_confint_agresti <- function(x, n, level, bins, freq) {
  rval <- add4ci(x, n, level)$conf.int * n
  if (!freq) rval <- rval / (n / bins)
  rval
}


## copy of `add4ci` package from package `PropCIs` by Ralph Scherer (licensed under GPL-2/GPL-3)
add4ci <- function(x, n, conf.level) {
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
