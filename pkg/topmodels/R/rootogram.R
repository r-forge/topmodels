## Rootogram
##
## - Oabserved y in-sample or out-of-sample (n x 1)
## - Breaks for observations (m x 1)
## - Predicted probabilities at breaks F_y(br1), ..., F_y(brm) (n x m)
##
## - Cut observations at breaks -> observed frequencies for (m-1) groups
## - Aggregate probabilities -> expected frequencies for (m-1) groups
## - Can be drawn in different style (standing vs. hanging / raw vs. sqrt)

## Functions:
## - rootogram() generic plus default method
## - Return object of class "rootogram" that is plotted by default
## - But has plot=FALSE so that suitable methods can be added afterwards
## - Methods: plot(), autoplot(), c()/rbind(), +

## TODO: The code below is simply copied from countreg.
## - No adaptations yet.
## - No 'newdata' argument yet.
## - Needs to be rewritten using procast().
## - Better defaults for breaks:
##   if(is.null(breaks)) breaks <- "Sturges"
##   if(length(breaks) == 1L) breaks <- hist(y, breaks = breaks)$breaks

rootogram <- function(object, ...) {
  UseMethod("rootogram")
}


rootogram.default <- function(object, 
                              newdata = NULL,
                              #na.action = na.pass, # FIXME: (ML) na.action must be na.omit, see newresponse()
                              plot = TRUE,
                              flavor = NULL,
                              style = c("hanging", "standing", "suspended"),
                              scale = c("sqrt", "raw"), 
                              breaks = NULL,
                              width = NULL, 
                              response_type = NULL,
                              xlab = NULL, 
                              ylab = NULL,
                              main = NULL, 
                              ...) {
  ## sanity checks
  ## `object` and `newdata` w/i `newrepsone()`
  ## `breaks` w/i `hist()`
  stopifnot(is.logical(plot))
  if (!is.null(flavor)) flavor <- try(match.arg(flavor, c("base", "tidyverse")), silent = TRUE)
  stopifnot(
    "`flavor` must either be NULL, or match the arguments 'base' or 'tidyverse'" = 
    is.null(flavor) || !inherits(flavor, "try-error")
  )
  stopifnot(is.null(breaks) || (is.numeric(breaks) && is.null(dim(breaks))))
  stopifnot(is.null(width) || (is.numeric(width) && length(width) == 1))
  stopifnot(length(xlab) == 1 | length(xlab) == 0)
  stopifnot(length(ylab) == 1 | length(ylab) == 0)
  stopifnot(length(main) == 1 | length(main) == 0)

  ## match arguments
  scale <- match.arg(scale)
  style <- match.arg(style)

  ## set flavor
  if (is.null(flavor) && "ggplot2" %in% (.packages()) && any(c("dplyr", "tibble") %in% (.packages()))) {
    flavor <- "tidyverse"
  } else if (is.null(flavor)) {
    flavor <- "base"
  } 

  ## default annotation
  if (is.null(xlab)) {
    xlab <- as.character(attr(terms(object), "variables"))[2L]  
  }
  if (is.null(ylab)) {
    ylab <- if (scale == "raw") "Frequency" else "sqrt(Frequency)"
  }
  if (is.null(main)) {
    main <- deparse(substitute(object))
  }

  ## data and weights
  y <- newresponse(object, newdata = newdata, na.action = na.omit)
  w <- attr(y, "weights")
  if(is.null(response_type)) response_type <- attr(y, "response_type")
  response_type <- match.arg(response_type, c("discrete", "logseries", "continuous"))

  ## set breaks and midpoints
  if (is.null(breaks) && response_type == "discrete") {
    breaks <- -1L:max(y[w > 0]) + 0.5
  } else if (is.null(breaks) && response_type == "logseries") {
    breaks <- 0L:max(y[w > 0]) + 0.5
  } else if (is.null(breaks)) {
    breaks <- "Sturges"
  }
  breaks <- hist(y[w > 0], plot = FALSE, breaks = breaks)$breaks
  x <- (head(breaks, -1L) + tail(breaks, -1L)) / 2

  
  ## set widths
  if (is.null(width) && (response_type == "discrete" || response_type == "logseries")) {
    width <- 0.9
  } else if (is.null(width)) {
    width <- 1
  }

  obsrvd <- as.vector(xtabs(w ~ cut(y, breaks, include.lowest = TRUE)))

  ## expected frequencies
  p <- matrix(NA, nrow = length(y), ncol = length(breaks) - 1L)
  for (i in 1L:ncol(p)) {
    p[, i] <- 
      procast(object, newdata = newdata, na.action = na.omit, type = "probability", 
        at = breaks[i + 1L], drop = TRUE) -
      procast(object, newdata = newdata, na.action = na.omit, type = "probability", 
        at = breaks[i], drop = TRUE) 
  }
  expctd <- colSums(p * w)
  ## FIXME: (ML) Do we need terms here? No info in `newresponse()` or `procast()` output? What did I mean?
  ## FIXME: (ML) 
  ## * This code will fail if NAs in response, newresponse() omits cases, procast() not.
  ## * Same error in countreg::procast.glm()
  ## * Should we allow other na.actions? What if NAs in covariates: newresponse will not match to forecast.

  ## raw vs. sqrt scale
  ## FIXME: (ML) Move scale in plot fun: Let it there (Z)
  if(scale == "sqrt") {
    y <- if(style == "hanging") sqrt(expctd) - sqrt(obsrvd) else 0
    height <- if(style == "suspended") sqrt(expctd) - sqrt(obsrvd) else sqrt(obsrvd)
  } else {
    y <- if(style == "hanging") expctd - obsrvd else 0
    height <- if(style == "suspended") expctd - obsrvd else obsrvd
  }

  ## collect everything as data.frame
  rval <- data.frame(
    observed = as.vector(obsrvd), 
    expected = as.vector(expctd),
    x = x, 
    y = y, 
    width = diff(breaks) * width, 
    height = height,
    line = if(scale == "sqrt") sqrt(expctd) else expctd
  )

  ## attributes for graphical display
  attr(rval, "style") <- style
  attr(rval, "scale") <- scale
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main

  if (flavor == "base") {
    class(rval) <- c("rootogram", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("rootogram", class(rval))
  }

  ## also plot by default
  if (plot & flavor == "tidyverse") {
    try(print(ggplot2::autoplot(rval, ...)))
  } else if (plot) {
    try(plot(rval, ...))
  }
  
  ## return invisibly
  invisible(rval)
 
}

c.rootogram <- rbind.rootogram <- function(...) {

  ## list of rootograms
  rval <- list(...)

  ## set flavor to tidyverse if any rval is a tibble
  if (any(do.call("c", lapply(rval, class)) %in% "tbl")) {
    flavor <- "tidyverse"
  } else {
    flavor <- "base"
  }

  ## convert always to data.frame
  if (flavor == "tidyverse") {
    rval <- lapply(rval, as.data.frame)
  }
  
  ## group sizes
  for(i in seq_along(rval)) {
    if(is.null(rval[[i]]$group)) rval[[i]]$group <- 1L
  }
  n <- lapply(rval, function(r) table(r$group))

  ## labels
  style <- unlist(lapply(rval, function(r) attr(r, "style")))
  scale <- unlist(lapply(rval, function(r) attr(r, "scale")))
  xlab <- unlist(lapply(rval, function(r) attr(r, "xlab")))
  ylab <- unlist(lapply(rval, function(r) attr(r, "ylab")))
  nam <- names(rval)
  main <- if(is.null(nam)) {
    as.vector(sapply(rval, function(r) attr(r, "main")))
  } else {
    make.unique(rep.int(nam, sapply(n, length)))
  }
  n <- unlist(n)

  ## combine
  rval <- do.call("rbind.data.frame", rval)
  rval$group <- if(length(n) < 2L) NULL else rep.int(seq_along(n), n)
  attr(rval, "style") <- style
  attr(rval, "scale") <- scale
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
 
  ## set class according to flavor
  if (flavor == "base") {
    class(rval) <- c("rootogram", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("rootogram", class(rval))
  }

  ## return
  return(rval)
}


plot.rootogram <- function(x,
                           ref = TRUE,
                           xlim = NULL, 
                           ylim = NULL, 
                           xlab = NULL, 
                           ylab = NULL, 
                           main = NULL,
                           border = "black", 
                           fill = adjustcolor("black", alpha.f = 0.2),
                           col = 2,
                           lwd = 2, 
                           pch = 19, 
                           lty = 1, 
                           type = NULL, 
                           axes = TRUE, 
                           box = FALSE,
                           ...) {

  ## sanity checks
  ## lengths of all arguments are checked by recycling
  ## `xlim`, `ylim`, `xlab`, `ylab`, `main` and `...` w/i `plot()`
  ## `border` and `fill` w/i `rect()`
  ##  `col`, `lwd`, `pch`, `lty` and `type` w/i `lines()`
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if(is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## arguments for plotting
  if(is.null(type)) type <- ifelse(any(table(x$group) > 20L), "l", "b")

  ## recycle arguments for plotting to match the number of groups
  if (is.null(xlim)) xlim <- c(NA, NA)
  if (is.null(ylim)) ylim <- c(NA, NA)
  if (is.list(xlim)) xlim <- as.data.frame(do.call("rbind", xlim))
  if (is.list(ylim)) ylim <- as.data.frame(do.call("rbind", ylim))
  plot_arg <- data.frame(1:n, ref,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]],
    border, fill, col, lwd, pch, lty, type, axes, box
  )[, -1]

  ## annotation
  if(is.null(xlab)) xlab <- TRUE
  if(is.null(ylab)) ylab <- TRUE
  if(is.null(main)) main <- TRUE
  xlab <- rep(xlab, length.out = n)
  ylab <- rep(ylab, length.out = n)
  main <- rep(main, length.out = n)
  if(is.logical(xlab)) xlab <- ifelse(xlab, attr(x, "xlab"), "")
  if(is.logical(ylab)) ylab <- ifelse(ylab, attr(x, "ylab"), "")
  if(is.logical(main)) main <- ifelse(main, attr(x, "main"), "")

  ## plotting function
  rootogram_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## rect elements
    xleft <- d$x - d$width/2
    xright <- d$x + d$width/2
    ybottom <- d$y
    ytop <- d$y + d$height
    
    ## get xlim and ylim
    if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) 
      plot_arg[j, c("xlim1", "xlim2")] <- range(c(xleft, xright))
    if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) 
      plot_arg[j, c("ylim1", "ylim2")] <- range(c(ybottom, ytop, d$line))

    ## trigger plot
    plot(0, 0, type = "n", xlim = c(plot_arg$xlim1[j], plot_arg$xlim2[j]),
      ylim = c(plot_arg$ylim1[j], plot_arg$ylim2[j]), 
      xlab = xlab[j], ylab = ylab[j], main = main[j], axes = FALSE, ...)
    if(plot_arg$axes[j]) {
      axis(1)
      axis(2)
    }
    if(plot_arg$box[j]) {
      box()
    }

    ## plot rootogram
    rect(xleft, ybottom, xright, ytop, border = plot_arg$border[j], col = plot_arg$fill[j])
  
    ## plot ref line 
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
      abline(h = 0, col = plot_arg$ref[j])
    }

    lines(d$x, d$line,
      col = plot_arg$col[j], pch = plot_arg$pch[j], type = plot_arg$type[j], 
      lty = plot_arg$lty[j], lwd =plot_arg$lwd[j]
    )
  }
   
  ## draw plots
  if(n > 1L) par(mfrow = n2mfrow(n))
  for(i in 1L:n) {
    rootogram_plot(x[x$group == i, ], ...)
  }
}


autoplot.rootogram <- function(object,
                               ref = TRUE,
                               #xlim = NULL, 
                               #ylim = NULL, 
                               #xlab = NULL, 
                               #ylab = NULL, 
                               #main = NULL,
                               border = "black",  
                               fill = "darkgray",
                               colour = 2, #col
                               size = 1.2, #lwd 
                               shape = 19, #pch
                               linetype = 1, #lty
                               type = NULL, 
                               #axes = TRUE, 
                               #box = FALSE,
                               ...) {


  ## determine grouping
  class(object) <- "data.frame"
  if(is.null(object$group)) object$group <- 1L
  n <- max(object$group)
  object$group <- factor(object$group, levels = 1L:n, 
    labels = make.names(attr(object, "main"), unique = TRUE))

  ## arguments for plotting
  if(is.null(type)) type <- ifelse(any(table(object$group) > 20L), "l", "b")

  ## recycle arguments for plotting to match the number of groups
  if(is.logical(ref)) ref <- ifelse(ref, 1, NA)  # get colors right 
  plot_arg <- data.frame(1:n, fill, colour, size, shape, ref, linetype)[, -1]

  ## recycle arguments for plotting to match the number of groups
  plot_arg2 <- data.frame(1:n, border, size, type, colour)[, -1]
  plot_arg2 <- as.data.frame(lapply(plot_arg2, rep, each = nrow(object) / n))
  plot_arg2$type[plot_arg2$type == "l"] <- NA  # get colors right
  plot_arg2$type[!is.na(plot_arg2$type)] <- plot_arg2$colour[!is.na(plot_arg2$type)]  # get colors right

  ## rectangles and fitted lines
  rval <- ggplot2::ggplot(object, ggplot2::aes_string(xmin = "x - width/2", xmax = "x + width/2", 
      ymin = "y", ymax = "y + height", x = "x", y = "line")) +
    ggplot2::geom_rect(ggplot2::aes_string(fill = "group"), colour = plot_arg2$border, show.legend = FALSE) + 
    ggplot2::geom_line(ggplot2::aes_string(colour = "group", size = "group", linetype = "group"), show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0, colour = plot_arg$ref) +
    ggplot2::geom_point(ggplot2::aes_string(shape = "group"), 
      colour = plot_arg2$type, size = plot_arg2$size * 2.2, show.legend = FALSE)

  ## grouping (if any)
  if(n > 1L) rval <- rval + ggplot2::facet_grid(group ~ .)

  rval <- rval + 
            ggplot2::scale_colour_manual(values = plot_arg$colour) + 
            ggplot2::scale_fill_manual(values = plot_arg$fill) + 
            ggplot2::scale_size_manual(values = plot_arg$size) + 
            ggplot2::scale_shape_manual(values = plot_arg$shape) + 
            ggplot2::scale_linetype_manual(values = plot_arg$linetype)
  
  ## annotation
  rval <- rval + ggplot2::xlab(paste(unique(attr(object, "xlab")), collapse = "/")) +
    ggplot2::ylab(paste(unique(attr(object, "ylab")), collapse = "/"))

  ## return with annotation
  rval
}


## FIXME: (ML) Should this be implemented?
#"+.rootogram" <- function(e1, e2) {
#  style <- unique(c(attr(e1, "style"), attr(e2, "style")))
#  if(length(style) > 1L) {
#    warning(sprintf("different styles (%s != %s) had been used, result now uses style = %s",
#      style[1L], style[2L], style[1L]))
#    style <- style[1L]
#  }
#  scale <- unique(c(attr(e1, "scale"), attr(e2, "scale")))
#  if(length(scale) > 1L) {
#    warning(sprintf("different scales (%s != %s) had been used, result now uses scale = %s",
#      scale[1L], scale[2L], scale[1L]))
#    scale <- scale[1L]
#  }
#
#  ylab <- attr(e1, "ylab")
#  xlab <- paste(unique(c(attr(e1, "xlab"), attr(e2, "xlab"))), collapse = " / ")
#  main <- paste(unique(c(attr(e1, "main"), attr(e2, "main"))), collapse = " / ")
#  e1 <- as.data.frame(e1)
#  e2 <- as.data.frame(e2)
#  e <- e1[e1$x %in% e2$x, ] + e2[e2$x %in% e1$x, ]
#  rootogram.default(structure(e$observed, .Names = e$x/2), e$expected,
#    style = style, scale = scale,
#    main = main, xlab = xlab, ylab = ylab, plot = FALSE)
#}
