## Wormplot
## shows the difference between the empirical quantile and the unit normal quantile 


wormplot <- function(object, ...) {
  UseMethod("wormplot")
}


wormplot.default <- function(object, 
                            newdata = NULL, 
                            plot = TRUE,
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

  ## sanity checks
  ## `object`, `newdata`, `delta and `prob` w/i `qresiduals()`
  ## `confint` w/i `polygon()`
  ## `delta` w/i `qresiduals()`
  ## `...` in `plot()`
  stopifnot(is.logical(plot))
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

  ## compute quantile residuals
  qres <- qresiduals(object, newdata = newdata, trafo = trafo, type = "random", nsim = nsim, delta = delta)
  if (is.null(dim(qres))) qres <- matrix(qres, ncol = 1L)

  ## compute corresponding quantiles on the transformed scale (default: normal)
  if (is.null(trafo)) trafo <- identity
  q2q <- function(y) trafo(ppoints(length(y)))[order(order(y))]
  qthe <- apply(qres, 2L, q2q)

  ## compute ci interval
  ## FIXME: (ML) Implement exact method if exists (see "inst/misc/2021_04_16_errorsearch_qqrplot.Rmd")
  if (!identical(confint, FALSE)) {
    set.seed(confint_seed)
    tmp <- qresiduals(object, newdata = newdata, trafo = trafo, type = "random", nsim = confint_nsim, 
      delta = delta)
    confint_prob <- (1 - confint_level) / 2
    confint_prob <- c(confint_prob, 1 - confint_prob)
    qres_ci_lwr <- apply(apply(tmp, 2, sort), 1, quantile, prob = confint_prob[1], na.rm = TRUE)
    qres_ci_upr <- apply(apply(tmp, 2, sort), 1, quantile, prob = confint_prob[2], na.rm = TRUE)
    qthe_ci_lwr <- q2q(qres_ci_lwr)
    qthe_ci_upr <- q2q(qres_ci_upr)

    ## FIXME: (ML) Dirty hack to get CI only for discrete values 
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

  ## setup function to calculate ref lines (ci interval according to Van Burren and Fredriks (2001) p. 1276)
  ## FIXME: (ML) Allow flexible trafos for `pnorm()`, `dnorm()` and `qnorm()`
  ciplot <- function(n, level = 0.95, lz = -2.75, hz = 2.75, dz = 0.25) {
    stopifnot(is.numeric(n), length(n) == 1)
    stopifnot(is.numeric(level), length(level) == 1, level >= 0, level <= 1)
    stopifnot(is.numeric(lz), length(lz) == 1)
    stopifnot(is.numeric(hz), length(hz) == 1)
    stopifnot(is.numeric(dz), length(dz) == 1)

    z <- seq(lz, hz, dz)
    p <- pnorm(z)
    se <- (1 / dnorm(z)) * (sqrt(p * (1 - p) / n))
    low <- qnorm((1 - level) / 2 ) * se
    high <- -low
    return(list(z = z, low = as.numeric(low), high = as.numeric(high)))
  }

  ## labels
  if (is.null(main)) main <- deparse(substitute(object))

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
  attr(rval, "ref_fun") <- list(ciplot)
  attr(rval, "confint_level") <- ifelse(confint, confint_level, NA)
  class(rval) <- c("wormplot", "data.frame")

  ## also plot by default
  if(plot){ 
    try(plot(rval, confint = confint, ...))
  }
  
  ## return coordinates invisibly
  invisible(rval)
}

## Combine several wormplots
c.wormplot <- rbind.wormplot <- function(...) {

  ## list of wormplots
  rval <- list(...)
  
  ## remove temporary the class
  for(i in 1:length(rval)) class(rval[[i]]) <- "data.frame"  # TODO: (ML) How does that work with lapply?

  ## group sizes
  for (i in seq_along(rval)) {
    if (is.null(rval[[i]]$group)) rval[[i]]$group <- 1L
  }
  n <- lapply(rval, function(r) table(r$group))

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

  ## combine and return
  all_names <- unique(unlist(lapply(rval, names)))
  if (any(grepl("x_1", all_names)) & any(grepl("^x$", all_names))) {
    for (i in 1:length(rval)) {
      names(rval[[i]])[grepl("^x$", names(rval[[i]]))] <- "x_1";
      names(rval[[i]])[grepl("^y$", names(rval[[i]]))] <- "y_1"
    }
  all_names <- unique(unlist(lapply(rval, names)))
  }

  rval <- do.call("rbind.data.frame", 
            c(lapply(
              rval,
              function(x) data.frame(c(x, sapply(setdiff(all_names, names(x)), function(y) NA)))),
              make.row.names = FALSE)
          )

  rval$group <- if (length(n) < 2L) NULL else rep.int(seq_along(n), n)
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "confint_level") <- confint_level
  attr(rval, "ref_fun") <- ref_fun
  class(rval) <- c("wormplot", "data.frame")
  return(rval)
}


## actual drawing
plot.wormplot <- function(x,
                         single_graph = FALSE,
                         confint = TRUE, 
                         ref = TRUE,
                         xlim = NULL, 
                         ylim = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         main = NULL,
                         col = "black",
                         fill = adjustcolor("black", alpha.f = 0.2), 
                         alpha_min = 0.2, # single or n values  
                         pch = 19,
                         axes = TRUE,
                         box = TRUE,
                         ...) {

  ## sanity checks
  ## lengths of all arguments are checked by recycling; `ref` w/i `abline()`
  ## `xlim`, `ylim`, `xlab`, `ylab`, `main` and `....` w/i `plot()`
  ## `col`, `pch` w/i `lines()`
  ## `confint`, `fill` in `polygon()`
  ## `alpha_min` w/i colorspace fun 
  stopifnot(is.logical(single_graph))
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  if (is.null(xlim)) xlim <- c(NA, NA)
  if (is.null(ylim)) ylim <- c(NA, NA)
  plot_arg <- data.frame(1:n, confint, ref,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]],
    col, fill, alpha_min, pch, axes, box 
  )[, -1]

  ## annotation
  if (single_graph) {
    if (is.null(xlab)) xlab <- "Theoretical quantiles"
    if (is.null(ylab)) ylab <- "Deviation"
    if (is.null(main)) main <- "Worm plot" # FIXME: (ML) Achim prefers other title
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

  ## plotting function
  wormplot_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get xlim and ylim conditional on confint
    if(
      !identical(plot_arg$confint[j], FALSE) && 
      !is.na(attr(d, "confint_level")[j])
    ) {
      if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) 
        xlim <- range(as.matrix(d[grepl("x", names(d))]), finite = TRUE)
      if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) 
        ylim <- range(as.matrix(d[grepl("y", names(d))]), finite = TRUE)
    } else {
      if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) 
        xlim <- range(as.matrix(d[grepl("^x$|x_[0-9]", names(d))]), finite = TRUE)
      if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) 
        ylim <- range(as.matrix(d[grepl("^y$|y_[0-9]", names(d))]), finite = TRUE)
    }

    ## calculate ref lines
    ## FIXME: (ML) All data should be prepared in `wormplot.default()`, but length does not fit to rval.
    ## FIXME: (ML) Adapt to other trafos (density and quantile functions).
    #level <- 0.95
    #ref_x <- seq(xlim[1] - 10, xlim[2] + 10, 0.25)
    #p <- pnorm(ref_x)
    #se <- (1 / dnorm(ref_x)) * (sqrt(p * (1 - p) / NROW(d)))
    #ref_low <- qnorm((1 - level) / 2) * se
    #ref_high <- -ref_low

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
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min =  plot_arg$alpha_min[j]),
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

    ## plot reference diagonal
    if (j == 1 || (!single_graph && j > 1)) {
      if (!identical(plot_arg$ref[j], FALSE)) {
        if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
        ref_lines <- attr(d, "ref_fun")[[j]](
          n = length(d$y), 
          level = 0.95, 
          lz = xlim[1] - 10, 
          hz = xlim[2] + 10, 
          dz = 0.25
        )
        abline(h = 0, col = plot_arg$ref[j], lty = 2)
        lines(ref_lines$z, ref_lines$low, col = plot_arg$ref[j], lty = 2)
        lines(ref_lines$z, ref_lines$high, col = plot_arg$ref[j], lty = 2)
      }
    }

  }

  ## draw plots
  if (n > 1L) par(mfrow = n2mfrow(n))
  for (i in 1L:n) wormplot_plot(x[x$group == i, ], ...)
}


points.wormplot <- function(x,
                           confint = FALSE, 
                           ref = FALSE,
                           col = "black",
                           fill = adjustcolor("black", alpha.f = 0.2), 
                           alpha_min = 0.2,
                           pch = 19,
                           ...) {

  ## sanity checks
  ## lengths of all arguments are checked by recycling; `ref` w/i `abline()`
  ## `col`, `pch` w/i `lines()`
  ## `confint`, `fill` in `polygon()`
  ## `alpha_min` w/i colorspace fun 

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(1:n, confint, ref, col, fill, alpha_min, pch )[, -1]

  ## plotting function
  wormplot_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get xlim and ylim conditional on confint
    if(
      !identical(plot_arg$confint[j], FALSE) && 
      !is.na(attr(d, "confint_level")[j])
    ) {
      if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) 
        xlim <- range(as.matrix(d[grepl("x", names(d))]), finite = TRUE)
      if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) 
        ylim <- range(as.matrix(d[grepl("y", names(d))]), finite = TRUE)
    } else {
      if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) 
        xlim <- range(as.matrix(d[grepl("^x$|x_[0-9]", names(d))]), finite = TRUE)
      if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) 
        ylim <- range(as.matrix(d[grepl("^y$|y_[0-9]", names(d))]), finite = TRUE)
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
      points(
        d[grepl("x", names(d))][, i], 
        d[grepl("y", names(d))][, i], 
        col = plot_arg$col[j], pch = plot_arg$pch[j], ...
      )
    }

    ## plot reference diagonal
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
      abline(h = 0, col = plot_arg$ref[j], lty = 2)
    }

  }

  ## draw plots
  for (i in 1L:n) {
    wormplot_plot(x[x$group == i, ], ...)
  }
}

## ggplot2 interface
autoplot.wormplot <- function(object, ...) {
  NULL
}

