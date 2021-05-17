## Wormplot
## shows the difference between the empirical quantile and the unit normal quantile 


wormplot <- function(object, ...) {
  UseMethod("wormplot")
}


wormplot.default <- function(object, 
                             newdata = NULL, 
                             plot = TRUE,
                             flavor = NULL,
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
  if (!is.null(flavor)) flavor <- try(match.arg(flavor, c("base", "tidyverse")), silent = TRUE)
  stopifnot(
    "`flavor` must either be NULL, or match the arguments 'base' or 'tidyverse'" =
    is.null(flavor) || !inherits(flavor, "try-error")
  )
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

  ## guess flavor
  if (is.null(flavor) && "ggplot2" %in% (.packages()) && any(c("dplyr", "tibble") %in% (.packages()))) {
    flavor <- "tidyverse"
  } else if (is.null(flavor)) {
    flavor <- "base"
  }

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

  ## setup function to calculate ref lines (ci interval according to Van Buuren and Fredriks (2001) p. 1276)
  ## FIXME: (ML) All data should be prepared in `wormplot.default()`, but length does not fit to rval.
  ## FIXME: (ML) Adapt to other trafos (density and quantile functions).
  ref_fun <- function(x, n, level = 0.95, which = c("lower", "upper")) {
    stopifnot(is.numeric(n), length(n) == 1)
    stopifnot(is.numeric(level), length(level) == 1, level >= 0, level <= 1)
    which <- match.arg(which)

    p <- pnorm(x)
    se <- (1 / dnorm(x)) * (sqrt(p * (1 - p) / n))
    rval <- as.numeric(qnorm((1 - level) / 2 ) * se)

    if (which == "lower") {
      rval
    } else {
      -rval
    }
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
  attr(rval, "ref_fun") <- list(ref_fun)
  attr(rval, "confint_level") <- ifelse(confint, confint_level, NA)

  if (flavor == "base") {
    class(rval) <- c("wormplot", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("wormplot", class(rval))
  }

  ## plot by default
  if (plot & flavor == "tidyverse") {
    try(print(ggplot2::autoplot(rval, confint = confint, ...)))
  } else if (plot) {
    try(plot(rval, confint = confint, ...))
  }
  
  ## return coordinates invisibly
  invisible(rval)
}

## Combine several wormplots
c.wormplot <- rbind.wormplot <- function(...) {

  ## list of wormplots
  rval <- list(...)

  ## set flavor to tidyverse if any rval is a tibble
  if (any(do.call("c", lapply(rval, class)) %in% "tbl")) {
    flavor <- "tidyverse"
  } else {
    flavor <- "base"
  }
  
  ## remove temporary the class
  ## TODO: (ML) How does that work with lapply?
  ## FIXME: (ML) Maybe this can be nicer?
  for(i in 1:length(rval)) class(rval[[i]]) <- class(rval[[i]])[!class(rval[[i]]) %in% "wormplot"]

  ## convert always to data.frame
  if (flavor == "tidyverse") {
    rval <- lapply(rval, as.data.frame)
  }

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

  ## set class according to flavor
  if (flavor == "base") {
    class(rval) <- c("wormplot", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("wormplot", class(rval))
  }

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
                         col = adjustcolor("black", alpha.f = 0.4), 
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

    ## plot ref lines
    if (j == 1 || (!single_graph && j > 1)) {
      if (!identical(plot_arg$ref[j], FALSE)) {
        if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
        curve(
          attr(d, "ref_fun")[[j]](
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
          attr(d, "ref_fun")[[j]](
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

  ## draw plots
  if (n > 1L) par(mfrow = n2mfrow(n))
  for (i in 1L:n) wormplot_plot(x[x$group == i, ], ...)
}


points.wormplot <- function(x,
                           confint = FALSE, 
                           ref = FALSE,
                           col = "black",
                           fill = adjustcolor("black", alpha.f = 0.4), 
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
      abline(h = 0, col = plot_arg$ref[j], lty = 2, lwd = 1.25)
    }

  }

  ## draw plots
  for (i in 1L:n) {
    wormplot_plot(x[x$group == i, ], ...)
  }
}

## ggplot2 interface
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

  ## sanity checks
  stopifnot(is.logical(single_graph))

  ## determine grouping
  class(object) <- "data.frame"
  if (is.null(object$group)) object$group <- 1L
  n <- max(object$group)

  ## get annotations in the right lengths
  if(is.null(xlab)) xlab <- attr(object, "xlab")
  xlab <- paste(unique(xlab), collapse = "/")
  if(is.null(ylab)) ylab <- attr(object, "ylab")
  ylab <- paste(unique(ylab), collapse = "/")
  if(is.null(main)) main <- attr(object, "main")
  main <- make.names(rep_len(main, n), unique = TRUE)

  ## prepare grouping
  object$group <- factor(object$group, levels = 1L:n, labels = main)

  ## FIXME: (ML) This must be done in base and somehow nicer!
  ref_fun <- attr(object, "ref_fun")[[1]]
  object <- tidyr::pivot_longer(object, cols = names(object)[grepl("^x$|x_[0-9]", names(object))],
                            names_to = "x_sim", values_to = "x")
  object <- tidyr::pivot_longer(object, cols = names(object)[grepl("^y$|y_[0-9]", names(object))],
                            names_to = "y_sim", values_to = "y")
  object <- object[which(gsub("x", "", object$x_sim) == gsub("y", "", object$y_sim)), ]
  object$y_sim <- NULL

  ## get x and y limit
  if (is.null(xlim)) xlim <- c(NA_real_, NA_real_)
  if (is.null(ylim)) ylim <- c(NA_real_, NA_real_)

  ## stat helper function to get left/right points from respective mid points
  calc_confint_polygon <- ggplot2::ggproto("calc_confint_polygon", ggplot2::Stat,

    # Required as we operate on groups (facetting)
    compute_group = function(data, scales) {
      ## Manipulate object  #TODO: (ML) Could maybe be improved?
      nd <- data.frame(
        x = c(data$x_ci_lwr, rev(data$x_ci_upr)),
        y = c(data$y_ci_lwr, rev(data$y_ci_upr))
      )
      nd
    },

    # Tells us what we need
    required_aes = c("x_ci_lwr", "x_ci_upr", "y_ci_lwr", "y_ci_upr")
  )

  ## recycle arguments for plotting to match the number of groups (for geom w/o aes)
  plot_arg <- data.frame(1:n,
    fill, colour, size, ref, linetype, confint, alpha_min, shape
  )[, -1]

  ## recycle arguments for plotting to match the object rows (for geom w/ aes)
  ## FIXME: (ML) Why does it need to be equal to the length of the object and not to the number of groups?
  plot_arg2 <- data.frame(1:n, size, colour, ref, confint, fill, linetype)[, -1]
  plot_arg2 <- as.data.frame(lapply(plot_arg2, rep, each = nrow(object) / n))

  ## set alpha for polygon
  plot_arg$fill <- sapply(seq_along(plot_arg$fill), function(idx)
    set_minimum_transparency(plot_arg$fill[idx], alpha_min = plot_arg$alpha_min[idx]))

  ## actual plotting
  rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y")) +
    ggplot2::geom_point(ggplot2::aes_string(colour = "group", shape = "group", size = "group"),
      show.legend = FALSE) 

  ## add conf
  if (!identical(confint, FALSE) && all(c("x_ci_lwr", "x_ci_upr", "y_ci_lwr", "y_ci_upr") %in% names(object))) {
    rval <- rval + 
      ggplot2::geom_polygon(ggplot2::aes_string(x_ci_lwr = "x_ci_lwr", x_ci_upr = "x_ci_upr", 
        y_ci_lwr = "y_ci_lwr", y_ci_upr = "y_ci_upr", fill = "group"), 
        stat = calc_confint_polygon, show.legend = FALSE) 
  }

  ## add ref lines
  if (!identical(ref, FALSE)) {
    if (isTRUE(ref)) plot_arg$ref <- "black"
    rval <- rval + 
      ggplot2::geom_function(
        fun = ref_fun, 
        args = list(n = NROW(object), level = 0.95, which = "lower"),
        linetype = 2, col = plot_arg$ref
      ) +  
      ggplot2::geom_function(
        fun = ref_fun,
        args = list(n = NROW(object), level = 0.95, which = "upper"),
        linetype = 2, col = plot_arg$ref
      ) + 
      ggplot2::geom_hline(yintercept = 0, col = plot_arg$ref, linetype = 2)

    ## do not include ref lines in ylim
    ylim[1] <- ifelse(is.na(ylim)[1], min(object$y), ylim[1])
    ylim[2] <- ifelse(is.na(ylim)[2], max(object$y), ylim[2])
  }

  ## set the colors, shapes, etc.
  rval <- rval +
    ggplot2::scale_colour_manual(values = plot_arg$colour) +
    ggplot2::scale_fill_manual(values = plot_arg$fill) +
    ggplot2::scale_shape_manual(values = plot_arg$shape) +
    ggplot2::scale_size_manual(values = plot_arg$size) +
    ggplot2::scale_linetype_manual(values = plot_arg$linetype)

  ## add legend
  if (legend) {
    rval <- rval + ggplot2::labs(colour = "Model") +
      ggplot2::guides(colour = "legend", size = "none", linetype = "none")
  } else {
    rval <- rval + ggplot2::guides(colour = "none", size = "none", linetype = "none")
  }

  ## set x and y limits 
  rval <- rval + ggplot2::coord_cartesian(xlim= xlim, ylim = ylim, expand = 0.01)

  ## grouping (if any)
  if (!single_graph && n > 1L) {
    rval <- rval + ggplot2::facet_grid(group ~ .) 
  }

  ## annotation
  rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

  ## return ggplot object
  rval
}
