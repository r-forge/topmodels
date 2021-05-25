# Reliagram (reliability diagram)
##
## - Oabserved y in-sample or out-of-sample (n x 1)
## - Thresholds in y (k x 1)
## - Predicted probabilities F_y(thresh) at thresholds (n x k)
## - Breaks for predicted probabilities in [0, 1] (m x 1)
##
## - Cut probabilities at breaks -> (m-1) groups
## - Cut y at thresholds, aggregate by groups -> (m-1) x k proportions

## Functions:
## - reliagram() generic plus default method
## - Return object of class "reliagram" that is plotted by default
## - But has plot=FALSE so that suitable methods can be added afterwards
## - At least methods: plot(), autoplot(), lines()

reliagram <- function(object, ...) {
  UseMethod("reliagram")
}


reliagram.default <- function(object,
                              newdata = NULL,
                              plot = TRUE,
                              class = NULL,
                              breaks = seq(0, 1, by = 0.1),
                              quantiles = 0.5,
                              thresholds = NULL,
                              confint = TRUE,
                              confint_level = 0.95,
                              confint_nboot = 250,
                              confint_seed = 1,
                              single_graph = FALSE,
                              xlab = "Forecast probability",
                              ylab = "Observed relative frequency",
                              main = NULL,
                              ...) {

  ## sanity checks
  ## `object` and `newdata` w/i `newrepsone()`; `breaks w/i `cut()`, `...` w/i `plot()`;
  ## `confint` in `polygon()`
  stopifnot(is.numeric(quantiles), is.null(dim(quantiles)))
  stopifnot(is.null(thresholds) || (is.numeric(thresholds) && is.null(dim(thresholds))))
  stopifnot(
    is.numeric(confint_level),
    length(confint_level) == 1,
    confint_level >= 0,
    confint_level <= 1
  )
  stopifnot(
    is.numeric(confint_nboot),
    length(confint_nboot) == 1,
    confint_nboot >= 0
  )
  stopifnot(is.logical(single_graph))
  stopifnot(length(xlab) == 1 || length(xlab) == length(quantiles))
  stopifnot(length(ylab) == 1 || length(ylab) == length(quantiles))
  stopifnot(is.null(main) || (length(main) == 1 || length(main) == length(quantiles)))
  stopifnot(is.numeric(breaks) || is.function(breaks))

  ## guess plotting flavor
  if (isFALSE(plot)) {
    plot <- "none"
  } else if (isTRUE(plot)) {
    plot <- if("package:ggplot2" %in% search()) "ggplot2" else "base"
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
    class <- if("tibble" %in% loadedNamespaces()) "tibble" else "data.frame"
  }
  class <- try(match.arg(class, c("tibble", "data.frame")))
  stopifnot(
    "`class` must either be NULL, or match the arguments 'tibble', or 'data.frame'" =
    !inherits(class, "try-error")
  )

  ## data + thresholds
  y <- newresponse(object, newdata = newdata)
  if (is.null(thresholds)){
    thresholds <- quantile(y, probs = quantiles, na.rm = TRUE)
    thresholds <- as.numeric(thresholds)
    thresholds_text <- sprintf("q_%.2f", signif(quantiles, 2))
  } else {
    thresholds_text <- as.character(signif(thresholds, 2))
  }
  
  ## fix length of annotations
  if (length(xlab) < length(quantiles)) xlab <- rep(xlab, length.out = length(quantiles))
  if (length(ylab) < length(quantiles)) ylab <- rep(ylab, length.out = length(quantiles))
  if (is.null(main)) {
    main <- deparse(substitute(object))
    main <- sprintf("%s (threshold = %s)", main, thresholds_text)
  }

  ## predicted probabilities  # FIXME: (ML) Check format and dim of thresholds and pred
  pred <- procast(object,
    newdata = newdata, type = "probability", at = matrix(thresholds, nrow = 1L),
    drop = FALSE
  )

  ## make sure lengths match (can't really happen after no arg `y` exists anymore)
  stopifnot(NROW(pred) == length(y)) 

  ## get and prepare observations
  y <- sapply(thresholds, function(x) y <= x)

  ## define convenience variables
  N <- NROW(y)

  ## loop over all quantiles
  rval <- vector(mode = "list", length = NCOL(y))
  for (idx in 1:NCOL(y)) {

    ## calculate breaks
    if (is.function(breaks)) {
      try(breaks <- as.numeric(breaks(pred[, idx])))
      stopifnot("`breaks` function must produce a numeric valid to be used w/i `cut()`" = is.numeric(breaks))
    }

    ## compute number of prediction and idx for minimum number of prediction per probability subset
    n_pred <- aggregate(
      pred[, idx],
      by = list(prob = cut(pred[, idx], breaks, include.lowest = TRUE)),
      FUN = length,
      drop = FALSE
    )[, "x"]

    ## compute observed relative frequencies of positive examples (obs_rf)
    obs_rf <- as.numeric(
      aggregate(
        y[, idx],
        by = list(prob = cut(pred[, idx], breaks, include.lowest = TRUE)),
        FUN = mean,
        drop = FALSE
      )[, "x"]
    )

    ## compute mean predicted probability (mean_pr)
    tmp_mean_pr <- aggregate(
      pred[, idx],
      by = list(prob = cut(pred[, idx], breaks, include.lowest = TRUE)),
      FUN = mean,
      drop = FALSE
    )
    mean_pr <- as.numeric(tmp_mean_pr[, "x"])

    ## calculate pred with mean values (for calulating bs) 
    lookup <- as.numeric(cut(pred[, idx], breaks, include.lowest = TRUE))
    pred_bin <- tmp_mean_pr$x[lookup]

    ## consistency resampling from Broecker (2007)
    if (!identical(confint, FALSE)) {
      set.seed(confint_seed)
      obs_rf_boot <- vector("list", length = N)
      for (i in 1:confint_nboot) {

        ## take bootstrap sample from predictions (surrogate forecasts)
        pred_hat <- sample(pred[, idx], replace = TRUE)

        ## surrogate observations that are reliable by construction
        yhat <- runif(N) < pred_hat

        ## compute observed relative frequencies of the surrogate observations
        obs_rf_boot[[i]] <- as.numeric(
          aggregate(
            yhat,
            by = list(prob = cut(pred_hat, breaks, include.lowest = TRUE)),
            FUN = mean,
            drop = FALSE
          )[, "x"]
        )
      }
      obs_rf_boot <- do.call("rbind", obs_rf_boot)

      ## compute lower and upper limits for reliable forecasts
      confint_prob <- (1 - confint_level) / 2
      confint_prob <- c(confint_prob, 1 - confint_prob)
      ci_lwr <- apply(obs_rf_boot, 2, quantile, prob = confint_prob[1], na.rm = TRUE)
      ci_upr <- apply(obs_rf_boot, 2, quantile, prob = confint_prob[2], na.rm = TRUE)
    } else {
      ci_lwr <- NA
      ci_upr <- NA
      confint_level <- NA
    }

    ## collect everything as data.frame
    rval_i <- data.frame(
      x = mean_pr,
      y = obs_rf,
      bin_lwr = breaks[-length(breaks)],
      bin_upr = breaks[-1] ,
      n_pred,
      ci_lwr,
      ci_upr
    )

    ## attributes for graphical display
    attr(rval_i, "xlab") <- xlab[idx]
    attr(rval_i, "ylab") <- ylab[idx]
    attr(rval_i, "main") <- main[idx]
    attr(rval_i, "threshold") <- thresholds_text[idx]
    attr(rval_i, "confint_level") <- confint_level

    ## add bs, rel, res, and unc
    ## TODO: (ML) Check what to do with NAs exactly
    ## NOTE: (ML) Here the unique forecasts equal the mean forecasts per bin (as in `verification` pkg):
    ##  * Hence, BS is not independent to the bins
    ##  * Hence, BS should vary conditional on the minimum
    attr(rval_i, "bs") <- mean((pred_bin - y[, idx])^2)
    attr(rval_i, "rel") <- sum(n_pred * (mean_pr - obs_rf)^2, na.rm = TRUE) / sum(n_pred, na.rm = TRUE)
    attr(rval_i, "res") <- sum(n_pred * (obs_rf - mean(y[, idx]))^2, na.rm = TRUE) / sum(n_pred, na.rm = TRUE)
    attr(rval_i, "unc") <- mean(y[, idx]) * (1 - mean(y[, idx]))

    ## add class
    if (class == "data.frame") {
      class(rval_i) <- c("reliagram", "data.frame")
    } else {
      rval_i <- tibble::as_tibble(rval_i)
      class(rval_i) <- c("reliagram", class(rval_i))
    }

    rval[[idx]] <- rval_i
  }

  ## combine different groups
  rval <- do.call(c, rval)

  ## plot by default
  if (plot == "ggplot2") {
    try(print(ggplot2::autoplot(rval,
      confint = confint, single_graph = single_graph, ...)
    ))
  } else if (plot == "base") {
    try(plot(rval,
      confint = confint, single_graph = single_graph, ...
    ))
  }

  ## return invisibly
  invisible(rval)
}


plot.reliagram <- function(x,
                           single_graph = FALSE, # logical
                           minimum = 0, # single or n values
                           confint = TRUE, # single or n values, logical or color 
                           ref = TRUE, # single or n values, logical or color
                           xlim = c(0, 1), # single vector of lenght 2, or list w/ 2 vectors of length n
                           ylim = c(0, 1), # single vector of lenght 2, or list w/ 2 vectors of length n
                           xlab = NULL, # single or n values
                           ylab = NULL, # single or n values
                           main = NULL, # single or n values
                           col = "black", # single or n colors
                           fill = adjustcolor("black", alpha.f = 0.2), # single or n colors (convenience for confint)
                           alpha_min = 0.2, # single or n values  
                           lwd = 2, # single or n values
                           pch = 19, # single or n values
                           lty = 1, # single or n values
                           type = NULL, # single or n values
                           add_hist = TRUE,
                           add_info = TRUE, # single or n values
                           add_rug = TRUE,
                           add_min = TRUE, 
                           axes = TRUE,
                           box = TRUE,
                           ...) {
  ## sanity checks
  ## lengths of all arguments are checked by recycling; `ref` w/i `abline()`; `xlim`, `ylim`,
  ## `xlab`, `ylab`, `main`, `col`, `fill`, `lwd`, `pch`, `lty`, `type` and `...` w/i `plot()`
  ## `confint` w/i `polygon()`
  ## `alpha_min` w/i colorspace fun 
  stopifnot(is.logical(single_graph))
  stopifnot(is.numeric(minimum), all(minimum >= 0))
  stopifnot(is.logical(add_info))
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## intern for single plot, single_graph is set to FALSE
  if (n == 1) {
    single_graph <- FALSE
  }

  ## determine if points should be plotted
  if (is.null(type)) type <- ifelse(table(x$group) > 20L, "l", "b")

  ## recycle arguments for plotting to match the number of groups
  if (is.list(xlim)) xlim <- as.data.frame(do.call("rbind", xlim))
  if (is.list(ylim)) ylim <- as.data.frame(do.call("rbind", ylim))
  plot_arg <- data.frame(1:n, minimum, confint, ref,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]], 
    col, fill, alpha_min, lwd, pch, lty, type, add_hist, add_info, add_rug, add_min,
    axes, box
  )[, -1]

  ## annotation
  if (single_graph) {
    if (is.null(xlab)) xlab <- "Forecast probability"
    if (is.null(ylab)) ylab <- "Observed relative frequency"
    if (is.null(main)) main <- "Reliability diagram"
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

  ## function to trigger figure and plot confint 
  reliagram_trigger <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get bins with sufficient observations
    min_idx <- which(d$n_pred >= plot_arg$minimum[j])
    if (length(min_idx) == 0) {
      stop(sprintf("no bin has sufficent cases for defined minimum = %s\n", plot_arg$minimum[j]))
    }

    ## get minimum and maximum breaks to decide if extend confint polygon to the corners
    extend_left <- min(d$bin_lwr) == min(d[min_idx, "bin_lwr"])
    extend_right <- max(d$bin_upr) == max(d[min_idx, "bin_upr"])

    ## modify main using subscript for quantiles
    if (grepl("threshold = q_[0-9]+\\.?([0-9]+)?)$", main[j])){
      tmp_quantile <- regmatches(main[j], regexpr("q_[0-9]+\\.?([0-9]+)?", main[j]))
      tmp_quantile <- regmatches(tmp_quantile, regexpr("[0-9]+\\.?([0-9]+)?", tmp_quantile))
      tmp_text <- sub('q_[0-9]+\\.?([0-9]+)?)', '', main[j])
      main[j] <- as.expression(bquote(bold(.(tmp_text)*q[.(tmp_quantile)]*")")))
    }

    ## trigger plot
    if (j == 1 || (!single_graph && j > 1)) {
      plot(0, 0,
        type = "n", xlim = c(plot_arg$xlim1[j], plot_arg$xlim2[j]),
        ylim = c(plot_arg$ylim1[j], plot_arg$ylim2[j]), xlab = xlab[j],
        ylab = ylab[j], main = main[j],
        xaxs = "i", yaxs = "i", axes = FALSE, ...
      )
      if(plot_arg$axes[j]) {
        axis(1)
        axis(2)
      }
      if(plot_arg$box[j]) {
        box()
      }
    }

    ## plot reference line
    if (j == 1 || (!single_graph && j > 1)) {
      if (!identical(plot_arg$ref[j], FALSE)) {
        if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
        abline(0, 1, col = plot_arg$ref[j], lty = 2, lwd = 1.25)
      }
    } 

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]
      polygon(
        na.omit(c(
          ifelse(extend_left, 0, NA),
          d[min_idx, "x"],
          ifelse(extend_right, 1, NA),
          rev(d[min_idx, "x"]),
          ifelse(extend_left, 0, NA)
        )),
        na.omit(c(
          ifelse(extend_left, 0, NA),
          d[min_idx, "ci_lwr"],
          ifelse(extend_right, 1, NA),
          rev(d[min_idx, "ci_upr"]),
          ifelse(extend_left, 0, NA)
        )),
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min = plot_arg$alpha_min[j]), 
        border = NA
      )
    }
  }

  ## plotting function
  reliagram_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get lines with sufficient observations
    min_idx <- which(d$n_pred >= plot_arg$minimum[j])

    ## plot reliability line
    lines(y ~ x, d[min_idx, ],
      type = plot_arg$type[j], lwd = plot_arg$lwd[j],
      pch = plot_arg$pch[j], lty = plot_arg$lty[j], col = plot_arg$col[j], ...
    )

    ## plot points below minimum 
    if (!identical(plot_arg$add_min[j], FALSE)) {
      if (isTRUE(plot_arg$add_min[j])) plot_arg$add_min[j] <- 4
      points(y ~ x, d[-min_idx, ], pch = plot_arg$add_min[j], col = plot_arg$col[j], ...)
    }

    ## add rugs
    if (!single_graph && !identical(plot_arg$add_rug[j], FALSE)) {
      if (isTRUE(plot_arg$add_rug[j])) plot_arg$add_rug[j] <- plot_arg$col[j]
      tmp_rugs <- c(d$bin_lwr, d$bin_upr[NROW(d)])
      tmp_rugs <- tmp_rugs[tmp_rugs >= plot_arg$xlim1[j] & tmp_rugs <= plot_arg$xlim2[j]]
      rug(tmp_rugs, lwd = 1, ticksize = 0.02, col = plot_arg$add_rug[j])
    }

     
    ## add hist
    if (!single_graph && !identical(plot_arg$add_hist[j], FALSE)) {
      if (isTRUE(plot_arg$add_hist[j])) plot_arg$add_hist[j] <- "lightgray"
      tmp_x <- par("pin")[1]
      tmp_y <- par("pin")[2]
      tmp_height <- 0.3  * diff(c(plot_arg$ylim1[j], plot_arg$ylim2[j]))
      tmp_width <- (0.3 * diff(c(plot_arg$xlim1[j], plot_arg$xlim2[j]))) * tmp_y / tmp_x

      add_hist_reliagram(
        d$n_pred, 
        c(d$bin_lwr, d$bin_upr[NROW(d)]), 
        plot_arg$minimum[j],
        xpos = 0.05 * diff(c(plot_arg$xlim1[j], plot_arg$xlim2[j])) + plot_arg$xlim1[j], 
        ypos = 0.925 * diff(c(plot_arg$ylim1[j], plot_arg$ylim2[j])) - tmp_height + plot_arg$ylim1[j], 
        width = tmp_width,
        height = tmp_height,
        col = plot_arg$add_hist[j]
      )
    }

    ## print info
    if (!single_graph && plot_arg$add_info[j]) {
      legend(
        "bottomright",
        c(
          "BS", signif(attr(d, "bs")[j], 3), 
          "REL", signif(attr(d, "rel")[j], 3), 
          "RES", signif(attr(d, "res")[j], 3), 
          "UNC", signif(attr(d, "unc")[j], 3)
        ),
        cex = 0.8,
        ncol = 4,
        bty = "n",
        inset = c(0.01, 0.01)
      )
    }
  }

  ## draw plots
  if (!single_graph && n > 1L) par(mfrow = n2mfrow(n))

  ## draw polygons first
  if (single_graph) {
    for (i in 1L:n) reliagram_trigger(x[x$group == i, ], ...)
    for (i in 1L:n) reliagram_plot(x[x$group == i, ], ...)
  } else {
    for (i in 1L:n) {
      reliagram_trigger(x[x$group == i, ], ...)
      reliagram_plot(x[x$group == i, ], ...)
    }
  }
}


lines.reliagram <- function(x,
                            minimum = 0,
                            confint = FALSE,
                            ref = FALSE,
                            col = "black",
                            fill = adjustcolor("black", alpha.f = 0.2), # single or n colors (convenience for confint)
                            alpha_min = 0.2, # single or n values  
                            lwd = 2,
                            pch = 19,
                            lty = 1,
                            type = "b",
                            ...) {
  ## sanity checks
  ## lengths of all arguments are checked by recycling,
  ## `col`, `lwd`, `pch`, `lty`, `type` and `...` w/i `plot()`
  ## `ref` w/i `abline()`; `confint` w/i `polygon()`
  ## `alpha_min` w/i colorspace fun 
  stopifnot(is.numeric(minimum), all(minimum >= 0))

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(
    1:n, minimum, confint, ref, col, fill, alpha_min, lwd, pch, lty, type
  )[, -1]

  ## plotting function
  reliagramplot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get lines with sufficient observations
    min_idx <- which(d$n_pred >= plot_arg$minimum[j])

    ## get minimum and maximum breaks to decide if extend confint polygon to the corners
    extend_left <- min(d$bin_lwr) == min(d[min_idx, "bin_lwr"])
    extend_right <- max(d$bin_upr) == max(d[min_idx, "bin_upr"])

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]
      polygon(
        na.omit(c(
          ifelse(extend_left, 0, NA),
          d[min_idx, "x"],
          ifelse(extend_right, 1, NA),
          rev(d[min_idx, "x"]),
          ifelse(extend_left, 0, NA)
        )),
        na.omit(c(
          ifelse(extend_left, 0, NA),
          d[min_idx, "ci_lwr"],
          ifelse(extend_right, 1, NA),
          rev(d[min_idx, "ci_upr"]),
          ifelse(extend_left, 0, NA)
        )),
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min = plot_arg$alpha_min[j]), 
        border = NA
      )
    }

    ## plot reference line
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
      abline(0, 1, col = plot_arg$ref[j], lty = 2, lwd = 1.25)
    }

    ## plot reliability line
    lines(y ~ x, d[min_idx, ],
      type = plot_arg$type[j], lwd = plot_arg$lwd[j],
      pch = plot_arg$pch[j], lty = plot_arg$lty[j], col = plot_arg$col[j], ...
    )
  }

  ## draw plots
  for (i in 1L:n) {
    reliagramplot(x[x$group == i, ], ...)
  }
}


c.reliagram <- rbind.reliagram <- function(...) {

  ## list of reliagrams
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

  ## labels
  xlab <- unlist(lapply(rval, function(r) attr(r, "xlab")))
  ylab <- unlist(lapply(rval, function(r) attr(r, "ylab")))
  prob <- unlist(lapply(rval, function(r) attr(r, "prob")))
  confint_level <- unlist(lapply(rval, function(r) attr(r, "confint_level")))
  bs <- unlist(lapply(rval, function(r) attr(r, "bs")))
  rel <- unlist(lapply(rval, function(r) attr(r, "rel")))
  res <- unlist(lapply(rval, function(r) attr(r, "res")))
  unc <- unlist(lapply(rval, function(r) attr(r, "unc")))
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
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "prob") <- prob
  attr(rval, "confint_level") <- confint_level
  attr(rval, "bs") <- bs
  attr(rval, "rel") <- rel
  attr(rval, "res") <- res
  attr(rval, "unc") <- unc

  ## set class to data.frame or tibble
  if (class == "data.frame") {
    class(rval) <- c("reliagram", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("reliagram", class(rval))
  }

  return(rval)
}



autoplot.reliagram <- function(object, 
                               single_graph = FALSE, 
                               minimum = 0,
                               confint = TRUE, 
                               ref = TRUE, 
                               xlim = c(0, 1),
                               ylim = c(0, 1),
                               xlab = NULL, 
                               ylab = NULL, 
                               main = NULL, 
                               colour = "black",
                               fill = adjustcolor("black", alpha.f = 0.2), 
                               alpha_min = 0.2, 
                               size = 1, 
                               shape = 19, 
                               linetype = 1, 
                               type = NULL,
                               add_info = TRUE,
                               add_min = TRUE,
                               legend = FALSE,
                               ...) {

  ## convert always to data.frame
  object <- as.data.frame(object)

  ## determine grouping
  if (is.null(object$group)) object$group <- 1L
  n <- max(object$group)

  ## get annotations in the right lengths
  if(is.null(xlab)) xlab <- attr(object, "xlab")
  xlab <- paste(unique(xlab), collapse = "/")
  if(is.null(ylab)) ylab <- attr(object, "ylab")
  ylab <- paste(unique(ylab), collapse = "/")
  if(is.null(main)) main <- attr(object, "main")
  main <- make.unique(rep_len(main, n))

  ## prepare grouping
  object$group <- factor(object$group, levels = 1L:n, labels = main)

  ## get base style arguments
  add_arg <- list(...)
  if (!is.null(add_arg$pch)) shape <- add_arg$pch
  if (!is.null(add_arg$lwd)) size <- add_arg$lwd
  if (!is.null(add_arg$lty)) linetype <- add_arg$lty

  ## get x and y limit in correct format
  if (is.null(xlim)) xlim <- c(NA_real_, NA_real_)
  if (is.null(ylim)) ylim <- c(NA_real_, NA_real_)

  ## determine if points should be plotted
  if (is.null(type)) type <- ifelse(table(object$group) > 20L, "l", "b")

  ## set color to NA for not plotting
  type <- ifelse(type == "l", 0, 1)
  if (is.logical(ref)) ref <- ifelse(ref, 1, NA)
  if (is.logical(add_min)) add_min <- ifelse(add_min, 4, NA)

  ## get min and max breaks to decide if extend confint polygon to the corners
  min_break <- min(object$bin_lwr)
  max_break <- max(object$bin_upr)

  ## stat helper function to get left/right points from respective mid points
  calc_confint_polygon <- ggplot2::ggproto("calc_confint_polygon", ggplot2::Stat,

    # Required as we operate on groups (facetting)
    compute_group = function(data, scales) {
      ## Manipulate object  #TODO: (ML) Could possibly be improved?
      if (min(data$bin_lwr) == min_break & max(data$bin_upr) == max_break) {
        nd <- data.frame(
          x = c(0, data$x, 1, rev(data$x), 0),
          y = c(0, data$ci_lwr, 1, rev(data$ci_upr), 0)
        )
      } else if (min(data$bin_lwr) == min_break) {
        nd <- data.frame(
          x = c(0, data$x, rev(data$x), 0),
          y = c(0, data$ci_lwr, rev(data$ci_upr), 0)
        )
      } else if (max(data$bin_upr) == max_break) {
        nd <- data.frame(
          x = c(data$x, 1, rev(data$x)),
          y = c(data$ci_lwr, 1, rev(data$ci_upr))
        )
      } else { 
        nd <- data.frame(
          x = c(data$x, rev(data$x)),
          y = c(data$ci_lwr, rev(data$ci_upr))
        )
      }
      nd
    },

    # Tells us what we need
    required_aes = c("x", "ci_lwr", "ci_upr", "bin_lwr", "bin_upr")
  )

  ## recycle arguments for plotting to match the number of groups (for `scale_<...>_manual()`)
  plot_arg <- data.frame(1:n,
    colour, fill, size, linetype, confint, alpha_min, minimum
  )[, -1]

  ## prepare fill color for confint (must be done on vector to match args) 
  if (is.logical(plot_arg$confint)) {

    ## use fill and set alpha
    plot_arg$fill <- sapply(seq_along(plot_arg$fill), function(idx)
      set_minimum_transparency(plot_arg$fill[idx], alpha_min = plot_arg$alpha_min[idx]))

    ## set color to NA for not plotting
    plot_arg$fill[!plot_arg$confint] <- NA

  } else {

    ## use confint and set alpha 
    plot_arg$fill <- sapply(seq_along(plot_arg$confint), function(idx)
      set_minimum_transparency(plot_arg$confint[idx], alpha_min = plot_arg$alpha_min[idx]))
  }

  ## recycle arguments for plotting to match the length (rows) of the object (for geom w/ aes)
  plot_arg2 <- data.frame(1:n, ref, type, size, shape, minimum, add_min)[, -1]
  plot_arg2 <- as.data.frame(lapply(plot_arg2, rep, each = nrow(object) / n))

  ## set points to NA with no sufficent number of predictions
  idx_min <- which(object$n_pred < plot_arg2$minimum)
  object2 <- object[idx_min, ]
  object[idx_min, c("x", "y")] <- NA

  ## actual plotting
  rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y")) +
    ggplot2::geom_abline(ggplot2::aes_string(intercept = 0, slope = 1), 
      linetype = 2, colour = plot_arg2$ref) + 
    ggplot2::geom_polygon(ggplot2::aes_string(
        ci_lwr = "ci_lwr", ci_upr = "ci_upr", 
        bin_lwr = "bin_lwr", bin_upr = "bin_upr", fill = "group"
      ), 
      stat = calc_confint_polygon, show.legend = FALSE, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes_string(colour = "group", size = "group", linetype = "group"), 
      na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes_string(colour = "group"),
      alpha = plot_arg2$type, size = plot_arg2$size * 2, shape = plot_arg2$shape, 
      show.legend = FALSE, na.rm = TRUE) 

  rval <- rval + 
    ggplot2::geom_point(ggplot2::aes_string(x = "x", y = "y"), data = object2,
      alpha = plot_arg2$type[idx_min], size = plot_arg2$size[idx_min] * 2, shape = plot_arg2$add_min[idx_min],
      show.legend = FALSE, na.rm = TRUE)

  ## set the colors, shapes, etc. for the groups
  rval <- rval +
    ggplot2::scale_colour_manual(values = plot_arg$colour) +
    ggplot2::scale_fill_manual(values = plot_arg$fill) +
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
  rval <- rval + ggplot2::scale_x_continuous(limits = xlim, expand = c(0.01, 0.01))
  rval <- rval + ggplot2::scale_y_continuous(limits = ylim, expand = c(0.01, 0.01))

  ## grouping (if any)
  if (!single_graph && n > 1L) {
    rval <- rval + ggplot2::facet_grid(group ~ .) 
  }

  ## annotation
  rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

  ## return ggplot object
  rval
}


add_hist_reliagram <- function(n, 
                               breaks, 
                               minimum,
                               xpos, 
                               ypos, 
                               width = .2,
                               height = .2, 
                               col = "lightgray", 
                               main = NULL) {

  idx_min <- which(n < minimum)
  n[is.na(n)] <- 0 
  #n[is.na(n) | n < minimum] <- 0 
  max_n <- max(n)
  col <- rep(col, max_n)
  col[n < minimum] <- 0
  for (i in seq_along(n)) {
    x <- xpos + breaks[c(i, i + 1)] * width
    y <- ypos + c(0, n[i] / max_n * height)
    rect(x[1L], y[1L], x[2L], y[2L], col = col[i])
  }
  #points(
  #  na.omit(filter(breaks, c(0.5, 0.5), sides = 2))[idx_min] * width + xpos, 
  #  rep(0, length(n))[idx_min] + ypos, 
  #  pch = 4
  #)

  ytick <- pretty(c(0, max_n * .8), 4)
  ytick <- ytick[ytick > 0 & ytick < max_n]
  text(xpos + 1.05 * width, ypos + ytick / max_n * height, ytick, cex = .8, adj = c(0.0,0.5))
  segments(
    x0 = rep(xpos, length(ytick)), 
    x1 = rep(xpos + 1 * width, length(ytick)), 
    y0 = ypos + ytick / max_n * height,
    col = "darkgray", lwd = .5, lty = 2)
  segments(
    x0 = rep(xpos + 0.975 * width, length(ytick)), 
    x1 = rep(xpos + 1.025 * width, length(ytick)), 
    y0 = ypos + ytick / max_n * height,
    lwd = .5)

  if (length(idx_min) > 0) {
    text(xpos + -0.05 * width, ypos + minimum / max_n * height, "Min.", cex = .8, adj = c(1, 0.5))
    segments(
      x0 = rep(xpos, length(ytick)), 
      x1 = rep(xpos + 1 * width, length(ytick)), 
      y0 = ypos + minimum / max_n * height,
      col = "darkgray", lwd = .5, lty = 1)
    segments(
      x0 = rep(xpos - 0.025 * width, length(ytick)), 
      x1 = rep(xpos + 0.025 * width, length(ytick)), 
      y0 = ypos + minimum / max_n * height,
      lwd = .5)
  }
  if (is.null(main)) {
    main <- sprintf("N=%d", sum(n[n > minimum]))
  }
  text(xpos + width / 2, ypos + 1.1 * height, font = 2, cex = 0.8, main)
}


### methods might only add custom y and thresholds  ## FIXME: (ML) What is that for?
# reliagram.crch <- function(object, newdata = NULL,
#  breaks = seq(0, 1, by = 0.1), thresholds = NULL, ...)
# {
#  y <- if((missing(newdata) || is.null(newdata)) && !is.null(object$y)) object$y else newresponse(object, newdata = newdata)
#  if(is.null(thresholds)) thresholds <- object$left
#  reliagram(object, breaks = breaks, thresholds = thresholds, y = y, ...)
# }
