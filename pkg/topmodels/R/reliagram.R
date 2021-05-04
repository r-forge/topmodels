## Reliagram (reliability diagram)
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
                              breaks = seq(0, 1, by = 0.1),
                              quantiles = 0.5,
                              thresholds = NULL,
                              confint = TRUE,
                              confint_level = 0.95,
                              confint_nboot = 250,
                              single_graph = FALSE,
                              xlab = "Forecast probability",
                              ylab = "Observed relative frequency",
                              main = NULL,
                              ...) {

  ## sanity checks
  ## `object` and `newdata` w/i `newrepsone()`; `breaks w/i `cut()`, `...` w/i `plot()`;
  ## `confint` in `polygon()`
  stopifnot(is.logical(plot))
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

  ## data + thresholds
  y <- newresponse(object, newdata = newdata)
  if (is.null(thresholds)){
    thresholds <- quantile(y, probs = quantiles, na.rm = TRUE)
    thresholds <- as.numeric(thresholds)
    thresholds_text <- sprintf("q_%.2f", signif(quantiles, 2))
  } else {
    thresholds
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
    mean_pr <- as.numeric(
      aggregate(
        pred[, idx],
        by = list(prob = cut(pred[, idx], breaks, include.lowest = TRUE)),
        FUN = mean,
        drop = FALSE
      )[, "x"]
    )

    ## consistency resampling from Broecker (2007)
    if (!identical(confint, FALSE)) {
      obs_rf_boot <- vector("list", length = N)
      for (i in 1:confint_nboot) {

        ## take bootstrap sample from predictions (surrogate forecasts)
        xhat <- sample(pred[, idx], replace = TRUE)

        ## surrogate observations that are reliable by construction
        yhat <- runif(N) < xhat

        ## compute observed relative frequencies of the surrogate observations
        obs_rf_boot[[i]] <- as.numeric(
          aggregate(
            yhat,
            by = list(prob = cut(xhat, breaks, include.lowest = TRUE)),
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
      obs_rf,
      mean_pr,
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
    attr(rval_i, "bs") <- mean(y[, idx] - pred[, idx])^2 # FIXME: (ML) Does not match for minimum != 0
    class(rval_i) <- c("reliagram", "data.frame")

    rval[[idx]] <- rval_i
  }

  ## combine different groups
  rval <- do.call(c, rval)

  ## plot by default
  if (plot) {
    try(plot(rval, confint = confint, single_graph = single_graph, ...))
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
                           type = "b", # single or n values
                           add_info = TRUE, # single or n values
                           extend_left = NULL, # either null or logical of length 1 / n
                           extend_right = NULL, # either null or logical of length 1 / n
                           ...) {
  ## sanity checks
  ## lengths of all arguments are checked by recycling; `ref` w/i `abline()`; `xlim`, `ylim`,
  ## `xlab`, `ylab`, `main`, `col`, `fill`, `lwd`, `pch`, `lty`, `type` and `...` w/i `plot()`
  ## `confint` w/i `polygon()`
  ## `alpha_min` w/i colorspace fun 
  stopifnot(is.logical(single_graph))
  stopifnot(is.numeric(minimum), all(minimum >= 0))
  stopifnot(is.logical(add_info))
  stopifnot(is.null(extend_left) || is.logical(extend_left))
  stopifnot(is.null(extend_right) || is.logical(extend_right))

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  if (is.null(extend_left)) extend_left <- NA
  if (is.null(extend_right)) extend_right <- NA
  plot_arg <- data.frame(1:n, minimum, confint, ref,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]], 
    col, fill, alpha_min, lwd, pch, lty, type, add_info, extend_left, extend_right
  )[, -1]

  ## annotation
  if (single_graph) {
    if (is.null(xlab)) xlab <- "Forecast probability"
    if (is.null(ylab)) ylab <- "Observed relative frequency"
    if (is.null(main)) main <- "Reliability Diagram"
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

    ## get lines with sufficient observations
    min_idx <- which(d$n_pred >= plot_arg$minimum[j])

    ## check where confint should be extended
    if (is.na(plot_arg$extend_left[j])) plot_arg$extend_left[j] <- 1 %in% min_idx
    if (is.na(plot_arg$extend_right[j])) plot_arg$extend_right[j] <- nrow(d) %in% min_idx

    ## modify main using subscript for quantiles
    if (grepl("threshold = q_[0-9]+\\.?([0-9]+)?)$", main[j])){
      tmp_quantile <- regmatches(main[j], regexpr("q_[0-9]+\\.?([0-9]+)?", main[j]))
      tmp_quantile <- regmatches(tmp_quantile, regexpr("[0-9]+\\.?([0-9]+)?", tmp_quantile))
      tmp_text <- sub('q_[0-9]+\\.?([0-9]+)?)', '', main[j])
      main[j] <- as.expression(bquote(.(tmp_text)*q[.(tmp_quantile)]*")"))
    }

    ## trigger plot
    if (j == 1 || (!single_graph && j > 1)) {
      plot(0, 0,
        type = "n", xlim = c(plot_arg$xlim1[j], plot_arg$xlim2[j]),
        ylim = c(plot_arg$ylim1[j], plot_arg$ylim2[j]), xlab = xlab[j],
        ylab = ylab[j], main = main[j],
        xaxs = "i", yaxs = "i", ...
      )
      box()
    }

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]
      polygon(
        na.omit(c(
          ifelse(plot_arg$extend_left[j], 0, NA),
          d[min_idx, "mean_pr"],
          ifelse(plot_arg$extend_right[j], 1, NA),
          rev(d[min_idx, "mean_pr"]),
          ifelse(plot_arg$extend_left[j], 0, NA)
        )),
        na.omit(c(
          ifelse(plot_arg$extend_left[j], 0, NA),
          d[min_idx, "ci_lwr"],
          ifelse(plot_arg$extend_right[j], 1, NA),
          rev(d[min_idx, "ci_upr"]),
          ifelse(plot_arg$extend_left[j], 0, NA)
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

    ## plot reference line
    if (j == 1 || (!single_graph && j > 1)) {
      if (!identical(plot_arg$ref[j], FALSE)) {
        if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
        abline(0, 1, col = plot_arg$ref[j], lty = 2)
      }
    } 

    ## plot reliability line
    lines(obs_rf ~ mean_pr, d[min_idx, ],
      type = plot_arg$type[j], lwd = plot_arg$lwd[j],
      pch = plot_arg$pch[j], lty = plot_arg$lty[j], col = plot_arg$col[j], ...
    )

    ## print info
    if (single_graph && j == n && plot_arg$add_info[j]) {
      legend(
        "bottomright",
        sprintf("%s = %.4f", unlist(attr(d, "main")), signif(attr(d, "bs"), 4)),
        pch = plot_arg$pch,
        col = plot_arg$col,
        bty = "n",
        title = "Brier Score",
        y.intersp = 0.9
      )
    } else if (!single_graph && plot_arg$add_info[j]) {
      legend(
        "bottomright",
        sprintf("BS = %.4f", signif(attr(d, "bs")[j], 4)),
        bty = "n"
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
                            extend_left = NULL, # either null or logical of length 1 / n
                            extend_right = NULL, # either null or logical of length 1 / n
                            ...) {
  ## sanity checks
  ## lengths of all arguments are checked by recycling,
  ## `col`, `lwd`, `pch`, `lty`, `type` and `...` w/i `plot()`
  ## `ref` w/i `abline()`; `confint` w/i `polygon()`
  ## `alpha_min` w/i colorspace fun 
  stopifnot(is.numeric(minimum), all(minimum >= 0))
  stopifnot(is.null(extend_left) || is.logical(extend_left))
  stopifnot(is.null(extend_right) || is.logical(extend_right))

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  if (is.null(extend_left)) extend_left <- NA
  if (is.null(extend_right)) extend_right <- NA
  plot_arg <- data.frame(
    1:n, minimum, confint, ref, col, fill, alpha_min, lwd, pch, lty, type, extend_left, extend_right
  )[, -1]

  ## plotting function
  reliagramplot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get lines with sufficient observations
    min_idx <- which(d$n_pred >= plot_arg$minimum[j])

    ## check where confint should be extended
    if (is.na(plot_arg$extend_left[j])) plot_arg$extend_left[j] <- 1 %in% min_idx
    if (is.na(plot_arg$extend_right[j])) plot_arg$extend_right[j] <- nrow(d) %in% min_idx

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]
      polygon(
        na.omit(c(
          ifelse(plot_arg$extend_left[j], 0, NA),
          d[min_idx, "mean_pr"],
          ifelse(plot_arg$extend_right[j], 1, NA),
          rev(d[min_idx, "mean_pr"]),
          ifelse(plot_arg$extend_left[j], 0, NA)
        )),
        na.omit(c(
          ifelse(plot_arg$extend_left[j], 0, NA),
          d[min_idx, "ci_lwr"],
          ifelse(plot_arg$extend_right[j], 1, NA),
          rev(d[min_idx, "ci_upr"]),
          ifelse(plot_arg$extend_left[j], 0, NA)
        )),
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min = plot_arg$alpha_min[j]), 
        border = NA
      )
    }

    ## plot reference line
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
      abline(0, 1, col = plot_arg$ref[j], lty = 2)
    }

    ## plot reliability line
    lines(obs_rf ~ mean_pr, d[min_idx, ],
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

  ## list of pithists
  rval <- list(...)

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
  class(rval) <- c("reliagram", "data.frame")
  return(rval)
}



autoplot.reliagram <- function(object, ...) {
  NULL
}


### methods might only add custom y and thresholds  ## FIXME: (ML) What is that for?
# reliagram.crch <- function(object, newdata = NULL,
#  breaks = seq(0, 1, by = 0.1), thresholds = NULL, ...)
# {
#  y <- if((missing(newdata) || is.null(newdata)) && !is.null(object$y)) object$y else newresponse(object, newdata = newdata)
#  if(is.null(thresholds)) thresholds <- object$left
#  reliagram(object, breaks = breaks, thresholds = thresholds, y = y, ...)
# }
