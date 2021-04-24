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
                              probs = 0.5, 
                              thresholds = NULL, 
                              y = NULL, 
                              confint = TRUE, 
                              confint_level = 0.95,
                              confint_nboot = 250, 
                              xlab = "Forecast probability", 
                              ylab = "Observed relative frequency",
                              main = NULL,
                              ...) {

  ## sanity checks (breaks is checked within `cut()`)
  stopifnot(is.logical(plot))
  stopifnot(is.numeric(probs) && is.null(dim(probs)))
  stopifnot(is.null(thresholds) || (is.numeric(thresholds) && is.null(dim(thresholds))))
  stopifnot(is.null(y) || (is.numeric(y) && is.null(dim(y))))
  stopifnot(is.logical(confint))
  stopifnot(
    is.numeric(confint_level) &&
    is.null(dim(confint_level)) &&
    length(confint_level) == 1 &&
    confint_level >= 0 &&
    confint_level <= 1
  )
  stopifnot(
    is.numeric(confint_nboot) &&
    is.null(dim(confint_nboot)) &&
    length(confint_level) == 1 &&
    confint_level >= 0
  )

  ## fix length of annotations
  if(length(xlab) < length(probs)) xlab <- rep(xlab, length.out = length(probs))
  if(length(ylab) < length(probs)) ylab <- rep(ylab, length.out = length(probs))
  if(is.null(main)) {
    main <- deparse(substitute(object))
    main <- sprintf("%s (prob = %.2f)", main, probs)
  }

  ## data + thresholds
  y <- if(is.null(y)) newresponse(object, newdata = newdata) else as.numeric(y)
  ## FIXME: (ML) This can lead to an error
  ## R> x <- 1
  ## R> x <- if(is.null(x)) 2
  ## R> x
  ## NULL 
  thresholds <- if(is.null(thresholds)) quantile(y, probs = probs, na.rm = TRUE) 
  thresholds <- as.numeric(thresholds)

  ## predicted probabilities
  pred <- procast(object, newdata = newdata, type = "probability", at = matrix(thresholds, nrow = 1L),
    drop = FALSE)

  ## get and prepare observations
  y <- sapply(thresholds, function(x) y <= x)

  ## define convenience variables
  N <- NROW(y)

  ## loop over all probs
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
    if (confint) {
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

    ## return value
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
    attr(rval_i, "prob") <- probs[idx]
    attr(rval_i, "confint_level") <- confint_level
    attr(rval_i, "bs") <- mean(y[, idx] - pred[, idx])^2  # FIXME: (ML) Does not match for minimum != 0
    class(rval_i) <- c("reliagram", "data.frame")

    rval[[idx]] <- rval_i
  }

  ## combine different groups
  rval <- do.call(c, rval)

  ## plot by default
  if (plot) {
    plot(rval, confint = confint, ...)
  }

  ## return invisibly
  invisible(rval)
}


## actual drawing
plot.reliagram <- function(x, 
                           minimum = 0,
                           confint = TRUE, 
                           abline = TRUE,
                           xlim = c(0, 1), 
                           ylim = c(0, 1), 
                           xlab = NULL,
                           ylab = NULL,
                           main = NULL,
                           col = "black",
                           fill = adjustcolor(1, alpha.f = 0.2),
                           lwd = 2, 
                           pch = 19,
                           type = "b", 
                           print_info = FALSE, 
                           extend_left = NULL, 
                           extend_right = NULL, 
                           ...) {

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## annotation
  if (is.null(xlab)) xlab <- TRUE
  if (is.null(ylab)) ylab <- TRUE
  if (is.null(main)) main <- TRUE
  xlab <- rep(xlab, length.out = n)
  ylab <- rep(ylab, length.out = n)
  main <- rep(main, length.out = n)
  if (is.logical(xlab)) xlab <- ifelse(xlab, attr(x, "xlab"), "")
  if (is.logical(ylab)) ylab <- ifelse(ylab, attr(x, "ylab"), "")
  if (is.logical(main)) main <- ifelse(main, attr(x, "main"), "")

  ## plotting function
  reliagramplot1 <- function(d, ...)  {

    ## get lines with sufficient observations
    min_idx <- which(d$n_pred >= minimum)

    ## check where ci should be extended
    if(is.null(extend_left)) extend_left <- 1 %in% min_idx
    if(is.null(extend_right)) extend_right <- nrow(d) %in% min_idx

    ## get group index
    j <- unique(d$group)

    ## trigger plot
    plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab[j], ylab = ylab[j], main = main[j], 
      xaxs = "i", yaxs = "i", ...)
    box()
  
    ## plot ci
    if (confint & !is.na(attr(d, "confint_level"))) {
      polygon(
        na.omit(c(
          ifelse(extend_left, 0, NA), 
          d[min_idx, "mean_pr"], 
          ifelse(extend_right, 1, NA), 
          rev(d[min_idx, "mean_pr"]), 
          ifelse(extend_left, 0, NA)
        )), 
        na.omit(c(
          ifelse(extend_left, 0, NA), 
          d[min_idx, "ci_lwr"], 
          ifelse(extend_right, 1, NA), 
          rev(d[min_idx, "ci_upr"]), 
          ifelse(extend_left, 0, NA)
        )), 
        col = fill, border = NA
      )
    }

    ## plot perfect prediction
    if(abline) abline(0, 1, lty = 2)

    ## plot reliability line
    lines(obs_rf ~ mean_pr, d[min_idx, ], type = type, lwd = lwd, pch = pch, col = col, ...)

    ## print info
    if(print_info) text(0, 1, paste0("BS = ", signif(attr(x, "bs"), 4)), adj = c(0,1))
  }

  ## draw plots
  if (n > 1L) par(mfrow = n2mfrow(n))
  for (i in 1L:n) reliagramplot1(x[x$group == i, ], ...)
}


lines.reliagram <- function(x, ...) {
  NULL
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
  attr(rval, "bs") <- bs
  attr(rval, "prob") <- prob
  class(rval) <- c("reliagram", "data.frame")
  return(rval)
}

## ggplot2 interface
autoplot.reliagram <- function(object, ...) {
  NULL
}


### methods might only add custom y and thresholds  ## FIXME: (ML) What is that for?
#reliagram.crch <- function(object, newdata = NULL,
#  breaks = seq(0, 1, by = 0.1), thresholds = NULL, ...)
#{
#  y <- if((missing(newdata) || is.null(newdata)) && !is.null(object$y)) object$y else newresponse(object, newdata = newdata)
#  if(is.null(thresholds)) thresholds <- object$left
#  reliagram(object, breaks = breaks, thresholds = thresholds, y = y, ...)
#}
