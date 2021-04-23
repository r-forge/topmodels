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
                              minimum = 0,
                              confint = TRUE, 
                              confint_level = 0.95,
                              confint_nboot = 250, 
                              extend_left = NULL,
                              extend_right = NULL,
                              xlab = NULL, 
                              ylab = NULL, 
                              main = NULL,
                              ...) {

  ## TODO: (ML) argument validy check

  ## TODO: (ML) annotations
  ## default annotation
  #if(is.null(xlab)) {
  #  xlab <- if(is.null(names(dimnames(object)))) {
  #    deparse(substitute(object))
  #  } else {
  #    names(dimnames(object))[1L]
  #  }
  #}
  #if(is.null(ylab)) {
  #  ylab <- if(scale == "raw") "Frequency" else "sqrt(Frequency)"
  #}
  if(is.null(main)) main <- deparse(substitute(object))

  ## data + thresholds
  y <- if(is.null(y)) newresponse(object, newdata = newdata) else as.numeric(y)
  thresholds <- if(is.null(thresholds)) quantile(y, probs = probs, na.rm = TRUE) 
  thresholds <- as.numeric(thresholds)

  ## predicted probabilities
  pred <- procast(object, newdata = newdata, type = "probability", at = matrix(thresholds, nrow = 1L),
    drop = FALSE)

  ## TODO: (Z) compute quantities to be plotted
  #rval <- data.frame("<cut prob>", "<aggregate y by prob and thresh>")

  ## get and prepare observations
  #y <- y <= thresholds  #FIXME: (ML) Dirty hack to tryout function
  y <- sapply(thresholds, function(x) y <= x)

  ## define convenience variables
  N <- NROW(y)

  rval <- vector(mode = "list", length = NCOL(y))
  for (idx in 1:NCOL(y)) {
    ## compute number of prediction and idx for minimum number of prediction per probability subset
    n_pred <- aggregate(
      pred[, idx], 
      by = list(prob = cut(pred[, idx], breaks, include.lowest = TRUE)), 
      FUN = length, 
      drop = FALSE
    )[, "x"]
    min_idx <- which(n_pred >= minimum)

    ## check if ci should be extended
    if(is.null(extend_left)) extend_left <- 1 %in% min_idx
    if(is.null(extend_right)) extend_right <- (length(breaks) - 1) %in% min_idx

    ## compute observed relative frequencies of positive examples (obs_rf)
    obs_rf <- as.numeric(
      aggregate(
        y[, idx], 
        by = list(prob = cut(pred[, idx], breaks, include.lowest = TRUE)), 
        FUN = mean, 
        drop = FALSE
        )[, "x"]
    )
    obs_rf <- obs_rf[min_idx]

    ## compute mean predicted probability (mean_pr)
    mean_pr <- as.numeric(
      aggregate(
        pred[, idx], 
        by = list(prob = cut(pred[, idx], breaks, include.lowest = TRUE)), 
        FUN = mean, 
        drop = FALSE
      )[, "x"]
    )
    mean_pr <- mean_pr[min_idx]

    ## consistency resampling from Broecker (2007) 
    obs_rf_boot <- vector("list", length = N)
    for (i in 1:confint_nboot) {
      ## take bootstrap sample from predictions (surrogate forecasts)
      xhat <- sample(pred[, idx], replace = TRUE)

      ## surrogate observations that are reliable by construction
      yhat <- runif(N) < xhat

      ## compute observed relative frequencies of the surrogate observations
      obs_rf_boot[[i]] <- as.numeric(
        aggregate(yhat, by = list(prob = cut(xhat, breaks, include.lowest = TRUE)), mean, drop = FALSE)[, "x"]
      )
    }
    obs_rf_boot <- do.call("rbind", obs_rf_boot)

    ## compute lower and upper limits for reliable forecasts
    confint_prob <- (1 - confint_level) / 2
    confint_prob <- c(confint_prob, 1 - confint_prob)
    lowerlimit <- apply(obs_rf_boot, 2, quantile, prob = confint_prob[1], na.rm = TRUE)[min_idx]
    upperlimit <- apply(obs_rf_boot, 2, quantile, prob = confint_prob[2], na.rm = TRUE)[min_idx]

    ## return value
    rval_i <- data.frame(
      obs_rf, 
      mean_pr, 
      lowerlimit, 
      upperlimit
    )

    ## attributes for graphical display
    attr(rval_i, "xlab") <- xlab
    attr(rval_i, "ylab") <- ylab
    attr(rval_i, "main") <- main
    attr(rval_i, "prob") <- probs[idx]
    attr(rval_i, "bs") <- mean(y[, idx] - pred[, idx])^2
    class(rval_i) <- c("reliagram", "data.frame")

    rval[[idx]] <- rval_i
  }

  ## Combine different groups
  rval <- do.call(c, rval)

  ## plot by default
  if (plot) {
    plot(rval, 
    confint = confint, extend_left = extend_left, extend_right = extend_right, ...)
  }

  ## return invisibly
  invisible(rval)
}


## actual drawing
plot.reliagram <- function(x, 
                           confint = TRUE, 
                           abline = TRUE,
                           xlim = c(0, 1), 
                           ylim = c(0, 1), 
                           xlab = "Forecast probability", 
                           ylab = "Observed relative frequency",
                           main = NULL,
                           lwd = 2, 
                           type = "b", 
                           print_info = FALSE, 
                           col = "black", 
                           fill = adjustcolor(1, alpha.f = 0.1),
                           extend_left = TRUE, 
                           extend_right = TRUE, 
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
  if (is.logical(main)) main <- ifelse(main, sprintf("%s (prob = %.2f)", attr(x, "main"), attr(x, "prob")), "")

  ## plotting function
  reliagramplot1 <- function(d, ...)  {

    ##as.vector(sapply(rval, function(r) sprintf("%s (prob = %.2f)", attr(r, "main"), attr(r, "prob"))))

    ## get group index
    j <- unique(d$group)

    plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab[j], ylab = ylab[j], main = main[j], 
      lwd = lwd, col = col, xaxs = "i", yaxs = "i", ...)
    box()

    polygon(
      na.omit(c(
        ifelse(extend_left, 0, NA), 
        d$mean_pr, 
        ifelse(extend_right, 1, NA), 
        rev(d$mean_pr), 
        ifelse(extend_left, 0, NA)
      )), 
      na.omit(c(
        ifelse(extend_left, 0, NA), 
        d$lowerlimit, 
        ifelse(extend_right, 1, NA), 
        rev(d$upperlimit), 
        ifelse(extend_left, 0, NA)
      )), 
      col = fill, border = NA
    )

    if(abline) abline(0, 1, lty = 2)

    lines(obs_rf ~ mean_pr, d, type = type, lwd = lwd, col = col, ...)

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
