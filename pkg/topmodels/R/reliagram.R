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
                              minimum = 10, 
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


  ## data + thresholds
  y <- if(is.null(y)) newresponse(object, newdata = newdata) else as.numeric(y)
  thresholds <- if(is.null(thresholds)) quantile(y, probs = probs, na.rm = TRUE) 
  thresholds <- as.numeric(thresholds)

  ## predicted probabilities
  pred <- procast(object, newdata = newdata, type = "probability", at = matrix(thresholds, nrow = 1L),
    drop = TRUE)

  ## TODO: (Z) compute quantities to be plotted
  #rval <- data.frame("<cut prob>", "<aggregate y by prob and thresh>")

  ## get and prepare observations
  y <- y <= thresholds  #FIXME: (ML) Dirty hack to tryout function

  ## define convenience variables
  N <- NROW(y)

  ## compute number of prediction and idx for minimum number of prediction per probability subset
  n_pred <- aggregate(
    pred, 
    by = list(prob = cut(pred, breaks, include.lowest = TRUE)), 
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
      y, 
      by = list(prob = cut(pred, breaks, include.lowest = TRUE)), 
      FUN = mean, 
      drop = FALSE
      )[, "x"]
  )
  obs_rf <- obs_rf[min_idx]

  ## compute mean predicted probability (mean_pr)
  mean_pr <- as.numeric(
    aggregate(
      pred, 
      by = list(prob = cut(pred, breaks, include.lowest = TRUE)), 
      FUN = mean, 
      drop = FALSE
    )[, "x"]
  )
  mean_pr <- mean_pr[min_idx]

  ## consistency resampling from Broecker (2007) 
  obs_rf_boot <- vector("list", length = N)
  for (i in 1:confint_nboot) {
    ## take bootstrap sample from predictions (surrogate forecasts)
    xhat <- sample(pred, replace = TRUE)

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
  rval <- data.frame(
    obs_rf, 
    mean_pr, 
    lowerlimit, 
    upperlimit
  )

  ## attributes for graphical display
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "bs") <- mean(y - pred)^2
  class(rval) <- c("reliagram", "data.frame")

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
                           lwd = 2, 
                           type = "b", 
                           print_info = FALSE, 
                           col = "black", 
                           fill = adjustcolor(1, alpha.f = 0.1),
                           extend_left = TRUE, 
                           extend_right = TRUE, 
                           ...) {

  plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, lwd = lwd, col = col, 
    xaxs = "i", yaxs = "i", ...)
  box()

  polygon(
    na.omit(c(
      ifelse(extend_left, 0, NA), 
      x$mean_pr, 
      ifelse(extend_right, 1, NA), 
      rev(x$mean_pr), 
      ifelse(extend_left, 0, NA)
    )), 
    na.omit(c(
      ifelse(extend_left, 0, NA), 
      x$lowerlimit, 
      ifelse(extend_right, 1, NA), 
      rev(x$upperlimit), 
      ifelse(extend_left, 0, NA)
    )), 
    col = fill, border = NA
  )

  if(abline) abline(0, 1, lty = 2)

  lines(obs_rf ~ mean_pr, x, type = type, lwd = lwd, col = col, ...)

  if(print_info) text(0, 1, paste0("BS = ", signif(attr(x, "bs"), 4)), adj = c(0,1))
}


lines.reliagram <- function(x, ...) {
  NULL
}

c.reliagram <- function(x, ...) {
  NULL
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
