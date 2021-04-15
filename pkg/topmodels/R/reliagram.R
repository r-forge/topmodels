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
                              breaks = seq(0, 1, by = 0.1), 
                              thresholds = NULL, 
                              y = NULL, 
                              plot = TRUE,
                              xlab = NULL, 
                              ylab = NULL, 
                              bins = TRUE, 
                              prob = c(0.05, 0.95), 
                              nboot = 250, 
                              add = FALSE, 
                              shaded = FALSE, 
                              print.bs = FALSE, 
                              confint = TRUE, 
                              minimum = 10, 
                              ...) {

  require("verification")

  ## convert prob into right format
  if(length(prob) == 1) prob <- c((0.5 - prob/2), (0.5 + prob/2))
  if(max(prob) >= 1 | min(prob) <= 0) stop("confidence interval probabilities outside [0,1]")

  ## data + thresholds
  y <- if(is.null(y)) newresponse(object, newdata = newdata) else as.numeric(y)
  thresholds <- if(is.null(thresholds)) median(y, na.rm = TRUE) else as.numeric(thresholds)

  ## predicted probabilities
  pred <- procast(object, newdata = newdata, type = "probability", at = matrix(thresholds, nrow = 1L),
    drop = TRUE)

  ## TODO: (Z) compute quantities to be plotted
  #rval <- data.frame("<cut prob>", "<aggregate y by prob and thresh>")

  ## get/prepare observations
  #obs <- y <= thresholds  #FIXME: (ML) Dirty hack to tryout function
  obs <- y

  ## define convenience variables
  N <- NROW(obs)

  ## compute observed relative frequencies (obar.i) and bin-centers (y.i)
  temp <- brier(obs = obs, pred = pred, thresholds = breaks, bins = bins)

  ## FIXME: (ML) Do we want to use `brier()`
  ## * Returns returns midpoints of intervals
  ## * Also removes NAs in case no observations are within breaks
  minimumIndex <- which(temp$prob.y * N >= minimum)
  y.i <- temp$y.i[minimumIndex]  
  obar.i <- temp$obar.i[minimumIndex]

  ## Jakob's style for y.i
  #y.i <- NULL
  #for (i in 1:(length(breaks) - 1)) {
  #  y.i <- c(y.i, mean(pred[pred >= breaks[i] & pred < breaks[i + 1]]))
  #}
  ## Moritz' style for y.i -> problem here, when no predictions within intervals
  #y.i <- as.numeric(aggregate(pred, by = list(cut(pred, breaks, include.lowest = TRUE)), mean)[,2])

  ## Moritz' style for y.i
  #obar.i <- as.numeric(
  #  aggregate(obs, by = list(cut(pred, breaks, include.lowest = TRUE)), mean)[, 2]
  #)


  ## consistency resampling from Broecker (2007) 
  obar.i.boot <- NULL
  for (i in 1:nboot) {
    ## take bootstrap sample from predictions (surrogate forecasts)
    xhat <- sample(pred, replace = TRUE)

    ## surrogate observations that are reliable by construction
    yhat <- runif(N) < xhat

    ### compute observed relative frequencies of the surrogate observations
    temp <- brier(obs = yhat, pred = xhat, thresholds = breaks, bins = bins)
    temp_obar.i <- temp$obar.i
    temp_obar.i[is.na(temp_obar.i)] <- 0  # FIXME: (ML) Included as I got NAs (see above) 
    obar.i.boot <- rbind(obar.i.boot, temp_obar.i)
  }

  ## compute lower and upper limits for reliable forecasts
  lowerlimit <- apply(obar.i.boot, 2, quantile, prob = prob[1], na.rm = TRUE)[minimumIndex]
  upperlimit <- apply(obar.i.boot, 2, quantile, prob = prob[2], na.rm = TRUE)[minimumIndex]

  ## return value
  rval <- data.frame(
    obar.i, 
    y.i, 
    lowerlimit, 
    upperlimit
    ## FIXME: (ML) Still needed? If yes subset with minimumIndex and use not temp of resampling
    #prob.y = N * temp$prob.y  
  )

  ## attributes for graphical display
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "bs") <- mean(obs - pred)^2
  class(rval) <- c("reliagram", "data.frame")

  ## plot by default
  if (plot) {
    plot(rval, shaded = shaded, confint = confint, print.bs = print.bs, ...)
  }

  ## return invisibly
  invisible(rval)
}


## actual drawing
plot.reliagram <- function(object, 
                           shaded = !confint, 
                           confint = FALSE, 
                           abline = TRUE,
                           xlim = c(0, 1), 
                           ylim = c(0, 1), 
                           xlab = "Forecast probability", 
                           ylab = "Observed relative frequency",
                           lwd = 2, 
                           type = "b", 
                           print.bs = FALSE, 
                           col = "black", 
                           fill = "lightgray",
                           extend = TRUE, 
                           ...) {

  ## auxiliary variable for polygon function (shading)
  shading <- na.omit(cbind(c(object$y.i, rev(object$y.i)), c(object$upperlimit, rev(object$lowerlimit))))

  plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, lwd = lwd, col = col, 
    xaxs = "i", yaxs = "i", ...)

  if(!identical(shaded, FALSE)) {
    if(identical(shaded, TRUE)) shaded <- rgb(col2rgb(col)[1]/255, col2rgb(col)[2]/255, col2rgb(col)[3]/255, alpha = 0.1)
    with(object, polygon(shading[,1], shading[,2], col = shaded, border = NA))
  }

  if(confint) {
    if(extend) { 
      polygon(
        c(0, object$y.i, 1, rev(object$y.i), 0), 
        c(0, object$lowerlimit, 1, rev(object$upperlimit), 0), 
        col = fill, border = fill
      )
    } else {
      polygon(
        c(object$y.i, rev(object$y.i)), 
        c(object$lowerlimit, rev(object$upperlimit)), 
        col = fill, border = fill
      )
    }
    box()
  }

  if(abline) abline(0, 1, lty = 2)
  lines(obar.i ~ y.i, object, type = type, lwd = lwd, col = col, ...)
  if(print.bs) text(0,1, paste0("BS = ", signif(attr(object, "bs"), 4)), adj = c(0,1))
}


lines.reliagram <- function(x, ...) {
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
