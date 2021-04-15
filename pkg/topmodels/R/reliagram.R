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

  ### data + thresholds
  #y <- if(is.null(y)) newresponse(object, newdata = newdata) else as.numeric(y)
  #thresholds <- if(is.null(thresholds)) median(y, na.rm = TRUE) else as.numeric(thresholds)

  ### predicted probabilities
  #pred <- procast(object, newdata = newdata, type = "probability", at = matrix(thresholds, nrow = 1L),
  #  drop = TRUE)

  ## compute quantities to be plotted
  #rval <- data.frame("<cut prob>", "<aggregate y by prob and thresh>")

  ## get/prepare observations
  #obs <- y <= thresholds  #FIXME: (ML) Dirty hack to tryout function
  obs <- y

  ## define convenience variables
  N <- NROW(obs)

  ## use function brier from verification package to compute 
  ## observed relative frequencies (obar.i) and bin-centers (y.i)
  temp   <- brier(obs = obs, pred = pred, thresholds = breaks, bins = bins)

  minimumIndex <- which(temp$prob.y * N >= minimum)  # FIXME: (ML) What is that for?

  ## mean values of pred in the intervals
  # y.i <- temp$y.i  # FIXME: (ML) Why not used instead of for loop?
  y.i <- NULL
  for (i in 1:(length(breaks) - 1)) {
    y.i <- c(y.i, mean(pred[pred >= breaks[i] & pred < breaks[i + 1]]))
  }

  y.i <- y.i[minimumIndex]
  obar.i <- temp$obar.i[minimumIndex]
  obar.i[is.na(obar.i)] <- 0  ## FIXME: (ML) Include by me as I got NAs, how could that happen? miminimIndex takes care of that! 
 
  ## FIXME: (ML) My version to not use `brier()` - works only with `bins = TRUE`
  #obar.i <- as.numeric(
  #  aggregate(obs, by = list(cut(pred, breaks, include.lowest = TRUE)), mean)[, 2]
  #)
  #y.i <- na.omit(as.numeric(
  #  filter(breaks, filter = c(0.5, 0.5))
  #))

  ### FIXME: (ML) Included by Reto, what is it for?
  #meanpred.i <- as.numeric(
  #  aggregate(pred, by = list(cut(pred, breaks, include.lowest = TRUE)), mean, drop = FALSE)[, 2]
  #)
  #meanpred.i[is.na(meanpred.i)] <- 0 ## FIXME: (ML) I got NAs, why can that happen?

  ## consistency resampling from Broecker (2007) 
  obar.i.boot <- NULL
  for (i in 1:nboot) {
    ## take bootstrap sample from predictions (surrogate forecasts)
    xhat <- sample(pred, replace = TRUE)

    ## surrogate observations that are reliable by construction
    yhat <- runif(N) < xhat

    ### compute observed relative frequencies of the surrogate observations ## can lead to Error
    temp <- brier(obs = yhat, pred = xhat, thresholds = breaks, bins = bins)
    temp_obar.i <- temp$obar.i
    temp_obar.i[is.na(temp_obar.i)] <- 0  ## FIXME: (ML) Include by me as I got NAs, how could that happen? 
    obar.i.boot <- rbind(obar.i.boot, temp_obar.i)
  }

  ## compute lower and upper limits for reliable forecasts
  lowerlimit <- apply(obar.i.boot, 2, quantile, prob = prob[1], na.rm = TRUE)[minimumIndex]
  upperlimit <- apply(obar.i.boot, 2, quantile, prob = prob[2], na.rm = TRUE)[minimumIndex]

  ## auxiliary variable for polygon function (shading)
  shadedarea <- na.omit(cbind(c(y.i, rev(y.i)), c(upperlimit, rev(lowerlimit))))

  ## return value
  rval <- list(
    diagram = data.frame(
      obar.i, 
      y.i, 
      lowerlimit, 
      upperlimit
    ), 
    shading = shadedarea, 
    bs = mean((obs - pred)^2), 
    prob.y = N * temp$prob.y
  )

  ## attributes for graphical display
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  class(rval) <- c("reliagram", "list")  ## FIXME: (ML) Change to data.frame

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
                           ...) {

  plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, lwd = lwd, col = col, ...)

  if(!identical(shaded, FALSE)) {
    if(identical(shaded, TRUE)) shaded <- rgb(col2rgb(col)[1]/255, col2rgb(col)[2]/255, col2rgb(col)[3]/255, alpha = 0.1)
    with(object, polygon(shading[,1], shading[,2], col = shaded, border = NA))
  }

  if(!identical(confint, FALSE)) {
    if(identical(confint, TRUE)) confint <- "darkgray"
    lines(lowerlimit ~ y.i, data = object$diagram, col = confint)
    lines(upperlimit ~ y.i, data = object$diagram, col = confint)
  }

  if(abline) abline(0, 1, lty = 3)
  lines(obar.i ~ y.i, object$diagram, type = type, lwd = lwd, col = col, ...)
  if(print.bs) text(0,1, paste0("BS = ", signif(object$bs, 4)), adj = c(0,1))
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
