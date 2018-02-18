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

reliagram.default <- function(object, newdata = NULL,
  breaks = seq(0, 1, by = 0.1), thresholds = NULL, y = NULL, plot = TRUE,
  xlab = NULL, ylab = NULL, ...)
{
  ## data + thresholds
  y <- if(is.null(y)) newresponse(object, newdata = newdata) else as.numeric(y)
  thresholds <- if(is.null(thresholds)) median(y, na.rm = TRUE) else as.numeric(thresholds)

  ## predicted probabilities
  prob <- p4(object, newdata = newdata, type = "probability", at = matrix(thresholds, nrow = 1L))

  ## compute quantities to be plotted
  rval <- data.frame("<cut prob>", "<aggregate y by prob and thresh>")

  ## attributes for graphical display
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  ## ...
  class(rval) <- c("reliagram", "data.frame")

  ## plot and/or return
  if(plot) {
    plot(rval, ...)
    invisible(rval)
  } else {
    return(rval)
  }
}

## methods might only add custom y and thresholds
reliagram.crch <- function(object, newdata = NULL,
  breaks = seq(0, 1, by = 0.1), thresholds = NULL, ...)
{
  y <- if((missing(newdata) || is.null(newdata)) && !is.null(object$y)) object$y else newresponse(object, newdata = newdata)
  if(is.null(thresholds)) thresholds <- object$left
  reliagram(object, breaks = breaks, thresholds = thresholds, y = y, ...)
}


## actual drawing
plot.reliagram <- function(x, ...) {
  NULL
}

lines.reliagram <- function(x, ...) {
  NULL
}

## ggplot2 interface
autoplot.reliagram <- function(object, ...) {
  NULL
}

