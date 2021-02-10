## PIT histogram
##
## - Observed y in-sample or out-of-sample (n x 1)
## - Predicted probabilities F_y(y - eps) and F_y(y) (n x 2)
## - Two columns can be essentially equal -> continuous
##   or different -> (partially) discrete
## - Breaks for predicted probabilities in [0, 1] (m x 1)
##
## - Cut probabilities at breaks -> (m-1) groups and draw histogram
## - TODO: In case of point masses either use a random draw
##   or distribute evenly across relevant intervalst
## - TODO: Random draws could be drawn by hist() (current solution)
##   but proportional distribution requires drawing rectangles by hand
## - TODO: add confidence interval as well.
## - TODO: Instead of shaded rectangles plus reference line and CI lines
##   support shaded CI plus step lines

## Functions:
## - pithist() generic plus default method
## - Return object of class "pithist" that is plotted by default
## - But has plot=FALSE so that suitable methods can be added afterwards
## - At least methods: plot(), autoplot(), lines(), possibly c(), +

pithist <- function(object, ...) {
  UseMethod("pithist")
}


pithist.default <- function(object, newdata = NULL, type = c("random", "proportional"), nsim = 1L,
  breaks = NULL, xlim = c(0, 1), ylim = NULL,
  xlab = "PIT", ylab = "Density", main = NULL,
  border = "black", fill = "lightgray", col = "#B61A51",
  lwd = 2, lty = 1, freq = FALSE, ...)
{
  ## either compute proportion exactly (to do...) or approximate by simulation
  type <- match.arg(type, c("random", "proportional"))
  if(type == "proportional") {
    stop("not yet implemented")
  } else {
    #p <- qresiduals.default(object, newdata = newdata, trafo = NULL, nsim = nsim)
    # TODO: (ML) What is the default fun for? Compare comment for `qresiduals.crch()'.
    p <- qresiduals(object, newdata = newdata, trafo = NULL, nsim = nsim, distribute_cens = "quantile")
  }

  ## breaks
  if(is.null(breaks)) breaks <- c(4, 10, 20, 25)[cut(NROW(p), c(0, 50, 5000, 1000000, Inf))]
  if(length(breaks) == 1L) breaks <- seq(xlim[1L], xlim[2L], length.out = breaks + 1L)

  ## labels
  if(is.null(main)) main <- deparse(substitute(object))

  ## histogram
  rval <- hist(p, breaks = breaks, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main,
    col = fill, border = border, freq = freq, ...)
  abline(h = 1, col = col, lty = lty, lwd = lwd)
  invisible(rval)
}

## actual drawing
plot.pithist <- function(x, ...) {
  NULL
}

lines.pithist <- function(x, ...) {
  NULL
}

## ggplot2 interface
autoplot.pithist <- function(object, ...) {
  NULL
}

