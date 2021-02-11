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
  mass_redist = c("quantile", "random", "none"),
  breaks = NULL, plot = TRUE, xlim = c(0, 1), ylim = NULL,
  xlab = "PIT", ylab = "Density", main = NULL,
  border = "black", fill = "lightgray", col = "#B61A51",
  lwd = 2, lty = 1, freq = FALSE, ...)
{

  ## point mass redistribution
  mass_redist <- match.arg(mass_redist, c("quantile", "random", "none"))

  ## either compute proportion exactly (to do...) or approximate by simulation
  type <- match.arg(type, c("random", "proportional"))
  if(type == "proportional") {
    stop("not yet implemented")
  } else {
    #p <- qresiduals.default(object, newdata = newdata, trafo = NULL, nsim = nsim)
    # TODO: (ML) What is the default fun for? Compare comment for `qresiduals.crch()'.
    p <- qresiduals(object, newdata = newdata, trafo = NULL, nsim = nsim, 
      mass_redist = mass_redist)
  }

  ## breaks
  if(is.null(breaks)) breaks <- c(4, 10, 20, 25)[cut(NROW(p), c(0, 50, 5000, 1000000, Inf))]
  if(length(breaks) == 1L) breaks <- seq(xlim[1L], xlim[2L], length.out = breaks + 1L)

  ## labels
  if(is.null(main)) main <- deparse(substitute(object))

  ## collect everything as list $ TODO: (ML) Should it be prepared as `data.frame' to be consistent?
  rval <- hist(p, breaks = breaks, plot = FALSE, ...)
  class(rval) <- c("pithist", "list")

  ## also plot by default
  if (plot) {
    plot(rval, freq = freq, fill = fill, col = col, lwd = lwd, border = border, main = main, 
      xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
  }

  ## return invisibly
  invisible(rval)
}

## actual drawing
# TODO: (ML) currently copy of `graphics:::plot.histogram()'
plot.pithist <- function(x,  
  xlab = "PIT", ylab = "Density", main = NULL,
  border = "black", fill = "lightgray", col = "#B61A51",
  lwd = 2, lty = 1, freq = FALSE, 
  
  density = NULL, angle = 45, 
  sub = NULL, xlim = range(x$breaks), 
  ylim = NULL, axes = TRUE, labels = FALSE, 
  add = FALSE, ann = TRUE, ...) {

  equidist <- if (is.logical(x$equidist)) 
      x$equidist
  else {
      h <- diff(x$breaks)
      diff(range(h)) < 1e-07 * mean(h)
  }
  if (freq && !equidist) 
      warning("the AREAS in the plot are wrong -- rather use 'freq = FALSE'")
  y <- if (freq) 
      x$counts
  else x$density
  nB <- length(x$breaks)
  if (is.null(y) || 0L == nB) 
      stop("'x' is wrongly structured")
  dev.hold()
  on.exit(dev.flush())
  if (!add) {
      if (is.null(ylim)) 
          ylim <- range(y, 0)
      if (missing(ylab)) 
          ylab <- if (!freq) 
              "Density"
          else "Frequency"
      plot.new()
      plot.window(xlim, ylim, "", ...)
      if (ann) 
          title(main = main, sub = sub, xlab = xlab, ylab = ylab, 
              ...)
      if (axes) {
          axis(1, ...)
          axis(2, ...)
      }
  }
  rect(x$breaks[-nB], 0, x$breaks[-1L], y, col = fill, border = border, 
      angle = angle, density = density, lty = lty)
  if ((logl <- is.logical(labels) && labels) || is.character(labels)) 
      text(x$mids, y, labels = if (logl) {
          if (freq) 
              x$counts
          else round(x$density, 3)
      }
      else labels, adj = c(0.5, -0.5))
  abline(h = 1, col = col, lty = lty, lwd = lwd)
}

lines.pithist <- function(x, ...) {
  NULL
}

## ggplot2 interface
autoplot.pithist <- function(object, ...) {
  NULL
}
