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
    # TODO: (ML) What is the default fun for? 
    p <- qresiduals.default(object, newdata = newdata, trafo = NULL, nsim = nsim, 
      mass_redist = mass_redist)
  }

  ## breaks
  if(is.null(breaks)) breaks <- c(4, 10, 20, 25)[cut(NROW(p), c(0, 50, 5000, 1000000, Inf))]
  if(length(breaks) == 1L) breaks <- seq(xlim[1L], xlim[2L], length.out = breaks + 1L)

  ## labels
  if(is.null(main)) main <- deparse(substitute(object))

  ## collect everything as data.frame
  ## TODO: (ML) Should it really be prepared as a `data.frame` and 
  ## return name of hist() should be renamed
  ## TODO: (ML) Maybe get rid of `hist()`
  tmp_hist <- hist(p, breaks = breaks, plot = FALSE, ...)
  rval <- data.frame(x = tmp_hist$mids, xleft = tmp_hist$breaks[-length(tmp_hist$breaks)], 
                     xright = tmp_hist$breaks[-1L], y = tmp_hist$density) 
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  class(rval) <- c("pithist", "data.frame")

  ## also plot by default
  if (plot) {
    plot(rval, ...)
  }

  ## return invisibly
  invisible(rval)
}


c.pithist <- rbind.pithist <- function(...)
{
  ## list of pithists
  rval <- list(...)

  ## group sizes
  for(i in seq_along(rval)) {
    if(is.null(rval[[i]]$group)) rval[[i]]$group <- 1L
  } 
  n <- lapply(rval, function(r) table(r$group))

  ## labels
  xlab <- unlist(lapply(rval, function(r) attr(r, "xlab")))
  ylab <- unlist(lapply(rval, function(r) attr(r, "ylab")))
  nam <- names(rval)
  main <- if(is.null(nam)) {
    as.vector(sapply(rval, function(r) attr(r, "main")))
  } else {
    make.unique(rep.int(nam, sapply(n, length)))
  } 
  n <- unlist(n)

  ## combine and return
  rval <- do.call("rbind.data.frame", rval)
  rval$group <- if(length(n) < 2L) NULL else rep.int(seq_along(n), n)
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  class(rval) <- c("pithist", "data.frame")
  return(rval)
}


plot.pithist <- function(x,
  xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, main = NULL,
  border = "black", fill = "lightgray", col = "#B61A51",
  lwd = 2, pch = 19, lty = 1, axes = TRUE, ...)
{

  ## handling of groups
  if(is.null(x$group)) x$group <- 1L
  n <- max(x$group) 
  
  ## annotation
  if(is.null(xlab)) xlab <- TRUE
  if(is.null(ylab)) ylab <- TRUE
  if(is.null(main)) main <- TRUE 
  xlab <- rep(xlab, length.out = n)
  ylab <- rep(ylab, length.out = n)
  main <- rep(main, length.out = n)
  if(is.logical(xlab)) xlab <- ifelse(xlab, attr(x, "xlab"), "")
  if(is.logical(ylab)) ylab <- ifelse(ylab, attr(x, "ylab"), "")
  if(is.logical(main)) main <- ifelse(main, attr(x, "main"), "")

  ## plotting function
  pithist1 <- function(d, ...) {
    ## rect elements
    xleft <- d$xleft
    xright <- d$xright
    y <- d$y 
    j <- unique(d$group)

    ## defaults
    if(is.null(xlim)) xlim <- range(c(xleft, xright))
    if(is.null(ylim)) ylim <- range(c(0, y))

    ## draw pithist 
    plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
      xlab = xlab[j], ylab = ylab[j], main = main[j], axes = FALSE, ...)
    if(axes) {
      axis(1)
      axis(2)
    }
    rect(xleft, 0, xright, y, border = border, col = fill)
    abline(h = 1, col = col, lty = lty, lwd = lwd)
   }

   ## draw plots
   if(n > 1L) par(mfrow = n2mfrow(n))
   for(i in 1L:n) pithist1(x[x$group == i, ], ...)
}


lines.pithist <- function(x, ...) {
  NULL
}


## ggplot2 interface
autoplot.pithist <- function(object, ...) {
  NULL
}
