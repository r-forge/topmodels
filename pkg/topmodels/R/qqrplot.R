## (Randomized) Q-Q residuals plot
##
## - Oabserved y in-sample or out-of-sample (n x 1)
## - Predicted probabilities F_y(y - eps) and F_y(y) (n x 2)
## - Two columns can be essentially equal -> continuous
##   or different -> (partially) discrete
## - Potentially transform uniform scale to different
##   distribution (default: Gaussian, via qnorm()).
##
## - Plot ordered empirical quantile residuals against
##   theoretical quantiles (from same distribution)
## - FIXME: To deal with point masses, draw either multiple random
##   draws (TODO: enable alpha blending by default) or
##   shade quantiles (FIXME: current implementation not correct)

## Functions:
## - qqrplot() generic plus default method
## - Return object of class "qqrplot" that is plotted by default
## - But has plot=FALSE so that suitable methods can be added afterwards
## - At least methods: plot(), autoplot()


## - qqrplot.default + plot.qqrplot + autoplot.qqrplot
##   trafo = qnorm vs. NULL vs. ...

qqrplot <- function(object, ...) {
  UseMethod("qqrplot")
}


qqrplot.default <- function(object, 
                            newdata = NULL, 
                            trafo = qnorm, 
                            type = c("random", "quantile"),
                            nsim = 1L, 
                            prob = 0.5, 
                            range = FALSE, 
                            diag = TRUE,
                            col = "black", 
                            fill = "lightgray", 
                            xlim = NULL, 
                            ylim = NULL,
                            main = "Q-Q residuals plot", 
                            xlab = "Theoretical quantiles", 
                            ylab = "Quantile residuals",
                            ...) {
  ## compute quantile residuals
  qres <- qresiduals(object, newdata = newdata, trafo = trafo, type = type, nsim = nsim, prob = prob)
  if(is.null(dim(qres))) qres <- matrix(qres, ncol = 1L)

  ## corresponding quantiles on the transformed scale (default: normal)
  if(is.null(trafo)) trafo <- identity
  q2q <- function(y) trafo(ppoints(length(y)))[order(order(y))]
  qthe <- apply(qres, 2L, q2q)
    
  ## default plotting ranges
  if(is.null(xlim)) xlim <- range(qthe)
  if(is.null(ylim)) ylim <- range(qres)

  ## set up coordinates
  plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)

  ## polygon for range:
  ## FIXME currently something goes wrong here - either in qresiduals() or here
  if(!identical(range, FALSE)) {
    if(isTRUE(range)) range <- c(0.01, 0.99)
    rg <- qresiduals(object, newdata = newdata, type = "quantile", prob = range)
    y1 <- sort(rg[,1])
    y2 <- sort(rg[,2])
    x <- c(q2q(y1), rev(q2q(y2)))
    y <- c(y1, rev(y2))
    y[!is.finite(y)] <- 100 * sign(y[!is.finite(y)])
    x[!is.finite(x)] <- 100 * sign(x[!is.finite(x)])
    polygon(x, y, col = fill, border = fill)
    box()
  }

  ## add Q-Q plot(s)
  for(i in 1L:ncol(qres)) points(qthe[, i], qres[, i], col = col, ...)
  
  ## reference diagonal
  if(!identical(diag, FALSE)) {
    if(isTRUE(diag)) diag <- "black"
    abline(0, 1, col = diag, lty = 2)
  }
  
  ## return coordinates invisibly
  invisible(list(theoretical = qthe, residuals = qres))
}

## actual drawing
c.qqrplot <- rbind.qqrplot <- function(x, ...) {
  NULL
}

## actual drawing
plot.qqrplot <- function(x, ...) {
  NULL
}

## ggplot2 interface
autoplot.qqrplot <- function(object, ...) {
  NULL
}

