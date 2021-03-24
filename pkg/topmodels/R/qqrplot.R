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
                            plot = TRUE,
                            range = FALSE, 
                            diag = TRUE,
                            main = "Q-Q residuals plot", 
                            xlab = "Theoretical quantiles", 
                            ylab = "Quantile residuals",
                            ...) {

  ## compute quantile residuals
  qres <- qresiduals(object, newdata = newdata, trafo = trafo, type = type, nsim = nsim, prob = prob)
  if(is.null(dim(qres))) qres <- matrix(qres, ncol = 1L)

  ## corresponding quantiles on the transformed scale (default: normal)
  if(is.null(trafo)) trafo <- identity
  q2q <- function(y) trafo(ppoints(length(y)))[order(order(y))]  ## FIXME: (ML) Why 2 times `order()`
  qthe <- apply(qres, 2L, q2q)

  ## collect everything as data.frame
  rval <- data.frame(theoretical = qthe, residuals = qres)
  names(rval) <- gsub("(\\.r|\\.q)", "", names(rval))

  ## polygon for range:
  ## FIXME: (Z) currently something goes wrong here - either in qresiduals() or here
  ## TODO: (ML) Does it really go wrong, or just previously applied to object w/o randomized residuals?
  ## TODO: (ML) Why first if sentence? 
  ## TODO: (ML) Is the theoretical lower and upper not always the same?
  ## TODO: (ML) Instead of range, confidence interval possible?
  if(!identical(range, FALSE)) {
    if(isTRUE(range)) range_vec <- c(0.01, 0.99)
    rg <- qresiduals(object, newdata = newdata, type = "quantile", prob = range_vec)
    
    rval$residuals_range_lower <- rg[, 1]
    rval$residuals_range_upper <- rg[, 2]
    rval$theoretical_range_lower <- q2q(rg[, 1])
    rval$theoretical_range_upper <- q2q(rg[, 2])
  }

  attr(rval, "main") <- main
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  class(rval) <- c("qqrplot", "data.frame")

  ## also plot by default
  if(plot){ 
    plot(rval, range = range, diag = diag, ...)
  }
  
  ## return coordinates invisibly
  invisible(rval)
}

## Combine several qqrplots
c.qqrplot <- rbind.qqrplot <- function(x, ...) {
  NULL
}

## actual drawing
plot.qqrplot <- function(x, 
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

  ## polygon for range:
  ## FIXME: (Z) Currently something goes wrong here - either in qresiduals() or here
  ## TODO: (ML) Does it really go wrong, or just not randomized residuals?
  ## TODO: (ML) Why first if sentence? 
  ## TODO: (ML) With increased difference in qresiduals, now check for equality does not work any more.
  if(
    !identical(range, FALSE) && 
    isTRUE(range) && 
    sum(grepl("range_lower|range_upper", names(x))) == 4 && 
    !(isTRUE(all.equal(x$theoretical_range_upper, x$theoretical_range_lower, tol =  .Machine$double.eps^0.4)) & 
      isTRUE(all.equal(x$residuals_range_upper, x$residuals_range_lower, tol =  .Machine$double.eps^0.4)))
    ) {

    ## default plotting ranges (as.matrix to get `finite = TRUE` working)
    if(is.null(xlim)) xlim <- range(as.matrix(x[grepl("theoretical", names(x))]), finite = TRUE)
    if(is.null(ylim)) ylim <- range(as.matrix(x[grepl("residuals", names(x))]), finite = TRUE)

    ## set up coordinates
    plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)
    ## plot polygon
    x_pol <- c(sort(x$theoretical_range_lower), sort(x$theoretical_range_upper, decreasing = TRUE))
    y_pol <- c(sort(x$residuals_range_lower), sort(x$residuals_range_upper, decreasing = TRUE))
    x_pol[!is.finite(x_pol)] <- 100 * sign(x_pol[!is.finite(x_pol)])
    y_pol[!is.finite(y_pol)] <- 100 * sign(y_pol[!is.finite(y_pol)])
    polygon(x_pol, y_pol, col = fill, border = fill)
    box()
  } else {
 
    ## set range to null 
    x$theoretical_range_lower <- x$theoretical_range_upper <- NULL
    x$residuals_range_lower <- x$residuals_range_upper <- NULL

    ## default plotting ranges
    if(is.null(xlim)) xlim <- range(x[grepl("theoretical", names(x))])
    if(is.null(ylim)) ylim <- range(x[grepl("residuals", names(x))])

    ## set up coordinates
    plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)

  }

  ## add Q-Q plot(s)
  for (i in 1L:ncol(x[grepl("residuals_[0-9]", names(x))])) {
    points(
      x[grepl("theoretical", names(x))][, i], 
      x[grepl("residuals", names(x))][, i], 
      col = col, ...
    )
  }
  
  ## reference diagonal
  if(!identical(diag, FALSE)) {
    if(isTRUE(diag)) diag <- "black"
    abline(0, 1, col = diag, lty = 2)
  }
  
}

## ggplot2 interface
autoplot.qqrplot <- function(object, ...) {
  NULL
}

