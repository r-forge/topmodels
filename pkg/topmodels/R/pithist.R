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
##   or distribute evenly across relevant intervalst TODO: (ML) hot to do the latter?
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


pithist.default <- function(object, 
                            newdata = NULL, 
                            plot = TRUE, 
                            style = c("histogram", "lines"),
                            type = c("random", "proportional"), 
                            nsim = 1L, 
                            freq = FALSE, 
                            breaks = NULL, 
                            xlim = c(0, 1), 
                            ylim = NULL, 
                            main = NULL, 
                            xlab = "PIT", 
                            ylab = if(freq) "Frequency" else "Density",
                            confint = TRUE, 
                            confint_level = 0.95, 
                            confint_type = c("exact", "approximation"), 
                            ...) {

  ## match arguments
  style <- match.arg(style, c("histogram", "lines"))
  confint_type <- match.arg(confint_type, c("exact", "approximation"))

  ## either compute proportion exactly (to do...) or approximate by simulation
  type <- match.arg(type, c("random", "proportional"))
  if(type == "proportional") {
    stop("not yet implemented")
  } else {
    # TODO: (ML) What is the default fun for? 
    p <- qresiduals.default(object, newdata = newdata, trafo = NULL, type = "random", nsim = nsim)
  }

  ## breaks
  if(is.null(breaks)) breaks <- c(4, 10, 20, 25)[cut(NROW(p), c(0, 50, 5000, 1000000, Inf))]
  if(length(breaks) == 1L) breaks <- seq(xlim[1L], xlim[2L], length.out = breaks + 1L)

  ## ci interval
  if (confint_type == "exact") {
    ci <- get_confint(NROW(p), length(breaks) - 1, confint_level, freq)
  } else {
    ci <- get_confint_agresti(
      NROW(p) / (length(breaks) - 1), 
      NROW(p), 
      confint_level, 
      length(breaks) - 1, 
      freq)
  }

  ## perfect prediction
  pp <- ifelse(freq, NROW(p) / (length(breaks) - 1), 1)

  ## labels
  if(is.null(main)) main <- deparse(substitute(object))

  ## collect everything as data.frame
  ## TODO: (ML) Should it really be prepared as a `data.frame` and 
  ## return name of hist() should be renamed
  ## TODO: (ML) Maybe get rid of `hist()`
  tmp_hist <- hist(p, breaks = breaks, plot = FALSE)
  if (freq) {
    rval <- data.frame(x = tmp_hist$mids, 
                       xleft = tmp_hist$breaks[-length(tmp_hist$breaks)], 
                       xright = tmp_hist$breaks[-1L], 
                       counts = tmp_hist$counts, 
                       ci_lower = ci[1], 
                       ci_upper = ci[2], 
                       pp = pp) 
  } else {
    rval <- data.frame(x = tmp_hist$mids, 
                       xleft = tmp_hist$breaks[-length(tmp_hist$breaks)], 
                       xright = tmp_hist$breaks[-1L], 
                       density = tmp_hist$density, 
                       ci_lower = ci[1], 
                       ci_upper = ci[2], 
                       pp = pp) 
  }
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  class(rval) <- c("pithist", "data.frame")

  ## also plot by default
  if (plot) {
    plot(rval, style = style, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
      confint = confint, ...)
  }

  ## return invisibly
  invisible(rval)
}


c.pithist <- rbind.pithist <- function(...) {

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
                         style = c("histogram", "lines"), 
                         confint = TRUE,  
                         xlim = NULL, 
                         ylim = NULL, 
                         xlab = NULL, 
                         ylab = NULL, 
                         main = NULL,
                         border = "black", 
                         fill = adjustcolor(1, alpha.f = 0.2), 
                         col = c(1, 1),
                         lwd = 2, 
                         lty = c(1, 2), 
                         axes = TRUE, 
                         ...) {

  style <- match.arg(style, c("histogram", "lines"))

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

    if (length(col) < 2) col <- rep(col[1], 2)
    if (length(lty) < 2) lty <- rep(lty[1], 2)

    ## rect elements
    xleft <- d$xleft
    xright <- d$xright
    y <- if (!is.null(d$density)) d$density else d$counts
    j <- unique(d$group)
    ci_lower <- ifelse(confint, unique(d$ci_lower), NULL)
    ci_upper <- ifelse(confint, unique(d$ci_upper), NULL)
    pp <-  unique(d$pp)

    ## defaults
    if(is.null(xlim)) xlim <- range(c(xleft, xright))
    if(is.null(ylim)) ylim <- range(c(0, y, ci_lower, ci_upper))

    ## draw pithist 
    plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
      xlab = xlab[j], ylab = ylab[j], main = main[j], axes = FALSE, ...)
    if(axes) {
      axis(1)
      axis(2)
    }
    rect(xleft, 0, xright, y, border = border, col = fill)
    abline(h = pp, col = col[1], lty = lty[1], lwd = lwd)
    if(confint) abline(h = ci_lower, col = col[2], lty = lty[2], lwd = lwd)
    if(confint) abline(h = ci_upper, col = col[2], lty = lty[2], lwd = lwd)
  }

  ## plotting function
  pithist2 <- function(d, ...) {
    ## rect elements
    x <- c(d$xleft, d$xright[length(d$xright)])
    y <- if (!is.null(d$density)) d$density else d$counts
    y <- c(y, y[length(y)])
    j <- unique(d$group)
    ci_lower <- if (confint) unique(d$ci_lower) else NULL
    ci_upper <- if (confint) unique(d$ci_upper) else NULL
    pp <-  unique(d$pp)

    ## defaults
    if(is.null(xlim)) xlim <- range(x)
    if(is.null(ylim)) ylim <- range(c(0, y, ci_lower, ci_upper))
    ## draw pithist 
    plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
      xlab = xlab[j], ylab = ylab[j], main = main[j], axes = FALSE, ...)
    if(axes) {
      axis(1)
      axis(2)
    }

    if(confint) polygon(c(0, 1, 1, 0), c(ci_lower, ci_lower, ci_upper, ci_upper),
      col = fill, border = NA)
    #abline(h = pp, col = col, lty = lty[1], lwd = lwd)
    lines(y ~ x, type = "s", lwd = lwd, lty = lty[1], col = col[1])
  }

  ## draw plots
  if(n > 1L) par(mfrow = n2mfrow(n))
  for (i in 1L:n) {
    if (style == "histogram") {
      pithist1(x[x$group == i, ], ...) 
    } else {
      pithist2(x[x$group == i, ], ...)
    }
  }
}


lines.pithist <- function(x, 
                          confint = FALSE,  
                          fill = adjustcolor(1, alpha.f = 0.2), 
                          col = 1,
                          lwd = 2, 
                          lty = 1,  
                          ...) {

  ## handling of groups
  if(is.null(x$group)) x$group <- 1L
  n <- max(x$group) 

  ## Get right length of `col` and `lty`
  if (length(col) < n) col <- rep(col[1], n)
  if (length(lty) < n) lty <- rep(lty[1], n)
  if (length(lwd) < n) lwd <- rep(lwd[1], n)
  if (length(fill) < n) fill <- rep(fill[1], n)
  
  ## plotting function
  pithist2 <- function(d, col, lty, lwd, fill, ...) {
    ## rect elements
    x <- c(d$xleft, d$xright[length(d$xright)])
    y <- if (!is.null(d$density)) d$density else d$counts
    y <- c(y, y[length(y)])
    j <- unique(d$group)

    ci_lower <- if (confint) unique(d$ci_lower) else NULL
    ci_upper <- if (confint) unique(d$ci_upper) else NULL

    ## draw pithist 
    if(confint) polygon(c(0, 1, 1, 0), c(ci_lower, ci_lower, ci_upper, ci_upper),
      col = fill, border = NA)
    lines(y ~ x, type = "s", col = col, lty = lty, lwd = lwd, ...)
  }

  ## draw plots
  for (i in 1L:n) {
    pithist2(x[x$group == i, ], col = col[i], lty = lty[i], lwd = lwd[i], fill = fill[i], ...)
  }
}


## ggplot2 interface
autoplot.pithist <- function(object,
                             colour = c("black", "#B61A51"), 
                             fill = "darkgray", 
                             size = 1.2, 
                             ...) {

  ## determine grouping
  class(object) <- "data.frame"
  if(is.null(object$group)) object$group <- 1L
  n <- max(object$group)
  object$group <- factor(object$group, levels = 1L:n, 
    labels = make.names(attr(object, "main"), unique = TRUE))
  
  ## rectangles and fitted lines
  rval <- ggplot2::ggplot(object, ggplot2::aes_string(xmin = "xleft", xmax = "xright", ymin = 0, ymax = "y")) + 
    ggplot2::geom_rect(colour = colour[1L], fill = fill) + 
    ggplot2::geom_hline(yintercept = 1, colour = colour[2L], size = size)
  
  ## grouping (if any)
  if(n > 1L) rval <- rval + ggplot2::facet_grid(group ~ .)

  ## annotation
  rval <- rval + ggplot2::xlab(paste(unique(attr(object, "xlab")), collapse = "/")) +
    ggplot2::ylab(paste(unique(attr(object, "ylab")), collapse = "/"))

  ## return with annotation
  rval
}


## helper function to calculate CI employing `qbinom()`
get_confint <- function(n, bins, level, freq) {
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  rval <- qbinom(a, size = n, prob = 1 / bins)
  if(!freq) rval <- rval / (n / bins) 
  rval
}


## helper function to calculate an approximated CI according to Agresti & Coull (1998)
## doi=10.1080/00031305.1998.10480550
get_confint_agresti <- function(x, n, level, bins, freq) {
  rval <- add4ci(x, n, level)$conf.int * n
  if(!freq) rval <- rval / (n / bins) 
  rval
}


## copy of `add4ci` package from package `PropCIs` by Ralph Scherer (licensed under GPL-2/GPL-3)
add4ci <- function (x, n, conf.level) {
  ptilde = (x + 2)/(n + 4)
  z = abs(qnorm((1 - conf.level)/2))
  stderr = sqrt(ptilde * (1 - ptilde)/(n + 4))
  ul = ptilde + z * stderr
  ll = ptilde - z * stderr
  if (ll < 0) 
      ll = 0
  if (ul > 1) 
      ul = 1
  cint <- c(ll, ul)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint, estimate = ptilde)
  class(rval) <- "htest"
  return(rval)
}

