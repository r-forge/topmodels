# -------------------------------------------------------------------
# Programming outline: PIT histogram
# -------------------------------------------------------------------
# - Observed y in-sample or out-of-sample (n x 1)
# - Predicted probabilities F_y(y - eps) and F_y(y) (n x 2)
# - Two columns can be essentially equal -> continuous
#   or different -> (partially) discrete
# - Breaks for predicted probabilities in [0, 1] (m x 1)
#
# - Cut probabilities at breaks -> (m-1) groups and draw histogram
# - In case of point masses either use a random draw
#   or distribute evenly across relevant intervals
# - Random draws could be drawn by hist() (current solution)
#   but proportional distribution requires drawing rectangles by hand
# - Add confidence interval as well.
# - Instead of shaded rectangles plus reference line and CI lines
#   support shaded CI plus step lines

# Functions:
# - pithist() generic plus default method
# - Return object of class "pithist" that is plotted by default
# - But has plot=FALSE so that suitable methods can be added afterwards
# - At least methods: plot(), autoplot(), lines(), possibly c(), +
# -------------------------------------------------------------------


#' PIT Histograms for Assessing Goodness of Fit of Probability Models
#' 
#' PIT histograms graphically compare empirical probabilities from fitted models
#' with a uniform distribution. If \code{plot = TRUE}, the resulting object of
#' class \code{"pithist"} is plotted by \code{\link{plot.pithist}} or
#' \code{\link{autoplot.pithist}} before it is returned, depending on whether the
#' package \code{ggplot2} is loaded.
#' 
#' PIT histograms graphically evaluate the probability integral transform (PIT),
#' i.e., the value that the predictive CDF attains at the observation, with a
#' uniform distribution. For a well calibrated model fit, the observation will be
#' drawn from the predictive distribution and the PIT will have a standard uniform
#' distribution. For computation, \code{\link{pithist}} leverages the function
#' \code{\link{qresiduals}} employing the \code{\link{procast}} generic and then
#' essentially draws a \code{\link[graphics]{hist}}. 
#' 
#' In case of discrete distributions the PIT can be either drawn randomly from the
#' corresponding interval or distributed proportionally in the histogram, whereby
#' the latter is not yet supported.
#'
#' In addition to the \code{plot} and \code{\link[ggplot2]{autoplot}} method for
#' pithist objects, it is also possible to combine two (or more) PIT histograms by
#' \code{c}/\code{rbind}, which creates a set of PIT histograms that can then be
#' plotted in one go. 
#' 
#' @aliases pithist pithist.default c.pithist rbind.pithist
#' @param object an object from which probability integral transforms can be
#' extracted using the generic function \code{\link{procast}}.
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param plot Should the \code{plot} or \code{autoplot} method be called to
#' draw the computed extended reliability diagram? Either set \code{plot}
#' expicitly to \code{"base"} vs. \code{"ggplot2"} to choose the type of plot, or for a
#' logical \code{plot} argument it's chosen conditional if the package
#' \code{ggplot2} is loaded.
#' @param class Should the invisible return value be either a \code{data.frame}
#' or a \code{tibble}. Either set \code{class} expicitly to \code{"data.frame"} vs.
#' \code{"tibble"}, or for \code{NULL} it's chosen automatically conditional if the package
#' \code{tibble} is loaded.
#' @param trafo function for tranforming residuals from probability scale to a
#' different distribution scale.
#' @param style character specifying the style of pithist. For \code{style = "histogram"}
#' a traditional PIT hisogram is drawn, for \code{style = "lines"} solely the upper border 
#' line is plotted. For \code{single_graph = TRUE}, always line-style PIT histograms are 
#' drawn.
#' @param type character. In case of discrete distributions should the PITs be
#' drawn randomly from the corresponding interval or distributed
#' proportionally? This argument is not fully supported yet, please keep to the default 
#' \code{"random"} for now. Additionally, expected (nonnormal) PIT histograms according to
#' Czado et al. (2009) can be chosen by setting the type to \code{"expected"}.
#' @param nsim integer. If \code{type} is \code{"random"} how many simulated
#' PITs should be drawn?
#' @param delta numeric. The minimal difference to compute the range of
#' proabilities corresponding to each observation according to get (randomized)
#' quantile residuals. For \code{NULL}, the minimal observed difference in the
#' resonse divided by \code{5e-6} is used.
#' @param freq logical. If \code{TRUE}, the PIT histogram is represented by
#' frequencies, the \code{counts} component of the result; if \code{FALSE},
#' probability densities, component \code{density}, are plotted (so that the
#' histogram has a total area of one).
#' @param breaks numeric. Breaks for the histogram intervals.
#' @param confint logical. Should confident intervals be drawn?
#' @param confint_level numeric. The confidence level required.
#' @param confint_type character. Which type of confidence interval should be plotted. According
#' to Agresti and Coull (1998), for interval estimation of binomial proportions
#' an approximation can be better than exact.
#' @param single_graph logical. Should all computed extended reliability
#' diagrams be plotted in a single graph? If yes, \code{style} must be set to \code{"lines"}.
#' @param xlim,ylim,xlab,ylab,main graphical parameters passed to
#' \code{\link{plot.pithist}} or \code{\link{autoplot.pithist}}.  
#' @param \dots further graphical parameters.
#' @return An object of class \code{"pithist"} inheriting from
#' \code{"data.frame"} or \code{"tibble"} conditional on the argument \code{class} 
#' with the following variables: \item{x}{histogram
#' interval midpoints on the x-axis,} \item{y}{bottom coordinate of the
#' histogram bars,} \item{width}{widths of the histogram bars,}
#' \item{ci_lwr, ci_upr}{lower and upper confidence interval bound,} \item{ref}{y-coordinates of
#' the reference curve.} Additionally, \code{freq}, \code{xlab},
#' \code{ylab}, \code{main}, and \code{confint_level} are stored as attributes.
#' @seealso \code{\link{plot.pithist}}, \code{\link{qresiduals}}, \code{\link{procast}}
#' @references 
#' Agresti A, Coull AB (1998). \dQuote{Approximate is Better than ``Exact''
#' for Interval Estimation of Binomial Proportions.} \emph{The American
#' Statistician}, \bold{52}(2), 119--126. \doi{10.1080/00031305.1998.10480550}
#'
#' Czado C, Gneiting T, Held L (2009). \dQuote{Predictive Model
#' Assessment for Count Data.} \emph{Biometrics}, \bold{65}(4), 1254--1261. 
#' \doi{10.2307/2981683}
#'  
#' Dawid AP (1984). \dQuote{Present Position and Potential Developments: Some
#' Personal Views: Statistical Theory: The Prequential Approach}, \emph{Journal of
#' the Royal Statistical Society: Series A (General)}, \bold{147}(2), 278--292.
#' \doi{10.2307/2981683}
#'
#' Diebold FX, Gunther TA, Tay AS (1998). \dQuote{Evaluating Density Forecasts
#' with Applications to Financial Risk Management}. \emph{International Economic
#' Review}, \bold{39}(4), 863--883. \doi{10.2307/2527342}
#' 
#' Gneiting T, Balabdaoui F, Raftery AE (2007). \dQuote{Probabilistic Forecasts,
#' Calibration and Sharpness}.  \emph{Journal of the Royal Statistical Society:
#' Series B (Methodological)}. \bold{69}(2), 243--268.
#' \doi{10.1111/j.1467-9868.2007.00587.x}
#' 
#' @keywords hplot
#' @examples
#' ## speed and stopping distances of cars
#' m1_lm <- lm(dist ~ speed, data = cars)
#' 
#' ## compute and plot pithist
#' pithist(m1_lm)
#' 
#' #-------------------------------------------------------------------------------
#' ## determinants for male satellites to nesting horseshoe crabs
#' data("CrabSatellites", package = "countreg")
#' 
#' ## linear poisson model
#' m1_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#' m2_pois <- glm(satellites ~ color, data = CrabSatellites, family = poisson)
#' 
#' ## compute and plot pithist as base graphic
#' p1 <- pithist(m1_pois, plot = FALSE)
#' p2 <- pithist(m2_pois, plot = FALSE)
#' 
#' ## plot combined pithist as "ggplot2" graphic
#' ggplot2::autoplot(c(p1, p2), single_graph = TRUE, style = "lines", col = c(1, 2))
#' 
#' @export
pithist <- function(object, ...) {
  UseMethod("pithist")
}


#' @rdname pithist
#' @method pithist default
#' @export
pithist.default <- function(object,
                            newdata = NULL,
                            plot = TRUE,
                            class = NULL,
                            trafo = NULL,
                            style = c("histogram", "lines"),
                            type = c("random", "expected", "proportional"), # FIXME: (ML) not yet implemented
                            nsim = 1L,
                            delta = NULL,
                            freq = FALSE,
                            breaks = NULL,
                            confint = TRUE,
                            confint_level = 0.95,
                            confint_type = c("exact", "approximation"),
                            single_graph = FALSE,
                            xlim = c(NA, NA),
                            ylim = c(0, NA),
                            xlab = "PIT",
                            ylab = if (freq) "Frequency" else "Density",
                            main = NULL,
                            ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * `object` and `newdata` w/i `newrepsone()`
  ## * `delta w/i `qresiduals()`
  ## * `confint` w/i `abline()`
  ## * `...` w/i `plot()` and `autoplot()`
  stopifnot(is.numeric(nsim), length(nsim) == 1)
  stopifnot(is.logical(freq))
  stopifnot(is.null(breaks) || (is.numeric(breaks) && is.null(dim(breaks))))
  stopifnot(
    is.numeric(confint_level),
    length(confint_level) == 1,
    confint_level >= 0,
    confint_level <= 1
  )
  stopifnot(is.logical(single_graph))
  stopifnot(length(xlim) == 2 && (all(is.na(xlim)) || is.numeric(xlim)))
  stopifnot(length(ylim) == 2 && (all(is.na(ylim)) || is.numeric(ylim)))
  stopifnot(length(xlab) == 1)
  stopifnot(length(ylab) == 1)
  stopifnot(length(main) == 1 | length(main) == 0)

  ## match arguments
  style <- match.arg(style)
  confint_type <- match.arg(confint_type)
  type <- match.arg(type)

  ## guess plotting flavor
  if (isFALSE(plot)) {
    plot <- "none"
  } else if (isTRUE(plot)) {
    plot <- if ("ggplot2" %in% .packages()) "ggplot2" else "base"
  } else if (!is.character(plot)) {
    plot <- "base"
  }
  plot <- try(match.arg(plot, c("none", "base", "ggplot2")))
  stopifnot(
    "`plot` must either be logical, or match the arguments 'none', 'base' or 'ggplot2'" =
      !inherits(plot, "try-error")
  )

  ## guess output class
  if (is.null(class)) {
    class <- if ("tibble" %in% .packages()) "tibble" else "data.frame"
  }
  class <- try(match.arg(class, c("tibble", "data.frame")))
  stopifnot(
    "`class` must either be NULL, or match the arguments 'tibble', or 'data.frame'" =
      !inherits(class, "try-error")
  )

  # -------------------------------------------------------------------
  # COMPUTATION OF PIT
  # -------------------------------------------------------------------
  ## get breaks: part 1 (due to "type = expected" must be done before)
  n <- NROW(newresponse(object, newdata = newdata))  # solely to get n to compute number of breaks
  if (is.null(breaks)) breaks <- c(4, 10, 20, 25)[cut(n, c(0, 50, 5000, 1000000, Inf))]

  ## either compute proportion exactly (to do...), "expected" (nonrandom) according to Czado et al. (2009) 
  ## or approximate by simulation
  if (type == "proportional") {
    ## FIXME: (ML)
    ## * implement proportional over the inteverals (e.g., below censoring point)
    ## * confusing naming, as `type` in `qresiduals()` must be `random` or `quantile`
    stop("not yet implemented")
  } else if (type == "random") {
    p <- qresiduals(object,
      newdata = newdata, trafo = trafo, delta = delta,
      type = "random", nsim = nsim
    )

    ## get breaks: part 2
    ## FIXME: (ML) maybe use xlim instead or `0` and `1`
    if (is.null(trafo)) {
      if (length(breaks) == 1L) breaks <- seq(0, 1, length.out = breaks + 1L)
    } else {
      tmp_range <- max(abs(p))
      if (length(breaks) == 1L) breaks <- seq(-tmp_range, tmp_range, length.out = breaks + 1L)
    }

    ## collect everything as data.frame
    tmp_hist <- hist(p, breaks = breaks, plot = FALSE)
    ## TODO: (ML) Maybe get rid of `hist()`

  } else if (type == "expected") {
    ## compare "nonrandom" in Czado et al. (2009) 

    ## minimum and maximum PIT for each observation (P_x-1 and P_x)
    p <- qresiduals(object, 
      newdata = newdata, trafo = trafo, delta = delta,
      type = "quantile", prob = c(0, 1)
    )

    ## equation 2: CDF for each PIT (continuous vs. discrete)
    F <- if (all(abs(p[, 2L] - p[, 1L]) < sqrt(.Machine$double.eps))) {
      function(u) as.numeric(u >= p[, 1L]) 
      ## FIXME: (Z) check inequality sign to cover include.lowest/right options
    } else {
      function(u) punif(u, min = p[, 1L], max = p[, 2L]) ## pmin(1, pmax(0, (u - p[, 1L]) / (p[, 2L] - p[, 1L])))
    }

    ## get breaks: part 2
    ## FIXME: (ML) maybe use xlim instead or `0` and `1`
    if (is.null(trafo)) {
      if (length(breaks) == 1L) breaks <- seq(0, 1, length.out = breaks + 1L)
    } else {
      tmp_range <- max(abs(p))
      if (length(breaks) == 1L) breaks <- seq(-tmp_range, tmp_range, length.out = breaks + 1L)
    }

    ## equation 3 and computation of probability for each interval (f_j)
    f <- diff(colMeans(sapply(breaks, F)))

    ## collect everything as data.frame
    tmp_hist <- list(
      breaks = breaks,
      counts = f * n,
      density  = f / diff(breaks),
      mids = 0.5 * (breaks[-1L] + breaks[-length(breaks)])
    )
  }

  ## compute ci interval
  ## FIXME: (Z) confidence interval really only valid for equidistant breaks
  if (!is.null(trafo)) {
    ci <- c(NA_real_, NA_real_)  # must be numeric for ggplot2
    warning("confint is not yet implemented when employing a `trafo`")
  } else if (confint_type == "approximation") {
    if (length(unique(round(diff(breaks), 10))) > 1) {
      warning("confint is not yet implemented for non equidistant breaks for argument `type = 'approximation'`")
      ci <- c(NA_real_, NA_real_)  # must be numeric for ggplot2
    } else {
      ci <- get_confint_agresti(
        NROW(p) / (length(breaks) - 1),
        NROW(p),
        confint_level,
        length(breaks) - 1,
        freq
      )
    }
  } else {
    ci <- get_confint(NROW(p), breaks, confint_level, freq)
  }

  ## perfect prediction
  #pp <- ifelse(freq, NROW(p) / (length(breaks) - 1), 1)
  ## FIXME: (ML) 
  ## * Is this correct?
  ## * This is not the same as uncommented above and for small n qbinom(0.5, ...) is not equal 1 
  if (!is.null(trafo)) {
    pp <- NA_real_  # must be numeric for ggplot2
    warning("ref is not yet implemented when employing a `trafo`")
  } else {
    pp <- get_pp(NROW(p), breaks, freq)
  }

  ## labels
  if (is.null(main)) main <- deparse(substitute(object))

  # -------------------------------------------------------------------
  # OUTPUT AND OPTIONAL PLOTTING
  # -------------------------------------------------------------------
  if (freq) {
    rval <- data.frame(
      x = tmp_hist$mids,
      y = tmp_hist$counts / nsim,  ## FIXME: (ML) Double check if correct even for point masses
      width = diff(tmp_hist$breaks),
      ci_lwr = ci[[1]],
      ci_upr = ci[[2]],
      ref = pp
    )
  } else {
    rval <- data.frame(
      x = tmp_hist$mids,
      y = tmp_hist$density,
      width = diff(tmp_hist$breaks),
      ci_lwr = ci[[1]],
      ci_upr = ci[[2]],
      ref = pp
    )
  }

  ## attributes for graphical display
  attr(rval, "freq") <- freq
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "confint_level") <- ifelse(confint, confint_level, NA)

  ## add class
  if (class == "data.frame") {
    class(rval) <- c("pithist", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("pithist", class(rval))
  }

  ## plot by default
  if (plot == "ggplot2") {
    try(print(ggplot2::autoplot(rval,
      style = style, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, confint = confint, ...
    )))
  } else if (plot == "base") {
    try(plot(rval,
      style = style, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, confint = confint, ...
    ))
  }

  ## return invisibly
  invisible(rval)
}


#' @export
c.pithist <- function(...) {
  # -------------------------------------------------------------------
  # GET DATA
  # -------------------------------------------------------------------
  ## list of pithists
  rval <- list(...)

  ## set class to tibble if any rval is a tibble
  if (any(do.call("c", lapply(rval, class)) %in% "tbl")) {
    class <- "tibble"
  } else {
    class <- "data.frame"
  }

  ## convert always to data.frame
  if (class == "tibble") {
    rval <- lapply(rval, as.data.frame)
  }

  ## group sizes
  for (i in seq_along(rval)) {
    if (is.null(rval[[i]]$group)) rval[[i]]$group <- 1L
  }
  n <- lapply(rval, function(r) table(r$group))

  # -------------------------------------------------------------------
  # PREPARE DATA
  # -------------------------------------------------------------------
  ## check if all of same `freq`
  freq <- unlist(lapply(rval, function(r) attr(r, "freq")))
  stopifnot(length(unique(freq)) == 1)

  ## labels
  xlab <- unlist(lapply(rval, function(r) attr(r, "xlab")))
  ylab <- unlist(lapply(rval, function(r) attr(r, "ylab")))
  confint_level <- unlist(lapply(rval, function(r) attr(r, "confint_level")))
  nam <- names(rval)
  main <- if (is.null(nam)) {
    as.vector(sapply(rval, function(r) attr(r, "main")))
  } else {
    make.unique(rep.int(nam, sapply(n, length)))
  }
  n <- unlist(n)

  # -------------------------------------------------------------------
  # RETURN DATA
  # -------------------------------------------------------------------
  ## combine and return
  rval <- do.call("rbind.data.frame", rval)
  rval$group <- if (length(n) < 2L) NULL else rep.int(seq_along(n), n)
  attr(rval, "freq") <- freq
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "confint_level") <- confint_level

  ## set class to data.frame or tibble
  if (class == "data.frame") {
    class(rval) <- c("pithist", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("pithist", class(rval))
  }

  ## return
  return(rval)
}


#' @export
rbind.pithist <- c.pithist


#' S3 Methods for Plotting PIT Histograms
#' 
#' Generic plotting functions for PIT histograms of the class \code{"pithist"}
#' computed by \code{link{pithist}}. 
#' 
#' PIT histograms graphically evaluate the probability integral transform (PIT),
#' i.e., the value that the predictive CDF attains at the observation, with a
#' uniform distribution. For a well calibrated model fit, the observation will be
#' drawn from the predictive distribution and the PIT will have a standard uniform
#' distribution. 
#'
#' PIT histograms can be rendered as \code{ggplot2} or base R graphics by using
#' the generics \code{\link[ggplot2]{autoplot}} or \code{\link[graphics]{plot}}. 
#' For a single base R graphically panel, \code{\link{lines}} adds an additional PIT histogram.
#' 
#' @aliases plot.pithist lines.pithist autoplot.pithist
#' @param object,x an object of class \code{\link{pithist}}.
#' @param single_graph logical. Should all computed extended reliability
#' diagrams be plotted in a single graph?
#' @param style character specifying the style of pithist. For \code{style = "histogram"}
#' a traditional PIT hisogram is drawn, for \code{style = "lines"} solely the upper border 
#' line is plotted. For \code{single_graph = TRUE}, always line-style PIT histograms are 
#' plotted.
#' @param confint logical. Should confident intervals be drawn?
#' @param xlim,ylim graphical parameters. These may pertain either to the whole
#' plot or just the histogram or just the fitted line.
#' @param xlab,ylab,main graphical parameters.
#' @param \dots further graphical parameters.
#' @param ref,col,fill,border,alpha_min,lwd,lty,axes,box additional graphical
#' parameters for base plots, whereby \code{x} is a object of class \code{pithist}.
#' @param colour,size,linetype,legend graphical parameters passed for 
#' \code{ggplot2} style plots, whereby \code{object} is a object of class \code{pithist}.
#' @seealso \code{\link{pithist}}, \code{\link{procast}}, \code{\link[graphics]{hist}}
#' @references 
#' Agresti A, Coull AB (1998). \dQuote{Approximate is Better than ``Exact''
#' for Interval Estimation of Binomial Proportions.} \emph{The American
#' Statistician}, \bold{52}(2), 119--126. \doi{10.1080/00031305.1998.10480550}
#'
#' Czado C, Gneiting T, Held L (2009). \dQuote{Predictive Model
#' Assessment for Count Data.} \emph{Biometrics}, \bold{65}(4), 1254--1261. 
#' \doi{10.2307/2981683}
#'  
#' Dawid AP (1984). \dQuote{Present Position and Potential Developments: Some
#' Personal Views: Statistical Theory: The Prequential Approach}, \emph{Journal of
#' the Royal Statistical Society: Series A (General)}, \bold{147}(2), 278--292.
#' \doi{10.2307/2981683}
#'
#' Diebold FX, Gunther TA, Tay AS (1998). \dQuote{Evaluating Density Forecasts
#' with Applications to Financial Risk Management}. \emph{International Economic
#' Review}, \bold{39}(4), 863--883. \doi{10.2307/2527342}
#' 
#' Gneiting T, Balabdaoui F, Raftery AE (2007). \dQuote{Probabilistic Forecasts,
#' Calibration and Sharpness}.  \emph{Journal of the Royal Statistical Society:
#' Series B (Methodological)}. \bold{69}(2), 243--268.
#' \doi{10.1111/j.1467-9868.2007.00587.x}
#' 
#' @keywords hplot
#' @examples
#' 
#' ## speed and stopping distances of cars
#' m1_lm <- lm(dist ~ speed, data = cars)
#' 
#' ## compute and plot pithist
#' pithist(m1_lm)
#' 
#' ## customize colors and style
#' pithist(m1_lm, ref = "blue", lty = 2, pch = 20, style = "lines")
#' 
#' ## add separate model
#' if (require("crch", quietly = TRUE)) {
#'   m1_crch <- crch(dist ~ speed | speed, data = cars)
#'   lines(pithist(m1_crch, plot = FALSE), col = 2, lty = 2, confint = 2)
#' }
#' 
#' #-------------------------------------------------------------------------------
#' if (require("crch")) {
#' 
#'   ## precipitation observations and forecasts for Innsbruck
#'   data("RainIbk", package = "crch")
#'   RainIbk <- sqrt(RainIbk)
#'   RainIbk$ensmean <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, mean)
#'   RainIbk$enssd <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, sd)
#'   RainIbk <- subset(RainIbk, enssd > 0)
#' 
#'   ## linear model w/ constant variance estimation
#'   m2_lm <- lm(rain ~ ensmean, data = RainIbk)
#' 
#'   ## logistic censored model 
#'   m2_crch <- crch(rain ~ ensmean | log(enssd), data = RainIbk, left = 0, dist = "logistic")
#' 
#'   ## compute pithists
#'   pit2_lm <- pithist(m2_lm, plot = FALSE)
#'   pit2_crch <- pithist(m2_crch, plot = FALSE)
#' 
#'   ## plot in single graph with style "lines"
#'   plot(c(pit2_lm, pit2_crch), col = c(1, 2), confint = c(1, 2), ref = 3, 
#'     style = "lines", single_graph = TRUE)
#' }
#' 
#' #-------------------------------------------------------------------------------
#' ## determinants for male satellites to nesting horseshoe crabs
#' data("CrabSatellites", package = "countreg")
#' 
#' ## linear poisson model
#' m3_pois  <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#' 
#' ## compute and plot pithist as "ggplot2" graphic
#' pithist(m3_pois, plot = "ggplot2")
#'
#' @export
plot.pithist <- function(x,
                         single_graph = FALSE,
                         style = c("histogram", "lines"),
                         confint = TRUE,
                         ref = TRUE,
                         xlim = c(0, 1),
                         ylim = c(0, NA),
                         xlab = NULL,
                         ylab = NULL,
                         main = NULL,
                         col = "black",
                         fill = adjustcolor("black", alpha.f = 0.2),
                         border = "black",
                         alpha_min = 0.2,
                         lwd = NULL,
                         lty = 1,
                         axes = TRUE,
                         box = TRUE,
                         ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * lengths of all arguments are checked by recycling
  ## * `ref` and `confint` w/i `abline()`
  ## * `xlab`, `ylab`, `main`, `col`, `fill`, `lwd`, `lty` and `...` w/i `plot()`
  ## * `alpha_min` w/i `set_minimum_transparency()`
  stopifnot(is.logical(single_graph))
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))
  if (single_graph) {
    stopifnot(
      "for `single_graph` all `freq` in attr of object `x` must be of the same type" =
        length(unique(attr(x, "freq"))) == 1
    )
  }
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## set style
  style <- match.arg(style)
  if (n > 1 && single_graph && style == "histogram") {
    message(" * For several histograms in a single graph solely line style histograms can be plotted. \n * For proper usage, set `style` = 'lines' when numbers of histograms greater one and `single_graph` = TRUE.")
    style <- "lines"
  }

  ## set lwd
  if (is.null(lwd)) lwd <- if (style == "histogram") 1.5 else 2

  ## recycle arguments for plotting to match the number of groups
  if (is.list(xlim)) xlim <- as.data.frame(do.call("rbind", xlim))
  if (is.list(ylim)) ylim <- as.data.frame(do.call("rbind", ylim))
  plot_arg <- data.frame(1:n, confint, ref,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]],
    border, col, fill, alpha_min, lwd, lty, axes, box
  )[, -1]

  ## annotation
  if (single_graph) {
    if (is.null(xlab)) xlab <- "PIT"
    if (is.null(ylab)) ylab <- if (all(attr(x, "freq"))) "Frequency" else "Density"
    if (is.null(main)) main <- "PIT histogram"
  } else {
    if (is.null(xlab)) xlab <- TRUE
    if (is.null(ylab)) ylab <- TRUE
    if (is.null(main)) main <- TRUE
    xlab <- rep(xlab, length.out = n)
    ylab <- rep(ylab, length.out = n)
    main <- rep(main, length.out = n)
    if (is.logical(xlab)) xlab <- ifelse(xlab, attr(x, "xlab"), "")
    if (is.logical(ylab)) ylab <- ifelse(ylab, attr(x, "ylab"), "")
    if (is.logical(main)) main <- ifelse(main, attr(x, "main"), "")
  }

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR 'HISTOGRAM-STYLE PITHIST'
  # -------------------------------------------------------------------
  pithist_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## rect elements
    xleft <- d$x - d$width / 2
    xright <- d$x + d$width / 2
    y <- d$y

    ## get xlim and ylim
    ylim_idx <- c(is.na(plot_arg$ylim1[j]), is.na(plot_arg$ylim2[j]))
    xlim_idx <- c(is.na(plot_arg$xlim1[j]), is.na(plot_arg$xlim2[j]))
    if (any(xlim_idx)) {
      plot_arg[j, c("xlim1", "xlim2")[xlim_idx]] <- range(c(xleft, xright))[xlim_idx]
    }
    if (any(ylim_idx)) {
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <- range(c(0, y, d$ci_lwr, d$ci_upr), na.rm = TRUE)[ylim_idx]
    }

    ## trigger plot
    if (j == 1 || (!single_graph && j > 1)) {
      plot(0, 0,
        type = "n", xlim = c(plot_arg$xlim1[j], plot_arg$xlim2[j]),
        ylim = c(plot_arg$ylim1[j], plot_arg$ylim2[j]),
        xlab = xlab[j], ylab = ylab[j], main = main[j], axes = FALSE, ...
      )
      if (plot_arg$axes[j]) {
        axis(1)
        axis(2)
      }
      if (plot_arg$box[j]) {
        box()
      }
    }

    if (plot_arg$fill[j] == adjustcolor("black", alpha.f = 0.2) &&
      plot_arg$col[j] != "black") {
      message(" * As the argument `col` is set but no argument `fill` is specified, \n   the former is used for colorizing the PIT histogram. \n * For proper usage, solely provide `fill` for histogram style plots.")
      plot_arg$fill[j] <- plot_arg$col[j]
    }

    ## plot pithist
    rect(xleft, 0, xright, y,
      border = plot_arg$border[j], col = plot_arg$fill[j],
      lty = plot_arg$lty[j]
    )

    ## plot ref line
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- 2  # red

      ref_z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
      ref_y <- c(d$ref, d$ref[NROW(d)])
      lines(ref_y ~ ref_z, type = "s", col = plot_arg$ref[j], lty = 1, lwd = plot_arg$lwd[j])
    }

    ## plot confint lines
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- 2  # red

      ## lower confint line
      ci_lwr_z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
      ci_lwr_y <- c(d$ci_lwr, d$ci_lwr[NROW(d)])
      lines(ci_lwr_y ~ ci_lwr_z, type = "s", col = plot_arg$confint[j], lty = 2, lwd = plot_arg$lwd[j])

      ## upper confint line
      ci_upr_z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
      ci_upr_y <- c(d$ci_upr, d$ci_upr[NROW(d)])
      lines(ci_upr_y ~ ci_upr_z, type = "s", col = plot_arg$confint[j], lty = 2, lwd = plot_arg$lwd[j])
    }
  }


  # -------------------------------------------------------------------
  # FUNCTION TO TRIGGER FIGURE AND PLOT CONFINT FOR 'LINE-STYLE PITHIST'
  # -------------------------------------------------------------------
  pitlines_trigger <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## step elements
    z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
    y <- c(d$y, d$y[NROW(d)])

    ## get xlim and ylim
    ylim_idx <- c(is.na(plot_arg$ylim1[j]), is.na(plot_arg$ylim2[j]))
    xlim_idx <- c(is.na(plot_arg$xlim1[j]), is.na(plot_arg$xlim2[j]))
    if (any(xlim_idx)) {
      plot_arg[j, c("xlim1", "xlim2")[xlim_idx]] <- range(z)[xlim_idx]
    }
    if (any(ylim_idx) && !single_graph) {
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <- range(c(y, d$ci_lwr, d$ci_upr), na.rm = TRUE)[ylim_idx]
    }
    if (any(ylim_idx) && single_graph) {
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <- range(c(x$y, x$ci_lwr, x$ci_upr), na.rm = TRUE)[ylim_idx]
    }

    ## trigger plot
    if (j == 1 || (!single_graph && j > 1)) {
      plot(0, 0,
        type = "n", xlim = c(plot_arg$xlim1[j], plot_arg$xlim2[j]),
        ylim = c(plot_arg$ylim1[j], plot_arg$ylim2[j]),
        xlab = xlab[j], ylab = ylab[j], xaxs = "i", main = main[j], axes = FALSE, ...
      )
      if (plot_arg$axes[j]) {
        axis(1)
        axis(2)
      }
      if (plot_arg$box[j]) {
        box()
      }
    }

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]
      ref_z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
      polygon(
        c(
          rep(ref_z, each = 2)[-c(1, length(ref_z) * 2)],
          rev(rep(ref_z, each = 2)[-c(1, length(ref_z) * 2)])
        ),
        c(
          rep(d$ci_lwr, each = 2),
          rev(rep(d$ci_upr, each = 2))
        ), 
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min = plot_arg$alpha_min[j]),
        border = NA
      )
    }
  }

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR 'LINE_STYLE PITHIST'
  # -------------------------------------------------------------------
  pitlines_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## step elements
    z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
    y <- c(d$y, d$y[NROW(d)])

    ## plot ref line
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"

      ref_z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
      ref_y <- c(d$ref, d$ref[NROW(d)])
      lines(ref_y ~ ref_z, type = "s", col = plot_arg$ref[j], lty = 2, lwd = 1)
    }

    ## plot stepfun
    lines(y ~ z, type = "s", lwd = plot_arg$lwd[j], lty = plot_arg$lty[j], col = plot_arg$col[j])
  }

  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  ## set up necessary panels
  if (!single_graph && n > 1L) {
    old_pars <- par(mfrow = n2mfrow(n))
    on.exit(par(old_pars), add = TRUE)
  }

  ## draw polygons first
  if (single_graph || n == 1) {
    if (style == "histogram") {
      for (i in 1L:n) pithist_plot(x[x$group == i, ], ...)
    } else {
      for (i in 1L:n) pitlines_trigger(x[x$group == i, ], ...)
      for (i in 1L:n) pitlines_plot(x[x$group == i, ], ...)
    }
  } else {
    if (style == "histogram") {
      for (i in 1L:n) pithist_plot(x[x$group == i, ], ...)
    } else {
      for (i in 1L:n) {
        pitlines_trigger(x[x$group == i, ], ...)
        pitlines_plot(x[x$group == i, ], ...)
      }
    }
  }
}


#' @rdname plot.pithist
#' @method lines pithist
#' @export
lines.pithist <- function(x,
                          confint = FALSE,
                          ref = FALSE,
                          col = "black",
                          fill = adjustcolor("black", alpha.f = 0.2),
                          alpha_min = 0.2,
                          lwd = 2,
                          lty = 1,
                          ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * lengths of all arguments are checked by recycling
  ## * `ref` and `confint` w/i `abline()`
  ## * `col`, `fill`, `lwd`, `lty` and `...` w/i `lines()`
  ## * `alpha_min` w/i `set_minimum_transparency()`
  stopifnot(
    "all `freq` in attr of object `x` must be of the same type" =
      length(unique(attr(x, "freq"))) == 1
  )

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(
    1:n, confint, ref, col, fill, alpha_min, lwd, lty
  )[, -1]

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR LINES
  # -------------------------------------------------------------------
  pitlines_lines <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## step elements
    z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
    y <- c(d$y, d$y[NROW(d)])

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]
     ref_z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
      polygon( 
        c(
          rep(ref_z, each = 2)[-c(1, length(ref_z) * 2)],
          rev(rep(ref_z, each = 2)[-c(1, length(ref_z) * 2)])
        ),
        c(
          rep(d$ci_lwr, each = 2),
          rev(rep(d$ci_upr, each = 2))
        ), 
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min = plot_arg$alpha_min[j]),
        border = NA
      )
    }

    ## plot ref line
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"

      ref_z <- c(d$x - d$width / 2, d$x[NROW(d)] + d$width[NROW(d)] / 2)
      ref_y <- c(d$ref, d$ref[NROW(d)])
      lines(ref_y ~ ref_z, type = "s", col = plot_arg$ref[j], lty = 2, lwd = 1)
    }

    ## plot stepfun
    lines.default(y ~ z, type = "s", lwd = plot_arg$lwd[j], lty = plot_arg$lty[j], col = plot_arg$col[j])
  }

  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  for (i in 1L:n) {
    pitlines_lines(x[x$group == i, ], ...)
  }
}


##autoplot.pithist <- function(object,
##                             single_graph = FALSE,
##                             style = c("histogram", "lines"),
##                             confint = TRUE,
##                             ref = TRUE,
##                             xlim = c(0, 1),
##                             ylim = c(0, NA),
##                             xlab = NULL,
##                             ylab = NULL,
##                             main = NULL,
##                             colour = "black",
##                             fill = "darkgray",
##                             border = "black",
##                             alpha_min = 0.2,
##                             size = NULL,
##                             linetype = 1,
##                             legend = FALSE,
##                             ...) {
##  # -------------------------------------------------------------------
##  # SET UP PRELIMINARIES
##  # -------------------------------------------------------------------
##  ## get base style arguments
##  add_arg <- list(...)
##  if (!is.null(add_arg$lwd)) size <- add_arg$lwd
##  if (!is.null(add_arg$lty)) linetype <- add_arg$lty
##
##  ## sanity checks
##  stopifnot(is.logical(single_graph))
##  if (single_graph) {
##    stopifnot(
##      "for `single_graph` all `freq` in attr of `object` must be of the same type" =
##        length(unique(attr(object, "freq"))) == 1
##    )
##  }
##  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
##  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))
##
##  ## convert data always to data.frame
##  object <- as.data.frame(object)
##
##  ## determine grouping
##  if (is.null(object$group)) object$group <- 1L
##  n <- max(object$group)
##
##  ## get title
##  if (!is.null(main)) {
##    title <- main[1]
##    object$title <- factor(title)
##  }
##
##  ## get annotations in the right lengths
##  if (is.null(xlab)) xlab <- attr(object, "xlab")
##  xlab <- paste(unique(xlab), collapse = "/")
##  if (is.null(ylab)) ylab <- attr(object, "ylab")
##  ylab <- paste(unique(ylab), collapse = "/")
##  if (is.null(main)) main <- attr(object, "main")
##  main <- make.names(rep_len(main, n), unique = TRUE)
##
##  ## prepare grouping
##  object$group <- factor(object$group, levels = 1L:n, labels = main)
##
##  # -------------------------------------------------------------------
##  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING
##  # -------------------------------------------------------------------
##  ## determine which style should be plotted
##  style <- match.arg(style)
##  if (n > 1 && single_graph && style == "histogram") {
##    message(" * For several histograms in a single graph solely line style histograms can be plotted. \n * For proper usage, set `style` = 'lines' when numbers of histograms greater one and `single_graph` = TRUE.")
##    style <- "lines"
##  }
##
##  ## set size
##  if (is.null(size)) size <- if (style == "histogram") 0.7 else 1
##
##  ## set color to 2 (red) or NA for not plotting
##  if (is.logical(ref) & style == "histogram") {
##    ref <- ifelse(ref, 2, NA)
##  } else if (is.logical(ref) & style == "lines") {
##    ref <- ifelse(ref, 1, NA)
##  }
##
##  ## only needed for `style == "lines"`
##  ## stat helper function to get left/right points from respective mid points
##  calc_pit_points <- ggplot2::ggproto("calc_pit_points", ggplot2::Stat,
##
##    # required as we operate on groups (facetting)
##    compute_group = function(data, scales) {
##      ## manipulate object
##      nd <- data.frame(
##        x = c(data$x - data$width / 2, data$x[NROW(data)] + data$width[NROW(data)] / 2),
##        y = c(data$y, data$y[NROW(data)])
##      )
##      nd
##    },
##
##    # tells us what we need
##    required_aes = c("x", "y")
##  )
##
##  ## helper function to get right number of indices after using `stat = calc_pit_points`
##  ## TODO: (ML) This must (!) be done smoother (UPDATE: at the moment not needed)
##  #calc_pit_points_index <- function(group) {
##  #  idx <- cumsum(rle(as.numeric(group))$lengths)
##  #  idx2 <- lapply(1:length(idx), function(i) c(seq(1, idx[i]), idx[i]))
##  #  if (length(idx2) > 1) {
##  #    idx3 <- unlist(c(idx2[[1]], sapply(2:length(idx2), function(j) idx2[[j]][idx2[[j]] > max(idx2[[j-1]])])))
##  #  } else {
##  #    idx3 <- idx2[[1]]
##  #  }
##  #  return(idx3)
##  #}
##
##  ## recycle arguments for plotting to match the number of groups (for `scale_<...>_manual()`)
##  plot_arg <- data.frame(
##    1:n,
##    fill, colour, size, linetype, confint, alpha_min
##  )[, -1]
##
##  ## prepare fill and confint color depending on style
##  if (style == "histogram") {
##
##    ## set color to 2 (red) or NA for not plotting
##    if (is.logical(confint)) confint <- ifelse(confint, 2, NA)
##
##    ## check if colour and no fill is set
##    if (all(fill == "darkgray") && any(colour != "black")) {
##      message(" * As the argument `colour` is set but no argument `fill` is specified, \n   the former is used for colorizing the PIT histogram. \n * For proper usage, solely provide `fill` for histogram style plots.")
##      plot_arg$fill <- plot_arg$colour
##    }
##  } else {
##    if (is.logical(plot_arg$confint)) {
##
##      ## use fill and set alpha
##      plot_arg$fill <- sapply(seq_along(plot_arg$fill), function(idx) {
##        set_minimum_transparency(plot_arg$fill[idx], alpha_min = plot_arg$alpha_min[idx])
##      })
##
##      ## set color to NA for not plotting
##      plot_arg$fill[!plot_arg$confint] <- NA
##    } else {
##      ## use confint and set alpha
##      plot_arg$fill <- sapply(seq_along(plot_arg$confint), function(idx) {
##        set_minimum_transparency(plot_arg$confint[idx], alpha_min = plot_arg$alpha_min[idx])
##      })
##    }
##  }
##
##  ## FIXME: (ML) Decide if plot_arg2 or plot_arg3 should be used
##  ## recycle arguments for plotting to match the length (rows) of the object (for geom w/ aes)
##  plot_arg2 <- data.frame(1:n, border, colour, ref, confint)[, -1]
##  plot_arg2 <- as.data.frame(lapply(plot_arg2, rep, table(object$group)))
##
##  ## new approach, maybe plot (ref, confint) same for all plots?!
##  plot_arg3 <- list("ref" = ref, "confint" = confint)
##
##  # -------------------------------------------------------------------
##  # MAIN PLOTTING
##  # -------------------------------------------------------------------
##  if (style == "histogram") {
##    ## actual plotting
##    rval <- ggplot2::ggplot(
##      object,
##      ggplot2::aes_string(x = "x", y = "y / 2", width = "width", height = "y")
##    ) +
##      ggplot2::geom_tile(ggplot2::aes_string(fill = "group"),
##        colour = plot_arg2$border
##      ) 
##
##    ## FIXME: (ML) Does not work for object w/ and w/o "ref" (NA)
##    if (!all(is.na(object$ref))) {
##      rval <- rval + 
##        ggplot2::geom_step(ggplot2::aes_string(x = "x", y = "ref", size = "group"),
##          colour = plot_arg3$ref, #plot_arg2$ref[calc_pit_points_index(object$group)],
##          stat = calc_pit_points,
##          na.rm = TRUE
##        ) 
##    }
##
##    ## FIXME: (ML) Does not work for object w/ and w/o "ci_lwr" (NA)
##    if (!all(is.na(object$ci_lwr))) {
##      rval <- rval +
##        ggplot2::geom_step(ggplot2::aes_string(x = "x", y = "ci_lwr", size = "group"),
##          colour = plot_arg3$confint, #plot_arg2$confint[calc_pit_points_index(object$group)],
##          linetype = 2,
##          stat = calc_pit_points,
##          na.rm = TRUE
##        )
##    }
##
##    ## FIXME: (ML) Does not work for object w/ and w/o "ci_upr" (NA)
##    if (!all(is.na(object$ci_upr))) {
##      rval <- rval + 
##        ggplot2::geom_step(ggplot2::aes_string(x = "x", y = "ci_upr", size = "group"),
##          colour = plot_arg3$confint, #plot_arg2$confint[calc_pit_points_index(object$group)],
##          linetype = 2,
##          stat = calc_pit_points,
##          na.rm = TRUE
##        ) 
##    }
##
##    ## set the colors, shapes, etc. for the groups
##    rval <- rval +
##      ggplot2::scale_fill_manual(values = plot_arg$fill) +
##      ggplot2::scale_size_manual(values = plot_arg$size)
##
##    ## annotation
##    rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
##
##    ## add legend
##    if (legend) {
##      rval <- rval + ggplot2::labs(fill = "Model") +
##        ggplot2::guides(fill = "legend", size = "none")
##    } else {
##      rval <- rval + ggplot2::guides(fill = "none", size = "none")
##    }
##
##    ## set x and y limits
##    rval <- rval + ggplot2::scale_x_continuous(limits = xlim, expand = c(0.01, 0.01))
##    rval <- rval + ggplot2::scale_y_continuous(limits = ylim, expand = c(0.01, 0.01))
##  } else {
##
##    ## actual plotting
##    rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y", width = "width")) +
##      ggplot2::geom_tile(
##        ggplot2::aes_string(
##          x = "x", 
##          y = "(ci_upr + ci_lwr) / 2", 
##          width = "width", 
##          height = "ci_upr - ci_lwr", 
##          fill = "group"
##        ),
##        colour = NA, show.legend = FALSE, na.rm = TRUE
##      ) +
##      ggplot2::geom_step(ggplot2::aes_string(colour = "group", size = "group", linetype = "group"),
##        stat = calc_pit_points
##      )
##
##    ## FIXME: (ML) Does not work for object w/ and w/o "ref" (NA)
##    if (!all(is.na(object$ref))) {
##      rval <- rval + 
##        ggplot2::geom_step(ggplot2::aes_string(x = "x", y = "ref"),
##          colour = plot_arg3$ref, #plot_arg2$ref[calc_pit_points_index(object$group)],
##          linetype = 2,
##          size = 1,
##          stat = calc_pit_points,
##          na.rm = TRUE
##        ) 
##    }
##
##    ## set the colors, shapes, etc.
##    rval <- rval +
##      ggplot2::scale_colour_manual(values = plot_arg$colour) +
##      ggplot2::scale_fill_manual(values = plot_arg$fill) +
##      ggplot2::scale_size_manual(values = plot_arg$size) +
##      ggplot2::scale_linetype_manual(values = plot_arg$linetype)
##
##    ## add annotation
##    rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
##
##    ## add legend
##    if (legend) {
##      rval <- rval + ggplot2::labs(colour = "Model") +
##        ggplot2::guides(colour = "legend", size = "none", linetype = "none")
##    } else {
##      rval <- rval + ggplot2::guides(colour = "none", size = "none", linetype = "none")
##    }
##
##    ## set x and y limits
##    rval <- rval + ggplot2::scale_x_continuous(limits = xlim, expand = c(0.01, 0.01))
##    rval <- rval + ggplot2::scale_y_continuous(limits = ylim, expand = c(0.01, 0.01))
##  }
##
##  # -------------------------------------------------------------------
##  # GROUPING (IF ANY) AND RETURN PLOT
##  # -------------------------------------------------------------------
##  ## grouping
##  if (!single_graph && n > 1L) {
##    rval <- rval + ggplot2::facet_grid(group ~ .)
##  } else if (!is.null(object$title)) {
##    rval <- rval + ggplot2::facet_wrap(title ~ .)
##  }
##
##  ## return ggplot object
##  rval
##}

#' @rdname plot.pithist
#' @method autoplot pithist
#' @exportS3Method ggplot2::autoplot
autoplot.pithist <- function(object,
                             single_graph = FALSE,
                             style = c("histogram", "lines"),
                             confint = TRUE,
                             ref = TRUE,
                             xlim = c(0, 1),
                             ylim = c(0, NA),
                             xlab = NULL,
                             ylab = NULL,
                             main = NULL,
                             colour = "black",
                             fill = "darkgray",
                             border = "black",
                             alpha_min = 0.2,
                             size = NULL,
                             linetype = 1,
                             legend = FALSE,
                             ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## get base style arguments
  add_arg <- list(...)
  if (!is.null(add_arg$lwd)) size <- add_arg$lwd
  if (!is.null(add_arg$lty)) linetype <- add_arg$lty

  ## sanity checks
  stopifnot(is.logical(single_graph))
  if (single_graph) {
    stopifnot(
      "for `single_graph` all `freq` in attr of `object` must be of the same type" =
        length(unique(attr(object, "freq"))) == 1
    )
  }
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))

  ## convert data always to data.frame
  object <- as.data.frame(object)

  ## determine grouping
  if (is.null(object$group)) object$group <- 1L
  n <- max(object$group)

  ## get title
  if (!is.null(main)) {
    title <- main[1]
    object$title <- factor(title)
  }

  ## get annotations in the right lengths
  if (is.null(xlab)) xlab <- attr(object, "xlab")
  xlab <- paste(unique(xlab), collapse = "/")
  if (is.null(ylab)) ylab <- attr(object, "ylab")
  ylab <- paste(unique(ylab), collapse = "/")
  if (is.null(main)) main <- attr(object, "main")
  main <- make.names(rep_len(main, n), unique = TRUE)

  ## prepare grouping
  object$group <- factor(object$group, levels = 1L:n, labels = main)

  # -------------------------------------------------------------------
  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING
  # -------------------------------------------------------------------
  ## determine which style should be plotted
  style <- match.arg(style)
  if (n > 1 && single_graph && style == "histogram") {
    message(" * For several histograms in a single graph solely line style histograms can be plotted. \n * For proper usage, set `style` = 'lines' when numbers of histograms greater one and `single_graph` = TRUE.")
    style <- "lines"
  }

  ## set size
  if (is.null(size)) size <- if (style == "histogram") 0.7 else 1

  ## set color to 2 (red) or NA for not plotting
  if (is.logical(ref) & style == "histogram") {
    ref <- ifelse(ref, 2, NA)
  } else if (is.logical(ref) & style == "lines") {
    ref <- ifelse(ref, 1, NA)
  }

  ## only needed for `style == "lines"`
  ## stat helper function to get left/right points from respective mid points
  calc_pit_points <- ggplot2::ggproto("calc_pit_points", ggplot2::Stat,

    # required as we operate on groups (facetting)
    compute_group = function(data, scales) {
      ## manipulate object
      nd <- data.frame(
        x = c(data$x - data$width / 2, data$x[NROW(data)] + data$width[NROW(data)] / 2),
        y = c(data$y, data$y[NROW(data)])
      )
      nd
    },

    # tells us what we need
    required_aes = c("x", "y")
  )

  ## helper function to get right number of indices after using `stat = calc_pit_points`
  ## TODO: (ML) This must (!) be done smoother (UPDATE: at the moment not needed)
  #calc_pit_points_index <- function(group) {
  #  idx <- cumsum(rle(as.numeric(group))$lengths)
  #  idx2 <- lapply(1:length(idx), function(i) c(seq(1, idx[i]), idx[i]))
  #  if (length(idx2) > 1) {
  #    idx3 <- unlist(c(idx2[[1]], sapply(2:length(idx2), function(j) idx2[[j]][idx2[[j]] > max(idx2[[j-1]])])))
  #  } else {
  #    idx3 <- idx2[[1]]
  #  }
  #  return(idx3)
  #}

  ## recycle arguments for plotting to match the number of groups (for `scale_<...>_manual()`)
  plot_arg <- data.frame(
    1:n,
    fill, colour, size, linetype, confint, alpha_min
  )[, -1]

  ## prepare fill and confint color depending on style
  if (style == "histogram") {

    ## set color to 2 (red) or NA for not plotting
    if (is.logical(confint)) confint <- ifelse(confint, 2, NA)

    ## check if colour and no fill is set
    if (all(fill == "darkgray") && any(colour != "black")) {
      message(" * As the argument `colour` is set but no argument `fill` is specified, \n   the former is used for colorizing the PIT histogram. \n * For proper usage, solely provide `fill` for histogram style plots.")
      plot_arg$fill <- plot_arg$colour
    }
  } else {
    if (is.logical(plot_arg$confint)) {

      ## use fill and set alpha
      plot_arg$fill <- sapply(seq_along(plot_arg$fill), function(idx) {
        set_minimum_transparency(plot_arg$fill[idx], alpha_min = plot_arg$alpha_min[idx])
      })

      ## set color to NA for not plotting
      plot_arg$fill[!plot_arg$confint] <- NA
    } else {
      ## use confint and set alpha
      plot_arg$fill <- sapply(seq_along(plot_arg$confint), function(idx) {
        set_minimum_transparency(plot_arg$confint[idx], alpha_min = plot_arg$alpha_min[idx])
      })
    }
  }

  ## FIXME: (ML) Decide if plot_arg2 or plot_arg3 should be used
  ## recycle arguments for plotting to match the length (rows) of the object (for geom w/ aes)
  plot_arg2 <- data.frame(1:n, border, colour, ref, confint)[, -1]
  plot_arg2 <- as.data.frame(lapply(plot_arg2, rep, table(object$group)))

  ## new approach, maybe plot (ref, confint) same for all plots?!
  plot_arg3 <- list("ref" = ref, "confint" = confint)

  # -------------------------------------------------------------------
  # MAIN PLOTTING
  # -------------------------------------------------------------------
  if (style == "histogram") {
    ## actual plotting
    rval <- ggplot2::ggplot(
      object,
      ggplot2::aes_string(x = "x", y = "y / 2", width = "width", height = "y")
    ) +
      ggplot2::geom_tile(ggplot2::aes_string(fill = "group"),
        colour = plot_arg2$border
      ) 

    ## FIXME: (ML) Does not work for object w/ and w/o "ref" (NA)
    if (!all(is.na(object$ref))) {
      rval <- rval + 
        ggplot2::geom_step(ggplot2::aes_string(x = "x", y = "ref", size = "group"),
          colour = plot_arg3$ref, #plot_arg2$ref[calc_pit_points_index(object$group)],
          stat = calc_pit_points,
          na.rm = TRUE
        ) 
    }

    ## FIXME: (ML) Does not work for object w/ and w/o "ci_lwr" (NA)
    if (!all(is.na(object$ci_lwr))) {
      rval <- rval +
        ggplot2::geom_step(ggplot2::aes_string(x = "x", y = "ci_lwr", size = "group"),
          colour = plot_arg3$confint, #plot_arg2$confint[calc_pit_points_index(object$group)],
          linetype = 2,
          stat = calc_pit_points,
          na.rm = TRUE
        )
    }

    ## FIXME: (ML) Does not work for object w/ and w/o "ci_upr" (NA)
    if (!all(is.na(object$ci_upr))) {
      rval <- rval + 
        ggplot2::geom_step(ggplot2::aes_string(x = "x", y = "ci_upr", size = "group"),
          colour = plot_arg3$confint, #plot_arg2$confint[calc_pit_points_index(object$group)],
          linetype = 2,
          stat = calc_pit_points,
          na.rm = TRUE
        ) 
    }

    ## set the colors, shapes, etc. for the groups
    rval <- rval +
      ggplot2::scale_fill_manual(values = plot_arg$fill) +
      ggplot2::scale_size_manual(values = plot_arg$size)

    ## annotation
    rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

    ## add legend
    if (legend) {
      rval <- rval + ggplot2::labs(fill = "Model") +
        ggplot2::guides(fill = "legend", size = "none")
    } else {
      rval <- rval + ggplot2::guides(fill = "none", size = "none")
    }

    ## set x and y limits
    rval <- rval + ggplot2::scale_x_continuous(limits = xlim, expand = c(0.01, 0.01))
    rval <- rval + ggplot2::scale_y_continuous(limits = ylim, expand = c(0.01, 0.01))
  } else {

    ## actual plotting
    rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y", width = "width")) +
      ggplot2::geom_tile(
        ggplot2::aes_string(
          x = "x", 
          y = "(ci_upr + ci_lwr) / 2", 
          width = "width", 
          height = "ci_upr - ci_lwr", 
          fill = "group"
        ),
        colour = NA, show.legend = FALSE, na.rm = TRUE
      ) +
      ggplot2::geom_step(ggplot2::aes_string(colour = "group", size = "group", linetype = "group"),
        stat = calc_pit_points
      )

    ## FIXME: (ML) Does not work for object w/ and w/o "ref" (NA)
    if (!all(is.na(object$ref))) {
      rval <- rval + 
        ggplot2::geom_step(ggplot2::aes_string(x = "x", y = "ref"),
          colour = plot_arg3$ref, #plot_arg2$ref[calc_pit_points_index(object$group)],
          linetype = 2,
          size = 1,
          stat = calc_pit_points,
          na.rm = TRUE
        ) 
    }

    ## set the colors, shapes, etc.
    rval <- rval +
      ggplot2::scale_colour_manual(values = plot_arg$colour) +
      ggplot2::scale_fill_manual(values = plot_arg$fill) +
      ggplot2::scale_size_manual(values = plot_arg$size) +
      ggplot2::scale_linetype_manual(values = plot_arg$linetype)

    ## add annotation
    rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

    ## add legend
    if (legend) {
      rval <- rval + ggplot2::labs(colour = "Model") +
        ggplot2::guides(colour = "legend", size = "none", linetype = "none")
    } else {
      rval <- rval + ggplot2::guides(colour = "none", size = "none", linetype = "none")
    }

    ## set x and y limits
    rval <- rval + ggplot2::scale_x_continuous(limits = xlim, expand = c(0.01, 0.01))
    rval <- rval + ggplot2::scale_y_continuous(limits = ylim, expand = c(0.01, 0.01))
  }

  # -------------------------------------------------------------------
  # GROUPING (IF ANY) AND RETURN PLOT
  # -------------------------------------------------------------------
  ## grouping
  if (!single_graph && n > 1L) {
    rval <- rval + ggplot2::facet_grid(group ~ .)
  } else if (!is.null(object$title)) {
    rval <- rval + ggplot2::facet_wrap(title ~ .)
  }

  ## return ggplot object
  rval
}

get_pp <- function(n, breaks, freq) {
  ## helper function to calculate perfect prediciton employing `qbinom()`

  
  ## calc bin specific confidence levels
  rval <- qbinom(0.5, size = n, prob = diff(breaks))

  ## transform counts to frequency
  if (!freq) rval <- rval / (n / (length(breaks) - 1))
  rval
}

get_confint <- function(n, breaks, level, freq) {
  ## helper function to calculate CI employing `qbinom()`

  ## get confidence level
  a <- (1 - level) / 2
  
  ## calc bin specific confidence levels
  rval <- list(
    qbinom(a, size = n, prob = diff(breaks)),
    qbinom(1 - a, size = n, prob = diff(breaks))
  )

  ## transform counts to frequency
  if (!freq) rval <- lapply(rval, function(x) x / (n / (length(breaks) - 1)))
  rval
}


get_confint_agresti <- function(x, n, level, bins, freq) {
  ## helper function to calculate an approximated CI according to Agresti & Coull (1998)
  ## doi=10.1080/00031305.1998.10480550
  rval <- add4ci(x, n, level)$conf.int * n
  if (!freq) rval <- rval / (n / bins)
  rval
}


add4ci <- function(x, n, conf.level) {
  ## copy of `add4ci` package from package `PropCIs` by Ralph Scherer (licensed under GPL-2/GPL-3)
  ptilde <- (x + 2) / (n + 4)
  z <- abs(qnorm((1 - conf.level) / 2))
  stderr <- sqrt(ptilde * (1 - ptilde) / (n + 4))
  ul <- ptilde + z * stderr
  ll <- ptilde - z * stderr
  if (ll < 0) {
    ll <- 0
  }
  if (ul > 1) {
    ul <- 1
  }
  cint <- c(ll, ul)
  attr(cint, "conf.level") <- conf.level
  rval <- list(conf.int = cint, estimate = ptilde)
  class(rval) <- "htest"
  return(rval)
}


#' \code{geom_*} and \code{stat_*} for Producing PIT Histograms with `ggplot2`
#' 
#' Various \code{geom_*} and \code{stat_*} used within
#' \code{\link[ggplot2]{autoplot}} for producing PIT histograms.
#' 
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#' @param style Fix description.
#' @param linejoin Fix description.
#' @examples
#' require("ggplot2")
#' ## Fit model
#' data("CrabSatellites", package = "countreg")
#' m1_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#' m2_pois <- glm(satellites ~ color, data = CrabSatellites, family = poisson)
#' 
#' ## Compute pithist
#' p1 <- pithist(m1_pois, plot = FALSE)
#' p2 <- pithist(m2_pois, plot = FALSE)
#' 
#' d <- c(p1, p2) 
#' 
#' ## Get label names
#' xlab <- unique(attr(d, "xlab"))
#' ylab <- unique(attr(d, "ylab"))
#' main <- attr(d, "main")
#' main <- make.names(main, unique = TRUE)
#' d$group <- factor(d$group, labels = main)
#' 
#' gg1 <- ggplot(data = d) + 
#'   geom_pit_line(aes(x = x, y = y, width = width, group = group)) + 
#'   geom_pit_confint(aes(x = x, ci_upr = ci_upr, ci_lwr = ci_lwr, width = width)) + 
#'   facet_grid(group~.)
#' gg1
#' @export
geom_pit_hist <- function(mapping = NULL, data = NULL, stat = "identity",
                          position = "identity", na.rm = FALSE,
                          show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    geom = GeomPitHist, mapping = mapping,
    data = data, stat = stat, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' @rdname geom_pit_hist
#' @format NULL
#' @usage NULL
#' @export
GeomPitHist <- ggplot2::ggproto("GeomPitHist", ggplot2::GeomTile,

  default_aes = ggplot2::aes(colour = "black", fill = "darkgray", size = 0.2, linetype = 1, 
    alpha = NA),

  required_aes = c("x", "y", "width"),

  setup_data = function(data, params) {
      data <- transform(data,
        y = y /2,
        height = y
      )
      ggplot2::GeomTile$setup_data(data, params)
  }
)


#' @rdname geom_pit_hist
#' @format NULL
#' @usage NULL
#' @export
geom_pit_line <- function(mapping = NULL, data = NULL, stat = "pit_line",
                            position = "identity", na.rm = FALSE,
                            show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    geom = GeomPitLine, mapping = mapping,
    data = data, stat = stat, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' @rdname geom_pit_hist
#' @format NULL
#' @usage NULL
#' @export
GeomPitLine <- ggplot2::ggproto("GeomPitLine", ggplot2::GeomStep,
  default_aes = ggplot2::aes(colour = "black", size = 1, linetype = 1,
  alpha = NA),
  required_aes = c("x", "y")
)


#' @rdname geom_pit_hist
#' @export
stat_pit_line <- function(mapping = NULL, data = NULL, geom = "pit_line",
                         position = "identity", na.rm = FALSE,
                         show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatPitLine,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}


#' @rdname geom_pit_hist
#' @format NULL
#' @usage NULL
#' @export
StatPitLine <- ggplot2::ggproto("StatPitLine", ggplot2::Stat,

  compute_group = function(data, scales) {
    nd <- data.frame(
      x = c(data$x - data$width / 2, data$x[NROW(data)] + data$width[NROW(data)] / 2),
      y = c(data$y, data$y[NROW(data)])
    )
    nd
  },

  required_aes = c("x", "y", "width")
)


#' @rdname geom_pit_hist
#' @format NULL
#' @usage NULL
#' @export
geom_pit_ref <- function(mapping = NULL, data = NULL, stat = "pit_line",
                            position = "identity", na.rm = FALSE,
                            show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    geom = GeomPitRef, mapping = mapping,
    data = data, stat = stat, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' @rdname geom_pit_hist
#' @format NULL
#' @usage NULL
#' @export
GeomPitRef <- ggplot2::ggproto("GeomPitRef", ggplot2::GeomStep,
  default_aes = ggplot2::aes(colour = 2, size = 0.75, linetype = 1,
  alpha = NA),
  required_aes = c("x", "y")
)


#' @rdname geom_pit_hist
#' @export
stat_pit_ref <- function(mapping = NULL, data = NULL, geom = "pit_ref",
                         position = "identity", na.rm = FALSE,
                         show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatPitLine,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}


#' @rdname geom_pit_hist
#' @export
stat_pit_confint <- function(mapping = NULL, data = NULL, geom = "pit_confint",
                         position = "identity", na.rm = FALSE,
                         show.legend = NA, inherit.aes = TRUE, 
                         style = c("polygon", "line"), ...) {

  style <- match.arg(style)

  ggplot2::layer(
    stat = StatPitConfint,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      style = style,
      na.rm = na.rm,
      ...
    )
  )
}


#' @rdname geom_pit_hist
#' @format NULL
#' @usage NULL
#' @export
StatPitConfint <- ggplot2::ggproto("StatPitConfint", ggplot2::Stat,

  compute_group = function(data, scales, style = "polygon") {

    if (style == "polygon") {
      nd <- data.frame(
        xmin = data$x - data$width / 2,
        xmax = data$x + data$width / 2,
        ymin = data$ci_lwr,
        ymax = data$ci_upr,
        x = NaN,
        ci_lwr = NaN, 
        ci_upr = NaN
      )
    } else {
      nd <- data.frame(
        x = c(data$x - data$width / 2, data$x[NROW(data)] + data$width[NROW(data)] / 2),
        ci_lwr = c(data$ci_lwr, data$ci_lwr[NROW(data)]),
        ci_upr = c(data$ci_upr, data$ci_upr[NROW(data)]),
        xmin = NaN,
        xmax = NaN, 
        ymin = NaN, 
        ymax = NaN
      )
    }
    nd
  },

  required_aes = c("x", "ci_upr", "ci_lwr", "width")
)


#' @rdname geom_pit_hist
#' @export
geom_pit_confint <- function(mapping = NULL, data = NULL, stat = "pit_confint",
                            position = "identity", na.rm = FALSE,
                            show.legend = NA, inherit.aes = TRUE,
                            linejoin = "mitre", style = c("polygon", "line"), ...) {
  style <- match.arg(style)

  ggplot2::layer(
    geom = GeomPitConfint,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      linejoin = linejoin,
      na.rm = na.rm,
      style = style,
      ...
    )
  )
}


#' @rdname geom_pit_hist
#' @export
GeomPitConfint <- ggplot2::ggproto("GeomPitConfint", ggplot2::Geom,

  # FIXME: (ML) Does not vary for style
  required_aes = c("x", "ci_lwr", "ci_upr", "xmin", "xmax", "ymin", "ymax"),

  extra_params = c("na.rm", "linejoin", "style"),

  # FIXME: (ML) Does not vary for style; this is a copy of `GeomPolygon$handle_na()`
  handle_na = function(data, params) {
    data
  },

  ## Setting up all defaults needed for `GeomPolygon` and `GeomStep`
  default_aes = ggplot2::aes(
    colour = NA,
    fill = NA,
    size = 0.5,
    linetype = NA,
    alpha = NA
  ),


  draw_panel = function(data, panel_params, coord, 
                        linejoin = "mitre", direction = "hv",
                        style = c("polygon", "line")) {

    style <- match.arg(style)

    ## Swap NAs in `default_aes` with own defaults 
    data <- my_modify_list(data, pit_default_aesthetics(style), force = FALSE)

    if (style == "polygon") {
      ggplot2::GeomRect$draw_panel(data, panel_params, coord, linejoin)

    } else { 
      ## Join two Grobs
      data1 <- transform(data, 
        y = ci_upr
      )
      data2 <- transform(data, 
        y = ci_lwr
      )
      grid::grobTree(
        ggplot2::GeomStep$draw_panel(data1, panel_params, coord, direction),
        ggplot2::GeomStep$draw_panel(data2, panel_params, coord, direction)
      )

    }
  },


  draw_key = function(data, params, size) {
    ## Swap NAs in `default_aes` with own defaults 
    data <- my_modify_list(data, pit_default_aesthetics(params$style), force = FALSE)
    if (params$style == "polygon") {
      draw_key_polygon(data, params, size)
    } else {
      draw_key_path(data, params, size)
    }
  }

)


## Helper function inspired by internal from `ggplot2` defined in `geom-sf.R`
pit_default_aesthetics <- function(style) {
  if (style == "line") {
    my_modify_list(ggplot2::GeomPath$default_aes, list(colour = 2, size = 0.75, linetype = 2, alpha = NA),
      force = TRUE)
  } else {
    my_modify_list(ggplot2::GeomPolygon$default_aes, list(colour = "NA", fill = "black", size = 0.5, 
      linetype = 1, alpha = 0.2, subgroup = NULL), force = TRUE)
  }
}
