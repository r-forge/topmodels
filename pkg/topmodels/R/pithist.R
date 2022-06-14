# -------------------------------------------------------------------
# Programming outline: PIT histogram
# -------------------------------------------------------------------
# - Observed y in-sample or out-of-sample (n x 1)
# - Predicted probabilities F_y(y - eps) and F_y(y) (n x 2)
#   eps is a tiny numeric value required for discrete/censored distributions
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
#
# Functions:
# - pithist() generic plus default method
# - Return object of class "pithist" that is plotted by default
# - But has plot=FALSE so that suitable methods can be added afterwards
# - At least methods: plot(), autoplot(), lines(), possibly c(), +
# -------------------------------------------------------------------


#' PIT Histograms for Assessing Goodness of Fit of Probability Models
#'
#' Probability integral transform (PIT) histograms graphically
#' compare empirical probabilities from fitted models
#' with a uniform distribution. If \code{plot = TRUE}, the resulting object of
#' class \code{"pithist"} is plotted by \code{\link{plot.pithist}} or
#' \code{\link{autoplot.pithist}} depending on whether the
#' package \code{ggplot2} is loaded, before the \code{"pithist"} object is returned.
#'
#' PIT histograms graphically evaluate the probability integral transform (PIT),
#' i.e., the value that the predictive CDF attains at the observation, with a
#' uniform distribution. For a well calibrated model fit, the PIT will have a
#' standard uniform distribution.
#' For computation, \code{\link{pithist}} leverages the function
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
#'
#' @param object an object from which probability integral transforms can be
#' extracted using the generic function \code{\link{procast}}.
#' @param newdata an optional data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param plot logical or character. Should the \code{plot} or \code{autoplot}
#' method be called to draw the computed extended reliability diagram? Logical
#' \code{FALSE} will suppress plotting, \code{TRUE} (default) will choose the
#' type of plot conditional if the package \code{ggplot2} is loaded.
#' Alternatively \code{"base"} or \code{"ggplot2"} can be specified to
#' explicitly choose the type of plot.
#' @param class should the invisible return value be either a \code{data.frame}
#' or a \code{tbl_df}. Can be set to \code{"data.frame"} or \code{"tibble"} to
#' explicitly specify the return class, or to \code{NULL} (default) in which
#' case the return class is conditional on whether the package \code{"tibble"}
#' is loaded.
#' @param scale controls the scale on which the PIT residuals are computed: on
#' the probability scale (\code{"uniform"}; default) or on the normal scale
#' (\code{"normal"}).
#' @param breaks \code{NULL} (default) or numeric to manually specify the breaks for
#' the rootogram intervals. A single numeric (larger \code{0}) specifies the number of breaks
#' to be automatically chosen, multiple numeric values are interpreted as manually specified breaks.
#' @param type character. In case of discrete distributions, should an expected
#' (non-normal) PIT histogram be computed according to Czado et al. (2009)
#' (\code{"expected"}; default) or should the PIT be drawn randomly from the corresponding
#' interval (\code{"random"})?
#' @param nsim positive integer, defaults to \code{1L}. Only used when
#' \code{type = "random"}; how many simulated PITs should be drawn?
#' @param delta \code{NULL} or numeric. The minimal difference to compute the range of
#' probabilities corresponding to each observation to get (randomized)
#' quantile residuals. For \code{NULL} (default), the minimal observed difference in the
#' response divided by \code{5e-6} is used.
#' @param simint \code{NULL} (default) or logical. In case of discrete
#' distributions, should the simulation (confidence) interval due to the
#' randomization be visualized?
#' @param simint_level numeric, defaults to \code{0.95}. The confidence level
#' required for calculating the simulation (confidence) interval due to the
#' randomization.
#' @param simint_nrep numeric, defaults to \code{250}. The repetition number of
#' simulated quantiles for calculating the simulation (confidence) interval due
#' to the randomization.
#' @param style character specifying plotting style. For \code{style = "bar"} (default)
#' a traditional PIT histogram is drawn, \code{style = "line"} solely plots the upper border
#' of the bars. If \code{single_graph = TRUE} is used (see \code{\link{plot.pithist}}),
#' line-style PIT histograms will be enforced.
#' @param freq logical. If \code{TRUE}, the PIT histogram is represented by
#' frequencies, the \code{counts} component of the result; if \code{FALSE},
#' probability densities, component \code{density}, are plotted (so that the
#' histogram has a total area of one).
#' @param expected logical. Should the expected values be plotted as reference?
#' @param confint logical. Should confident intervals be drawn?
#' @param xlab,ylab,main graphical parameters passed to
#' \code{\link{plot.pithist}} or \code{\link{autoplot.pithist}}.
#' @param \dots further graphical parameters forwarded to the plotting functions.
#'
#' @return An object of class \code{"pithist"} inheriting from
#' \code{data.frame} or \code{tbl_df} conditional on the argument \code{class}
#' including the following variables:
#' \item{x}{histogram interval midpoints on the x-axis,}
#' \item{y}{bottom coordinate of the histogram bars,}
#' \item{width}{widths of the histogram bars,}
#' \item{confint_lwr}{lower bound of the confidence interval,}
#' \item{confint_upr}{upper bound of the confidence interval,}
#' \item{expected}{y-coordinate of the expected curve.}
#' Additionally, \code{freq}, \code{xlab}, \code{ylab}, \code{main}, and
#' \code{confint_level} are stored as attributes.
#'
#' @seealso \code{\link{plot.pithist}}, \code{\link{qresiduals}}, \code{\link{procast}}
#'
#' @references
#' Agresti A, Coull AB (1998). \dQuote{Approximate is Better than ``Exact''
#' for Interval Estimation of Binomial Proportions.} \emph{The American
#' Statistician}, \bold{52}(2), 119--126. \doi{10.1080/00031305.1998.10480550}
#'
#' Czado C, Gneiting T, Held L (2009). \dQuote{Predictive Model
#' Assessment for Count Data.} \emph{Biometrics}, \bold{65}(4), 1254--1261.
#' \url{https://www.jstor.org/stable/20640646}.
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
#' Series B (Statistical Methodology)}. \bold{69}(2), 243--268.
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
#' ggplot2::autoplot(c(p1, p2), single_graph = TRUE, style = "line", col = c(1, 2))
#' @export
pithist <- function(object, ...) {
  UseMethod("pithist")
}


#' @rdname pithist
#' @method pithist default
#' @export
pithist.default <- function( 
                            ## computation arguments
                            object,
                            newdata = NULL,
                            plot = TRUE,
                            class = NULL,
                            scale = c("uniform", "normal"),
                            breaks = NULL,
                            type = c("expected", "random"),
                            nsim = 1L,
                            delta = NULL,
                            simint = NULL, # needed also for plotting
                            simint_level = 0.95,
                            simint_nrep = 250,

                            ## plotting arguments
                            style = c("bar", "line"),
                            freq = FALSE,
                            expected = TRUE,
                            confint = TRUE,
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
  ## * `expected`, `confint`, `...` w/i `plot()` and `autoplot()`
  stopifnot(is.null(breaks) || (is.numeric(breaks) && length(breaks) > 0 && is.null(dim(breaks))))
  if (length(breaks) == 1) stopifnot(breaks >= 1)
  stopifnot(is.null(simint) || isTRUE(simint) || isFALSE(simint))
  stopifnot(is.numeric(simint_level), length(simint_level) == 1, simint_level >= 0, simint_level <= 1)
  stopifnot(is.numeric(simint_nrep), length(simint_nrep) == 1, simint_nrep >= 1)
  stopifnot(isTRUE(freq) || isFALSE(freq))
  stopifnot(is.character(xlab), length(xlab) == 1)
  stopifnot(is.character(ylab), length(ylab) == 1)
  stopifnot(is.null(main) || (length(main) == 1 && is.character(main)))
  stopifnot(isTRUE(plot) || isFALSE(plot) || (is.character(plot) && length(plot) == 1L))
  stopifnot(is.null(class) || (is.character(class) && length(class) == 1L))
  stopifnot(isTRUE(expected) || isFALSE(expected))
  stopifnot(isTRUE(confint) || isFALSE(confint))

  ## match arguments
  style <- match.arg(style)
  type <- match.arg(type)
  scale <- match.arg(scale)

  ## determine other arguments conditional on `style` and `type`
  if (is.null(simint)) {
    simint <- if (style == "bar" && type == "random") TRUE else FALSE
  }

  ## default annotation
  if (is.null(main)) {
    main <- deparse(substitute(object))
  }

  ## guess plotting flavor
  if (is.logical(plot)) {
      plot <- ifelse(isFALSE(plot), "none", if ("ggplot2" %in% .packages()) "ggplot2" else "base")
  }
  plot <- try(match.arg(plot, c("none", "base", "ggplot2")), silent = TRUE)
  if (inherits(plot, "try-error"))
    stop("`plot` must be logical `TRUE`/`FALSE` or one of \"none\", \"base\", or \"ggplot2\"")

  ## guess output class
  if (is.null(class)) {
    class <- if ("tibble" %in% .packages()) "tibble" else "data.frame"
  }
  class <- try(match.arg(class, c("tibble", "data.frame")), silent = TRUE)
  if (inherits(class, "try-error"))
    stop("`class` must be `NULL` or one of \"tibble\", \"data.frame\"")

  # -------------------------------------------------------------------
  # COMPUTATION OF PIT
  # -------------------------------------------------------------------
  ## get breaks: part 1 (due to "type = expected" must be done before)
  n <- NROW(newresponse(object, newdata = newdata)) # solely to get n for computing breaks
  if (is.null(breaks)) breaks <- c(4, 10, 20, 25)[cut(n, c(0, 50, 5000, 1000000, Inf))]

  # -------------------------------------------------------------------
  if (type == "proportional") {
    # -------------------------------------------------------------------
    ## TODO: (ML)
    ## * implement proportional over the inteverals (e.g., below censoring point)
    ## * confusing naming, as `type` in `qresiduals()` must be `random` or `quantile`
    ## NOTE: (RS) Cannot be tested as we currently don't allow type = "proportional"
    ##       (see @param; match.arg(type)).
    stop("not yet implemented")

    # -------------------------------------------------------------------
  } else if (type == "random") {
    # -------------------------------------------------------------------
    p <- qresiduals(object,
      newdata = newdata, scale = scale, delta = delta,
      type = "random", nsim = nsim
    )

    ## get breaks: part 2
    ## TODO: (ML) maybe use xlim instead or `0` and `1`
    if (scale == "uniform") {
      if (length(breaks) == 1L) breaks <- seq(0, 1, length.out = breaks + 1L)
    } else {
      tmp_range <- range(p, finite = TRUE)
      if (length(breaks) == 1L) breaks <- seq(tmp_range[1], tmp_range[2], length.out = breaks + 1L)
    }

    ## collect everything as data.frame
    ## TODO: (ML) Maybe get rid of `hist()`
    tmp_rval <- hist(p, breaks = breaks, plot = FALSE)
    tmp_rval$counts <- tmp_rval$counts / nsim

    if (!isFALSE(simint)) {
      ## helper function for calculating simulation interval
      compute_simint <- function(object, newdata, scale, delta, nsim, breaks) {
        simint_p <- qresiduals(object,
          newdata = newdata, scale = scale, delta = delta,
          type = "random", nsim = nsim
        )

        ## TODO: (ML) Maybe nicer workaround to get not values outside breaks
        simint_p[simint_p < min(breaks)] <- min(breaks)
        simint_p[simint_p > max(breaks)] <- max(breaks)
        simint_hist <- hist(simint_p, breaks = breaks, plot = FALSE)

        return(simint_hist$counts / nsim)
      }

      ## compute simulation interval bases on quantiles
      simint_tmp <- replicate(
        simint_nrep,
        compute_simint(object, newdata, scale, delta, nsim, breaks)
      )

      simint_prob <- (1 - simint_level) / 2
      simint_prob <- c(simint_prob, 1 - simint_prob)
      simint_lwr <- apply(simint_tmp, 1, quantile, probs = simint_prob[1], na.rm = TRUE)
      simint_upr <- apply(simint_tmp, 1, quantile, probs = simint_prob[2], na.rm = TRUE)
    }

    # -------------------------------------------------------------------
  } else if (type == "expected") {
    # -------------------------------------------------------------------
    ## compare "nonrandom" in Czado et al. (2009)

    ## minimum and maximum PIT for each observation (P_x-1 and P_x)
    p <- qresiduals(object,
      newdata = newdata, scale = "uniform", delta = delta,
      type = "quantile", prob = c(0, 1)
    )

    ## TODO: (ML) Adapt by employing `distributions3`.
    if (scale == "uniform") {
      qFun <- identity
      pFun <- punif
    } else {
      qFun <- qnorm
      pFun <- pnorm
    }

    ## equation 2: CDF for each PIT (continuous vs. discrete)
    F <- if (all(abs(p[, 2L] - p[, 1L]) < sqrt(.Machine$double.eps))) {
      function(u) as.numeric(u >= p[, 1L])
      ## TODO: (Z) Check inequality sign to cover include.lowest/right options.
    } else {
      #function(u) punif(u, min = p[, 1L], max = p[, 2L]) # original 
      #function(u) pmin(1, pmax(0, (u - p[, 1L]) / (p[, 2L] - p[, 1L]))) # equal to original
      function(u) pmin(1, pmax(0, 1/(p[, 2L] - p[, 1L]) * (pFun(u) - p[, 1L]))) # new working w trafo
    }

    ## get breaks: part 2
    ## TODO: (ML) Improve this educated guess.
    tmp_range <- range(
      c(
        qFun(0), qFun(0.000001), qFun(0.0001), qFun(0.001), qFun(0.005), qFun(0.01), 
        qFun(0.99), qFun(0.995), qFun(0.999), qFun(0.9999), qFun(0.999999), qFun(1)
      ), 
      finite = TRUE
    )
    if (length(breaks) == 1L) breaks <- seq(tmp_range[1], tmp_range[2], length.out = breaks + 1L)

    ## equation 3 and computation of probability for each interval (f_j)
    ## TODO: (RS2ML) pure diff loses 'point mass' on 0; code adjusted to
    ##       account for possible point mass (probability on breaks[1])
    tmp <- colMeans(sapply(breaks, F))
    f   <- diff(tmp) + c(tmp[1], rep(0, length(tmp) - 2L)) #diff(colMeans(sapply(breaks, F)))
    rm(tmp)

    ## collect everything as data.frame
    tmp_rval <- list(
      breaks = breaks,
      counts = f * n,
      density = f / diff(breaks),
      mid = 0.5 * (breaks[-1L] + breaks[-length(breaks)])
    )
  }

  ## compute expected line (freq = TRUE, as scaling below)
  val_expected <- compute_pithist_expected(
    n = n,
    breaks = tmp_rval$breaks,
    freq = TRUE,
    scale = scale
  )

  # -------------------------------------------------------------------
  # OUTPUT AND OPTIONAL PLOTTING
  # -------------------------------------------------------------------
  ## collect everything as data.frame
  rval <- data.frame(
    observed = tmp_rval$counts,
    expected = val_expected,
    mid = tmp_rval$mid,
    width = diff(tmp_rval$breaks)
  )

  ## add simint
  if (!isFALSE(simint) && type == "random") {
    ## TODO: (ML) Check if simint is correct.
    ## TODO: (RS2ML) It is possible that the resulting simint_lwr
    ##       is getting negative; we have one test case in `test_pithist_04_simint.R`.
    ##       Just pmax(0, rval$simint_lwr)?
    rval$simint_lwr <- rval$observed - (simint_upr - simint_lwr) / 2
    rval$simint_upr <- rval$observed + (simint_upr - simint_lwr) / 2
  } else {
    rval$simint_lwr <- NA
    rval$simint_upr <- NA
  }

  ## not allow freq = TRUE for non-equidistant breaks
  ## TODO: (ML) Round values to get unique widths?
  if (freq && any(abs(diff(rval$width)) > sqrt(.Machine$double.eps))) {
    warning(
      "For non-equidistant breaks `freq = FALSE` must be used and has been set accordingly.", 
      call. = FALSE
    )
    freq <- FALSE
  }

  ## optional return counts
  if (!freq) {
    rval <- transform(rval,
      observed = rval$observed / (n * rval$width),
      expected = rval$expected / (n * rval$width),
      simint_lwr = rval$simint_lwr / (n * rval$width),
      simint_upr = rval$simint_upr / (n * rval$width)
    )
  }

  ## add attributes
  attr(rval, "scale") <- scale
  attr(rval, "type") <- type
  attr(rval, "simint") <- simint
  attr(rval, "style") <- style
  attr(rval, "freq") <- freq
  attr(rval, "expected") <- expected
  attr(rval, "confint") <- confint
  attr(rval, "counts") <- n
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main

  ## set class to data.frame or tibble
  if (class == "data.frame") {
    class(rval) <- c("pithist", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("pithist", class(rval))
  }

  ## plot by default
  if (plot == "ggplot2") {
    try(print(ggplot2::autoplot(rval, ...)))
  } else if (plot == "base") {
    try(plot(rval, ...))
  }

  ## return invisibly
  invisible(rval)
}


#' @export
#' @method c pithist
c.pithist <- function(...) {
  ## TODO: (RS2ML) The method is listed nowhere and not documented.
  ##               Should there be an entry for c/rbind pithist with
  ##               a description? Furthermore, no sanity checks are
  ##               carried out. Does it make sense to check the
  ##               rvals?
  ##               All allowed is (named or unnamed) series of
  ##               pithist objects, right?
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

  ## remove temporary the class (needed below for `c()`)
  rval <- lapply(rval, function(x) structure(x, class = class(x)[-1]))

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
  ## labels
  xlab <- prepare_arg_for_attributes(rval, "xlab")
  ylab <- prepare_arg_for_attributes(rval, "ylab")
  nam  <- names(rval)
  main <- if (is.null(nam)) {
    prepare_arg_for_attributes(rval, "main")
  } else {
    make.unique(rep.int(nam, sapply(n, length)))
  }

  ## parameters
  scale    <- prepare_arg_for_attributes(rval, "scale",    force_single = FALSE) # check/force below
  type     <- prepare_arg_for_attributes(rval, "type",     force_single = FALSE)
  simint   <- prepare_arg_for_attributes(rval, "simint")
  style    <- prepare_arg_for_attributes(rval, "style",    force_single = TRUE)
  freq     <- prepare_arg_for_attributes(rval, "freq",     force_single = TRUE)
  expected <- prepare_arg_for_attributes(rval, "expected")
  confint  <- prepare_arg_for_attributes(rval, "confint")
  counts   <- prepare_arg_for_attributes(rval, "counts")

  ## fix `ylabel` according to possible new `freq`
  if (isTRUE(freq)) {
    ylab[grepl("^Density$", ylab)] <- "Frequency"
  } else if (isFALSE(freq)) {
    ylab[grepl("^Frequency$", ylab)] <- "Density"
  }

  # -------------------------------------------------------------------
  # CHECK FOR COMPATIBILITY
  # -------------------------------------------------------------------
  if (length(unique(scale)) > 1) {
    stop("Can't combine pit histograms which are on different scales.")
  } else {
    scale <- scale[[1]]
  }

  # -------------------------------------------------------------------
  # RETURN DATA
  # -------------------------------------------------------------------
  ## get all names needed if extended object should be computed and for combination
  all_names <- unique(unlist(lapply(rval, names)))

  ## get both objects on the same scale
  if (any(grepl("confint_lwr|confint_upr|expected", all_names))) {
    rval <- lapply(rval, summary.pithist, freq = freq, extend = TRUE)
  } else {
    ## TODO: (RS2ML): Is this case even possible?
    rval <- lapply(rval, summary.pithist, freq = freq, extend = FALSE)
  }
  rval <- lapply(rval, as.data.frame) # remove inner class

  ## combine and return (fill up missing variables with NAs)
  rval <- do.call(
    "rbind.data.frame",
    c(lapply(
      rval,
      function(x) {
        data.frame(c(x, sapply(setdiff(all_names, names(x)), function(y) NA)))
      }
    ),
    make.row.names = FALSE
    )
  )

  ## add group
  n <- unlist(n)
  rval$group <- if (length(n) < 2L) NULL else rep.int(seq_along(n), n)

  ## re-order, make sure 'group' is the last variable such that
  ## we get identical object when c-binding objects
  ## returned by pithist() or from summary(pithist())
  rval <- rval[, c(names(rval)[names(rval) != "group"], "group")]

  ## add attributes
  attr(rval, "scale")    <- scale
  attr(rval, "type")     <- type
  attr(rval, "simint")   <- simint
  attr(rval, "style")    <- style
  attr(rval, "freq")     <- freq
  attr(rval, "expected") <- expected
  attr(rval, "confint")  <- confint
  attr(rval, "counts")   <- counts
  attr(rval, "xlab")     <- xlab
  attr(rval, "ylab")     <- ylab
  attr(rval, "main")     <- main

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
#' @method rbind pithist
rbind.pithist <- c.pithist


#' S3 Methods for Plotting PIT Histograms
#'
#' Generic plotting functions for probability integral transform (PIT)
#' histograms of the class \code{"pithist"} computed by \code{link{pithist}}.
#'
#' PIT histograms evaluates the predictive cumulative distribution
#' function (CDF) evaluated at the observation and compares the resulting
#' values to a uniform distribution or normal distribution.
#' For a well calibrated model fit, the distribution of the PIT residuals
#' will show a standard uniform distrbution or normal distribution depending
#' on the scale selected by the user.
#'
#' PIT histograms can be rendered as \code{ggplot2} or base R graphics by using
#' the generics \code{\link[ggplot2]{autoplot}} or \code{\link[graphics]{plot}}.
#' For a single base R graphically panel, \code{\link{lines}} adds an additional PIT histogram.
#'
#' @aliases plot.pithist lines.pithist autoplot.pithist
#'
#' @param object,x an object of class \code{\link{pithist}}.
#' @param single_graph logical. Should all computed extended reliability
#' diagrams be plotted in a single graph? If yes, \code{style} must be set to \code{"line"}.
#' @param style \code{NULL} or character specifying the style of pithist. For
#' \code{style = "bar"} a traditional PIT hisogram is drawn, for \code{style =
#' "line"} solely the upper border line is plotted.  \code{single_graph = TRUE}
#' always results in a combined line-style PIT histogram.
#' @param freq \code{NULL} or logical.
#' \code{TRUE} will enforce the PIT to be represented by frequencies (counts) while
#' \code{FALSE} will enforce densities.
#' @param expected logical. Should the expected values be plotted as reference?
#' @param confint \code{NULL} or logical. Should confident intervals be drawn? Either logical or as 
#' @param confint_level numeric in \code{[0, 1]}. The confidence level to be shown.
#' @param confint_type character. Which type of confidence interval should be
#' plotted: `"exact"` or `"approximation"`. According to Agresti and Coull
#' (1998), for interval estimation of binomial proportions an approximation can
#' be better than exact.
#' @param simint \code{NULL} or logical. In case of discrete distributions, should the simulation
#' (confidence) interval due to the randomization be visualized?
#' character string defining one of `"polygon"`, `"line"` or `"none"`.
#' If \code{freq = NULL} it is taken from the \code{object}.
#' @param xlim,ylim,xlab,ylab,main,axes,box graphical parameters.
#' @param col,border,lwd,lty,alpha_min graphical parameters for the main part of the base plot.
#' @param colour,fill,size,linetype,alpha graphical parameters for the histogram style part in the \code{autoplot}.
#' @param legend logical. Should a legend be added in the \code{ggplot2} style
#' graphic?
#' @param theme Which `ggplot2` theme should be used. If not set, \code{\link[ggplot2]{theme_bw}} is employed.
#' @param simint_col,simint_lty,simint_lwd,confint_col,confint_lty,confint_lwd,confint_alpha,expected_col,expected_lty,expected_lwd Further graphical parameters for the `confint` and `simint` line/polygon in the base plot.
#' @param simint_colour,simint_size,simint_linetype,simint_alpha,confint_colour,confint_fill,confint_size,confint_linetype,expected_colour,expected_size,expected_linetype,expected_alpha Further graphical parameters for the `confint` and `simint` line/polygon using \code{\link[ggplot2]{autoplot}}.
#' @param \dots further graphical parameters passed to the plotting function.
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
#' pithist(m1_lm, expected_col = "blue", lty = 2, pch = 20, style = "line")
#'
#' ## add separate model
#' if (require("crch", quietly = TRUE)) {
#'   m1_crch <- crch(dist ~ speed | speed, data = cars)
#'   #lines(pithist(m1_crch, plot = FALSE), col = 2, lty = 2, confint_col = 2) #FIXME
#' }
#'
#' #-------------------------------------------------------------------------------
#' if (require("crch")) {
#'
#'   ## precipitation observations and forecasts for Innsbruck
#'   data("RainIbk", package = "crch")
#'   RainIbk <- sqrt(RainIbk)
#'   RainIbk$ensmean <- apply(RainIbk[, grep("^rainfc", names(RainIbk))], 1, mean)
#'   RainIbk$enssd <- apply(RainIbk[, grep("^rainfc", names(RainIbk))], 1, sd)
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
#'   ## plot in single graph with style "line"
#'   plot(c(pit2_lm, pit2_crch),
#'     col = c(1, 2), confint_col = c(1, 2), expected_col = 3,
#'     style = "line", single_graph = TRUE
#'   )
#' }
#'
#' #-------------------------------------------------------------------------------
#' ## determinants for male satellites to nesting horseshoe crabs
#' data("CrabSatellites", package = "countreg")
#'
#' ## linear poisson model
#' m3_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#'
#' ## compute and plot pithist as "ggplot2" graphic
#' pithist(m3_pois, plot = "ggplot2")
#' @export
plot.pithist <- function(x,
                         single_graph = FALSE,
                         style = NULL,
                         freq = NULL,
                         expected = TRUE,
                         confint = NULL,
                         confint_level = 0.95,
                         confint_type = c("exact", "approximation"),
                         simint = NULL,
                         xlim = c(NA, NA),
                         ylim = c(0, NA),
                         xlab = NULL,
                         ylab = NULL,
                         main = NULL,
                         axes = TRUE,
                         box = TRUE,
                         col = "black",
                         border = "black",
                         lwd = NULL,
                         lty = 1,
                         alpha_min = 0.2,
                         expected_col = NULL,
                         expected_lty = NULL,
                         expected_lwd = 1.75,
                         confint_col = NULL,
                         confint_lty = 2,
                         confint_lwd = 1.75,
                         confint_alpha = NULL,
                         simint_col = "black",
                         simint_lty = 1,
                         simint_lwd = 1.75,
                         ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## check if ylab is defined
  ylab_missing <- missing(ylab)

  ## get default arguments
  type     <- use_arg_from_attributes(x, "type",     default = "expected", force_single = FALSE)
  style    <- use_arg_from_attributes(x, "style",    default = "bar",      force_single = TRUE)
  freq     <- use_arg_from_attributes(x, "freq",     default = FALSE,      force_single = TRUE)
  expected <- use_arg_from_attributes(x, "expected", default = TRUE,       force_single = TRUE)
  confint  <- use_arg_from_attributes(x, "confint",  default = TRUE,       force_single = TRUE)
  simint   <- use_arg_from_attributes(x, "simint",   default = NULL,       force_single = TRUE)

  ## sanity checks
  ## * lengths of most arguments are checked by recycling
  ## * graphical parameters are checked w/i function calls for plotting
  stopifnot(isTRUE(single_graph) || isFALSE(single_graph))
  stopifnot(isTRUE(freq) || isFALSE(freq))
  stopifnot(isTRUE(expected) || isFALSE(expected))
  stopifnot(isTRUE(confint) || isFALSE(confint) || confint %in% c("polygon", "line", "none"))
  stopifnot(
    is.numeric(confint_level),
    length(confint_level) == 1,
    confint_level >= 0,
    confint_level <= 1
  )
  stopifnot(is.null(simint) || isTRUE(simint) || isFALSE(simint))
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(is.null(xlab) || is.character(xlab))
  stopifnot(is.null(ylab) || is.character(ylab))
  stopifnot(is.null(main) || is.character(main))
  stopifnot(isTRUE(axes) || isFALSE(axes))
  stopifnot(isTRUE(box) || isFALSE(box))

  ## match arguments
  style <- match.arg(style, c("bar", "line"))
  confint_type <- match.arg(confint_type)

  ## TODO: (RS2ML) Do we have to test col, border, lwd, lty, alpha_min,
  ##       expected_col, expected_lty, confint_col, confint_lty,
  ##       confint_lwd, confint_alpha, simint_col, simint_lty, simint_lwd


  ## extend input object on correct scale (compute ci, ...)
  x <- summary(
    x, 
    freq = freq,
    confint_level = confint_level, 
    confint_type = confint_type
  ) 

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  # -------------------------------------------------------------------
  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING
  # -------------------------------------------------------------------
  ## determine which style should be plotted
  if (n > 1 && single_graph && style == "bar") {
    message("For several histograms in a single graph solely line style histograms can be plotted: \n * `style` = 'line' and `single_graph` = TRUE has been set accordingly.")
    style <- "line"
  }

  ## determine which confint should be plotted
  if (is.logical(confint)) {
    confint <- ifelse(confint,
      if (style == "bar") "line" else "polygon", 
      "none"
    )
  }
  stopifnot(all(confint %in% c("polygon", "line", "none")))

  ## determine other arguments conditional on `style` and `type`
  if (is.null(lwd)) lwd <- if (style == "bar") 1 else 2
  if (is.null(confint_col)) confint_col <- if (style == "bar") 2 else "black"
  # FIXME: (ML) Check this again May 2nd.
  if (is.null(confint_alpha)) confint_alpha <- if (single_graph && any(confint == "line")) 1 else 0.2 / n
  if (is.null(expected_col)) expected_col <- if (style == "bar") 2 else "black"
  if (is.null(expected_lty)) expected_lty <- if (style == "bar") 1 else 2

  if (is.null(simint)) {
    simint <- if (style == "bar" && any(type == "random")) TRUE else FALSE
  }

  ## prepare xlim and ylim
  if (is.list(xlim)) xlim <- as.data.frame(do.call("rbind", xlim))
  if (is.list(ylim)) ylim <- as.data.frame(do.call("rbind", ylim))

  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(1:n, expected, confint, simint,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]],
    axes, box,
    col, border, lwd, lty, alpha_min,
    expected_col, expected_lty, expected_lwd,
    confint_col, confint_lty, confint_lwd, confint_alpha,
    simint_col, simint_lty, simint_lwd
  )[, -1]

  ## prepare annotation
  if (single_graph) {
    xlab <- use_arg_from_attributes(x, "xlab", default = "PIT", force_single = TRUE)
    ylab <- use_arg_from_attributes(x, "ylab",
      default = if (freq) "Frequency" else "Density",
      force_single = TRUE
    )
    if (is.null(main)) main <- "PIT histogram"
  } else {
    xlab <- use_arg_from_attributes(x, "xlab", default = "PIT", force_single = FALSE)
    ylab <- use_arg_from_attributes(x, "ylab",
      default = if (freq) "Frequency" else "Density",
      force_single = FALSE
    )
    main <- use_arg_from_attributes(x, "main", default = "model", force_single = FALSE)
  }

  ## fix `ylabel` according to possible new `freq`
  if (ylab_missing) {
    ylab[(!freq & ylab == "Frequency")] <- "Density"
    ylab[(freq & ylab == "Density")]    <- "Frequency"
  }

  # -------------------------------------------------------------------
  # PREPARE DATA FOR PLOTTING
  # -------------------------------------------------------------------

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR 'HISTOGRAM-STYLE PITHIST'
  # -------------------------------------------------------------------
  pitbar_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## rect elements
    xleft  <- d$mid - d$width / 2
    xright <- d$mid + d$width / 2

    ## step elements (only needed for expected, confint)
    z <- compute_breaks(d$mid, d$width)

    ## get xlim and ylim
    ylim_idx <- c(is.na(plot_arg$ylim1[j]), is.na(plot_arg$ylim2[j]))
    xlim_idx <- c(is.na(plot_arg$xlim1[j]), is.na(plot_arg$xlim2[j]))
    if (any(xlim_idx)) {
      plot_arg[j, c("xlim1", "xlim2")[xlim_idx]] <- range(c(xleft, xright))[xlim_idx]
    }
    if (any(ylim_idx)) {
      ylim_use <- c(
        TRUE,
        rep(
          c(TRUE, plot_arg$expected[j], plot_arg$confint[j] != "none", plot_arg$confint[j] != "none",
            plot_arg$simint[j], plot_arg$simint[j]),
          each = NROW(d)
        )
      )
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <- range(
        c(0, d$observed, d$expected, d$confint_lwr, d$confint_upr, d$simint_lwr, d$simint_upr)[ylim_use],
        na.rm = TRUE
      )[ylim_idx]
    }

    ## trigger plot
    if (j == 1 || (!single_graph && j > 1)) {
      plot(NA,
           xlim = c(plot_arg$xlim1[j], plot_arg$xlim2[j]),
           ylim = c(plot_arg$ylim1[j], plot_arg$ylim2[j]),
           xlab = xlab[j], ylab = ylab[j], main = main[j], axes = FALSE, ...)
      if (plot_arg$axes[j]) sapply(1:2, axis)
      if (plot_arg$box[j])  box()
    }

    ## plot pithist
    rect(xleft, 0, xright, d$observed,
         border = plot_arg$border[j],
         col = set_minimum_transparency(plot_arg$col[j], alpha_min = plot_arg$alpha_min[j]),
         lty = plot_arg$lty[j],
         lwd = plot_arg$lwd[j])

    ## plot sim lines (vertical 'bars/lines')
    if (plot_arg$simint[j]) {
      segments(x0  = d$mid,
               y0  = d$simint_lwr,
               y1  = d$simint_upr,
               col = plot_arg$simint_col,
               lty = plot_arg$simint_lty,
               lwd = plot_arg$simint_lwd)
    }

    ## add expected line
    if (plot_arg$expected[j]) {
      expected_y <- duplicate_last_value(d$expected)
      lines(expected_y ~ z,
            type = "s",
            col  = plot_arg$expected_col[j],
            lty  = plot_arg$expected_lty[j],
            lwd  = plot_arg$expected_lwd[j])
    }

    ## plot confint lines
    if (plot_arg$confint[j] == "line") {
      ## lower confint line
      confint_lwr_y <- duplicate_last_value(d$confint_lwr)
      lines(confint_lwr_y ~ z,
            type = "s",
            col  = plot_arg$confint_col[j],
            lty  = plot_arg$confint_lty[j],
            lwd  = plot_arg$confint_lwd[j])

      ## upper confint line
      confint_upr_y <- duplicate_last_value(d$confint_upr)
      lines(confint_upr_y ~ z,
            type = "s",
            col  = plot_arg$confint_col[j],
            lty  = plot_arg$confint_lty[j],
            lwd  = plot_arg$confint_lwd[j])

    ## plot confint polygons
    } else if (plot_arg$confint[j] == "polygon") {
      polygon(x   = c(rep(z, each = 2)[-c(1, length(z) * 2)],
                      rev(rep(z, each = 2)[-c(1, length(z) * 2)])),
              y   = c(rep(d$confint_lwr, each = 2),
                      rev(rep(d$confint_upr, each = 2))),
              col = set_minimum_transparency(plot_arg$confint_col[j],
                                             alpha_min = plot_arg$confint_alpha[j]),
              border = NA)
    }
  }


  # -------------------------------------------------------------------
  # FUNCTION TO TRIGGER FIGURE AND PLOT CONFINT FOR 'LINE-STYLE PITHIST'
  # -------------------------------------------------------------------
  pitlines_trigger <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## step elements
    z <- compute_breaks(d$mid, d$width)
    y <- duplicate_last_value(d$observed)

    ## get xlim and ylim (needs data for all groups)
    ylim_idx <- c(is.na(plot_arg$ylim1[j]), is.na(plot_arg$ylim2[j]))
    xlim_idx <- c(is.na(plot_arg$xlim1[j]), is.na(plot_arg$xlim2[j]))

    if (any(xlim_idx)) {
      plot_arg[j, c("xlim1", "xlim2")[xlim_idx]] <- range(z)[xlim_idx]
    }

    ## 
    if (any(ylim_idx) && !single_graph) {
      ylim_use <- rep(
        c(TRUE, plot_arg$expected[j], plot_arg$confint[j] != "none", plot_arg$confint[j] != "none", 
          plot_arg$simint[j], plot_arg$simint[j]), 
        each = NROW(d)
      )
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <- range(
        c(d$observed, d$expected, d$confint_lwr, d$confint_upr, d$simint_lwr, d$simint_upr)[ylim_use],
        na.rm = TRUE
      )[ylim_idx]
    } else if (any(ylim_idx) && single_graph) {
      ylim_use <- rep(
        c(TRUE, plot_arg$expected[j], plot_arg$confint[j] != "none", plot_arg$confint[j] != "none", 
          plot_arg$simint[j], plot_arg$simint[j]), 
        each = NROW(x)
      )
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <- range(
        c(x$observed, x$expected, x$confint_lwr, x$confint_upr, x$simint_lwr, x$simint_upr)[ylim_use],
        na.rm = TRUE
      )[ylim_idx]
    }

    ## trigger plot
    if (j == 1 || (!single_graph && j > 1)) {
      plot(NA,
           xlim = c(plot_arg$xlim1[j], plot_arg$xlim2[j]),
           ylim = c(plot_arg$ylim1[j], plot_arg$ylim2[j]),
           xlab = xlab[j], ylab = ylab[j], xaxs = "i", main = main[j], axes = FALSE, ...)
      if (plot_arg$axes[j]) sapply(1:2, axis)
      if (plot_arg$box[j])  box()
    }

  }

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR 'LINE_STYLE PITHIST'
  # -------------------------------------------------------------------
  pitlines_plot <- function(d, ...) {
    ## get group index
    j <- unique(d$group)

    ## step elements
    z <- compute_breaks(d$mid, d$width)
    y <- duplicate_last_value(d$observed)

    ## plot confint lines
    if (plot_arg$confint[j] == "line") {

      ## lower confint line
      confint_lwr_y <- duplicate_last_value(d$confint_lwr)
      lines(confint_lwr_y ~ z,
            type = "s",
            col  = plot_arg$confint_col[j],
            lty  = plot_arg$confint_lty[j],
            lwd  = plot_arg$confint_lwd[j])

      ## upper confint line
      confint_upr_y <- duplicate_last_value(d$confint_upr)
      lines(confint_upr_y ~ z,
            type = "s",
            col  = plot_arg$confint_col[j],
            lty  = plot_arg$confint_lty[j],
            lwd  = plot_arg$confint_lwd[j])

    ## plot confint polygons
    } else if (plot_arg$confint[j] == "polygon") {

      polygon(x = c(rep(z, each = 2)[-c(1, length(z) * 2)],
                    rev(rep(z, each = 2)[-c(1, length(z) * 2)])),
              y = c(rep(d$confint_lwr, each = 2),
                    rev(rep(d$confint_upr, each = 2))),
              col = set_minimum_transparency(plot_arg$confint_col[j], alpha_min = plot_arg$confint_alpha[j]),
              border = NA)
    }

    ## plot expected line
    if (plot_arg$expected[j]) {

      expected_y <- duplicate_last_value(d$expected)
      lines(expected_y ~ z,
            type = "s",
            col  = plot_arg$expected_col[j],
            lty  = plot_arg$expected_lty[j],
            lwd  = plot_arg$expected_lwd[j])
    }

    ## plot sim lines
    if (plot_arg$simint[j]) {
      segments(x0  = d$mid,
               y0  = d$simint_lwr,
               y1  = d$simint_upr,
               col = plot_arg$simint_col,
               lty = plot_arg$simint_lty,
               lwd = plot_arg$simint_lwd)
    }

    ## plot stepfun
    lines(y ~ z,
          type = "s",
          lwd  = plot_arg$lwd[j],
          lty  = plot_arg$lty[j],
          col  = plot_arg$col[j])
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
    if (style == "bar") {
      for (i in 1L:n) pitbar_plot(x[x$group == i, ], ...)
    } else {
      for (i in 1L:n) pitlines_trigger(x[x$group == i, ], ...)
      for (i in 1L:n) pitlines_plot(x[x$group == i, ], ...)
    }
  } else {
    if (style == "bar") {
      for (i in 1L:n) pitbar_plot(x[x$group == i, ], ...)
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
                          freq = NULL,
                          expected = FALSE,
                          confint = FALSE,
                          confint_level = 0.95,
                          confint_type = c("exact", "approximation"),
                          simint = FALSE,
                          col = "black",
                          lwd = 2,
                          lty = 1,
                          expected_col = "black",
                          expected_lty = 2,
                          expected_lwd = 1.75,
                          confint_col = "black",
                          confint_lty = 1,
                          confint_lwd = 1.75,
                          confint_alpha = 1,
                          simint_col = "black",
                          simint_lty = 1,
                          simint_lwd = 1.75,
                          ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## get default arguments
  freq <- use_arg_from_attributes(x, "freq", default = FALSE, force_single = TRUE)
  expected <- use_arg_from_attributes(x, "expected", default = FALSE, force_single = FALSE)
  confint <- use_arg_from_attributes(x, "confint", default = FALSE, force_single = FALSE)
  simint <- use_arg_from_attributes(x, "simint", default = FALSE, force_single = FALSE)

  ## sanity checks
  ## * lengths of most arguments are checked by recycling
  ## * graphical parameters are checked w/i function calls for plotting
  stopifnot(is.logical(freq))
  stopifnot(is.logical(expected))
  stopifnot(is.logical(confint) || confint %in% c("polygon", "line", "none"))
  stopifnot(
    is.numeric(confint_level),
    length(confint_level) == 1,
    confint_level >= 0,
    confint_level <= 1
  )
  stopifnot(is.logical(simint))

  ## match arguments
  confint_type <- match.arg(confint_type)

  ## extend input object on correct scale (compute ci, ...)
  x <- summary(
    x,
    freq = freq,
    confint_level = confint_level,
    confint_type = confint_type
  )

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  # -------------------------------------------------------------------
  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING
  # -------------------------------------------------------------------
  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(1:n, 
    confint, expected, simint,
    col, lwd, lty,
    expected_col, expected_lty, expected_lwd,
    confint_col, confint_lty, confint_lwd, confint_alpha,
    simint_col, simint_lty, simint_lwd
  )[, -1]

  # -------------------------------------------------------------------
  # PREPARE DATA FOR PLOTTING
  # -------------------------------------------------------------------

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR LINES
  # -------------------------------------------------------------------
  pitlines_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## step elements
    z <- compute_breaks(d$mid, d$width)
    y <- duplicate_last_value(d$observed)

    ## plot confint lines
    if (plot_arg$confint[j] == "line") {

      ## lower confint line
      confint_lwr_y <- duplicate_last_value(d$confint_lwr)
      lines(
        confint_lwr_y ~ z,
        type = "s",
        col = plot_arg$confint_col[j],
        lty = plot_arg$confint_lty[j],
        lwd = plot_arg$confint_lwd[j]
      )

      ## upper confint line
      confint_upr_y <- duplicate_last_value(d$confint_upr)
      lines(
        confint_upr_y ~ z,
        type = "s",
        col = plot_arg$confint_col[j],
        lty = plot_arg$confint_lty[j],
        lwd = plot_arg$confint_lwd[j]
      )

    } else if (plot_arg$confint[j] == "polygon") {

      polygon(
        c(
          rep(z, each = 2)[-c(1, length(z) * 2)],
          rev(rep(z, each = 2)[-c(1, length(z) * 2)])
        ),
        c(
          rep(d$confint_lwr, each = 2),
          rev(rep(d$confint_upr, each = 2))
        ),
        col = set_minimum_transparency(plot_arg$confint_col[j], alpha_min = plot_arg$confint_alpha[j]),
        border = NA
      )
    }

    ## plot expected line
    if (plot_arg$expected[j]) {

      expected_y <- duplicate_last_value(d$expected)
      lines(
        expected_y ~ z,
        type = "s",
        col = plot_arg$expected_col[j],
        lty = plot_arg$expected_lty[j],
        lwd = plot_arg$expected_lwd[j]
      )
    }

    ## plot sim lines
    if (plot_arg$simint[j]) {
      segments(
        x0 = d$mid,
        y0 = d$simint_lwr,
        y1 = d$simint_upr,
        col = plot_arg$simint_col,
        lty = plot_arg$simint_lty,
        lwd = plot_arg$simint_lwd)
    }


    ## plot stepfun
    lines(
      y ~ z,
      type = "s",
      lwd = plot_arg$lwd[j],
      lty = plot_arg$lty[j],
      col = plot_arg$col[j]
    )
  }


  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  for (i in 1L:n) {
    pitlines_plot(x[x$group == i, ], ...)
  }
}


#' @rdname plot.pithist
#' @method autoplot pithist
#' @exportS3Method ggplot2::autoplot pithist
autoplot.pithist <- function(object,
                             single_graph = FALSE,
                             style = NULL,
                             freq = NULL,
                             expected = NULL,
                             confint = NULL,
                             confint_level = 0.95,
                             confint_type = c("exact", "approximation"),
                             simint = NULL,
                             xlim = c(NA, NA),
                             ylim = c(0, NA),
                             xlab = NULL,
                             ylab = NULL,
                             main = NULL,
                             legend = FALSE,
                             theme = NULL,
                             colour = NULL,
                             fill = NULL,
                             size = NULL,
                             linetype = NULL,
                             alpha = NULL,
                             expected_colour = NULL,
                             expected_size = 0.75,
                             expected_linetype = NULL,
                             expected_alpha = NA,
                             confint_colour = NULL,
                             confint_fill = NULL,
                             confint_size = 0.75,
                             confint_linetype = NULL,
                             confint_alpha = NULL,
                             simint_colour = "black",
                             simint_size = 0.5,
                             simint_linetype = 1,
                             simint_alpha = NA,
                             ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## check if ylab is defined
  ylab_missing <- missing(ylab)

  ## get default arguments
  type <- use_arg_from_attributes(object, "type", default = NULL, force_single = FALSE)
  style <- use_arg_from_attributes(object, "style", default = "bar", force_single = TRUE)
  freq <- use_arg_from_attributes(object, "freq", default = FALSE, force_single = TRUE)
  scale <- use_arg_from_attributes(object, "scale", default = NULL, force_single = TRUE)
  expected <- use_arg_from_attributes(object, "expected", default = NULL, force_single = TRUE)
  confint <- use_arg_from_attributes(object, "confint", default = TRUE, force_single = TRUE)
  simint <- use_arg_from_attributes(object, "simint", default = TRUE, force_single = TRUE)
  xlab <- use_arg_from_attributes(object, "xlab", default = "PIT", force_single = TRUE)
  ylab <- use_arg_from_attributes(object, "ylab",
    default = if (freq) "Frequency" else "Density",
    force_single = TRUE
  )

  ## get base style arguments
  add_arg <- list(...)
  if (!is.null(add_arg$lwd)) size <- add_arg$lwd
  if (!is.null(add_arg$lty)) linetype <- add_arg$lty

  ## fix `ylab` according to possible new `freq`
  if (ylab_missing) {
    if (freq && ylab == "Density") ylab <- "Frequency"
    if (!freq && ylab == "Frequency") ylab <- "Density"
  }

  ## sanity checks
  stopifnot(is.logical(single_graph))
  stopifnot(is.logical(freq))
  stopifnot(is.null(simint) || is.logical(simint))
  stopifnot(is.logical(confint) || confint %in% c("polygon", "line", "none"))
  stopifnot(
    is.numeric(confint_level),
    length(confint_level) == 1,
    confint_level >= 0,
    confint_level <= 1
  )
  stopifnot(is.null(expected) || is.logical(expected))
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(is.character(xlab), length(xlab) == 1)
  stopifnot(is.character(ylab), length(ylab) == 1)
  stopifnot(is.logical(legend))

  ## match arguments
  style <- match.arg(style, c("bar", "line"))
  confint_type <- match.arg(confint_type)
  scale <- match.arg(scale, c("uniform", "normal"))

  ## set all aesthetics equal NULL to NA
  alpha <- if (is.null(alpha)) NA else alpha
  colour <- if (is.null(colour)) NA else colour
  fill <- if (is.null(fill)) NA else fill
  size <- if (is.null(size)) NA else size
  linetype <- if (is.null(linetype)) NA else linetype

  ## transform to `freq = TRUE` scale and check if allowed (not for non-equidistant breaks)
  ## TODO: (ML) Improve and maybe move that into `geom`s.
  object <- summary.pithist(object, freq = TRUE, extend = FALSE, suppress_warnings = !freq)

  ## convert data always to data.frame
  object <- as.data.frame(object)

  ## determine grouping
  if (is.null(object$group)) object$group <- 1L
  n <- max(object$group)

  ## get title (must be done before handling of `main`)
  if (!is.null(main)) {
    title <- main[1]
    object$title <- factor(title)
  }

  ## get main and into the right length (must be done after handling of `title`)
  main <- use_arg_from_attributes(object, "main", default = "model", force_single = FALSE)
  stopifnot(is.character(main))
  main <- make.names(rep_len(main, n), unique = TRUE)

  ## prepare grouping
  object$group <- factor(object$group, levels = 1L:n, labels = main)

  # -------------------------------------------------------------------
  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING
  # -------------------------------------------------------------------
  ## determine which style should be plotted
  if (n > 1 && single_graph && style == "bar") {
    message("For several histograms in a single graph solely line style histograms can be plotted: \n * `style` = 'line' and `single_graph` = TRUE has been set accordingly.")
    style <- "line"
  }

  ## determine which confint should be plotted
  if (isFALSE(confint)) {
    confint <- "none"
  } else if (isTRUE(confint)) {
    confint <- if (style == "bar") "line" else "polygon"
  }
  confint <- match.arg(confint, c("polygon", "line", "none"))

  ## determine other arguments conditional on `style`
  if (is.null(simint)) {
    simint <- if (style == "bar" && any(type == "random")) TRUE else FALSE
  }

  if (is.null(expected_colour)){ 
    expected_colour <- if (style == "bar") 2 else "black" 
  }

  if (is.null(confint_colour)){ 
    confint_colour <- if (style == "bar") 2 else NA 
  }

  if (is.null(confint_alpha) && confint == "polygon"){ 
    confint_alpha <- if (style == "bar") NA else 0.2 / n
  }

  ## set plotting aes
  if (style == "line") {
    aes_expected_default <- list(colour = "black", linetype = 2)
  } else {
    aes_expected_default <- NULL
  }
  aes_expected <- set_aes_helper_geoms(
    GeomPithistExpected$default_aes,
    list(
      colour = expected_colour,
      size = expected_size,
      linetype = expected_linetype,
      alpha = expected_alpha
    ),
    aes_expected_default
  )

  aes_confint <- set_aes_helper_geoms(
    set_default_aes_pithist_confint(confint),
    list(
      colour = confint_colour,
      fill = confint_fill,
      size = confint_size,
      linetype = confint_linetype,
      alpha = confint_alpha
    )
  )

  aes_simint <- set_aes_helper_geoms(
    GeomPithistSimint$default_aes,
    list(
      colour = simint_colour,
      size = simint_size,
      linetype = simint_linetype,
      alpha = simint_alpha
    )
  )

  ## recycle arguments for plotting to match the number of groups (for `scale_<...>_manual()`)
  plot_arg <- data.frame(
    1:n,
    colour, fill, size, linetype, alpha
  )[, -1]

  # -------------------------------------------------------------------
  # MAIN PLOTTING
  # -------------------------------------------------------------------
  ## actual plotting
  rval <- ggplot2::ggplot(
    object,
    ggplot2::aes_string(
      x = "mid",
      y = "observed",
      width = "width",
      group = "group"
    )
  )

  ## add histogram
  rval <- rval +
    geom_pithist(
      ggplot2::aes_string(
        colour = "group",
        fill = "group",
        size = "group",
        linetype = "group",
        alpha = "group"
      ),
      freq = freq,
      style = style
    )

  ## add simint
  if (simint) {
    rval <- rval +
      geom_pithist_simint(
        ggplot2::aes_string(
          x = "mid",
          ymin = "simint_lwr",
          ymax = "simint_upr",
          group = "group"
        ),
        na.rm = TRUE,
        freq = freq,
        colour = aes_simint$colour,
        size = aes_simint$size,
        linetype = aes_simint$linetype,
        alpha = aes_simint$alpha
      )
  }

  ## add confint
  if (confint != "none") {
    rval <- rval +
      geom_pithist_confint(
        ggplot2::aes_string(
          x = "mid",
          y = "observed",
          width = "width"
        ),
        level = confint_level,
        type = confint_type,
        freq = freq,
        scale = scale,
        style = confint,
        colour = aes_confint$colour,
        fill = aes_confint$fill,
        size = aes_confint$size,
        linetype = aes_confint$linetype,
        alpha = aes_confint$alpha
      )
  }

  ## add expected
  if (expected) {
    rval <- rval +
      geom_pithist_expected(
        ggplot2::aes_string(
          x = "mid",
          y = "observed",
          width = "width"
        ),
        freq = freq,
        scale = scale,
        colour = aes_expected$colour,
        size = aes_expected$size,
        linetype = aes_expected$linetype,
        alpha = aes_expected$alpha
      )
  }

  ## set the colors, shapes, etc. for the groups
  ## w/i `set_default_aes_pithist()` by `my_modify_list()` all NAs or 0.999 values
  ## are replaced with intern defaults
  rval <- rval +
    ggplot2::scale_alpha_manual(values = plot_arg$alpha, na.value = NA) +
    ggplot2::scale_colour_manual(values = plot_arg$colour, na.value = NA) +
    ggplot2::scale_fill_manual(values = plot_arg$fill, na.value = NA) +
    ggplot2::scale_size_manual(values = plot_arg$size, na.value = 0.999) +
    ggplot2::scale_linetype_manual(values = plot_arg$linetype, na.value = NA)

  ## annotation
  rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

  ## add legend
  if (legend) {
    rval <- rval +
      ggplot2::labs(
        alpha = "Model",
        colour = "Model",
        fill = "Model",
        size = "Model",
        linetype = "Model"
      ) +
      ggplot2::guides(
        alpha = "legend",
        colour = "legend",
        fill = "legend",
        size = "legend",
        linetype = "legend"
      )
  } else {
    rval <- rval +
      ggplot2::guides(
        alpha = "none",
        colour = "none",
        fill = "none",
        size = "none",
        linetype = "none"
      )
  }

  ## set x and y limits
  rval <- rval + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE) +
    ggplot2::scale_x_continuous(expand = c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.01, 0.01))

  ## set ggplot2 theme
  if (is.character(theme)) {
    theme_tmp <- try(eval(parse(text = theme)), silent = TRUE)
    if (inherits(theme_tmp, "try-error") && !grepl("^ggplot2::", theme)) {
      theme_tmp <- try(eval(parse(text = paste0("ggplot2::", theme))), silent = TRUE)
    }
    theme <- theme_tmp
    if (!is.function(theme)) {
        warning("The argument `theme` must be a `ggplot2` theme, a theme-generating function or a valid 'character string'.", call. = FALSE)
      theme <- ggplot2::theme_bw()
    }
  }

  if (is.function(theme)) {
    theme <- try(theme(), silent = TRUE)
    if (inherits(theme, "try-error") || !inherits(theme, "theme")) {
        warning("The argument `theme` must be a `ggplot2` theme, a theme-generating function or a valid 'character string'.", call. = FALSE)
      theme <- ggplot2::theme_bw()
    }
  }

  if (inherits(theme, "theme")) {
    rval <- rval + theme
  } else if (isTRUE(all.equal(ggplot2::theme_get(), ggplot2::theme_gray()))) {
    rval <- rval + ggplot2::theme_bw()
  }

  # -------------------------------------------------------------------
  # ORDER LAYERS
  # -------------------------------------------------------------------
  # FIXME: (ML) Order layers conditional on plotting style.
  # if (style == "line") {
  #  rval$layer <- c(rval$layer[[3]], rval$layer[[4]], rval$layer[[1]], rval$layer[[2]])
  # }

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


# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_pithist()`
# -------------------------------------------------------------------

#' @rdname geom_pithist
#' @export
stat_pithist <- function(mapping = NULL,
                         data = NULL,
                         geom = "pithist",
                         position = "identity",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE,
                         freq = FALSE,
                         style = c("bar", "line"),
                         ...) {
  style <- match.arg(style)

  ggplot2::layer(
    stat = StatPithist,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      freq = freq,
      style = style,
      ...
    )
  )
}


#' @rdname geom_pithist
#' @format NULL
#' @usage NULL
#' @export
StatPithist <- ggplot2::ggproto("StatPithist", ggplot2::Stat,
  required_aes = c("x", "y", "width"),
  compute_group = function(data, scales, freq = FALSE, style = "bar") {

    ## transform observed to freq
    if (!freq) {
      data <- transform(data,
        y = y / (sum(y) * width)
      )
    }

    ## prepare data for different styles
    if (style == "bar") {
      transform(data,
        y = y / 2,
        height = y
      )
    } else { # "line" style
      data.frame(
        x = compute_breaks(data$x, data$width),
        y = duplicate_last_value(data$y)
      )
    }
  }
)


#' \code{geom_*} and \code{stat_*} for Producing PIT Histograms with `ggplot2`
#'
#' Various \code{geom_*} and \code{stat_*} used within
#' \code{\link[ggplot2]{autoplot}} for producing PIT histograms.
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#' @param style character specifying the style of pithist. For \code{style = "bar"}
#' a traditional PIT hisogram is drawn, for \code{style = "line"} solely the upper border
#' line is plotted.
#' @param type character. Which type of confidence interval should be plotted: 
#' `"exact"` or `"approximation"`. According
#' to Agresti and Coull (1998), for interval estimation of binomial proportions
#' an approximation can be better than exact.
#' @param level numeric. The confidence level required.
#' @param freq logical. If \code{TRUE}, the PIT histogram is represented by
#' frequencies, the \code{counts} component of the result; if \code{FALSE},
#' probability densities, component \code{density}, are plotted (so that the
#' histogram has a total area of one).
#' @param scale On which scale should the PIT residuals be computed: on the probability scale 
#' (\code{"uniform"}) or on the normal scale (\code{"normal"}).
#' @examples
#' if (require("ggplot2")) {
#'   ## Fit model
#'   data("CrabSatellites", package = "countreg")
#'   m1_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#'   m2_pois <- glm(satellites ~ color, data = CrabSatellites, family = poisson)
#'
#'   ## Compute pithist
#'   p1 <- pithist(m1_pois, type = "random", plot = FALSE)
#'   p2 <- pithist(m2_pois, type = "random", plot = FALSE)
#'
#'   d <- c(p1, p2)
#'
#'   ## Create factor
#'   main <- attr(d, "main")
#'   main <- make.names(main, unique = TRUE)
#'   d$group <- factor(d$group, labels = main)
#'
#'   ## Plot bar style PIT histogram
#'   gg1 <- ggplot(data = d) +
#'     geom_pithist(aes(x = mid, y = observed, width = width, group = group), freq = TRUE) +
#'     geom_pithist_simint(aes(x = mid, ymin = simint_lwr, ymax = simint_upr), freq = TRUE) +
#'     geom_pithist_confint(aes(x = mid, y = observed, width = width), style = "line", freq = TRUE) +
#'     geom_pithist_expected(aes(x = mid, y = observed, width = width), freq = TRUE) +
#'     facet_grid(group ~ .) +
#'     xlab("PIT") +
#'     ylab("Frequency")
#'   gg1
#'
#'   gg2 <- ggplot(data = d) +
#'     geom_pithist(aes(x = mid, y = observed, width = width, group = group), freq = FALSE) +
#'     geom_pithist_simint(aes(
#'       x = mid, ymin = simint_lwr, ymax = simint_upr, y = observed,
#'       width = width
#'     ), freq = FALSE) +
#'     geom_pithist_confint(aes(x = mid, y = observed, width = width), style = "line", freq = FALSE) +
#'     geom_pithist_expected(aes(x = mid, y = observed, width = width), freq = FALSE) +
#'     facet_grid(group ~ .) +
#'     xlab("PIT") +
#'     ylab("Density")
#'   gg2
#'
#'   ## Plot line style PIT histogram
#'   gg3 <- ggplot(data = d) +
#'     geom_pithist(aes(x = mid, y = observed, width = width, group = group), style = "line") +
#'     geom_pithist_confint(aes(x = mid, y = observed, width = width), style = "polygon") +
#'     facet_grid(group ~ .) +
#'     xlab("PIT") +
#'     ylab("Density")
#'   gg3
#' }
#' @export
geom_pithist <- function(mapping = NULL,
                         data = NULL,
                         stat = "pithist",
                         position = "identity",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE,
                         freq = FALSE, # only needed w/i stat_*
                         style = c("bar", "line"),
                         ...) {
  style <- match.arg(style)

  ggplot2::layer(
    geom = GeomPithist,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      freq = freq,
      style = style,
      ...
    )
  )
}


#' @rdname geom_pithist
#' @format NULL
#' @usage NULL
#' @export
GeomPithist <- ggplot2::ggproto("GeomPithist", ggplot2::Geom,
  default_aes = ggplot2::aes(
    colour = NA,
    fill = NA,
    size = 0.999, # TODO: (ML) resolve this hack (NA leads to error, even if it is modified w/i draw_panel() -> why?)
    linetype = NA,
    alpha = NA,
    subgroup = NULL
  ),
  required_aes = c("x", "y"),
  optional_aes = c("width", "height"), # TODO: (ML) "width" and "height" actually required for `style = "bar"`

  extra_params = c("na.rm", "style"),
  setup_data = function(data, params) {
    if (params$style == "bar") {
      ## uses `geom_tile()`, which uses center of tiles and its size (x, y, width, height)
      ggplot2::GeomTile$setup_data(data, params)
    } else if (params$style == "line") {
      ## uses `geom_step()`, which uses requires all x and y (including min/max, so one more than center points)
      data
    }
  },
  draw_panel = function(data, panel_params, coord, linejoin = "mitre", style = c("bar", "line")) {

    ## get style
    style <- match.arg(style)

    ## swap NAs in `default_aes` with own defaults
    data <- my_modify_list(data, set_default_aes_pithist(style), force = FALSE)

    ## create geom
    if (style == "bar") {
      ggplot2::GeomTile$draw_panel(
        data = data,
        panel_params = panel_params,
        coord = coord,
        linejoin = linejoin
      )
    } else { # "line" style
      ggplot2::GeomStep$draw_panel(
        data = data,
        panel_params = panel_params,
        coord = coord
      )
    }
  },
  draw_key = function(data, params, size) {
    ## Swap NAs in `default_aes` with own defaults
    data <- my_modify_list(data, set_default_aes_pithist(params$style), force = FALSE)

    if (params$style == "bar") {
      draw_key_polygon(data, params, size)
    } else { # "line" style
      draw_key_path(data, params, size)
    }
  }
)


## helper function inspired by internal from `ggplot2` defined in `geom-sf.r`
set_default_aes_pithist <- function(style) {
  if (style == "bar") {
    my_modify_list(
      ggplot2::GeomPolygon$default_aes, 
      list(colour = "black", fill = "darkgray", size = 0.5, linetype = 1, alpha = NA, subgroup = NULL),
      force = TRUE
    )
  } else { # "line" style
    my_modify_list(ggplot2::GeomPath$default_aes, list(colour = "black", size = 0.75, linetype = 1, alpha = NA),
      force = TRUE
    )
  }
}


# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_pithist_expected()`
# -------------------------------------------------------------------

#' @rdname geom_pithist
#' @export
stat_pithist_expected <- function(mapping = NULL,
                             data = NULL,
                             geom = "pithist_expected",
                             position = "identity",
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             scale = c("uniform", "normal"),
                             freq = FALSE,
                             ...) {

  scale <- match.arg(scale)

  ggplot2::layer(
    stat = StatPithistExpected,
    mapping = mapping,
    data = data,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      scale = scale,
      freq = freq,
      ...
    )
  )
}


#' @rdname geom_pithist
#' @format NULL
#' @usage NULL
#' @export
StatPithistExpected <- ggplot2::ggproto("StatPithistExpected", ggplot2::Stat,
  required_aes = c("x", "y", "width"),
  compute_group = function(data,
                           scales,
                           scale = "uniform",
                           freq = FALSE) {

    ## compute expected line
    expected <- compute_pithist_expected(
      n = sum(data$y),
      breaks = compute_breaks(data$x, data$width),
      freq = freq,
      scale = scale
    )

    data.frame(
      x = compute_breaks(data$x, data$width),
      y = duplicate_last_value(expected)
    )
  }
)

#' @rdname geom_pithist
#' @export
geom_pithist_expected <- function(mapping = NULL,
                             data = NULL,
                             stat = "pithist_expected",
                             position = "identity",
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             scale = c("uniform", "normal"),
                             freq = FALSE, # only needed w/i stat_*
                             ...) {

  scale <- match.arg(scale)

  ggplot2::layer(
    geom = GeomPithistExpected,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      scale = scale,
      freq = freq,
      ...
    )
  )
}


#' @rdname geom_pithist
#' @format NULL
#' @usage NULL
#' @export
GeomPithistExpected <- ggplot2::ggproto("GeomPithistExpected", ggplot2::GeomStep,
  default_aes = ggplot2::aes(
    colour = 2,
    size = 0.75,
    linetype = 1,
    alpha = NA
  ),
  required_aes = c("x", "y")
)


# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_pithist_confint()`
# -------------------------------------------------------------------

#' @rdname geom_pithist
#' @export
stat_pithist_confint <- function(mapping = NULL,
                                 data = NULL,
                                 geom = "pithist_confint",
                                 position = "identity",
                                 na.rm = FALSE,
                                 show.legend = NA,
                                 inherit.aes = TRUE,
                                 scale = c("uniform", "normal"),
                                 level = 0.95,
                                 type = "approximation",
                                 freq = FALSE,
                                 style = c("polygon", "line"),
                                 ...) {
  style <- match.arg(style)
  scale <- match.arg(scale)

  ggplot2::layer(
    stat = StatPithistConfint,
    mapping = mapping,
    data = data,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      scale = scale,
      level = level,
      type = type,
      freq = freq,
      style = style,
      ...
    )
  )
}


#' @rdname geom_pithist
#' @format NULL
#' @usage NULL
#' @export
StatPithistConfint <- ggplot2::ggproto("StatPithistConfint", ggplot2::Stat,
  required_aes = c("x", "y", "width"),
  compute_group = function(data,
                           scales,
                           scale = "uniform",
                           level = 0.95,
                           type = "approximation",
                           freq = FALSE,
                           style = "polygon") {
    ## compute ci interval
    ci <- compute_pithist_confint(
      n = sum(data$y),
      breaks = compute_breaks(data$x, data$width),
      level = level,
      type = type,
      freq = freq,
      scale = scale
    )

    ## return new data.frame condition on plotting `style`
    if (style == "polygon") {
      nd <- data.frame(
        xmin = data$x - data$width / 2,
        xmax = data$x + data$width / 2,
        ymin = ci[[1]],
        ymax = ci[[2]]
      )
    } else {
      nd <- data.frame(
        x = compute_breaks(data$x, data$width),
        ymin = duplicate_last_value(ci[[1]]),
        ymax = duplicate_last_value(ci[[2]])
      )
    }
    nd
  }
)


#' @rdname geom_pithist
#' @export
geom_pithist_confint <- function(mapping = NULL,
                                 data = NULL,
                                 stat = "pithist_confint",
                                 position = "identity",
                                 na.rm = FALSE,
                                 show.legend = NA,
                                 inherit.aes = TRUE,
                                 scale = c("uniform", "normal"),
                                 level = 0.95,
                                 type = "approximation",
                                 freq = FALSE, # only needed w/i stat_*
                                 style = c("polygon", "line"),
                                 ...) {
  style <- match.arg(style)
  scale <- match.arg(scale)

  ggplot2::layer(
    geom = GeomPithistConfint,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      scale = scale,
      level = level,
      type = type,
      freq = freq,
      style = style,
      ...
    )
  )
}


#' @rdname geom_pithist
#' @format NULL
#' @usage NULL
#' @export
GeomPithistConfint <- ggplot2::ggproto("GeomPithistConfint", ggplot2::Geom,
  required_aes = c("x|xmin", "ymin", "ymax"),
  optional_aes = c("xmax"), # TODO: (ML) "xmax" actually required for `style = "polygon"`

  extra_params = c("na.rm", "style"),

  # TODO: (ML) Does not vary for style; this is a copy of `GeomPolygon$handle_na()`
  handle_na = function(data, params) {
    data
  },

  ## Setting up all defaults needed for `GeomPolygon` and `GeomStep`
  default_aes = ggplot2::aes(
    colour = NA,
    fill = NA,
    size = 0.999,
    linetype = NA,
    alpha = NA
  ),
  draw_panel = function(data, panel_params, coord,
                        linejoin = "mitre", direction = "hv",
                        style = c("polygon", "line")) {
    style <- match.arg(style)

    ## Swap NAs in `default_aes` with own defaults
    data <- my_modify_list(data, set_default_aes_pithist_confint(style), force = FALSE)

    if (style == "polygon") {
      ggplot2::GeomRect$draw_panel(data, panel_params, coord, linejoin)
    } else { # "line" style
      ## Join two Grobs
      data1 <- transform(data,
        y = ymin
      )
      data2 <- transform(data,
        y = ymax
      )
      grid::grobTree(
        ggplot2::GeomStep$draw_panel(data1, panel_params, coord, direction),
        ggplot2::GeomStep$draw_panel(data2, panel_params, coord, direction)
      )
    }
  },
  draw_key = function(data, params, size) {
    ## Swap NAs in `default_aes` with own defaults
    data <- my_modify_list(data, set_default_aes_pithist_confint(params$style), force = FALSE)
    if (params$style == "polygon") {
      draw_key_polygon(data, params, size)
    } else { # "line" style
      draw_key_path(data, params, size)
    }
  }
)


## helper function inspired by internal from `ggplot2` defined in `geom-sf.r`
set_default_aes_pithist_confint <- function(style) {
  if (style == "line") {
    my_modify_list(ggplot2::GeomPath$default_aes, list(colour = 2, fill = NA, size = 0.75, linetype = 2, alpha = NA),
      force = TRUE
    )
  } else { # "polygon" style
    my_modify_list(ggplot2::GeomPolygon$default_aes, list(
      colour = NA, fill = "black", size = 0.5,
      linetype = 1, alpha = 0.2, subgroup = NULL
    ), force = TRUE)
  }
}


# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_pithist_simint()`
# -------------------------------------------------------------------

#' @rdname geom_pithist
#' @export
stat_pithist_simint <- function(mapping = NULL,
                                data = NULL,
                                geom = "pithist_simint",
                                position = "identity",
                                na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE,
                                freq = FALSE,
                                ...) {
  style <- match.arg(style)

  ggplot2::layer(
    stat = StatPithistSimint,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      freq = freq,
      ...
    )
  )
}


#' @rdname geom_pithist
#' @format NULL
#' @usage NULL
#' @export
StatPithistSimint <- ggplot2::ggproto("StatPithistSimint", ggplot2::Stat,
  required_aes = c("x", "ymin", "ymax"),
  optional_aes = c("y", "width"), # TODO: (ML) "y' and "width" actually required for `freq = FALSE`

  compute_group = function(data, scales, freq = FALSE) {

    ## transform counts to freq
    if (!freq) {
      if (is.null(data$width)) stop("for arg `freq = FALSE`, aesthetics `width` must be provided.")
      if (is.null(data$y)) stop("for arg `freq = FALSE`, aesthetics `y` must be provided.")
      data <- transform(data,
        ymin = ymin / (sum(y) * width),
        ymax = ymax / (sum(y) * width),
        y = NULL
      )
    }
  }
)


#' @rdname geom_pithist
#' @export
geom_pithist_simint <- function(mapping = NULL,
                                data = NULL,
                                stat = "pithist_simint",
                                position = "identity",
                                na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE,
                                freq = FALSE, # only needed w/i stat_*
                                ...) {
  ggplot2::layer(
    geom = GeomPithistSimint,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      freq = freq,
      ...
    )
  )
}


#' @rdname geom_pithist
#' @format NULL
#' @usage NULL
#' @export
GeomPithistSimint <- ggplot2::ggproto("GeomPithistSimint", ggplot2::GeomLinerange,
  required_aes = c("x", "ymin", "ymax"),
  default_aes = ggplot2::aes(
    colour = "black",
    size = 0.5,
    linetype = 1,
    alpha = NA
  )
)


# -------------------------------------------------------------------
# HELPER FUNCTIONS FOR PIT HISTOGRAMS
# -------------------------------------------------------------------
#' Compute Reference Line for PIT Histograms
#'
#' Helper function for computing expected lines showing perfect prediction for PIT histograms.
#'
#' @noRd
#' @param n number of observations.
#' @param breaks vector with breakpoints.
#' @param type Which kind of confindence interval should be computed? Type \code{exact} using
#' the quantile function of the binomial distribution or \code{approximation} uses an approximation
#' by Agresti and Coull (1998).
#' @param freq Should confidence intervals returned for reported frequencies or
#' for counts of observation.
compute_pithist_expected <- function(n, breaks, freq, scale) {

  ## get inverse trafo
  ## TODO: (ML) Must be extended using `distributions3`
  if (scale == "uniform") {
    pFun <- identity
  } else {
    pFun <- pnorm
  }

  ## TODO: (ML)
  ## * Is this correct?
  ## * For small n qbinom(0.5, ...) is not equal 1.

  ## get probs
  probs <- diff(pFun(breaks))

  ## calc bin specific expected line
  rval <- qbinom(0.5, size = n, prob = probs)

  ## transform counts to density
  if (!freq) {
    rval <- rval / (n * diff(breaks))
  }
  rval
}


#' Compute Confidence Interval for PIT Histograms
#'
#' Helper function for computing confidence intervals for PIT histograms.
#'
#' @noRd
#' @param n number of observations.
#' @param breaks vector with breakpoints.
#' @param level confidence level.
#' @param type Which kind of confindence interval should be computed? Type \code{exact} using
#' the quantile function of the binomial distribution or \code{approximation} uses an approximation
#' by Agresti and Coull (1998).
#' @param freq Should confidence intervals returned for reported frequencies or
#' for counts of observation.
compute_pithist_confint <- function(n, breaks, level, type = c("exact", "approximation"), freq, scale) {

  type <- match.arg(type)

  ## get confidence level
  a <- (1 - level) / 2

  ## get inverse trafo
  ## TODO: (ML) Must be extended using `distributions3`
  if (scale == "uniform") {
    pFun <- identity
  } else {
    pFun <- pnorm
  }

  ## get probs
  probs <- diff(pFun(breaks))

  if (type == "exact") {
    ## calc bin specific confidence levels
    rval <- data.frame(
      qbinom(a, size = n, prob = probs),
      qbinom(1 - a, size = n, prob = probs)
    )
  } else {
    rval <- data.frame(
      do.call(
        rbind,
        lapply(n * probs, function(x) PropCIs_add4ci(x, n, level)$conf.int * n)
      )
    )
  }

  ## transform counts to density
  if (!freq) {
    rval <- data.frame(sapply(rval, function(x) x / (n * diff(breaks)), simplify = FALSE))
  }
  colnames(rval) <- c("confint_lwr", "confint_upr")
  rval
}


#' Agresti-Coull add-4 CI for a binomial proportion
#'
#' Copy of \code{add4ci} package from package `PropCIs` by Ralph Scherer
#' (licensed under GPL-2/GPL-3).
#' Agresti-Coull add-4 CI for a binomial proportion, based on adding
#' 2 successes and 2 failures before computing the Wald CI. The CI is
#' truncated, when it overshoots the boundary.
#'
#' @param x number of successes.
#' @param n number of trials.
#' @param conf.level confidence coefficient.
#' @noRd
PropCIs_add4ci <- function(x, n, conf.level) {
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


#' @export
summary.pithist <- function(object,
                            freq = NULL,
                            confint_level = 0.95,
                            confint_type = c("exact", "approximation"),
                            extend = TRUE,
                            suppress_warnings = FALSE,
                            ...) {

  stopifnot(is.logical(extend))

  ## get arg `freq`
  freq_object <- attr(object, "freq")
  freq <- use_arg_from_attributes(object, "freq", default = FALSE, force_single = TRUE)
  counts <- use_arg_from_attributes(object, "counts", default = NULL, force_single = FALSE)

  stopifnot(is.logical(freq))
  stopifnot(!"group" %in% names(object) || length(counts) == length(unique(object$group)))

  ## get arg `scale`
  scale <-  attr(object, "scale")

  ## ensure column group
  if (!any(grepl("group", names(object)))) {
    object$group <- 1 
    return_group <- FALSE
  } else {
    return_group <- TRUE
  }

  ## loop over groups
  rval <- list()
  for (i in seq_along(counts)) { 

    ## not allow freq = TRUE for non-equidistant breaks
    ## TODO: (ML) Not exactly unique widths
    if (!suppress_warnings && freq && any(abs(diff(object[object$group == i, "width"])) > sqrt(.Machine$double.eps))) {
      warning(
        "For non-equidistant breaks `freq = TRUE` should not be used due to incorrect areas.",
        call. = FALSE  
      )
    }

    ## ensure values are counts
    if (!freq_object) {
    
      ## TODO: (ML) How else to pass checks for visible binding
      tmp <- object[object$group == i, ]

      tmp <- transform(tmp,
        observed = tmp$observed * (counts[i] * tmp$width),
        expected = tmp$expected * (counts[i] * tmp$width),
        simint_upr = tmp$simint_upr * (counts[i] * tmp$width),
        simint_lwr = tmp$simint_lwr * (counts[i] * tmp$width)
      )
    } else {
      tmp <- object[object$group == i, ]
    }

    if (extend) {

      ## compute confidence intervals for all groups (freq = TRUE, as scaling below)
      ci <- compute_pithist_confint(
        n = counts[i],
        breaks = compute_breaks(tmp$mid, tmp$width),
        level = confint_level,
        type = confint_type,
        freq = TRUE,
        scale = scale
      )

      if (freq) {
        rval[[i]] <- transform(tmp,
          observed = tmp$observed,
          expected = tmp$expected,
          simint_upr = tmp$simint_upr,
          simint_lwr = tmp$simint_lwr,
          confint_lwr = ci$confint_lwr,
          confint_upr = ci$confint_upr
        )
      } else {
        rval[[i]] <- transform(tmp,
          observed = tmp$observed / (counts[i] * tmp$width),
          expected = tmp$expected / (counts[i] * tmp$width),
          simint_upr = tmp$simint_upr / (counts[i] * tmp$width),
          simint_lwr = tmp$simint_lwr / (counts[i] * tmp$width),
          confint_lwr = ci$confint_lwr / (counts[i] * tmp$width),
          confint_upr = ci$confint_upr / (counts[i] * tmp$width)
        )
      }

    } else {

      if (freq) {
        rval[[i]] <- transform(tmp,
          observed = tmp$observed,
          expected = tmp$expected,
          simint_upr = tmp$simint_upr,
          simint_lwr = tmp$simint_lwr
        )
      } else {
        rval[[i]] <- transform(tmp,
          observed = tmp$observed / (counts[i] * tmp$width),
          expected = tmp$expected / (counts[i] * tmp$width),
          simint_upr = tmp$simint_upr / (counts[i] * tmp$width),
          simint_lwr = tmp$simint_lwr / (counts[i] * tmp$width)
        )
      }
    }
  }

  rval <- do.call("rbind", rval)

  if (!return_group) {
    rval$group <- NULL
  }

  ## set attributes
  attr(rval, "simint") <- attr(object, "simint")
  attr(rval, "confint") <- attr(object, "confint")
  attr(rval, "expected") <- attr(object, "expected")
  attr(rval, "xlab") <- attr(object, "xlab")
  attr(rval, "ylab") <- attr(object, "ylab")
  attr(rval, "main") <- attr(object, "main")
  attr(rval, "type") <- attr(object, "type")
  attr(rval, "scale") <- scale
  attr(rval, "style") <- attr(object, "style")
  attr(rval, "freq") <- attr(object, "freq")
  attr(rval, "counts") <- attr(object, "counts")

  ## return as `data.frame` or `tibble`
  if ("data.frame" %in% class(object)) {
    class(rval) <- c("pithist", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("pithist", class(rval))
  }

  rval
}


#' @export
print.pithist <- function(x, ...) {
  ## get arg `type`, `style` and `freq`
  type <- freq <- style <- NULL # needed for `use_arg_from_attributes()` # FIXME: (ML) Still needed?!
  style <- use_arg_from_attributes(x, "style", default = NULL, force_single = TRUE)
  type <- use_arg_from_attributes(x, "type", default = NULL, force_single = TRUE)
  freq <- use_arg_from_attributes(x, "freq", default = NULL, force_single = TRUE)

  ## return custom print statement
  if (is.null(type) || is.null(freq) || is.null(style)) {
    cat("A `pithist` object without mandatory attributes `type`, `freq` and `style`\n\n")
  } else if (all(c("confint_lwr", "confint_upr", "expected") %in% names(x))) {
    cat(
      paste0(
        sprintf(
          "A `pithist` object with `type = \"%s\"`, `freq = \"%s\"` and `style = \"%s\"`",
          type,
          freq,
          style
        ),
        " with columns: `confint_lwr`, `confint_upr` and `expected`\n\n"
      )
    )
  } else {
    cat(
      sprintf(
        "A `pithist` object with `type = \"%s\"`, `freq = \"%s\"` and `style = \"%s\"`\n\n",
        type,
        freq,
        style
      )
    )
  }

  ## call next print method
  NextMethod()
}
