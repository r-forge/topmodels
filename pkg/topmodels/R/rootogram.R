# -------------------------------------------------------------------
## Programming outline: Rootogram
# -------------------------------------------------------------------
# - Observed y in-sample or out-of-sample (n x 1)
# - Breaks for observations (m x 1)
# - Predicted probabilities at breaks F_y(br1), ..., F_y(brm) (n x m)
#
# - Cut observations at breaks -> observed frequencies for (m-1) groups
# - Aggregate probabilities -> expected frequencies for (m-1) groups
# - Can be drawn in different style (standing vs. hanging / raw vs. sqrt)

# Functions:
# - rootogram() generic plus default method
# - Return object of class "rootogram" that is plotted by default
# - But has plot=FALSE so that suitable methods can be added afterwards
# - Methods: plot(), autoplot(), c()/rbind(), +
# -------------------------------------------------------------------


#' Rootograms for Assessing Goodness of Fit of Probability Models
#'
#' Rootograms graphically compare (square roots) of empirical frequencies with
#' expected (fitted) frequencies from a probability model. If \code{plot = TRUE}, the
#' resulting object of class \code{"rootogram"} is plotted by
#' \code{\link{plot.rootogram}} or \code{\link{autoplot.rootogram}} before it is
#' returned, depending on whether the package \code{ggplot2} is loaded.
#'
#' Rootograms graphically compare frequencies of empirical distributions and
#' expected (fitted) probability models. For the observed distribution the histogram is
#' drawn on a square root scale (hence the name) and superimposed with a line
#' for the expected frequencies. The histogram can be \code{"standing"} on the
#' x-axis (as usual), or \code{"hanging"} from the expected curve, or a
#' \code{"suspended"} histogram of deviations can be drawn.
#'
#' The function \code{\link{rootogram}} leverages the \code{\link{procast}}
#' generic in order to compute all necessary coordinates based on observed and
#' expected (fitted) frequencies.
#'
#' In addition to the \code{plot} and \code{\link[ggplot2]{autoplot}} method for
#' rootogram objects, it is also possible to combine two (or more) rootograms by
#' \code{c}/\code{rbind}, which creates a set of rootograms that can then be
#' plotted in one go.
#'
#' @aliases rootogram rootogram.default c.rootogram rbind.rootogram
#' @param object an object from which an rootogram can be extracted with
#' \code{\link{procast}}.
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
#' @param response_type character. To set the default values for \code{breaks} and
#' \code{widths}.  Currently different defaults are available for \code{"discrete"}
#' and \code{"continous"} response distribution, as well as for the special case of a
#' \code{"logseries"} response.
#' @param breaks numeric. Breaks for the histogram intervals.
#' @param width numeric. Widths of the histogram bars.
#' @param style character specifying the syle of rootogram (see below).
#' @param scale character specifying whether raw frequencies or their square
#' roots (default) should be drawn.
#' @param expected Should the expected (fitted) frequencies be plotted? Either logical or as character string defining one of `"both"`, `"line"` or `"point"`.
#' @param confint logical. Should confident intervals be drawn?
#' @param ref logical. Should a reference line be plotted?
#' @param xlab,ylab,main graphical parameters.
#' @param \dots further graphical parameters passed to the plotting function.
#' @return An object of class \code{"rootogram"} inheriting from
#' \code{"data.frame"} or \code{"tibble"} conditional on the argument \code{class}
#' with the following variables: \item{observed}{observed
#' frequencies,} \item{expected}{expected (fitted) frequencies,} \item{x}{histogram
#' interval midpoints on the x-axis,} \item{y}{bottom coordinate of the
#' histogram bars,} \item{width}{widths of the histogram bars,}
#' \item{height}{height of the histogram bars,} \item{line}{y-coordinates of
#' the fitted curve,} \item{confint_lwr, confint_upr}{lower and upper confidence interval bound.} 
#' Additionally, \code{style}, \code{scale}, \code{xlab},
#' \code{ylab} and \code{main}, and \code{confint_level} are stored as attributes.
#' @note Note that there is also a \code{\link[vcd]{rootogram}} function in the
#' \pkg{vcd} package that is similar to the \code{numeric} method provided
#' here. However, it is much more limited in scope, hence a function has been
#' created here.
#' @seealso \code{\link{plot.rootogram}}, \code{\link{procast}}
#' @references Friendly M (2000), \emph{Visualizing Categorical Data}. SAS
#' Institute, Cary.
#'
#' Kleiber C, Zeileis A (2016).  \dQuote{Visualizing Count Data Regressions
#' Using Rootograms.} \emph{The American Statistician}, \bold{70}(3), 296--303.
#' \doi{10.1080/00031305.2016.1173590}.
#'
#' Tukey JW (1977). \emph{Exploratory Data Analysis}. Addison-Wesley, Reading.
#' @keywords hplot
#' @examples
#'
#' ## plots and output
#'
#' ## number of deaths by horsekicks in Prussian army (Von Bortkiewicz 1898)
#' deaths <- rep(0:4, c(109, 65, 22, 3, 1))
#'
#' ## fit glm model
#' m1_pois <- glm(deaths ~ 1, family = poisson)
#' rootogram(m1_pois)
#'
#' ## inspect output (without plotting)
#' r1 <- rootogram(m1_pois, plot = FALSE)
#' r1
#'
#' ## combine plots
#' plot(c(r1, r1), col = c(1, 2), expected_col = c(1, 2))
#'
#'
#' #-------------------------------------------------------------------------------
#' ## different styles
#'
#' ## artificial data from negative binomial (mu = 3, theta = 2)
#' ## and Poisson (mu = 3) distribution
#' set.seed(1090)
#' y <- rnbinom(100, mu = 3, size = 2)
#' x <- rpois(100, lambda = 3)
#'
#' ## glm method: fitted values via glm()
#' m2_pois <- glm(y ~ x, family = poisson)
#'
#' ## correctly specified Poisson model fit
#' par(mfrow = c(1, 3))
#' r1 <- rootogram(m2_pois, style = "standing", ylim = c(-2.2, 4.8), main = "Standing")
#' r2 <- rootogram(m2_pois, style = "hanging", ylim = c(-2.2, 4.8), main = "Hanging")
#' r3 <- rootogram(m2_pois, style = "suspended", ylim = c(-2.2, 4.8), main = "Suspended")
#' par(mfrow = c(1, 1))
#'
#' #-------------------------------------------------------------------------------
#' ## linear regression with normal/Gaussian response: anorexia data
#'
#' data("anorexia", package = "MASS")
#'
#' m3_gauss <- glm(Postwt ~ Prewt + Treat + offset(Prewt), family = gaussian, data = anorexia)
#'
#' ## plot rootogram as "ggplot2" graphic
#' rootogram(m3_gauss, plot = "ggplot2")
#' @export
rootogram <- function(object, ...) {
  UseMethod("rootogram")
}


#' @rdname rootogram
#' @method rootogram default
#' @export
rootogram.default <- function(
                              ## computation arguments
                              object,
                              newdata = NULL,
                              plot = TRUE,
                              class = NULL,
                              response_type = NULL,
                              breaks = NULL,
                              width = NULL,

                              ## plotting arguments
                              style = c("hanging", "standing", "suspended"),
                              scale = c("sqrt", "raw"),
                              expected = TRUE,
                              confint = TRUE,
                              ref = TRUE,
                              xlab = NULL,
                              ylab = NULL,
                              main = NULL,
                              ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * `object`, `newdata` w/i `newresponse()`
  ## * `expected`, `ref`, `confint`...` in `plot()` and `autoplot()`
  stopifnot(is.null(breaks) || (is.numeric(breaks) && is.null(dim(breaks))))
  stopifnot(is.null(width) || (is.numeric(width) && length(width) == 1))
  stopifnot(is.null(xlab) || (length(xlab) == 1 && is.character(xlab)))
  stopifnot(is.null(ylab) || (length(ylab) == 1 && is.character(ylab)))
  stopifnot(is.null(main) || (length(main) == 1 && is.character(main)))

  ## match arguments
  scale <- match.arg(scale)
  style <- match.arg(style)

  ## default annotation
  if (is.null(xlab)) {
    xlab <- as.character(attr(terms(object), "variables"))[2L]
  }
  if (is.null(ylab)) {
    ylab <- if (scale == "raw") "Frequency" else "sqrt(Frequency)"
  }
  if (is.null(main)) {
    main <- deparse(substitute(object))
  }

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
    "The argument `plot` must be logical or match the arguments 'none', 'base' or 'ggplot2'." =
      !inherits(plot, "try-error")
  )

  ## guess output class
  if (is.null(class)) {
    class <- if ("tibble" %in% .packages()) "tibble" else "data.frame"
  }
  class <- try(match.arg(class, c("tibble", "data.frame")))
  stopifnot(
    "The argument `class` must be NULL or match the arguments 'tibble' or 'data.frame'." =
      !inherits(class, "try-error")
  )

  # -------------------------------------------------------------------
  # PREPARE DATA
  # -------------------------------------------------------------------
  ## get data and weights
  y <- newresponse(object, newdata = newdata, na.action = na.pass)
  w <- attr(y, "weights")
  if (is.null(response_type)) response_type <- attr(y, "response_type")
  response_type <- match.arg(response_type, c("discrete", "logseries", "continuous"))

  ## set breaks and midpoints
  ## TODO: (ML) Extend breaks to the left, in case still expected frequency exists.
  ## TODO: (Z) Try to get rid of 'response_type'.
  if (is.null(breaks) && response_type == "discrete") {
    breaks <- -1L:max(y[w > 0]) + 0.5
  } else if (is.null(breaks) && response_type == "logseries") {
    breaks <- 0L:max(y[w > 0]) + 0.5
  } else if (is.null(breaks)) {
    breaks <- "Sturges"
  }

  breaks <- hist(y[w > 0], plot = FALSE, breaks = breaks)$breaks
  mid <- (head(breaks, -1L) + tail(breaks, -1L)) / 2

  ## fix pointmasses
  ## TODO: (ML) Check if that always works or could be improved.
  breaks[1] <- breaks[1] - 1e-12
  breaks[length(breaks)] <- breaks[length(breaks)] + 1e-12

  ## set widths
  if (is.null(width) && (response_type == "discrete" || response_type == "logseries")) {
    width <- 0.9
  } else if (is.null(width)) {
    width <- 1
  }

  # -------------------------------------------------------------------
  # COMPUTATION OF EXPECTED AND OBSERVED FREQUENCIES
  # -------------------------------------------------------------------
  ## expected frequencies (part1)
  p <- matrix(NA, nrow = length(y), ncol = length(breaks) - 1L)
  for (i in 1L:ncol(p)) {
    p[, i] <-
      procast(object,
        newdata = newdata, na.action = na.pass, type = "probability",
        at = breaks[i + 1L], drop = TRUE
      ) -
      procast(object,
        newdata = newdata, na.action = na.pass, type = "probability",
        at = breaks[i], drop = TRUE
      )
  }

  ## TODO: (ML) Sometime neg. expected occurs (see example for underdispersed model fit).
  p[abs(p) < sqrt(.Machine$double.eps)] <- 0

  ## handle NAs
  ## TODO: (ML) Maybe allow arg `na.action` in the future.
  idx_not_na <- as.logical(complete.cases(y) * complete.cases(p))
  y <- y[idx_not_na]
  p <- p[idx_not_na, ]
  w <- w[idx_not_na]

  ## observed frequencies
  val_observed <- as.vector(xtabs(w ~ cut(y, breaks, include.lowest = TRUE)))

  ## expected frequencies (part2)
  val_expected <- colSums(p * w)

  # -------------------------------------------------------------------
  # OUTPUT AND OPTIONAL PLOTTING
  # -------------------------------------------------------------------
  ## collect everything as data.frame
  rval <- data.frame(
    observed = if (scale == "raw") as.vector(val_observed) else sqrt(as.vector(val_observed)),
    expected = if (scale == "raw") as.vector(val_expected) else sqrt(as.vector(val_expected)),
    mid = mid,
    width = diff(breaks) * width
  )

  ## add attributes
  attr(rval, "style") <- style
  attr(rval, "scale") <- scale
  attr(rval, "expected") <- expected
  attr(rval, "confint") <- confint
  attr(rval, "ref") <- ref
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main

  ## set class to data.frame or tibble
  if (class == "data.frame") {
    class(rval) <- c("rootogram", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("rootogram", class(rval))
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
c.rootogram <- function(...) {
  # -------------------------------------------------------------------
  # GET DATA
  # -------------------------------------------------------------------
  ## list of rootograms
  rval <- list(...)

  ## set class to tibble if any rval is a tibble
  if (any(do.call("c", lapply(rval, class)) %in% "tbl")) {
    class <- "tibble"
  } else {
    class <- "data.frame"
  }

  ## remove temporary the class (needed below for `c()`)
  ## TODO: (ML) Rewrite by, e.g., employing `lapply()`.
  for (i in 1:length(rval)) class(rval[[i]]) <- class(rval[[i]])[!class(rval[[i]]) %in% "rootogram"]

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
  nam <- names(rval)
  main <- if (is.null(nam)) {
    prepare_arg_for_attributes(rval, "main")
  } else {
    make.unique(rep.int(nam, sapply(n, length)))
  }

  ## parameters
  style <- prepare_arg_for_attributes(rval, "style", force_single = TRUE)
  scale <- prepare_arg_for_attributes(rval, "scale", force_single = TRUE)
  expected <- prepare_arg_for_attributes(rval, "expected")
  ref <- prepare_arg_for_attributes(rval, "ref")
  confint <- prepare_arg_for_attributes(rval, "confint")

  ## fix `ylabel` according to possible new `scale`
  if (scale == "sqrt") {
    ylab[grepl("^Frequency$", ylab)] <- "sqrt(Frequency)"
  } else if (scale == "raw") {
    ylab[grepl("^sqrt\\(Frequency\\)$", ylab)] <- "Frequency"
  }

  # -------------------------------------------------------------------
  # RETURN DATA
  # -------------------------------------------------------------------
  ## get all names needed if extended object should be computed and for combination
  all_names <- unique(unlist(lapply(rval, names)))

  ## get both objects on the same scale
  if (any(grepl("ymin|ymax|cofint_lwr|confint_upr", all_names))) {
    rval <- lapply(rval, summary.rootogram, style = style, scale = scale, extend = TRUE)
  } else {
    rval <- lapply(rval, summary.rootogram, style = style, scale = scale, extend = FALSE)
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

  ## add attributes
  attr(rval, "style") <- style
  attr(rval, "scale") <- scale
  attr(rval, "expected") <- expected
  attr(rval, "confint") <- confint
  attr(rval, "ref") <- ref
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main

  ## set class to data.frame or tibble
  if (class == "data.frame") {
    class(rval) <- c("rootogram", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("rootogram", class(rval))
  }

  ## return
  return(rval)
}


#' @export
rbind.rootogram <- c.rootogram


#' S3 Methods for Plotting Rootograms
#'
#' Generic plotting functions for rootograms of the class \code{"rootogram"}
#' computed by \code{link{rootogram}}.
#'
#' Rootograms graphically compare (square roots) of empirical frequencies with
#' expected (fitted) frequencies from a probability model.
#'
#' Rootograms graphically compare frequencies of empirical distributions and
#' expected (fitted) probability models. For the observed distribution the histogram is
#' drawn on a square root scale (hence the name) and superimposed with a line
#' for the expected frequencies. The histogram can be \code{"standing"} on the
#' x-axis (as usual), or \code{"hanging"} from the expected (fitted) curve, or a
#' \code{"suspended"} histogram of deviations can be drawn.
#'
#' @aliases plot.rootogram autoplot.rootogram
#' @param x,object an object of class \code{\link{rootogram}}.
#' @param style character specifying the syle of rootogram.
#' @param scale character specifying whether raw frequencies or their square
#' roots (default) should be drawn.
#' @param expected Should the expected (fitted) frequencies be plotted?
#' @param ref logical. Should a reference line be plotted?
#' @param confint logical. Should confident intervals be drawn?
#' @param confint_level numeric. The confidence level required.
#' @param confint_type character. Should \code{"pointwise"} or \code{"simultaneous"} confidence intervals be visualized. 
#' @param confint_nrep numeric. The repetition number of simulation for computing the confidence intervals.
#' @param xlim,ylim,xlab,ylab,main,axes,box graphical parameters.
#' @param col,border,lwd,lty,alpha_min graphical parameters for the histogram style part of the base plot.
#' @param colour,fill,size,linetype,alpha graphical parameters for the histogram style part in the \code{autoplot}.
#' @param legend logical. Should a legend be added in the \code{ggplot2} style
#' graphic?
#' @param theme Which `ggplot2` theme should be used. If not set, \code{\link[ggplot2]{theme_bw}} is employed.
#' @param expected_col,expected_pch,expected_lty,expected_lwd,ref_col,ref_lty,ref_lwd,expected_colour,expected_size,expected_linetype,expected_alpha,expected_fill,expected_stroke,expected_shape,ref_colour,ref_size,ref_linetype,ref_alpha,confint_col,confint_lty,confint_lwd,confint_colour,confint_size,confint_linetype,confint_alpha Further graphical parameters for the `expected` and `ref` line using either \code{\link[ggplot2]{autoplot}} or \code{plot}.
#' @param \dots further graphical parameters passed to the plotting function.
#' @seealso \code{\link{rootogram}}, \code{\link{procast}}
#' @references Friendly M (2000), \emph{Visualizing Categorical Data}. SAS
#' Institute, Cary.
#'
#' Kleiber C, Zeileis A (2016).  \dQuote{Visualizing Count Data Regressions
#' Using Rootograms.} \emph{The American Statistician}, \bold{70}(3), 296--303.
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_expr_doi(\"#1\")}",
#' "10.1080/00031305.2016.1173590")\Sexpr{tools:::Rd_expr_doi("10.1080/00031305.2016.1173590")}.
#'
#' Tukey JW (1977). \emph{Exploratory Data Analysis}. Addison-Wesley, Reading.
#' @keywords hplot
#' @examples
#'
#' ## speed and stopping distances of cars
#' m1_lm <- lm(dist ~ speed, data = cars)
#'
#' ## compute and plot rootogram
#' rootogram(m1_lm)
#'
#' ## customize colors
#' rootogram(m1_lm, ref_col = "blue", lty = 2, pch = 20)
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
#'   ### compute rootograms FIXME
#'   #r2_lm <- rootogram(m2_lm, plot = FALSE)
#'   #r2_crch <- rootogram(m2_crch, plot = FALSE)
#'
#'   ### plot in single graph
#'   #plot(c(r2_lm, r2_crch), col = c(1, 2))
#' }
#'
#' #-------------------------------------------------------------------------------
#' ## determinants for male satellites to nesting horseshoe crabs
#' data("CrabSatellites", package = "countreg")
#'
#' ## linear poisson model
#' m3_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#'
#' ## compute and plot rootogram as "ggplot2" graphic
#' rootogram(m3_pois, plot = "ggplot2")
#'
#' #-------------------------------------------------------------------------------
#' ## artificial data from negative binomial (mu = 3, theta = 2)
#' ## and Poisson (mu = 3) distribution
#' set.seed(1090)
#' y <- rnbinom(100, mu = 3, size = 2)
#' x <- rpois(100, lambda = 3)
#'
#' ## glm method: fitted values via glm()
#' m4_pois <- glm(y ~ x, family = poisson)
#'
#' ## correctly specified Poisson model fit
#' par(mfrow = c(1, 3))
#' r4a_pois <- rootogram(m4_pois, style = "standing", ylim = c(-2.2, 4.8), main = "Standing")
#' r4b_pois <- rootogram(m4_pois, style = "hanging", ylim = c(-2.2, 4.8), main = "Hanging")
#' r4c_pois <- rootogram(m4_pois, style = "suspended", ylim = c(-2.2, 4.8), main = "Suspended")
#' par(mfrow = c(1, 1))
#' @export
plot.rootogram <- function(x,
                           style = NULL,
                           scale = NULL,
                           expected = NULL,
                           ref = NULL,
                           confint = NULL,
                           confint_level = 0.95,
                           confint_type = c("pointwise", "simultaneous"),
                           confint_nrep = 1000,
                           xlim = c(NA, NA),
                           ylim = c(NA, NA),
                           xlab = NULL,
                           ylab = NULL,
                           main = NULL,
                           axes = TRUE,
                           box = FALSE,
                           col = "darkgray",
                           border = "black",
                           lwd = 1,
                           lty = 1,
                           alpha_min = 0.8,
                           expected_col = 2,
                           expected_pch = 19,
                           expected_lty = 1,
                           expected_lwd = 2,
                           confint_col = "black",
                           confint_lty = 2,
                           confint_lwd = 1.75,
                           ref_col = "black",
                           ref_lty = 1,
                           ref_lwd = 1.25,
                           ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## check if ylab is defined
  ylab_missing <- missing(ylab)

  ## get default arguments
  style <- use_arg_from_attributes(x, "style", default = "hanging", force_single = TRUE)
  scale <- use_arg_from_attributes(x, "scale", default = "sqrt", force_single = TRUE)
  expected <- use_arg_from_attributes(x, "expected", default = TRUE, force_single = FALSE)
  ref <- use_arg_from_attributes(x, "ref", default = TRUE, force_single = FALSE)
  confint <- use_arg_from_attributes(x, "confint", default = TRUE, force_single = FALSE)

  ## sanity checks
  ## * lengths of most arguments are checked by recycling
  ## * graphical parameters are checked w/i function calls for plotting
  stopifnot(is.logical(ref))
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(is.logical(confint))
  stopifnot(
    is.numeric(confint_level),
    length(confint_level) == 1,
    confint_level >= 0,
    confint_level <= 1
  )
  stopifnot(
    is.numeric(confint_nrep),
    length(confint_nrep) == 1,
    confint_nrep >= 0
  )

  ## match arguments
  scale <- match.arg(scale, c("sqrt", "raw"))
  style <- match.arg(style, c("hanging", "standing", "suspended"))
  confint_type <- match.arg(confint_type)

  ## extend input object on correct scale (compute heights, ...)
  x <- summary(
    x, 
    scale = scale, 
    style = style,
    confint_level = confint_level,
    confint_type = confint_type,
    confint_nrep = confint_nrep
  )

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  # -------------------------------------------------------------------
  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING
  # -------------------------------------------------------------------
  ## determine in which style the `expected` should be plotted
  expected[expected == FALSE | expected == "FALSE"] <- "none"
  expected[expected == TRUE | expected == "TRUE"] <- "both"
  expected <- c("none", "b", "l", "p")[match(expected, c("none", "both", "line", "point"))]
  stopifnot(all(expected %in% c("none", "b", "l", "p")))

  ## prepare xlim and ylim
  if (is.list(xlim)) xlim <- as.data.frame(do.call("rbind", xlim))
  if (is.list(ylim)) ylim <- as.data.frame(do.call("rbind", ylim))

  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(1:n, expected, ref, confint,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]],
    axes, box,
    border, col, lwd, lty, alpha_min,
    expected_col, expected_pch, expected_lty, expected_lwd, 
    ref_col, ref_lty, ref_lwd,
    confint_col, confint_lty, confint_lwd
  )[, -1]

  ## prepare annotation
  xlab <- use_arg_from_attributes(x, "xlab", default = "Rootogram", force_single = FALSE)
  ylab <- use_arg_from_attributes(x, "ylab",
    default = if (scale == "raw") "Frequency" else "sqrt(Frequency)", force_single = FALSE
  )
  main <- use_arg_from_attributes(x, "main", default = "model", force_single = FALSE)

  ## fix `ylab` according to possible new `freq`
  if (ylab_missing && scale == "sqrt") {
    ylab[grepl("^Frequency$", ylab)] <- "sqrt(Frequency)"
  } else if (ylab_missing && scale == "raw") {
    ylab[grepl("^sqrt\\(Frequency\\)$", ylab)] <- "Frequency"
  }

  # -------------------------------------------------------------------
  # PREPARE DATA FOR PLOTTING
  # -------------------------------------------------------------------

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION
  # -------------------------------------------------------------------
  rootogram_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## rect elements
    ybottom <- d$ymin
    ytop <- d$ymax
    xleft <- d$mid - d$width / 2
    xright <- d$mid + d$width / 2

    ## step elements (only needed for expected, confint)
    z <- compute_breaks(d$mid, d$width, offset = TRUE)

    ## get xlim and ylim
    ylim_idx <- c(is.na(plot_arg$ylim1[j]), is.na(plot_arg$ylim2[j]))
    xlim_idx <- c(is.na(plot_arg$xlim1[j]), is.na(plot_arg$xlim2[j]))
    if (any(xlim_idx)) {
      plot_arg[j, c("xlim1", "xlim2")[xlim_idx]] <- range(c(xleft, xright))[xlim_idx]
    }
    if (any(ylim_idx)) {
      ylim_use <- rep(
        c(TRUE, TRUE, plot_arg$expected[j] != "none", plot_arg$confint[j], plot_arg$confint[j]),
          each = NROW(d)
      )
      plot_arg[j, c("ylim1", "ylim2")[ylim_idx]] <- range(
        c(d$ymin, d$ymax, d$expected, d$confint_lwr, d$confint_upr)[ylim_use],
        na.rm = TRUE
      )[ylim_idx]
    }

    ## trigger plot
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

    ## plot rootogram
    rect(xleft, ybottom, xright, ytop,
      col = set_minimum_transparency(plot_arg$col[j], alpha_min = plot_arg$alpha_min[j]),
      border = plot_arg$border[j],
      lty = plot_arg$lty[j],
      lwd = plot_arg$lwd[j]
    )

    ## add ref line
    if (plot_arg$ref[j]) {
      abline(
        h = 0,
        col = plot_arg$ref_col[j],
        lty = plot_arg$ref_lty[j],
        lwd = plot_arg$ref_lwd[j]
      )
    }

    ## add expected line
    if (plot_arg$expected[j] != "none") {
      lines(
        d$mid,
        d$expected,
        col =  plot_arg$expected_col[j],
        pch = plot_arg$expected_pch[j],
        type = plot_arg$expected[j],
        lty = plot_arg$expected_lty[j],
        lwd = plot_arg$expected_lwd[j]
      )
    }

    if (plot_arg$confint[j]) {
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
    }
  }

  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  ## set up necessary panels
  if (n > 1L) {
    old_pars <- par(mfrow = n2mfrow(n))
    on.exit(par(old_pars), add = TRUE)
  }

  ## draw rootograms
  for (i in 1L:n) rootogram_plot(x[x$group == i, ], ...)
}


#' @rdname plot.rootogram
#' @method autoplot rootogram
#' @exportS3Method ggplot2::autoplot
autoplot.rootogram <- function(object,
                               style = NULL,
                               scale = NULL,
                               expected = NULL,
                               ref = NULL,
                               confint = NULL,
                               confint_level = 0.95,
                               confint_type = c("pointwise", "simultaneous"),
                               confint_nrep = 1000,
                               xlim = c(NA, NA),
                               ylim = c(NA, NA),
                               xlab = NULL,
                               ylab = NULL,
                               main = NULL,
                               legend = FALSE,
                               theme = NULL,
                               colour = "black",
                               fill = "darkgray",
                               size = 0.5,
                               linetype = 1,
                               alpha = NA,
                               expected_colour = 2,
                               expected_size = 1,
                               expected_linetype = 1,
                               expected_alpha = 1,
                               expected_fill = NA,
                               expected_stroke = 0.5,
                               expected_shape = 19,
                               confint_colour = "black",
                               confint_size = 0.5,
                               confint_linetype = 2,
                               confint_alpha = NA,
                               ref_colour = "black",
                               ref_size = 0.5,
                               ref_linetype = 1,
                               ref_alpha = NA,
                               ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## check if ylab is defined
  ylab_missing <- missing(ylab)

  ## get default arguments
  style <- use_arg_from_attributes(object, "style", default = "hanging", force_single = TRUE)
  scale <- use_arg_from_attributes(object, "scale", default = "sqrt", force_single = TRUE)
  expected <- use_arg_from_attributes(object, "expected", default = TRUE, force_single = TRUE)
  ref <- use_arg_from_attributes(object, "ref", default = TRUE, force_single = TRUE)
  confint <- use_arg_from_attributes(object, "confint", default = TRUE, force_single = TRUE)
  xlab <- use_arg_from_attributes(object, "xlab", default = "Rootogram", force_single = TRUE)
  ylab <- use_arg_from_attributes(object, "ylab",
    default = if (scale == "raw") "Frequency" else "sqrt(Frequency)", force_single = TRUE
  )

  ## fix `ylabel` according to possible new `scale`
  if (ylab_missing) {
    if (scale == "sqrt" && grepl("^Frequency$", ylab)) ylab <- "sqrt(Frequency)"
    if (scale == "raw" && grepl("^sqrt\\(Frequency\\)$", ylab)) ylab <- "Frequency"
  }

  ## get base style arguments
  add_arg <- list(...)
  if (!is.null(add_arg$lwd)) size <- add_arg$lwd
  if (!is.null(add_arg$lty)) linetype <- add_arg$lty

  ## sanity checks
  stopifnot(is.logical(ref))
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(is.logical(confint))
  stopifnot(
    is.numeric(confint_level),
    length(confint_level) == 1,
    confint_level >= 0,
    confint_level <= 1
  )
  stopifnot(
    is.numeric(confint_nrep),
    length(confint_nrep) == 1,
    confint_nrep >= 0
  )

  ## match arguments
  confint_type <- match.arg(confint_type)

  ## get line style for expected
  if (isFALSE(expected)) {
    expected <- "none"
  } else if (isTRUE(expected)) {
    expected <- "both"
  }
  expected <- match.arg(expected, c("none", "both", "line", "point"))

  ## transform to correct scale
  object <- summary.rootogram(object, scale = scale, style = style, extend = FALSE)

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

  ## get main and transform to the right length (must be done after handling of `title`)
  main <- use_arg_from_attributes(object, "main", default = "model", force_single = FALSE)
  stopifnot(is.character(main))
  main <- make.names(rep_len(main, n), unique = TRUE)

  ## prepare grouping
  object$group <- factor(object$group, levels = 1L:n, labels = main)

  # -------------------------------------------------------------------
  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING
  # -------------------------------------------------------------------
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
      observed = "observed",
      expected = "expected",
      mid = "mid",
      width = "width",
      group = "group"
    )
  )

  ## add rootogram
  rval <- rval +
    geom_rootogram(
      ggplot2::aes_string(
        colour = "group",
        fill = "group",
        size = "group",
        linetype = "group",
        alpha = "group"
      ),
      scale = "raw",
      style = style
    )

  ## add expected line
  if (expected != "none") {
    rval <- rval +
      geom_rootogram_expected(
        scale = "raw",
        linestyle = expected,
        colour = expected_colour,
        size = expected_size,
        linetype = expected_linetype,
        alpha = expected_alpha,
        fill = expected_fill,
        stroke = expected_stroke,
        shape = expected_shape,
      )
  }

  ## add ref
  if (ref) {
    rval <- rval +
      geom_rootogram_ref(
        colour = ref_colour,
        size = ref_size,
        linetype = ref_linetype,
        alpha = ref_alpha,
      )
  }

  ## add confint
  if (confint) { 
    rval <- rval +
      geom_rootogram_confint(
        level = confint_level,
        type = confint_type,
        scale = scale,
        rootogram_style = style,
        colour = confint_colour,
        size = confint_size,
        linetype = confint_linetype,
        alpha = confint_alpha
      )
  }

  ## set the colors, shapes, etc.
  rval <- rval +
    ggplot2::scale_colour_manual(values = plot_arg$colour) +
    ggplot2::scale_fill_manual(values = plot_arg$fill) +
    ggplot2::scale_size_manual(values = plot_arg$size) +
    ggplot2::scale_linetype_manual(values = plot_arg$linetype) +
    ggplot2::scale_alpha_manual(values = plot_arg$alpha)

  ## annotation
  rval <- rval + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

  ## add legend
  if (legend) {
    rval <- rval +
      ggplot2::labs(
        colour = "Model",
        fill = "Model",
        size = "Model",
        linetype = "Model",
        alpha = "Model"
      ) +
      ggplot2::guides(
        colour = "legend",
        fill = "legend",
        size = "legend",
        linetype = "legend",
        alpha = "legend"
      )
  } else {
    rval <- rval +
      ggplot2::guides(
        colour = "none",
        fill = "none",
        size = "none",
        linetype = "none",
        alpha = "none"
      )
  }

  ## set x and y limits
  rval <- rval +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE) +
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
  # GROUPING (IF ANY) AND RETURN PLOT
  # -------------------------------------------------------------------
  ## grouping
  if (n > 1L) {
    rval <- rval + ggplot2::facet_grid(group ~ .)
  } else if (!is.null(object$title)) {
    rval <- rval + ggplot2::facet_wrap(title ~ .)
  }

  ## return ggplot object
  rval
}


# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_rootogram()`
# -------------------------------------------------------------------

#' @rdname geom_rootogram
#' @export
stat_rootogram <- function(mapping = NULL,
                           data = NULL,
                           geom = "rootogram",
                           position = "identity",
                           na.rm = FALSE,
                           show.legend = NA,
                           inherit.aes = TRUE,
                           scale = c("sqrt", "raw"),
                           style = c("hanging", "standing", "suspended"),
                           ...) {
  scale <- match.arg(scale)
  style <- match.arg(style)

  ggplot2::layer(
    stat = StatRootogram,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      scale = scale,
      style = style,
      ...
    )
  )
}


#' @rdname geom_rootogram
#' @format NULL
#' @usage NULL
#' @export
StatRootogram <- ggplot2::ggproto("StatRootogram", ggplot2::Stat,
  required_aes = c("observed", "expected", "mid", "width"),
  compute_group = function(data,
                           scales,
                           scale = c("sqrt", "raw"),
                           style = c("hanging", "standing", "suspended")) {
    scale <- match.arg(scale)
    style <- match.arg(style)

    tmp_heights <- compute_rootogram_heights(data$expected, data$observed, scale = scale, style = style)

    data <- transform(data,
      xmin = mid - width / 2,
      xmax = mid + width / 2,
      ymin = tmp_heights$ymin,
      ymax = tmp_heights$ymax,
      observed = NULL,
      expected = NULL,
      mid = NULL,
      width = NULL
    )
    data
  }
)


#' \code{geom_*} and \code{stat_*} for Producing PIT Histograms with `ggplot2`
#'
#' Various \code{geom_*} and \code{stat_*} used within
#' \code{\link[ggplot2]{autoplot}} for producing PIT histograms.
#'
#' @inheritParams ggplot2::layer
#' @inheritParams ggplot2::geom_point
#' @param style character specifying the syle of rootogram (see below).
#' @param scale character specifying whether values should be transformed to the square root scale (not checking for original scale, so maybe applied again).
#' @param linestyle Character string defining one of `"both"`, `"line"` or `"point"`.
#' @param level numeric. The confidence level required.
#' @param type character. Should \code{"pointwise"} or \code{"simultaneous"} confidence intervals be visualized. 
#' @param nrep numeric. The repetition number of simulation for computing the confidence intervals.
#' @param rootogram_style character specifying the syle of rootogram.
#' @examples
#' if (require("ggplot2")) {
#'   ## Fit model
#'   data("CrabSatellites", package = "countreg")
#'   m1_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#'   m2_pois <- glm(satellites ~ color, data = CrabSatellites, family = poisson)
#'
#'   ## Compute rootogram (on raw scale)
#'   p1 <- rootogram(m1_pois, scale = "raw", plot = FALSE)
#'   p2 <- rootogram(m2_pois, scale = "raw", plot = FALSE)
#'
#'   d <- c(p1, p2)
#'
#'   ## Get label names
#'   main <- attr(d, "main")
#'   main <- make.names(main, unique = TRUE)
#'   d$group <- factor(d$group, labels = main)
#'
#'   ## Plot rootograms w/ on default "sqrt" scale
#'   gg1 <- ggplot(data = d) +
#'     geom_rootogram(aes(
#'       observed = observed, expected = expected, mid = mid,
#'       width = width, group = group
#'     )) +
#'     geom_rootogram_expected(aes(expected = expected, mid = mid)) +
#'     geom_rootogram_ref() +
#'     facet_grid(group ~ .) + 
#'     xlab("satellites") +
#'     ylab("sqrt(Frequency)")
#'   gg1
#' }
#' @export
geom_rootogram <- function(mapping = NULL,
                           data = NULL,
                           stat = "rootogram",
                           position = "identity",
                           na.rm = FALSE,
                           show.legend = NA,
                           inherit.aes = TRUE,
                           scale = c("sqrt", "raw"),
                           style = c("hanging", "standing", "suspended"),
                           ...) {
  scale <- match.arg(scale)
  style <- match.arg(style)

  ggplot2::layer(
    geom = GeomRootogram,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      scale = scale,
      style = style,
      ...
    )
  )
}


#' @rdname geom_rootogram
#' @format NULL
#' @usage NULL
#' @export
GeomRootogram <- ggplot2::ggproto("GeomRootogram", ggplot2::GeomRect,
  default_aes = ggplot2::aes(
    colour = "black", fill = "darkgray", size = 0.5, linetype = 1,
    alpha = NA
  ),
  required_aes = c("xmin", "xmax", "ymin", "ymax")
)


# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_rootogram_expected()`
# -------------------------------------------------------------------

#' @rdname geom_rootogram
#' @export
stat_rootogram_expected <- function(mapping = NULL,
                                  data = NULL,
                                  geom = "rootogram_expected",
                                  position = "identity",
                                  na.rm = FALSE,
                                  show.legend = NA,
                                  inherit.aes = TRUE,
                                  scale = c("sqrt", "raw"),
                                  ...) {
  scale <- match.arg(scale)
  style <- match.arg(style)

  ggplot2::layer(
    stat = StatRootogramExpected,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      scale = scale,
      ...
    )
  )
}


#' @rdname geom_rootogram
#' @format NULL
#' @usage NULL
#' @export
StatRootogramExpected <- ggplot2::ggproto("StatRootogramExpected", ggplot2::Stat,
  required_aes = c("expected", "mid"),
  compute_group = function(data,
                           scales,
                           scale = c("sqrt", "raw")) {
    scale <- match.arg(scale)

    ## raw vs. sqrt scale
    if (scale == "sqrt") {
      data <- transform(data,
        x = mid,
        y = sqrt(expected),
        expected = NULL,
        mid = NULL
      )
    } else {
      data <- transform(data,
        x = mid,
        y = expected,
        expected = NULL,
        mid = NULL
      )
    }
    data
  }
)

#' @rdname geom_rootogram
#' @export
geom_rootogram_expected <- function(mapping = NULL,
                                  data = NULL,
                                  stat = "rootogram_expected",
                                  position = "identity",
                                  na.rm = FALSE,
                                  show.legend = NA,
                                  inherit.aes = TRUE,
                                  scale = c("sqrt", "raw"),
                                  linestyle = c("both", "line", "point"),
                                  ...) {
  scale <- match.arg(scale)
  linestyle <- match.arg(linestyle)

  ggplot2::layer(
    geom = GeomRootogramExpected,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      scale = scale,
      linestyle = linestyle,
      ...
    )
  )
}


#' @rdname geom_rootogram
#' @export
GeomRootogramExpected <- ggplot2::ggproto("GeomRootogramExpected", ggplot2::GeomPath,
  default_aes = ggplot2::aes(
    colour = 2, size = 1, linetype = 1,
    alpha = 1, fill = NA, stroke = 0.5, shape = 19
  ),
  draw_panel = function(data, panel_params, coord, arrow = NULL,
                        lineend = "butt", linejoin = "round", linemitre = 10,
                        na.rm = FALSE,
                        linestyle = c("both", "line", "point")) {
    linestyle <- match.arg(linestyle)

    if (linestyle == "both") {
      ## TODO: (ML) Do not copy data.
      data2 <- transform(data, size = size * 2)

      grid::grobTree(
        ggplot2::GeomPath$draw_panel(data, panel_params, coord, arrow = NULL,
          lineend = "butt", linejoin = "round", linemitre = 10,
          na.rm = FALSE
        ),
        ggplot2::GeomPoint$draw_panel(data2, panel_params, coord, na.rm = FALSE)
      )
    } else if (linestyle == "line") {
      ggplot2::GeomPath$draw_panel(data, panel_params, coord, arrow = NULL,
        lineend = "butt", linejoin = "round", linemitre = 10,
        na.rm = FALSE
      )
    } else {
      data <- transform(data, size = size * 2)
      ggplot2::GeomPoint$draw_panel(data, panel_params, coord, na.rm = FALSE)
    }
  }
)


# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_rootogram_ref()`
# -------------------------------------------------------------------

#' @rdname geom_rootogram
#' @export
geom_rootogram_ref <- function(mapping = NULL,
                               data = NULL,
                               stat = "identity",
                               position = "identity",
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = TRUE,
                               ...) {
  ggplot2::layer(
    geom = GeomRootogramRef,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}


#' @rdname geom_rootogram
#' @format NULL
#' @usage NULL
#' @export
GeomRootogramRef <- ggplot2::ggproto("GeomRootogramRef", ggplot2::GeomHline,
  default_aes = ggplot2::aes(
    colour = "black",
    size = 0.5,
    linetype = 1,
    alpha = NA,
    yintercept = 0
  ),
  required_aes = NULL
)

# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_rootogram_confint()`
# -------------------------------------------------------------------

#' @rdname geom_rootogram
#' @export
stat_rootogram_confint <- function(mapping = NULL,
                                   data = NULL,
                                   geom = "rootogram_confint",
                                   position = "identity",
                                   na.rm = FALSE,
                                   show.legend = NA,
                                   inherit.aes = TRUE,
                                   level = 0.95,
                                   nrep = 1000,
                                   type = c("pointwise", "simultaneous"),
                                   scale = c("sqrt", "raw"),
                                   rootogram_style =  c("hanging", "standing", "suspended"),
                                   ...) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  rootogram_style <- match.arg(rootogram_style)

  ggplot2::layer(
    stat = StatRootogramConfint,
    mapping = mapping,
    data = data,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      level = level,
      nrep = nrep,
      type = type,
      scale = scale,
      rootogram_style = rootogram_style,
      ...
    )
  )
}


#' @rdname geom_rootogram
#' @format NULL
#' @usage NULL
#' @export
StatRootogramConfint <- ggplot2::ggproto("StatRootogramConfint", ggplot2::Stat,
  required_aes = c("observed", "expected", "mid", "width"),
  compute_group = function(data,
                           scales,
                           level = 0.95,
                           nrep = 1000,
                           type = "pointwise",
                           scale = "sqrt",
                           rootogram_style = "hanging") {

    ## compute ci interval
    ci <- compute_rootogram_confint(
      observed = data$observed,
      expected = data$expected,
      mid = data$mid,
      width = data$width,
      level = level,
      nrep = nrep,
      type = type,
      scale = scale,
      style = rootogram_style
    )

    ## return new data.frame condition on plotting `style`
    nd <- data.frame(
      x = compute_breaks(data$mid, data$width, offset = TRUE),
      ymin = duplicate_last_value(ci[[1]]),
      ymax = duplicate_last_value(ci[[2]])
    )
    nd
  }
)


#' @rdname geom_rootogram
#' @export
geom_rootogram_confint <- function(mapping = NULL,
                                   data = NULL,
                                   stat = "rootogram_confint",
                                   position = "identity",
                                   na.rm = FALSE,
                                   show.legend = NA,
                                   inherit.aes = TRUE,
                                   level = 0.95,
                                   nrep = 1000,
                                   type = c("pointwise", "simultaneous"),
                                   scale = c("sqrt", "raw"),
                                   rootogram_style =  c("hanging", "standing", "suspended"),
                                   ...) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  rootogram_style <- match.arg(rootogram_style)

  ggplot2::layer(
    geom = GeomRootogramConfint,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      level = level,
      nrep = nrep,
      type = type,
      scale = scale,
      rootogram_style = rootogram_style,
      ...
    )
  )
}


#' @rdname geom_rootogram
#' @format NULL
#' @usage NULL
#' @export
GeomRootogramConfint <- ggplot2::ggproto("GeomRootogramConfint", ggplot2::Geom,
  required_aes = c("x", "ymin", "ymax"),

  extra_params = c("na.rm"),

  # TODO: (ML) Does not vary for style; this is a copy of `GeomPolygon$handle_na()`
  handle_na = function(data, params) {
    data
  },

  ## Setting up all defaults needed for `GeomPolygon` and `GeomStep`
  default_aes = ggplot2::aes(
    colour = "black",
    size = 0.5,
    linetype = 2,
    alpha = NA,
    yintercept = 0
  ),

  draw_panel = function(data, panel_params, coord,
                        linejoin = "mitre", direction = "hv") {

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
  },
  draw_key = function(data, params, size) {
    draw_key_path(data, params, size)
  }
)


# -------------------------------------------------------------------
# HELPER FUNCTIONS FOR GETTING AN EXTENDED ROOTOGRAM OBJECT
# -------------------------------------------------------------------
compute_rootogram_confint <- function(object,
                                      observed,
                                      expected,
                                      mid,
                                      width,
                                      level = 0.95,
                                      nrep = 1000,
                                      type = c("pointwise", "simultaneous"),
                                      scale = c("sqrt", "raw"),
                                      style = c("hanging", "standing", "suspended")) {

  ## checks
  scale <- match.arg(scale)
  type <- match.arg(type)
  style <- match.arg(style)
  stopifnot(is.numeric(level), length(level) == 1, level >= 0, level <= 1)
  stopifnot(is.numeric(nrep), length(nrep) == 1, nrep >= 0)

  ## transform back to "raw"
  if (!missing(object)) {
    object <- summary(object, scale = "raw", style = style, extend = FALSE)
  } else if (!missing(observed) & !missing(expected) & !missing(mid) & !missing(width)) {
    object <- if (scale == "raw") {
      data.frame(observed, expected, mid, width) 
    } else {
      data.frame(observed^2, expected^2, mid, width)
    }
  } else {
    stop("No appropriate input data is given.")
  }

  ## helper function to compute one observed table (on input scale)
  ytab <- function(rgram) {
    y <- sample(rgram$mid, sum(rgram$observed), prob = rgram$expected,
      replace = TRUE)
    table(factor(y, levels = rgram$mid))
  }
  ## repeat nrep times
  ytab <- replicate(nrep, ytab(object))

  ## compute CIs
  if (scale == "sqrt") {
    if (type == "pointwise") {
      yq <- apply(sqrt(object$expected) - sqrt(ytab), 1, quantile, c((1 - level) / 2, 1 - (1 - level) / 2))
    } else {
      yq <- rbind(
        rep.int(quantile(apply(sqrt(object$expected) - sqrt(ytab), 2, min), (1 - level) / 2), nrow(object)),
        rep.int(quantile(apply(sqrt(object$expected) - sqrt(ytab), 2, max), 1 - (1 - level) / 2), nrow(object))
      )
    }

    if (style == "standing") yq <- rep(sqrt(object$expected), each = 2) + yq
  } else {
    if (type == "pointwise") {
      yq <- apply(object$expected - ytab, 1, quantile, c((1 - level) / 2, 1 - (1 - level) / 2))
    } else {
      yq <- rbind(
        rep.int(quantile(apply(object$expected - ytab, 2, min), (1 - level) / 2), nrow(object)),
        rep.int(quantile(apply(object$expected - ytab, 2, max), 1 - (1 - level) / 2), nrow(object))
      )
    }
    if (style == "standing") yq <- rep(object$expected, each = 2) + yq
  }

  ## return
  df <- data.frame(
    confint_lwr = yq[1, ],
    confint_upr = yq[2, ]
  )

  df
  
}

compute_rootogram_heights <- function(expected,
                                      observed,
                                      scale = c("sqrt", "raw"),
                                      style = c("hanging", "standing", "suspended")) {
  scale <- match.arg(scale)
  style <- match.arg(style)

  ## raw vs. sqrt scale
  if (scale == "sqrt") {
    y <- if (style == "hanging") sqrt(expected) - sqrt(observed) else 0
    height <- if (style == "suspended") sqrt(expected) - sqrt(observed) else sqrt(observed)
  } else {
    y <- if (style == "hanging") expected - observed else 0
    height <- if (style == "suspended") expected - observed else observed
  }

  data.frame(
    ymin = y,
    ymax = y + height
  )
}


#' @export
summary.rootogram <- function(object,
                              scale = NULL,
                              style = NULL,
                              confint_level = 0.95,
                              confint_type = c("pointwise", "simultaneous"),
                              confint_nrep = 1000,
                              extend = TRUE,
                              ...) {


  stopifnot(is.logical(extend))

  ## get arg `style` and `scale`
  scale_object <- attr(object, "scale")
  scale <- use_arg_from_attributes(object, "scale", default = "sqrt", force_single = TRUE)
  style <- use_arg_from_attributes(object, "style", default = "hanging", force_single = TRUE)

  scale <- match.arg(scale, c("sqrt", "raw"))
  style <- match.arg(style, c("hanging", "standing", "suspended"))

  if (scale != scale_object && scale_object == "raw") {
    object <- transform(object,
      observed = sqrt(object$observed),
      expected = sqrt(object$expected)
    )
  } else if (scale != scale_object && scale_object == "sqrt") {
    object <- transform(object,
      observed = object$observed^2,
      expected = object$expected^2
    )
  } else if (scale != scale_object) {
    stop('attribute `scale` must be unique and must match one of c("sqrt", "raw")`')
  }

  if (extend) {
    ## compute heights (must always be `scale = "raw"` as already transformed)
    tmp <- compute_rootogram_heights(object$expected, object$observed, scale = "raw", style = style)

    ## compute confidence intervals (must be done per each group)
    if (!any(grepl("group", names(object)))) {
      tmp2 <- compute_rootogram_confint(
        object, 
        level = confint_level, 
        type = confint_type,
        nrep = confint_nrep, 
        scale = scale, 
        style = style
      )
    } else {
      tmp2 <- list()
      for (i in unique(object$group)) {
        tmp2[[i]] <- compute_rootogram_confint(
          object[object$group == i, ], 
          level = confint_level, 
          type = confint_type,
          nrep = confint_nrep, 
          scale = scale, 
          style = style
        )
      }
      tmp2 <- do.call("rbind", tmp2)
    }

    rval <- transform(object,
      ymin = tmp$ymin,
      ymax = tmp$ymax,
      confint_lwr = tmp2$confint_lwr,
      confint_upr = tmp2$confint_upr
    )

  } else {
    rval <- object
  }

  ## set attributes
  attr(rval, "style") <- style
  attr(rval, "scale") <- scale
  attr(rval, "xlab") <- attr(object, "xlab")
  attr(rval, "ylab") <- attr(object, "ylab")
  attr(rval, "main") <- attr(object, "main")

  ## return as `data.frame` or `tibble`
  if ("data.frame" %in% class(object)) {
    class(rval) <- c("rootogram", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("rootogram", class(rval))
  }

  rval
}

#' @export
print.rootogram <- function(x, ...) {

  ## get arg `style` and `scale`
  scale <- style <- NULL # needed for `use_arg_from_attributes()` # FIXME: (ML) Still needed?!
  style <- use_arg_from_attributes(x, "style", default = NULL, force_single = TRUE)
  scale <- use_arg_from_attributes(x, "scale", default = NULL, force_single = TRUE)

  ## return custom print statement
  if (is.null(scale) || is.null(style)) {
    cat("A `rootogram` object without mandatory attributes `scale` and `style`\n\n")
  } else if (all(c("ymin", "ymax") %in% names(x))) {
    cat(
      paste0(
        sprintf(
          "A `rootogram` object with `scale = \"%s\"` and `style = \"%s\"`",
          scale,
          style
        ),
        " with columns: `ymin` and `ymax`\n\n"
      )
    )
  } else {
    cat(
      sprintf(
        "A `rootogram` object with `scale = '%s'` and `style = '%s'`\n\n",
        scale,
        style
      )
    )
  }

  ## call next print method
  NextMethod()
}
