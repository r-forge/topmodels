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
#' fitted frequencies from a probability model. If \code{plot = TRUE}, the
#' resulting object of class \code{"rootogram"} is plotted by
#' \code{\link{plot.rootogram}} or \code{\link{autoplot.rootogram}} before it is
#' returned, depending on whether the package \code{ggplot2} is loaded.
#' 
#' Rootograms graphically compare frequencies of empirical distributions and
#' fitted probability models. For the observed distribution the histogram is
#' drawn on a square root scale (hence the name) and superimposed with a line
#' for the fitted frequencies. The histogram can be \code{"standing"} on the
#' x-axis (as usual), or \code{"hanging"} from the fitted curve, or a
#' \code{"suspended"} histogram of deviations can be drawn.
#' 
#' The function \code{\link{rootogram}} leverages the \code{\link{procast}}
#' generic in order to compute all necessary coordinates based on observed and
#' fitted frequencies.
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
#' @param style character specifying the syle of rootogram (see below).
#' @param scale character specifying whether raw frequencies or their square
#' roots (default) should be drawn.
#' @param breaks numeric. Breaks for the histogram intervals.
#' @param width numeric. Widths of the histogram bars.
#' @param response_type character. To set the default values for \code{breaks} and
#' \code{widths}.  Currently different defaults are available for \code{"discrete"}
#' and \code{"continous"} response distribution, as well as for the special case of a 
#' \code{"logseries"} response.
#' @param xlab,ylab,main graphical parameters.
#' @param \dots further graphical parameters passed to the plotting function.
#' @return An object of class \code{"rootogram"} inheriting from
#' \code{"data.frame"} or \code{"tibble"} conditional on the argument \code{class} 
#' with the following variables: \item{observed}{observed
#' frequencies,} \item{expected}{fitted frequencies,} \item{x}{histogram
#' interval midpoints on the x-axis,} \item{y}{bottom coordinate of the
#' histogram bars,} \item{width}{widths of the histogram bars,}
#' \item{height}{height of the histogram bars,} \item{line}{y-coordinates of
#' the fitted curve.} Additionally, \code{style}, \code{scale}, \code{xlab},
#' \code{ylab}, and \code{main} are stored as attributes.
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
#' plot(c(r1, r1), col = c(1, 2), ref = 4, lty = c(1, 2))
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
#' r1 <- rootogram(m2_pois, style = "standing",  ylim = c(-2.2, 4.8), main = "Standing")
#' r2 <- rootogram(m2_pois, style = "hanging",   ylim = c(-2.2, 4.8), main = "Hanging")
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
#' 
#' @export
rootogram <- function(object, ...) {
  UseMethod("rootogram")
}


#' @rdname rootogram
#' @method rootogram default
#' @export
rootogram.default <- function(object,
                              newdata = NULL,
                              plot = TRUE,
                              class = NULL,
                              style = c("hanging", "standing", "suspended"),
                              scale = c("sqrt", "raw"),
                              breaks = NULL,
                              width = NULL,
                              response_type = NULL,
                              xlab = NULL,
                              ylab = NULL,
                              main = NULL,
                              ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * `object`, `newdata` w/i `newresponse()`
  ## * `breaks` w/i `hist()`
  ## * `...` in `plot()` and `autoplot()`
  stopifnot(is.null(breaks) || (is.numeric(breaks) && is.null(dim(breaks))))
  stopifnot(is.null(width) || (is.numeric(width) && length(width) == 1))
  stopifnot(length(xlab) == 1 | length(xlab) == 0)
  stopifnot(length(ylab) == 1 | length(ylab) == 0)
  stopifnot(length(main) == 1 | length(main) == 0)

  ## match arguments
  scale <- match.arg(scale)
  style <- match.arg(style)

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

  # -------------------------------------------------------------------
  # PREPARE DATA
  # -------------------------------------------------------------------
  ## get data and weights
  y <- newresponse(object, newdata = newdata, na.action = na.pass)
  w <- attr(y, "weights")
  if (is.null(response_type)) response_type <- attr(y, "response_type")
  response_type <- match.arg(response_type, c("discrete", "logseries", "continuous"))

  ## set breaks and midpoints
  ## FIXME: (ML) Extend breaks to the left, in case still expected frequency exists
  if (is.null(breaks) && response_type == "discrete") {
    breaks <- -1L:max(y[w > 0]) + 0.5
  } else if (is.null(breaks) && response_type == "logseries") {
    breaks <- 0L:max(y[w > 0]) + 0.5
  } else if (is.null(breaks)) {
    breaks <- "Sturges"
  }

  breaks <- hist(y[w > 0], plot = FALSE, breaks = breaks)$breaks
  mids <- (head(breaks, -1L) + tail(breaks, -1L)) / 2

  ## fix pointmasses 
  ## FIXME: (ML) Check if that always works or could be improved
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

  ## FIXME: (ML) Sometime neg. expected occured (see example for underdispersed model fit)
  p[abs(p) < sqrt(.Machine$double.eps)] <- 0

  ## handle NAs
  ## TODO: (ML) Maybe allow arg `na.action` in the future
  idx_not_na <- as.logical(complete.cases(y) * complete.cases(p))
  y <- y[idx_not_na]
  p <- p[idx_not_na, ]
  w <- w[idx_not_na]

  ## observed frequencies
  observed <- as.vector(xtabs(w ~ cut(y, breaks, include.lowest = TRUE)))

  ## expected frequencies (part2)
  expected <- colSums(p * w)

  # -------------------------------------------------------------------
  # OUTPUT AND OPTIONAL PLOTTING
  # -------------------------------------------------------------------
  ## collect everything as data.frame
  rval <- data.frame(
    observed = as.vector(observed),
    expected = as.vector(expected),
    mids = mids,
    width = diff(breaks) * width
  )

  ## attributes for graphical display
  attr(rval, "style") <- style
  attr(rval, "scale") <- scale
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main

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
  style <- unlist(lapply(rval, function(r) attr(r, "style")))
  scale <- unlist(lapply(rval, function(r) attr(r, "scale")))
  xlab <- unlist(lapply(rval, function(r) attr(r, "xlab")))
  ylab <- unlist(lapply(rval, function(r) attr(r, "ylab")))
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
  attr(rval, "style") <- style
  attr(rval, "scale") <- scale
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

#' S3 Methods for Plotting Rootograms
#' 
#' Generic plotting functions for rootograms of the class \code{"rootogram"}
#' computed by \code{link{rootogram}}. 
#'
#' Rootograms graphically compare (square roots) of empirical frequencies with
#' fitted frequencies from a probability model.
#' 
#' Rootograms graphically compare frequencies of empirical distributions and
#' fitted probability models. For the observed distribution the histogram is
#' drawn on a square root scale (hence the name) and superimposed with a line
#' for the fitted frequencies. The histogram can be \code{"standing"} on the
#' x-axis (as usual), or \code{"hanging"} from the fitted curve, or a
#' \code{"suspended"} histogram of deviations can be drawn.
#' 
#' @aliases plot.rootogram autoplot.rootogram
#' @param x,object an object of class \code{\link{rootogram}}.
#' @param ref logical. Should a reference line be plotted?
#' @param xlab,ylab,main graphical parameters.
#' @param xlim,ylim,fill,col,lwd,lty,axes,box graphical
#' parameters. These may pertain either to the whole plot or just the histogram
#' or just the fitted line.
#' @param colour,size,fitted,alpha,linetype graphical parameters passed to
#' \code{geom_line} and \code{geom_point}, respectively.
#' @param legend logical. Should a legend be added in the \code{ggplot2} style
#' graphic?
#' @param style,scale,alpha_min Fix description.
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
#' rootogram(m1_lm, ref = "blue", lty = 2, pch = 20)
#' 
#' ## add separate model
#' if (require("crch", quietly = TRUE)) {
#'   m1_crch <- crch(dist ~ speed | speed, data = cars)
#'   points(rootogram(m1_crch, plot = FALSE), col = 2, lty = 2)
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
#'   ## compute rootograms
#'   r2_lm <- rootogram(m2_lm, plot = FALSE)
#'   r2_crch <- rootogram(m2_crch, plot = FALSE)
#' 
#'   ## plot in single graph
#'   plot(c(r2_lm, r2_crch), col = c(1, 2))
#' }
#' 
#' #-------------------------------------------------------------------------------
#' ## determinants for male satellites to nesting horseshoe crabs
#' data("CrabSatellites", package = "countreg")
#' 
#' ## linear poisson model
#' m3_pois  <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
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
#' 
#' @export
plot.rootogram <- function(x,
                           style = NULL,
                           scale = NULL,
                           fitted = TRUE,
                           ref = TRUE,
                           xlim = c(NA, NA),
                           ylim = c(NA, NA),
                           xlab = NULL,
                           ylab = NULL,
                           main = NULL,
                           col = "black",
                           fill = "darkgray",
                           lwd = 2,
                           lty = 1,
                           alpha_min = 0.8,
                           axes = TRUE,
                           box = FALSE,
                           ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## get default arguments
  style <- use_arg_from_attributes(x, "style", default = "hanging", force_single = TRUE)
  scale <- use_arg_from_attributes(x, "scale", default = "sqrt", force_single = TRUE)
  fitted <- use_arg_from_attributes(x, "fitted", default = TRUE, force_single = FALSE)
  ref <- use_arg_from_attributes(x, "ref", default = NULL, force_single = FALSE)

  ## sanity checks
  ## * lengths of all arguments are checked by recycling
  ## * `xlab`, `ylab`, `main` and `...` w/i `plot()`
  ## * `border` and `fill` w/i `rect()`
  ## * `col`, `lwd`, `pch`, `lty` and `type` w/i `lines()`
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))

 
  ## compute heights 
  x <- summary(x, scale = scale, style = style)

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  # -------------------------------------------------------------------
  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING
  # -------------------------------------------------------------------
  ## get line style
  if (isFALSE(fitted)) {
    fitted <- "none"
  } else if (isTRUE(fitted)) {
    fitted <- "b"
  } else {
    fitted <- match.arg(fitted, c("both", "line", "point"))
    ifelse(fitted == "both", "b", ifelse(fitted == "line", "l", "p"))
  }

  ## recycle arguments for plotting to match the number of groups
  if (is.list(xlim)) xlim <- as.data.frame(do.call("rbind", xlim))
  if (is.list(ylim)) ylim <- as.data.frame(do.call("rbind", ylim))
  plot_arg <- data.frame(1:n, ref,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]],
    fill, col, lwd, lty, alpha_min, fitted, axes, box 
  )[, -1]

  ## annotation
  xlab <- use_arg_from_attributes(x, "xlab", default = "Rootogram", force_single = FALSE)
  ylab <- use_arg_from_attributes(x, "ylab", 
    default = if (scale == "raw") "Frequency" else "sqrt(Frequency)", force_single = FALSE)
  main <- use_arg_from_attributes(x, "main", default = "model", force_single = FALSE)

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
    ybottom = d$ymin
    ytop = d$ymax
    xleft <- d$mids - d$width / 2
    xright <- d$mids + d$width / 2

    if (scale == "sqrt") {
      d$expected = sqrt(d$expected)
    }

    ## get xlim and ylim
    if (any(is.na(c(plot_arg$xlim1[j], plot_arg$xlim2[j])))) {
      plot_arg[j, c("xlim1", "xlim2")] <- range(c(xleft, xright))
    }
    if (any(is.na(c(plot_arg$ylim1[j], plot_arg$ylim2[j])))) {
      plot_arg[j, c("ylim1", "ylim2")] <- range(c(ybottom, ytop, d$line))
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
      col = set_minimum_transparency(plot_arg$fill[j], alpha_min = plot_arg$alpha_min[j]),
      border = plot_arg$col[j], 
      lty = plot_arg$lty[j],
      lwd = plot_arg$lwd[j]
    )

    ## add ref line
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
      abline(h = 0, col = plot_arg$ref[j], lty = 1, lwd = 1.25)
    }

    ## add fitted line
    if (plot_arg$fitted[j] != "none") {
      lines(d$mids, d$expected,
        col = 2, pch = 19, type = plot_arg$fitted[j],
        lty = 1, lwd = 2
      )
    }
  }

  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  ## set up necessary panels
  if (n > 1L){
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
                               ref = TRUE,
                               fitted = TRUE,
                               xlim = c(NA, NA),
                               ylim = c(NA, NA),
                               xlab = NULL,
                               ylab = NULL,
                               main = NULL,
                               colour = "black",
                               fill = "darkgray",
                               size = 0.5,
                               linetype = 1,
                               alpha = NA,
                               legend = FALSE,
                               ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## get default arguments
  style <- use_arg_from_attributes(object, "style", default = "hanging", force_single = TRUE)
  scale <- use_arg_from_attributes(object, "scale", default = "sqrt", force_single = TRUE)
  fitted <- use_arg_from_attributes(object, "fitted", default = TRUE, force_single = TRUE)
  ref <- use_arg_from_attributes(object, "ref", default = NULL, force_single = TRUE)
  xlab <- use_arg_from_attributes(object, "xlab", default = "Rootogram", force_single = TRUE)
  ylab <- use_arg_from_attributes(object, "ylab", 
    default = if (scale == "raw") "Frequency" else "sqrt(Frequency)", force_single = TRUE)

  ## get base style arguments
  add_arg <- list(...)
  #if (!is.null(add_arg$pch)) shape <- add_arg$pch
  if (!is.null(add_arg$lwd)) size <- add_arg$lwd
  if (!is.null(add_arg$lty)) linetype <- add_arg$lty

  ## sanity checks
  stopifnot(all(sapply(xlim, function(x) is.numeric(x) || is.na(x))))
  stopifnot(all(sapply(ylim, function(x) is.numeric(x) || is.na(x))))

  ## get line style
  if (isFALSE(fitted)) {
    fitted <- "none"
  } else if (isTRUE(fitted)) {
    fitted <- "both"
  }
  fitted <- match.arg(fitted, c("both", "line", "point"))

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

  ## get main and into the right length (must be done after handling of `title`)
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
      mids = "mids", 
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
        alpha = "group"), 
      scale = scale, 
      style = style
    )

  ## add fitted line
  rval <- rval + 
    geom_rootogram_fitted(scale = scale, style = style) 

  ## add ref
  rval <- rval + 
    geom_rootogram_ref(ggplot2::aes_string(yintercept = 0))

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


## FIXME: (ML) Should this be implemented?
# "+.rootogram" <- function(e1, e2) {
#   style <- unique(c(attr(e1, "style"), attr(e2, "style")))
#   if(length(style) > 1L) {
#     warning(sprintf("different styles (%s != %s) had been used, result now uses style = %s",
#       style[1L], style[2L], style[1L]))
#     style <- style[1L]
#   }
#   scale <- unique(c(attr(e1, "scale"), attr(e2, "scale")))
#   if(length(scale) > 1L) {
#     warning(sprintf("different scales (%s != %s) had been used, result now uses scale = %s",
#       scale[1L], scale[2L], scale[1L]))
#     scale <- scale[1L]
#   }
#
#   ylab <- attr(e1, "ylab")
#   xlab <- paste(unique(c(attr(e1, "xlab"), attr(e2, "xlab"))), collapse = " / ")
#   main <- paste(unique(c(attr(e1, "main"), attr(e2, "main"))), collapse = " / ")
#   e1 <- as.data.frame(e1)
#   e2 <- as.data.frame(e2)
#   e <- e1[e1$x %in% e2$x, ] + e2[e2$x %in% e1$x, ]
#   rootogram.default(structure(e$observed, .Names = e$x/2), e$expected,
#     style = style, scale = scale,
#     main = main, xlab = xlab, ylab = ylab, plot = FALSE)
# }


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

  required_aes = c("observed", "expected", "mids", "width"),

  compute_group = function(data,  
                           scales, 
                           scale = c("sqrt", "raw"), 
                           style = c("hanging", "standing", "suspended")) {

    scale <- match.arg(scale)
    style <- match.arg(style)

    tmp_heights <- compute_rootogram_heights(data$expected, data$observed, scale = scale, style = style)

    data <- transform(data,
      xmin = mids - width/2, 
      xmax = mids + width/2,
      ymin = tmp_heights$ymin,
      ymax = tmp_heights$ymax,
      observed = NULL,
      expected = NULL,
      mids = NULL,
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
#' @param style Fix description.
#' @param scale Fix description.
#' @param linestyle Fix description.
#' @examples
#' if (require("ggplot2")) {
#'   ## Fit model
#'   data("CrabSatellites", package = "countreg")
#'   m1_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#'   m2_pois <- glm(satellites ~ color, data = CrabSatellites, family = poisson)
#'   
#'   ## Compute rootogram
#'   p1 <- rootogram(m1_pois, plot = FALSE)
#'   p2 <- rootogram(m2_pois, plot = FALSE)
#'   
#'   d <- c(p1, p2) 
#'   
#'   ## Get label names
#'   xlab <- unique(attr(d, "xlab"))
#'   ylab <- unique(attr(d, "ylab"))
#'   main <- attr(d, "main")
#'   main <- make.names(main, unique = TRUE)
#'   d$group <- factor(d$group, labels = main)
#'   
#'   gg1 <- ggplot(data = d) + 
#'     geom_rootogram(aes(observed = observed, expected = expected, mids = mids, 
#'       width = width, group = group)) + 
#'     geom_rootogram_fitted(aes(expected = expected, mids = mids)) + 
#'     geom_rootogram_ref(yintercept = 0) + 
#'     facet_grid(group~.)
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

  default_aes = ggplot2::aes(colour = "black", fill = "darkgray", size = 0.5, linetype = 1,
    alpha = NA),

  required_aes = c("xmin", "xmax", "ymin", "ymax")

)


# -------------------------------------------------------------------
# GGPLOT2 IMPLEMENTATIONS FOR `geom_rootogram_fitted()`
# -------------------------------------------------------------------

#' @rdname geom_rootogram
#' @export
stat_rootogram_fitted <- function(mapping = NULL,
                           data = NULL,
                           geom = "rootogram_fitted",
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
    stat = StatRootogramFitted,
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
StatRootogramFitted <- ggplot2::ggproto("StatRootogramFitted", ggplot2::Stat,

  required_aes = c("expected", "mids"),

  compute_group = function(data, 
                           scales, 
                           scale = c("sqrt", "raw"),
                           style = c("hanging", "standing", "suspended")) {

    scale <- match.arg(scale)
    style <- match.arg(style)

    ## raw vs. sqrt scale
    if (scale == "sqrt") {
      data <- transform(data,
        x = mids,
        y = sqrt(expected),
        expected = NULL,
        mids = NULL
      )
    } else {
      data <- transform(data,
        x = mids,
        y = expected,
        expected = NULL,
        mids = NULL
      )
    }
    data
  }
)

#' @rdname geom_rootogram
#' @export
geom_rootogram_fitted <- function(mapping = NULL, 
                                  data = NULL, 
                                  stat = "rootogram_fitted",
                                  position = "identity", 
                                  na.rm = FALSE,
                                  show.legend = NA, 
                                  inherit.aes = TRUE, 
                                  scale = c("sqrt", "raw"),
                                  style = c("hanging", "standing", "suspended"),
                                  linestyle = c("both", "line", "point"),
                                  ...) {

  scale <- match.arg(scale)
  style <- match.arg(style)
  linestyle <- match.arg(linestyle)

  ggplot2::layer(
    geom = GeomRootogramFitted, 
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
      linestyle = linestyle,
      ...
    )
  )
}


#' @rdname geom_rootogram
#' @export
GeomRootogramFitted <- ggplot2::ggproto("GeomRootogramFitted", ggplot2::GeomPath,

  default_aes = ggplot2::aes(colour = 2, size = 1, linetype = 1,
    alpha = 1, fill = NA, stroke = 0.5, shape = 19),

  draw_panel = function(data, panel_params, coord, arrow = NULL,
                        lineend = "butt", linejoin = "round", linemitre = 10,
                        na.rm = FALSE, 
                        linestyle = c("both", "line", "point")) {

    linestyle <- match.arg(linestyle)

    if (linestyle == "both") {
      ## TODO: (ML) Do not copy data
      data2 <- transform(data, size = size * 2)

      grid::grobTree(
        ggplot2::GeomPath$draw_panel(data, panel_params, coord, arrow = NULL,
                            lineend = "butt", linejoin = "round", linemitre = 10,
                            na.rm = FALSE),
        ggplot2::GeomPoint$draw_panel(data2, panel_params, coord, na.rm = FALSE)
      )
    } else if (linestyle == "line") {
        ggplot2::GeomPath$draw_panel(data, panel_params, coord, arrow = NULL,
                            lineend = "butt", linejoin = "round", linemitre = 10,
                            na.rm = FALSE)
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
  default_aes = ggplot2::aes(colour = "black", size = 0.5, linetype = 1,
  alpha = NA)
)


# -------------------------------------------------------------------
# HELPER FUNCTIONS FOR GETTING AN EXTENDED ROOTOGRAM OBJECT
# -------------------------------------------------------------------
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
                              ...) { 

  ## get arg `style` and `scale`
  scale <- use_arg_from_attributes(object, "scale", default = "sqrt", force_single = TRUE)
  style <- use_arg_from_attributes(object, "style", default = "hanging", force_single = TRUE)

  ## compute heights
  tmp <- compute_rootogram_heights(object$expected, object$observed, scale = scale, style = style)

  rval <- transform(object,
    ymin = tmp$ymin, 
    ymax = tmp$ymax
  )

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
  scale <- style <- NULL  # needed for `use_arg_from_attributes()`
  style <- use_arg_from_attributes(x, "style", default = NULL, force_single = TRUE)
  scale <- use_arg_from_attributes(x, "scale", default = NULL, force_single = TRUE)

  ## return custom print statement
  if (is.null(scale) || is.null(style)) {
    cat("A rootogram object without mandatory attributes `scale` and `style` \n\n")
  } else if ("ymin" %in% names(x) && "ymax" %in% names(x)) {
    cat(
      sprintf("An extended rootogram object with `scale = '%s'` and `style = '%s'` \n\n", 
      scale, 
      style
      )
    )
  } else {
    cat(
      sprintf("A rootogram object with `scale = '%s'` and `style = '%s'` \n\n", 
      scale, 
      style
      )
    )
  }

  ## call next print method
  NextMethod()
}
