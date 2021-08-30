# -------------------------------------------------------------------
# Programming outline: Reliagram (reliability diagram)
# -------------------------------------------------------------------
# - Observed y in-sample or out-of-sample (n x 1)
# - Thresholds in y (k x 1)
# - Predicted probabilities F_y(thresh) at thresholds (n x k)
# - Breaks for predicted probabilities in [0, 1] (m x 1)
#
# - Cut probabilities at breaks -> (m-1) groups
# - Cut y at thresholds, aggregate by groups -> (m-1) x k proportions
#
# Functions:
# - reliagram() generic plus default method
# - Return object of class "reliagram" that is plotted by default
# - But has plot=FALSE so that suitable methods can be added afterwards
# - At least methods: plot(), autoplot(), lines()
# -------------------------------------------------------------------


#' Reliagram (Extended Reliability Diagram)
#' 
#' Reliagram (extended reliability diagram) assess the reliability of a fitted
#' probabilistic distributional forecast for a binary event. If \code{plot =
#' TRUE}, the resulting object of class \code{"reliagram"} is plotted by
#' \code{\link{plot.reliagram}} or \code{\link{autoplot.reliagram}} before it is
#' returned, depending on whether the package \code{ggplot2} is loaded.
#' 
#' Reliagrams evaluate if a probability model is calibrated (reliable) by first
#' partitioning the predicted probability for a binary event into a certain number
#' of bins and then plotting (within each bin) the averaged forecast probability
#' against the observered/empirical relative frequency.  For computation,
#' \code{\link{reliagram}} leverages the \code{\link{procast}} generic to 
#' forecast the respective predictive probabilities.
#' 
#' For continous probability forecasts, reliability diagrams can be computed either
#' for a pre-specified threshold or for a specific quantile probability of the
#' response values. Per default, reliagrams are computed for the 50\%-quantile
#' of the reponse. 
#' 
#' In addition to the \code{plot} and \code{\link[ggplot2]{autoplot}} method for
#' reliagram objects, it is also possible to combine two (or more) reliability
#' diagrams by \code{c}/\code{rbind}, which creates a set of reliability diagrams
#' that can then be plotted in one go. 
#' 
#' @aliases reliagram reliagram.default c.reliagram 
#' @param object an object from which an extended reliability diagram can be
#' extracted with \code{\link{procast}}.
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
#' @param breaks numeric vector passed on to \code{\link[base]{cut}} in order to bin
#' the observations and the predicted probabilities or a function applied to
#' the predicted probabilities to calculate a numeric value for
#' \code{\link[base]{cut}}. Typically quantiles to ensure equal number of predictions
#' per bin, e.g., by \code{breaks = function(x) quantile(x)}.
#' @param quantiles numeric vector of quantile probabilities with values in
#' [0,1] to calculate single or several thresholds. Only used if
#' \code{thresholds} is not specified. For binary responses typically the
#' 50\%-quantile is used.
#' @param thresholds numeric vector specifying both where to cut the
#' observations into binary values and at which values the predicted
#' probabilities should be calculated (\code{\link{procast}}).
#' @param confint logical. Should confident intervals be calculated and drawn?
#' @param confint_level numeric. The confidence level required.
#' @param confint_nboot numeric. The number of bootstrap steps.
#' @param confint_seed numeric. The seed to be set for the bootstrapping.
#' @param single_graph logical. Should all computed extended reliability
#' diagrams be plotted in a single graph?
#' @param xlab,ylab,main graphical parameters.
#' @param \dots further graphical parameters.
#' @return An object of class \code{"reliagram"} inheriting from
#' \code{"data.frame"} or \code{"tibble"} conditional on the argument \code{class}
#' with the following variables: \item{x}{forecast probabilities,}
#' \item{y}{observered/empirical relative frequencies,} \item{bin_lwr, bin_upr}{
#' lower and upper bound of the binned forecast probabilities,}
#' \item{n_pred}{number of predictions within the binned forecasts probabilites,}
#' \item{ci_lwr, ci_upr}{lower and upper confidence interval bound.} Additionally,
#' \code{xlab}, \code{ylab}, \code{main}, and \code{treshold},
#' \code{confint_level}, as well as the total and the decomposed Brier Score 
#' (\code{bs, rel, res, unc}) are stored as attributes.
#' @note Note that there is also a \code{\link[verification]{reliability.plot}} function in the
#' \pkg{verification} package. However, it only works for numeric
#' forecast probabilities and numeric observed relative frequencies, hence a function has been
#' created here.
#' @seealso \code{link{plot.reliagram}}, \code{\link{procast}}
#' @references Wilks DS (2011) \emph{Statistical Methods in the Atmospheric
#' Sciences}, 3rd ed., Academic Press, 704 pp.
#' @examples
#'
#' ## speed and stopping distances of cars
#' m1_lm <- lm(dist ~ speed, data = cars)
#' 
#' ## compute and plot reliagram
#' reliagram(m1_lm)
#' 
#' #-------------------------------------------------------------------------------
#' ## determinants for male satellites to nesting horseshoe crabs
#' data("CrabSatellites", package = "countreg")
#' 
#' ## linear poisson model
#' m1_pois <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#' m2_pois <- glm(satellites ~ color, data = CrabSatellites, family = poisson)
#' 
#' ## compute and plot reliagram as base graphic
#' r1 <- reliagram(m1_pois, plot = FALSE)
#' r2 <- reliagram(m2_pois, plot = FALSE)
#' 
#' ## plot combined reliagram as "ggplot2" graphic
#' ggplot2::autoplot(c(r1, r2), single_graph = TRUE, col = c(1, 2), fill = c(1, 2))
#'
#' @export
reliagram <- function(object, ...) {
  UseMethod("reliagram")
}


#' @rdname reliagram
#' @method reliagram default
#' @export
reliagram.default <- function(object,
                              newdata = NULL,
                              plot = TRUE,
                              class = NULL,
                              breaks = seq(0, 1, by = 0.1),
                              quantiles = 0.5,
                              thresholds = NULL,
                              confint = TRUE,
                              confint_level = 0.95,
                              confint_nboot = 250,
                              confint_seed = 1,
                              single_graph = FALSE,
                              xlab = "Forecast probability",
                              ylab = "Observed relative frequency",
                              main = NULL,
                              ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * `object` and `newdata` w/i `newrepsone()`
  ## * `breaks w/i `cut()`;
  ## * `confint` w/i `polygon()`
  ## * `...` w/i `plot()` and `autoplot()`
  stopifnot(is.numeric(quantiles), is.null(dim(quantiles)))
  stopifnot(is.null(thresholds) || (is.numeric(thresholds) && is.null(dim(thresholds))))
  stopifnot(
    is.numeric(confint_level),
    length(confint_level) == 1,
    confint_level >= 0,
    confint_level <= 1
  )
  stopifnot(
    is.numeric(confint_nboot),
    length(confint_nboot) == 1,
    confint_nboot >= 0
  )
  stopifnot(is.logical(single_graph))
  stopifnot(length(xlab) == 1 || length(xlab) == length(quantiles))
  stopifnot(length(ylab) == 1 || length(ylab) == length(quantiles))
  stopifnot(is.null(main) || (length(main) == 1 || length(main) == length(quantiles)))
  stopifnot(is.numeric(breaks) || is.function(breaks))

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

  ## arguments xlab/ylab/main needed in case of `single_graph = TRUE`
  if (single_graph) {
    arg_xlab <- xlab[1]
    arg_ylab <- ylab[1]
    arg_main <- main[1]
  } else {
    arg_xlab <- NULL
    arg_ylab <- NULL
    arg_main <- NULL
  }

  # -------------------------------------------------------------------
  # PREPARE DATA
  # -------------------------------------------------------------------
  ## get data and threshold(s)
  y <- newresponse(object, newdata = newdata)
  if (is.null(thresholds)) {
    thresholds <- quantile(y, probs = quantiles, na.rm = TRUE)
    thresholds <- as.numeric(thresholds)
    thresholds_text <- sprintf("q_%.2f", signif(quantiles, 2))
  } else {
    thresholds_text <- as.character(signif(thresholds, 2))
  }

  ## fix length of annotations
  if (length(xlab) < length(quantiles)) xlab <- rep(xlab, length.out = length(quantiles))
  if (length(ylab) < length(quantiles)) ylab <- rep(ylab, length.out = length(quantiles))
  if (is.null(main)) {
    main <- deparse(substitute(object))
    main <- sprintf("%s (threshold = %s)", main, thresholds_text)
  } else if (length(main) < length(quantiles)) {
    main <- rep(main, length.out = length(quantiles))
  }

  ## predicted probabilities
  pred <- procast(object,
    newdata = newdata, type = "probability", at = matrix(thresholds, nrow = 1L),
    drop = FALSE
  )

  ## make sure lengths match (can't really go wrong after no arg `y` exists anymore)
  stopifnot(NROW(pred) == length(y))

  ## get and prepare observations
  y <- sapply(thresholds, function(x) y <= x)

  ## define convenience variables
  N <- NROW(y)

  # -------------------------------------------------------------------
  # COMPUTATION OF RELIABILITY DIAGRAM W/ CONSISTENCY RESAMPLING
  # -------------------------------------------------------------------
  ## loop over all quantiles (several possible thresholds)
  rval <- vector(mode = "list", length = NCOL(y))
  for (idx in 1:NCOL(y)) {

    ## calculate breaks
    if (is.function(breaks)) {
      try(breaks <- as.numeric(breaks(pred[, idx])))
      stopifnot("`breaks` function must produce a numeric valid to be used w/i `cut()`" = is.numeric(breaks))
    }

    ## compute number of prediction and idx for minimum number of prediction per probability subset
    n_pred <- aggregate(
      pred[, idx],
      by = list(prob = cut(pred[, idx], breaks, include.lowest = TRUE)),
      FUN = length,
      drop = FALSE
    )[, "x"]

    ## compute observed relative frequencies of positive examples (obs_rf)
    obs_rf <- as.numeric(
      aggregate(
        y[, idx],
        by = list(prob = cut(pred[, idx], breaks, include.lowest = TRUE)),
        FUN = mean,
        drop = FALSE
      )[, "x"]
    )

    ## compute mean predicted probability (mean_pr)
    tmp_mean_pr <- aggregate(
      pred[, idx],
      by = list(prob = cut(pred[, idx], breaks, include.lowest = TRUE)),
      FUN = mean,
      drop = FALSE
    )
    mean_pr <- as.numeric(tmp_mean_pr[, "x"])

    ## calculate pred with mean values (for calulating bs)
    lookup <- as.numeric(cut(pred[, idx], breaks, include.lowest = TRUE))
    pred_bin <- tmp_mean_pr$x[lookup]

    ## consistency resampling from Broecker (2007)
    if (!identical(confint, FALSE)) {
      set.seed(confint_seed)
      obs_rf_boot <- vector("list", length = N)
      for (i in 1:confint_nboot) {

        ## take bootstrap sample from predictions (surrogate forecasts)
        pred_hat <- sample(pred[, idx], replace = TRUE)

        ## surrogate observations that are reliable by construction
        yhat <- runif(N) < pred_hat

        ## compute observed relative frequencies of the surrogate observations
        obs_rf_boot[[i]] <- as.numeric(
          aggregate(
            yhat,
            by = list(prob = cut(pred_hat, breaks, include.lowest = TRUE)),
            FUN = mean,
            drop = FALSE
          )[, "x"]
        )
      }
      obs_rf_boot <- do.call("rbind", obs_rf_boot)

      ## compute lower and upper limits for reliable forecasts
      confint_prob <- (1 - confint_level) / 2
      confint_prob <- c(confint_prob, 1 - confint_prob)
      ci_lwr <- apply(obs_rf_boot, 2, quantile, prob = confint_prob[1], na.rm = TRUE)
      ci_upr <- apply(obs_rf_boot, 2, quantile, prob = confint_prob[2], na.rm = TRUE)
    } else {
      ci_lwr <- NA
      ci_upr <- NA
      confint_level <- NA
    }

    ## collect everything as data.frame
    rval_i <- data.frame(
      x = mean_pr,
      y = obs_rf,
      bin_lwr = breaks[-length(breaks)],
      bin_upr = breaks[-1],
      n_pred,
      ci_lwr,
      ci_upr
    )

    ## attributes for graphical display
    attr(rval_i, "xlab") <- xlab[idx]
    attr(rval_i, "ylab") <- ylab[idx]
    attr(rval_i, "main") <- main[idx]
    attr(rval_i, "threshold") <- thresholds_text[idx]
    attr(rval_i, "confint_level") <- confint_level

    ## add bs, rel, res, and unc
    ## NOTE: (ML) Here the unique forecasts equal the mean forecasts per bin (as in `verification` pkg):
    ## * Hence, BS is not independent to the bins
    ## * Hence, BS should vary conditional on the minimum
    ## * na.rm = TRUE: Should this be changed?!
    attr(rval_i, "bs") <- mean((pred_bin - y[, idx])^2, na.rm = TRUE)
    attr(rval_i, "rel") <- sum(n_pred * (mean_pr - obs_rf)^2, na.rm = TRUE) / sum(n_pred, na.rm = TRUE)
    attr(rval_i, "res") <- sum(n_pred * (obs_rf - mean(y[, idx]))^2, na.rm = TRUE) / sum(n_pred, na.rm = TRUE)
    attr(rval_i, "unc") <- mean(y[, idx]) * (1 - mean(y[, idx]))

    ## add class
    if (class == "data.frame") {
      class(rval_i) <- c("reliagram", "data.frame")
    } else {
      rval_i <- tibble::as_tibble(rval_i)
      class(rval_i) <- c("reliagram", class(rval_i))
    }

    rval[[idx]] <- rval_i
  }

  # -------------------------------------------------------------------
  # OUTPUT AND OPTIONAL PLOTTING
  # -------------------------------------------------------------------
  ## combine different groups
  rval <- do.call(c, rval)

  ## plot by default
  if (plot == "ggplot2") {
    try(print(ggplot2::autoplot(rval,
      confint = confint, single_graph = single_graph, main = arg_main, ...
    )))
  } else if (plot == "base") {
    try(plot(rval,
      confint = confint, single_graph = single_graph, main = arg_main, xlab = arg_xlab, ylab = arg_ylab, ...
    ))
  }

  ## return invisibly
  invisible(rval)
}


#' @export
c.reliagram <- function(...) {
  # -------------------------------------------------------------------
  # GET DATA
  # -------------------------------------------------------------------
  ## list of reliagrams
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
  xlab <- unlist(lapply(rval, function(r) attr(r, "xlab")))
  ylab <- unlist(lapply(rval, function(r) attr(r, "ylab")))
  prob <- unlist(lapply(rval, function(r) attr(r, "prob")))
  confint_level <- unlist(lapply(rval, function(r) attr(r, "confint_level")))
  bs <- unlist(lapply(rval, function(r) attr(r, "bs")))
  rel <- unlist(lapply(rval, function(r) attr(r, "rel")))
  res <- unlist(lapply(rval, function(r) attr(r, "res")))
  unc <- unlist(lapply(rval, function(r) attr(r, "unc")))
  nam <- names(rval)
  main <- if (is.null(nam)) {
    as.vector(unlist(lapply(rval, function(r) attr(r, "main"))))
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
  attr(rval, "xlab") <- xlab
  attr(rval, "ylab") <- ylab
  attr(rval, "main") <- main
  attr(rval, "prob") <- prob
  attr(rval, "confint_level") <- confint_level
  attr(rval, "bs") <- bs
  attr(rval, "rel") <- rel
  attr(rval, "res") <- res
  attr(rval, "unc") <- unc

  ## set class to data.frame or tibble
  if (class == "data.frame") {
    class(rval) <- c("reliagram", "data.frame")
  } else {
    rval <- tibble::as_tibble(rval)
    class(rval) <- c("reliagram", class(rval))
  }

  ## return
  return(rval)
}


#' @export
rbind.reliagram <- c.reliagram


#' S3 Methods for a Reliagram (Extended Reliability Diagram)
#' 
#' Generic plotting functions for reliability diagrams of the class \code{"reliagram"}
#' computed by \code{link{reliagram}}. 
#' 
#' Reliagrams evaluate if a probability model is calibrated (reliable) by first
#' partitioning the forecast probability for a binary event into a certain number
#' of bins and then plotting (within each bin) the averaged forecast probability
#' against the observered/empirical relative frequency. 
#' 
#' For continous probability forecasts, reliability diagrams can be plotted either
#' for a pre-specified threshold or for a specific quantile probability of the
#' response values.
#' 
#' Reliagrams can be rendered as \code{ggplot2} or base R graphics by using
#' the generics \code{\link[ggplot2]{autoplot}} or \code{\link[graphics]{plot}}. 
#' For a single base R graphically panel, \code{\link{points}} adds an additional 
#' reliagram.
#' 
#' @aliases plot.reliagram lines.reliagram autoplot.reliagram
#' @param object,x an object of class \code{reliagram}.
#' @param single_graph logical. Should all computed extended reliability
#' diagrams be plotted in a single graph?
#' @param confint logical. Should confident intervals be calculated and drawn?
#' @param xlab,ylab,main graphical parameters.
#' @param \dots further graphical parameters.
#' @param minimum,ref,xlim,ylim,col,fill,alpha_min,lwd,pch,lty,type,add_hist,add_info,add_rug,add_min,axes,box additional graphical
#' parameters for base plots, whereby \code{x} is a object of class \code{reliagram}.
#' @param colour,size,shape,linetype,legend graphical parameters passed for 
#' \code{ggplot2} style plots, whereby \code{object} is a object of class \code{reliagram}.
#' @seealso \code{link{reliagram}}, \code{\link{procast}}
#' @references Wilks DS (2011) \emph{Statistical Methods in the Atmospheric
#' Sciences}, 3rd ed., Academic Press, 704 pp.
#' @examples
#' 
#' ## speed and stopping distances of cars
#' m1_lm <- lm(dist ~ speed, data = cars)
#' 
#' ## compute and plot reliagram
#' reliagram(m1_lm)
#' 
#' ## customize colors
#' reliagram(m1_lm, ref = "blue", lty = 2, pch = 20)
#' 
#' ## add separate model
#' if (require("crch", quietly = TRUE)) {
#'   m1_crch <- crch(dist ~ speed | speed, data = cars)
#'   lines(reliagram(m1_crch, plot = FALSE), col = 2, lty = 2, confint = 2)
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
#'   ## compute reliagrams
#'   rel2_lm <- reliagram(m2_lm, plot = FALSE)
#'   rel2_crch <- reliagram(m2_crch, plot = FALSE)
#' 
#'   ## plot in single graph
#'   plot(c(rel2_lm, rel2_crch), col = c(1, 2), confint = c(1, 2), ref = 3, single_graph = TRUE)
#' }
#' 
#' #-------------------------------------------------------------------------------
#' ## determinants for male satellites to nesting horseshoe crabs
#' data("CrabSatellites", package = "countreg")
#' 
#' ## linear poisson model
#' m3_pois  <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
#' 
#' ## compute and plot reliagram as "ggplot2" graphic
#' reliagram(m3_pois, plot = "ggplot2")
#'
#' @export
plot.reliagram <- function(x,
                           single_graph = FALSE,
                           minimum = 0,
                           confint = TRUE,
                           ref = TRUE,
                           xlim = c(0, 1),
                           ylim = c(0, 1),
                           xlab = NULL,
                           ylab = NULL,
                           main = NULL,
                           col = "black",
                           fill = adjustcolor("black", alpha.f = 0.2),
                           alpha_min = 0.2,
                           lwd = 2,
                           pch = 19,
                           lty = 1,
                           type = NULL,
                           add_hist = TRUE,
                           add_info = TRUE,
                           add_rug = TRUE,
                           add_min = TRUE,
                           axes = TRUE,
                           box = TRUE,
                           ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks:
  ## * lengths of all arguments are checked by recycling
  ## * `ref` w/i `abline()`; `xlim`, `ylim`, `xlab`, `ylab`, `main`, `col`, `fill`,
  ##     `lwd`, `pch`, `lty`, `type` and `...` w/i `plot()`
  ## * `confint` w/i `polygon()`
  ## * `alpha_min` w/i `set_minimum_transparency()`
  stopifnot(is.logical(single_graph))
  stopifnot(is.numeric(minimum), all(minimum >= 0))
  stopifnot(is.logical(add_info))
  stopifnot(is.logical(axes))
  stopifnot(is.logical(box))

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## set single_multigraph for single plot always to FALSE
  if (n == 1) {
    single_multigraph <- FALSE
  } else {
    single_multigraph <- single_graph
  }

  ## determine if points should be plotted
  if (is.null(type)) type <- ifelse(table(x$group) > 20L, "l", "b")

  ## recycle arguments for plotting to match the number of groups
  if (is.list(xlim)) xlim <- as.data.frame(do.call("rbind", xlim))
  if (is.list(ylim)) ylim <- as.data.frame(do.call("rbind", ylim))
  plot_arg <- data.frame(1:n, minimum, confint, ref,
    xlim1 = xlim[[1]], xlim2 = xlim[[2]], ylim1 = ylim[[1]], ylim2 = ylim[[2]],
    col, fill, alpha_min, lwd, pch, lty, type, add_hist, add_info, add_rug, add_min,
    axes, box
  )[, -1]

  ## annotation
  ## FIXME: (ML) main title must not be unique; hence also works w/ empty main. Different in `autoplot()`
  ##   e.g., `main <- make.unique(rep_len(main, n))`
  if (single_multigraph) {
    xlab <- if (is.null(xlab)) "Forecast probability" else xlab
    ylab <- if (is.null(ylab)) "Observed relative frequency" else ylab
    main <- if (is.null(main)) "Reliability diagram" else main
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
  # FUNCTION TO TRIGGER FIGURE AND PLOT CONFINT
  # -------------------------------------------------------------------
  reliagram_trigger <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get bins with sufficient observations
    min_idx <- which(d$n_pred >= plot_arg$minimum[j])
    if (length(min_idx) == 0) {
      stop(sprintf("no bin has sufficent cases for defined minimum = %s\n", plot_arg$minimum[j]))
    }

    ## get minimum and maximum breaks to decide if extend confint polygon to the corners
    extend_left <- min(d$bin_lwr) == min(d[min_idx, "bin_lwr"])
    extend_right <- max(d$bin_upr) == max(d[min_idx, "bin_upr"])

    ## modify main using subscript for quantiles
    if (grepl("threshold = q_[0-9]+\\.?([0-9]+)?)$", main[j])) {
      tmp_quantile <- regmatches(main[j], regexpr("q_[0-9]+\\.?([0-9]+)?", main[j]))
      tmp_quantile <- regmatches(tmp_quantile, regexpr("[0-9]+\\.?([0-9]+)?", tmp_quantile))
      tmp_text <- sub("q_[0-9]+\\.?([0-9]+)?)", "", main[j])
      main[j] <- as.expression(bquote(bold(.(tmp_text) * q[.(tmp_quantile)] * ")")))
    }

    ## trigger plot
    if (j == 1 || (!single_multigraph && j > 1)) {
      plot(0, 0,
        type = "n", xlim = c(plot_arg$xlim1[j], plot_arg$xlim2[j]),
        ylim = c(plot_arg$ylim1[j], plot_arg$ylim2[j]), xlab = xlab[j],
        ylab = ylab[j], main = main[j],
        xaxs = "i", yaxs = "i", axes = FALSE, ...
      )
      if (plot_arg$axes[j]) {
        axis(1)
        axis(2)
      }
      if (plot_arg$box[j]) {
        box()
      }
    }

    ## plot reference line
    if (j == 1 || (!single_multigraph && j > 1)) {
      if (!identical(plot_arg$ref[j], FALSE)) {
        if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
        abline(0, 1, col = plot_arg$ref[j], lty = 2, lwd = 1.25)
      }
    }

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]
      polygon(
        na.omit(c(
          ifelse(extend_left, 0, NA),
          d[min_idx, "x"],
          ifelse(extend_right, 1, NA),
          rev(d[min_idx, "x"]),
          ifelse(extend_left, 0, NA)
        )),
        na.omit(c(
          ifelse(extend_left, 0, NA),
          d[min_idx, "ci_lwr"],
          ifelse(extend_right, 1, NA),
          rev(d[min_idx, "ci_upr"]),
          ifelse(extend_left, 0, NA)
        )),
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min = plot_arg$alpha_min[j]),
        border = NA
      )
    }
  }

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION
  # -------------------------------------------------------------------
  reliagram_plot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get lines with sufficient observations
    min_idx <- which(d$n_pred >= plot_arg$minimum[j])

    ## plot reliability line
    lines(y ~ x, d[min_idx, ],
      type = plot_arg$type[j], lwd = plot_arg$lwd[j],
      pch = plot_arg$pch[j], lty = plot_arg$lty[j], col = plot_arg$col[j], ...
    )

    ## plot points below minimum
    if (!identical(plot_arg$add_min[j], FALSE)) {
      if (isTRUE(plot_arg$add_min[j])) plot_arg$add_min[j] <- 4
      points(y ~ x, d[-min_idx, ], pch = plot_arg$add_min[j], col = plot_arg$col[j], ...)
    }

    ## add rugs
    if (!identical(plot_arg$add_rug[j], FALSE)) {
      if (isTRUE(plot_arg$add_rug[j])) plot_arg$add_rug[j] <- plot_arg$col[j]
      tmp_rugs <- c(d$bin_lwr, d$bin_upr[NROW(d)])
      tmp_rugs <- tmp_rugs[tmp_rugs >= plot_arg$xlim1[j] & tmp_rugs <= plot_arg$xlim2[j]]
      rug(tmp_rugs, lwd = 1, ticksize = 0.02, col = plot_arg$add_rug[j])
    }

    ## add hist
    if (!single_multigraph && !identical(plot_arg$add_hist[j], FALSE)) {
      if (isTRUE(plot_arg$add_hist[j])) plot_arg$add_hist[j] <- "lightgray"
      tmp_x <- par("pin")[1]
      tmp_y <- par("pin")[2]
      tmp_height <- 0.3 * diff(c(plot_arg$ylim1[j], plot_arg$ylim2[j]))
      tmp_width <- (0.3 * diff(c(plot_arg$xlim1[j], plot_arg$xlim2[j]))) * tmp_y / tmp_x

      add_hist_reliagram(
        d$n_pred,
        c(d$bin_lwr, d$bin_upr[NROW(d)]),
        plot_arg$minimum[j],
        xpos = 0.05 * diff(c(plot_arg$xlim1[j], plot_arg$xlim2[j])) + plot_arg$xlim1[j],
        ypos = 0.925 * diff(c(plot_arg$ylim1[j], plot_arg$ylim2[j])) - tmp_height + plot_arg$ylim1[j],
        width = tmp_width,
        height = tmp_height,
        col = plot_arg$add_hist[j]
      )
    }

    ## print info
    if (!single_multigraph && plot_arg$add_info[j]) {
      legend(
        "bottomright",
        c(
          "BS",  sprintf("%.3f", signif(attr(d, "bs")[j], 3)),
          "REL", sprintf("%.3f", signif(attr(d, "rel")[j], 3)),
          "RES", sprintf("%.3f", signif(attr(d, "res")[j], 3)),
          "UNC", sprintf("%.3f", signif(attr(d, "unc")[j], 3))
        ),
        cex = 0.8,
        ncol = 4,
        bty = "n",
        inset = c(0.01, 0.01)
      )
    }
  }

  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  ## set up necessary panels
  if (!single_multigraph && n > 1L){
    old_pars <- par(mfrow = n2mfrow(n))
    on.exit(par(old_pars), add = TRUE)
  }

  ## draw polygons first
  if (single_multigraph) {
    for (i in 1L:n) reliagram_trigger(x[x$group == i, ], ...)
    for (i in 1L:n) reliagram_plot(x[x$group == i, ], ...)
  } else {
    for (i in 1L:n) {
      reliagram_trigger(x[x$group == i, ], ...)
      reliagram_plot(x[x$group == i, ], ...)
    }
  }
}


#' @rdname plot.reliagram
#' @method lines reliagram
#' @export
lines.reliagram <- function(x,
                            minimum = 0,
                            confint = FALSE,
                            ref = FALSE,
                            col = "black",
                            fill = adjustcolor("black", alpha.f = 0.2),
                            alpha_min = 0.2,
                            lwd = 2,
                            pch = 19,
                            lty = 1,
                            type = "b",
                            ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## sanity checks
  ## * lengths of all arguments are checked by recycling
  ## * `col`, `lwd`, `pch`, `lty`, `type` and `...` w/i `plot()`
  ## * `ref` w/i `abline()`; `confint` w/i `polygon()`
  ## * `alpha_min` w/i `set_minimum_transparency()`
  stopifnot(is.numeric(minimum), all(minimum >= 0))

  ## convert always to data.frame
  x <- as.data.frame(x)

  ## handling of groups
  if (is.null(x$group)) x$group <- 1L
  n <- max(x$group)

  ## recycle arguments for plotting to match the number of groups
  plot_arg <- data.frame(
    1:n, minimum, confint, ref, col, fill, alpha_min, lwd, pch, lty, type
  )[, -1]

  # -------------------------------------------------------------------
  # MAIN PLOTTING FUNCTION FOR LINES
  # -------------------------------------------------------------------
  reliagramplot <- function(d, ...) {

    ## get group index
    j <- unique(d$group)

    ## get lines with sufficient observations
    min_idx <- which(d$n_pred >= plot_arg$minimum[j])

    ## get minimum and maximum breaks to decide if extend confint polygon to the corners
    extend_left <- min(d$bin_lwr) == min(d[min_idx, "bin_lwr"])
    extend_right <- max(d$bin_upr) == max(d[min_idx, "bin_upr"])

    ## plot confint polygon
    if (!identical(plot_arg$confint[j], FALSE) && !is.na(attr(d, "confint_level")[j])) {
      if (isTRUE(plot_arg$confint[j])) plot_arg$confint[j] <- plot_arg$fill[j]
      polygon(
        na.omit(c(
          ifelse(extend_left, 0, NA),
          d[min_idx, "x"],
          ifelse(extend_right, 1, NA),
          rev(d[min_idx, "x"]),
          ifelse(extend_left, 0, NA)
        )),
        na.omit(c(
          ifelse(extend_left, 0, NA),
          d[min_idx, "ci_lwr"],
          ifelse(extend_right, 1, NA),
          rev(d[min_idx, "ci_upr"]),
          ifelse(extend_left, 0, NA)
        )),
        col = set_minimum_transparency(plot_arg$confint[j], alpha_min = plot_arg$alpha_min[j]),
        border = NA
      )
    }

    ## plot reference line
    if (!identical(plot_arg$ref[j], FALSE)) {
      if (isTRUE(plot_arg$ref[j])) plot_arg$ref[j] <- "black"
      abline(0, 1, col = plot_arg$ref[j], lty = 2, lwd = 1.25)
    }

    ## plot reliability line
    lines(y ~ x, d[min_idx, ],
      type = plot_arg$type[j], lwd = plot_arg$lwd[j],
      pch = plot_arg$pch[j], lty = plot_arg$lty[j], col = plot_arg$col[j], ...
    )
  }

  # -------------------------------------------------------------------
  # DRAW PLOTS
  # -------------------------------------------------------------------
  for (i in 1L:n) {
    reliagramplot(x[x$group == i, ], ...)
  }
}


#' @rdname plot.reliagram
#' @method autoplot reliagram
#' @exportS3Method ggplot2::autoplot
autoplot.reliagram <- function(object,
                               single_graph = FALSE,
                               minimum = 0,
                               confint = TRUE,
                               ref = TRUE,
                               xlim = c(0, 1),
                               ylim = c(0, 1),
                               xlab = NULL,
                               ylab = NULL,
                               main = NULL,
                               colour = "black",
                               fill = adjustcolor("black", alpha.f = 0.2),
                               alpha_min = 0.2,
                               size = 1,
                               shape = 19,
                               linetype = 1,
                               type = NULL,
                               add_hist = TRUE,
                               add_info = TRUE,
                               add_rug = TRUE,
                               add_min = TRUE,
                               legend = FALSE,
                               ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## get base style arguments
  add_arg <- list(...)
  if (!is.null(add_arg$pch)) shape <- add_arg$pch
  if (!is.null(add_arg$lwd)) size <- add_arg$lwd
  if (!is.null(add_arg$lty)) linetype <- add_arg$lty

  ## sanity checks
  stopifnot(is.logical(single_graph))
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
  main <- make.unique(rep_len(main, n))

  ## prepare grouping
  object$group <- factor(object$group, levels = 1L:n, labels = main)

  # -------------------------------------------------------------------
  # PREPARE AND DEFINE ARGUMENTS FOR PLOTTING
  # -------------------------------------------------------------------
  ## determine if points should be plotted
  if (is.null(type)) type <- ifelse(table(object$group) > 20L, "l", "b")

  ## set alpha to 0 or color to NA for not plotting
  type <- ifelse(type == "l", 0, 1)
  if (is.logical(ref)) ref <- ifelse(ref, 1, NA)
  if (is.logical(add_min)) add_min <- ifelse(add_min, 4, NA)
  if (is.logical(add_rug)) add_rug <- ifelse(add_rug, colour, NA)

  ## get min and max breaks to decide if extend confint polygon to the corners
  min_break <- min(object$bin_lwr)
  max_break <- max(object$bin_upr)

  ## stat helper function to get polygon 
  calc_confint_polygon <- ggplot2::ggproto("calc_confint_polygon", ggplot2::Stat,

    # required as we operate on groups (facetting)
    compute_group = function(data, scales) {
      ## manipulate object
      if (min(data$bin_lwr) == min_break & max(data$bin_upr) == max_break) {
        nd <- data.frame(
          x = c(0, data$x, 1, rev(data$x), 0),
          y = c(0, data$ci_lwr, 1, rev(data$ci_upr), 0)
        )
      } else if (min(data$bin_lwr) == min_break) {
        nd <- data.frame(
          x = c(0, data$x, rev(data$x), 0),
          y = c(0, data$ci_lwr, rev(data$ci_upr), 0)
        )
      } else if (max(data$bin_upr) == max_break) {
        nd <- data.frame(
          x = c(data$x, 1, rev(data$x)),
          y = c(data$ci_lwr, 1, rev(data$ci_upr))
        )
      } else {
        nd <- data.frame(
          x = c(data$x, rev(data$x)),
          y = c(data$ci_lwr, rev(data$ci_upr))
        )
      }
      nd
    },

    # tells us what we need
    required_aes = c("x", "ci_lwr", "ci_upr", "bin_lwr", "bin_upr")
  )

  ## recycle arguments for plotting to match the number of groups (for `scale_<...>_manual()`)
  plot_arg <- data.frame(
    1:n,
    colour, fill, size, linetype, confint, alpha_min, minimum, add_rug
  )[, -1]

  ## prepare fill color for confint (must be done on vector to match args)
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

  ## recycle arguments for plotting to match the length (rows) of the object (for geom w/ aes)
  plot_arg2 <- data.frame(1:n, ref, type, size, shape, minimum, add_min, add_rug)[, -1]
  plot_arg2 <- as.data.frame(lapply(plot_arg2, rep, table(object$group)))

  ## set points to NA with no sufficent number of predictions
  idx_min <- which(object$n_pred < plot_arg2$minimum)
  object2 <- object[idx_min, ]
  object[idx_min, c("x", "y")] <- NA

  # -------------------------------------------------------------------
  # MAIN PLOTTING
  # -------------------------------------------------------------------
  ## actual plotting
  rval <- ggplot2::ggplot(object, ggplot2::aes_string(x = "x", y = "y")) +
    ggplot2::geom_abline(ggplot2::aes_string(intercept = 0, slope = 1),
      linetype = 2, colour = plot_arg2$ref
    ) +
    ggplot2::geom_polygon(
      ggplot2::aes_string(
        ci_lwr = "ci_lwr", ci_upr = "ci_upr",
        bin_lwr = "bin_lwr", bin_upr = "bin_upr", fill = "group"
      ),
      stat = calc_confint_polygon, show.legend = FALSE, na.rm = TRUE
    ) +
    ggplot2::geom_line(ggplot2::aes_string(colour = "group", size = "group", linetype = "group"),
      na.rm = TRUE
    ) +
    ggplot2::geom_point(ggplot2::aes_string(colour = "group"),
      alpha = plot_arg2$type, size = plot_arg2$size * 2, shape = plot_arg2$shape,
      show.legend = FALSE, na.rm = TRUE
    )

  ## add points below minimum
  rval <- rval +
    ggplot2::geom_point(ggplot2::aes_string(x = "x", y = "y"),
      data = object2,
      alpha = plot_arg2$type[idx_min], size = plot_arg2$size[idx_min] * 2, shape = plot_arg2$add_min[idx_min],
      show.legend = FALSE, na.rm = TRUE
    )

  ## add rugs
  rval <- rval +
    ggplot2::geom_rug(ggplot2::aes_string(x = "x"),
      data = data.frame(x = object$bin_lwr[1], group = factor(1L:n, labels = main)),
      inherit.aes = FALSE, colour = plot_arg$add_rug
    ) +
    ggplot2::geom_rug(ggplot2::aes_string(x = "bin_upr"), y = NA, colour = plot_arg2$add_rug)

  ## set the colors, shapes, etc. for the groups
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

  # -------------------------------------------------------------------
  # ADD HISTOGRAM IF WANTED
  # -------------------------------------------------------------------
  if (!identical(add_hist, FALSE) & (n == 1 | !single_graph && n > 1L)) {
    if (isTRUE(add_hist)) add_hist <- "lightgray"

    get_inset <- function(df) {
      df$n_pred[is.na(df$n_pred)] <- 0

      ## draw histogram
      rval_inset <- ggplot2::ggplot(
        df,
        ggplot2::aes_string(
          x = "x", y = "y", xmin = "bin_lwr", xmax = "bin_upr", ymin = 0, ymax = "n_pred",
          fill = "above_min"
        )
      ) +
        ggplot2::geom_rect(show.legend = FALSE, colour = "black", size = 0.25) +
        ggplot2::scale_fill_manual(values = c(add_hist, "white")) +
        ggplot2::theme_void()

      ## add minimum line
      if (any(df$minimum > 0)) {
        rval_inset <- rval_inset +
          ggplot2::geom_segment(ggplot2::aes_string(x = "x", xend = "xend", y = "y", yend = "yend"),
            data = data.frame(x = 0, xend = 1, y = unique(df$minimum), yend = unique(df$minimum)),
            inherit.aes = FALSE
          ) +
          ggplot2::geom_text(ggplot2::aes_string(x = "x", y = "y", label = "label"),
            data = data.frame(x = 0, y = unique(df$minimum), label = "Min."), inherit.aes = FALSE,
            size = 3, hjust = 1, nudge_x = -0.01
          )
      }

      ## add simple y axis
      ytick <- pretty(c(0, max(df$n_pred)), 4)
      ytick <- ytick[ytick > 0 & ytick < max(df$n_pred)]

      rval_inset <- rval_inset +
        ggplot2::geom_segment(ggplot2::aes_string(x = "x", xend = "xend", y = "y", yend = "yend"),
          data = data.frame(x = 0.985, xend = 1.015, y = ytick, yend = ytick),
          inherit.aes = FALSE
        ) +
        ggplot2::geom_text(ggplot2::aes_string(x = "x", y = "y", label = "label"),
          data = data.frame(x = 1.015, y = ytick, label = ytick), inherit.aes = FALSE,
          size = 3, hjust = 0, nudge_x = 0.01
        )

      # add nobs
      rval_inset <- rval_inset +
        ggplot2::geom_text(ggplot2::aes_string(x = "x", y = "y", label = "label"),
          data = data.frame(
            x = mean(c(df$bin_lwr, df$bin_upr), na.rm = TRUE),
            y = max(df$n_pred) + max(df$n_pred) * 0.15,
            label = paste0("n = ", sum(df$n_pred[df$n_pred >= df$minimum], na.rm = TRUE))
          ),
          inherit.aes = FALSE, size = 4
        )

      ## return graph
      rval_inset <- rval_inset +
        ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
        ggplot2::scale_y_continuous(expand = c(0.1, 0.1))
    }

    ## add needed variables to df
    object$above_min <- factor(object$n_pred >= plot_arg2$minimum, levels = c(TRUE, FALSE))
    object$minimum <- plot_arg2$minimum

    ## loop over different groups
    insets <- lapply(split(object, object$group), function(x) {
      annotation_custom2(
        grob = ggplot2::ggplotGrob(get_inset(x)),
        data = data.frame(x),
        ymin = 0.7, ymax = 0.95, xmin = 0.05, xmax = 0.3
      )
    })
    rval <- rval + insets
  }

  # -------------------------------------------------------------------
  # ADD INFORMATION IF WANTED
  # -------------------------------------------------------------------
  if (add_info & (n == 1 | !single_graph && n > 1L)) {
    rval <- rval +
      ggplot2::geom_text(ggplot2::aes_string(x = "x", y = "y", label = "label"),
        data = data.frame(
          x = 1, y = 0.06,
          label = sprintf(
            "BS    REL   RES   UNC\n%.3f  %.3f  %.3f  %.3f",
            signif(attr(object, "bs"), 3),
            signif(attr(object, "rel"), 3),
            signif(attr(object, "res"), 3),
            signif(attr(object, "unc"), 3)
          ),
          group = factor(1L:n, labels = main)
        ),
        inherit.aes = FALSE, size = 3, hjust = 1
      )
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


add_hist_reliagram <- function(n,
                               breaks,
                               minimum,
                               xpos,
                               ypos,
                               width = .2,
                               height = .2,
                               col = "lightgray",
                               main = NULL) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  idx_min <- which(n < minimum)
  n[is.na(n)] <- 0
  max_n <- max(n)
  col <- rep(col, max_n)
  col[n < minimum] <- 0

  # -------------------------------------------------------------------
  # PLOT 
  # -------------------------------------------------------------------
  ## plot histogram
  for (i in seq_along(n)) {
    x <- xpos + breaks[c(i, i + 1)] * width
    y <- ypos + c(0, n[i] / max_n * height)
    rect(x[1L], y[1L], x[2L], y[2L], col = col[i])
  }

  ## plot y-axis
  ytick <- pretty(c(0, max_n * .8), 4)
  ytick <- ytick[ytick > 0 & ytick < max_n]
  text(xpos + 1.05 * width, ypos + ytick / max_n * height, ytick, cex = .8, adj = c(0.0, 0.5))
  segments(
    x0 = rep(xpos, length(ytick)),
    x1 = rep(xpos + 1 * width, length(ytick)),
    y0 = ypos + ytick / max_n * height,
    col = "darkgray", lwd = .5, lty = 2
  )
  segments(
    x0 = rep(xpos + 0.975 * width, length(ytick)),
    x1 = rep(xpos + 1.025 * width, length(ytick)),
    y0 = ypos + ytick / max_n * height,
    lwd = .5
  )

  ## plot minimum if exists
  if (minimum > 0) {
    text(xpos + -0.05 * width, ypos + minimum / max_n * height, "Min.", cex = .8, adj = c(1, 0.5))
    segments(
      x0 = rep(xpos, length(ytick)),
      x1 = rep(xpos + 1 * width, length(ytick)),
      y0 = ypos + minimum / max_n * height,
      col = "darkgray", lwd = .5, lty = 1
    )
    segments(
      x0 = rep(xpos - 0.025 * width, length(ytick)),
      x1 = rep(xpos + 0.025 * width, length(ytick)),
      y0 = ypos + minimum / max_n * height,
      lwd = .5
    )
  }

  ## add title 
  if (is.null(main)) {
    main <- sprintf("n=%d", sum(n[n > minimum]))
  }

  text(xpos + width / 2, ypos + 1.1 * height, font = 2, cex = 0.8, main)
}

## FIXME: (ML) Idea to setup up a crch specific reliagram:
##   * Set threshold to censor point
##   * Only programming outline
# reliagram.crch <- function(object,
#                           newdata = NULL,
#                           breaks = seq(0, 1, by = 0.1),
#                           thresholds = NULL,
#                           ...) {
#   if((missing(newdata) || is.null(newdata)) && !is.null(object$y)) {
#     y <- object$y
#   } else {
#     y <- newresponse(object, newdata = newdata)
#   }
#
#   if(is.null(thresholds)) thresholds <- object$left
#
#   ## call reliagram
#   reliagram(object, breaks = breaks, thresholds = thresholds, y = y, ...)
# }
