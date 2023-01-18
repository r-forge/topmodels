#' Procast: Probabilistic Forecasting
#' 
#' Generic function and methods for computing various kinds of probabilistic
#' forecasts from (regression) models.
#' 
#' The function \code{procast} provides a unified framework for probabilistic
#' forcasting (or procasting, for short) based on probabilistic (regression)
#' models, also known as distributional regression approaches. Typical types
#' of predictions include quantiles, probabilities, (conditional) expectations,
#' variances, and (log-)densities. Internally, \code{procast} methods typically
#' compute the predicted parameters for each observation and then compute the
#' desired outcome for the distributions with the respective parameters.
#'
#' Some quantities, e.g., the moments of the distribution (like mean or variance),
#' can be computed directly from the predicted parameters of the
#' distribution while others require an additional argument \code{at} which the
#' distribution is evaluated (e.g., the probability of a quantile or an
#' observation of the response).
#' 
#' The default \code{procast} method leverages the S3 classes and methods for
#' probability distributions from the \pkg{distributions3} package. In a first step
#' the predicted probability distribution object is obtained and, by default
#' (\code{type = "distribution"}), returned in order to reflect the distributional
#' nature of the forecast. For all other \code{type}s (e.g., \code{"mean"},
#' \code{"quantile"}, or \code{"density"}), the corresponding extractor methods
#' (e.g., \code{mean}, \code{quantile}, or \code{\link[distributions3]{pdf}}) are used to
#' compute the desired quantity from the distribution objects. The examples
#' provide some worked illustrations.
#'
#' Package authors or users, who want to enable \code{procast} for new types
#' of model objects, only need to provide a suitable \code{\link[distributions3]{prodist}}
#' extractor for the predicted probability distribution. Then the default \code{procast}
#' works out of the box. However, if the \pkg{distributions3} package does not support
#' the necessary probability distribution, then it may also be necessary to
#' implement a new distribution objects, see \code{\link[distributions3]{apply_dpqr}}.
#'
#' @aliases procast
#' @param object a fitted model object. For the \code{default} method this
#' needs to have a \code{\link[distributions3]{prodist}} method (or \code{object}
#' can inherit from \code{distribution} directly).
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param na.action function determining what should be done with missing
#' values in \code{newdata}.  The default is to employ \code{NA}.
#' @param type character specifying the type of probabilistic forecast to
#' compute. Note that \code{type = "probability"} corresponds to cumulative
#' probability as in \code{pnorm}, \code{pbinom}, etc.
#' @param at specification of values at which the forecasts should be
#' evaluated, typically a numeric vector but possibly also a matrix or data
#' frame.  Additionally, \code{at} can be the character string
#' \code{"function"} or \code{"list"}, see details below.
#' @param drop logical. Should forecasts be returned in a data frame (default)
#' or (if possible) dropped to a vector, see return value description below.
#' @param \dots further parameters passed to methods. In particular, this includes
#' the logical argument \code{elementwise = NULL}. Should each element of distribution only be evaluated at the
#' corresponding element of \code{at} (\code{elementwise = TRUE}) or at all elements
#' in \code{at} (\code{elementwise = FALSE}). Elementwise evaluation is only possible
#' if the number of observations is length of \code{at} are the same and in that case a vector of
#' the same length is returned. Otherwise a matrix is returned. The default is to use
#' \code{elementwise = TRUE} if possible, and otherwise \code{elementwise = FALSE}.
#' @param sigma character or numeric or \code{NULL}. Specification of the standard
#' deviation \code{sigma} to be used for the \code{\link[distributions3]{Normal}} distribution in the
#' \code{lm} method. The default \code{"ML"} (or equivalently \code{"MLE"} or \code{NULL})
#' uses the maximum likelihood estimate based on the residual sum of squares divided
#' by the number of observations, n. Alternatively, \code{sigma = "OLS"} uses the
#' least-squares estimate (divided by the residual degrees of freedom, n - k). Finally,
#' a concrete numeric value can also be specified in \code{sigma}.
#' @param dispersion character or numeric or \code{NULL}. Specification of the
#' dispersion parameter in the \code{glm} method. The default \code{NULL}
#' (or equivalently \code{"deviance"}) is to use the \code{\link[stats]{deviance}}
#' divided by the number of observations, n. Alternatively, \code{dispersion = "Chisquared"}
#' uses the Chi-squared statistic divided by the residual degrees of freedom, n - k.
#' Finally, a concrete numeric value can also be specified in \code{dispersion}.
#' @return Either a \code{data.frame} of predictions with the same number of rows
#' as the \code{newdata} (or the original observations if that is \code{NULL}).
#' If \code{drop = TRUE} predictions with just a single column are simplified
#' to a vector and predictions with multiple columns to a matrix.
#' @keywords regression
#' @examples
#' ## load packages
#' library("topmodels")
#' library("distributions3")
#' 
#' ## Poisson regression model for FIFA 2018 data:
#' ## number of goals scored by each team in each game, explained by
#' ## predicted ability difference of the competing teams
#' data("FIFA2018", package = "distributions3")
#' m <- glm(goals ~ difference, data = FIFA2018, family = poisson)
#' 
#' ## predicted probability distributions for all matches (in sample)
#' head(procast(m))
#' head(procast(m, drop = TRUE))
#' 
#' ## procasts for new data
#' ## much lower, equal, and much higher ability than opponent
#' nd <- data.frame(difference = c(-1, 0, 1))
#'
#' ## predicted goal distribution object
#' goals <- procast(m, newdata = nd, drop = TRUE)
#' goals
#'
#' ## predicted densities/probabilities for scoring 0, 1, ..., 5 goals
#' procast(m, newdata = nd, type = "density", at = 0:5)
#' ## by hand
#' pdf(goals, 0:5)
#'
#' ## means and medians
#' procast(m, newdata = nd, type = "mean")
#' procast(m, newdata = nd, type = "quantile", at = 0.5)
#' ## by hand
#' mean(goals)
#' quantile(goals, 0.5)
#'
#' ## evaluate procast elementwise or for all possible combinations
#' ## of distributions from 'nd' and observations in 'at'
#' procast(m, newdata = nd, type = "probability", at = 1:3, elementwise = TRUE)
#' procast(m, newdata = nd, type = "probability", at = 1:3, elementwise = FALSE)
#' 
#' ## compute in-sample log-likelihood sum via procast
#' sum(procast(m, type = "density", at = FIFA2018$goals, log = TRUE))
#' logLik(m)
#' 
#' @export 
procast <- function(object, newdata = NULL, na.action = na.pass, type = "distribution", at = 0.5, drop = FALSE, ...) {
  UseMethod("procast")
}

#' @rdname procast
#' @importFrom distributions3 prodist cdf pdf log_pdf variance skewness kurtosis
#' @export 
procast.default <- function(object, newdata = NULL, na.action = na.pass,
  type = c("distribution", "mean", "variance", "quantile", "probability", "density", "loglikelihood", "parameters", "kurtosis", "skewness"),
  at = 0.5, drop = FALSE, ...)
{
  ## match type
  type <- match.arg(type[1L], c(
    "quantile", "mean", "variance",
    "probability", "cdf",
    "density", "pdf", "pmf",
    "loglikelihood", "log_pdf",
    "distribution", "parameters",
    "kurtosis", "skewness"))
  if(type == "cdf") type <- "probability"
  if(type %in% c("pdf", "pmf")) type <- "density"
  if(type == "log_pdf") type <- "loglikelihood"
  label <- type
  
  ## TODO: one could also support type = function but it's not clear how easy that would really be for the users...
  ## fargs <- FALSE
  ## if(!is.function(type)) stop("'type' must either be a character string or a function")
  ## label <- "user-defined" ## FIXME
  ## f <- type
  ## fargs <- length(setdiff(names(formals(f)), "...")) > 1L
  ## type <- "custom"

  ## FIXME: how to handle 'size' in binomial family?
  ## extract probability distribution object
  pd <- if(inherits(object, "distribution")) {
    object
  } else if(is.null(newdata)) {
    distributions3::prodist(object)
  } else {
    distributions3::prodist(object, newdata = newdata, na.action = na.action)
  }
  
  ## evaluate type of procast
  pc <- switch(type,
    "distribution"  = pd,
    "quantile"      = quantile(pd, at, ...),
    "mean"          = mean(pd),
    "variance"      = distributions3::variance(pd),
    "probability"   = distributions3::cdf(pd, at, ...),
    "density"       = distributions3::pdf(pd, at, ...),
    "loglikelihood" = distributions3::log_pdf(pd, at, ...),
    "parameters"    = as.matrix(pd),
    "skewness"      = distributions3::skewness(pd, at, ...),
    "kurtosis"      = distributions3::kurtosis(pd, at, ...)
    ## "custom"        = if(fargs) f(pd, at, ...) else f(pd, ...)
  )
  
  ## convert to data frame if drop = FALSE
  if(drop) {
    if(!is.null(dim(pc)) && NCOL(pc) == 1L) pc <- drop(pc)
  } else {
    if(inherits(pc, "distribution")) {
      pc <- as.data.frame(pc)
      colnames(pc) <- label
    }
    if(is.null(dim(pc))) {
      pc <- as.matrix(pc)
      if(ncol(pc) == 1L) colnames(pc) <- label
    }
    if(!inherits(pc, "data.frame")) pc <- as.data.frame(pc)
  }
  
  return(pc)
}

#' @rdname procast
#' @method procast lm
#' @export 
procast.lm <- function(object, newdata = NULL, na.action = na.pass, type = "distribution", at = 0.5, drop = FALSE, ..., sigma = "ML") {
  object <- if(is.null(newdata)) {
    distributions3::prodist(object, sigma = sigma)
  } else {
    distributions3::prodist(object, newdata = newdata, na.action = na.action, sigma = sigma)
  }
  procast(object, newdata = newdata, na.action = na.action, type = type, at = at, drop = drop, ...)
}

#' @rdname procast
#' @method procast glm
#' @export 
procast.glm <- function(object, newdata = NULL, na.action = na.pass, type = "distribution", at = 0.5, drop = FALSE, ..., dispersion = NULL) {
  object <- if(is.null(newdata)) {
    distributions3::prodist(object, dispersion = dispersion)
  } else {
    distributions3::prodist(object, newdata = newdata, na.action = na.action, dispersion = dispersion)
  }
  procast(object, newdata = newdata, na.action = na.action, type = type, at = at, drop = drop, ...)
}

#' @rdname procast
#' @method procast bamlss
#' @export 
procast.bamlss <- function(object, newdata = NULL, na.action = na.pass, type = "distribution", at = 0.5, drop = FALSE, ..., distributions3 = FALSE) {
  object <- if(is.null(newdata)) {
    distributions3::prodist(object, distributions3 = distributions3)
  } else {
    distributions3::prodist(object, newdata = newdata, na.action = na.action, distributions3 = distributions3)
  }
  procast(object, newdata = newdata, na.action = na.action, type = type, at = at, drop = drop, ...)
}

