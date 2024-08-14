#' Residuals for Probabilistic Regression Models
#' 
#' Generic function and default method for (randomized) quantile residuals, PIT,
#' Pearson, and raw response residuals based on \pkg{distributions3} support.
#' 
#' The new generic function \code{proresiduals} comes with a powerful default
#' method that is based on the following idea: \code{\link[topmodels]{newresponse}}
#' and \code{\link[distributions3]{prodist}} can be used to extract the observed
#' response and expected distribution for it, respectively. For all model classes
#' that have methods for these two generic functions, \code{proresiduals} can
#' compute a range of different \code{type}s of residuals.
#' 
#' The simplest definition of residuals are the so-called \code{"response"} residuals
#' which simply compute the difference between the observations and the expected means.
#' The \code{"pearson"} residuals additionally standardize these residuals by the
#' square root of the expected variance. Thus, these residuals are based only on the
#' first and on the first two moments, respectively.
#' 
#' To assess the entire distribution and not just the first moments, there are also
#' residuals based on the probability integral transform (PIT).
#' For regression models with a continuous response distribution, \code{"pit"} residuals
#' (see Warton 2007) are simply the expected cumulative distribution (CDF) evaluated at the
#' observations (Dawid, 1984). For discrete distributions, a uniform random value is drawn
#' from the range of probabilities between the CDF at the observation and the supremum
#' of the CDF to the left of it. If the model fits well the PIT residuals should be uniformly
#' distributed.
#'
#' In order to obtain normally distributed residuals for well-fitting models (like often
#' desired in linear regression models), \code{"quantile"} residuals, proposed by Dunn and
#' Smyth (1996), additionally transform the PIT residuals by the standard normal quantile function.
#'
#' As quantile residuals and PIT residuals are subject to randomness for discrete distributions
#' (and also for mixed discrete-continuous distributions), it is sometimes
#' useful to explore the extent of the random variation. This can be done either by
#' obtaining multiple replications (via \code{nsim}) or by computing fixed quantiles
#' of each probability interval such as \code{prob = 0.5} (corresponding to mid-quantile residuals,
#' see Feng et al. 2020). Another common setting is \code{prob = c(0, 1)} yielding the range of
#' possible residuals.
#'
#' @param object an object for which a \code{\link[topmodels]{newresponse}} and a
#'   \code{\link[distributions3]{prodist}} method is available.
#' @param newdata optionally, a data frame in which to look for variables with
#'   which to predict. If omitted, the original observations are used.
#' @param type character indicating whether quantile (default), PIT, Pearson, or raw response
#'   residuals should be computed.
#' @param nsim integer. The number of randomly simulated residuals of \code{type = "quantile"}
#'   or \code{"pit"}. By default one simulation is returned.
#' @param prob numeric. Instead of simulating the probabilities (between 0 and 1) for
#'   \code{type = "quantile"} or \code{"pit"}, a vector of probabilities can be specified,
#'   e.g., \code{prob = 0.5} corresponding to mid-quantile residuals.
#' @param delta numeric. The minimal difference to compute the range of
#'   proabilities corresponding to each observation according to get (randomized)
#'   \code{"quantile"} or \code{"pit"} residuals.  For \code{NULL}, the minimal observed difference in the
#'   resonse divided by \code{5e-6} is used. Ignored for continuous distributions.
#' @param \dots further parameters passed to methods.
#'
#' @return A vector or matrix of residuals. A matrix of residuals is returned
#' if more than one replication of quantile or PIT residuals is computed, i.e., if either
#' \code{random > 1} or \code{random = FALSE} and \code{length(prob) > 1}.
#'
#' @seealso \code{\link[stats]{qnorm}}, \code{\link{qqrplot}}
#'
#' @references
#' Dawid AP (1984).
#' \dQuote{Present Position and Potential Developments: Some Personal Views:
#' Statistical Theory: The Prequential Approach.}
#' \emph{Journal of the Royal Statistical Society A}, \bold{147}(2), 278--292.
#' \doi{10.2307/2981683}.
#' 
#' Dunn KP, Smyth GK (1996).
#' \dQuote{Randomized Quantile Residuals.}
#' \emph{Journal of Computational and Graphical Statistics}, \bold{5}(3), 236--244.
#' \doi{10.2307/1390802}
#'
#' Feng C, Li L, Sadeghpour A (2020).
#' \dQuote{A Comparison of Residual Diagnosis Tools for Diagnosing Regression Models for Count Data}
#' \emph{BMC Medical Research Methodology}, \bold{20}(175), 1--21.
#' \doi{10.1186/s12874-020-01055-2}
#' 
#' Warton DI, Thibaut L, Wang YA (2017)
#' \dQuote{The PIT-Trap -- A \sQuote{Model-Free} Bootstrap Procedure for Inference
#' about Regression Models with Discrete, Multivariate Responses}.
#' \emph{PLOS ONE}, \bold{12}(7), 1--18.
#' \doi{10.1371/journal.pone.0181790}.
#' 
#' @keywords regression
#'
#' @examples
#' ## Poisson GLM for FIFA 2018 data
#' data("FIFA2018", package = "distributions3")
#' m <- glm(goals ~ difference, data = FIFA2018, family = poisson)
#'
#' ## random quantile residuals (on original data)
#' proresiduals(m)
#' 
#' ## various flavors of residuals on small new data
#' nd <- data.frame(goals = c(1, 1, 1), difference = c(-1, 0, 1))
#'
#' ## quantile residuals: random (1 sample), random (5 samples), mid-quantile (non-random)
#' proresiduals(m, newdata = nd, type = "quantile")
#' proresiduals(m, newdata = nd, type = "quantile", nsim = 5)
#' proresiduals(m, newdata = nd, type = "quantile", prob = 0.5)
#' 
#' ## PIT residuals (without transformation to normal): random vs. minimum/maximum quantile
#' proresiduals(m, newdata = nd, type = "pit", nsim = 5)
#' proresiduals(m, newdata = nd, type = "pit", prob = c(0, 1))
#' 
#' ## raw response residuals (observation - expected mean)
#' proresiduals(m, newdata = nd, type = "response")
#' 
#' ## standardized Pearson residuals (response residuals divided by standard deviation)
#' proresiduals(m, newdata = nd, type = "pearson")
#' 
#' ## compute residuals by manually obtaining distribution and response
#' ## proresiduals(procast(m, newdata = nd, drop = TRUE), nd$goals)

#' @export
proresiduals <- function(object, ...) {
  UseMethod("proresiduals")
}

#' @rdname proresiduals
#' @method proresiduals default
#' @importFrom distributions3 prodist cdf is_continuous
#' @export
proresiduals.default <- function(object, newdata = NULL, type = c("quantile", "pit", "pearson", "response"),
  nsim = NULL, prob = NULL, delta = NULL, ...) {

  ## type of residuals
  type <- match.arg(tolower(as.character(type))[1L], c("quantile", "pit", "pearson", "response", "raw"))
  if (type == "raw") type <- "response"

  ## extract observed response
  y <- newresponse(object, newdata = newdata, na.action = na.pass)
  y[["(weights)"]] <- NULL
  if (ncol(y) > 1L) stop("multivariate responses not supported yet")
  y <- y[[1L]]
  n <- NROW(y)

  ## predicted distribution
  pd <- if (is.null(newdata)) prodist(object, ...) else prodist(object, newdata = newdata, ...)

  ## number of random replications (if applicable) and
  ## quantile probabilities: sampled randomly vs. specified by prob
  if (type %in% c("pearson", "response")) {
    if (!is.null(nsim)) warning(sprintf("argument 'nsim' ignored for %s residuals", type))
    if (!is.null(prob)) warning(sprintf("argument 'prob' ignored for %s residuals", type))
  } else {
    if (random <- is.null(prob)) {
      if (is.null(nsim)) nsim <- 1L
      if (!is.numeric(nsim) || length(nsim) != 1L) stop("'nsim' must be an integer of length 1")
      nsim <- as.integer(nsim)
      nc <- nsim
      clab <- 1L:nc
    } else {
      if (!is.null(nsim)) warning(sprintf("argument 'nsim' ignored for 'prob'-based %s residuals", type))    
      prob <- as.numeric(prob)
      if (any(prob < 0) || any(prob > 1)) stop("'prob' for quantiles must be in [0, 1]")
      clab <- format(prob, digits = 3L, trim = TRUE, drop0trailing = TRUE)
      nc <- length(prob)
    }
  }

  ## minimal distance for observation to the left of y
  stopifnot(is.null(delta) || (is.numeric(delta) && length(delta) == 1 && delta > 0.0))
  if (is.null(delta)) {
    delta <- .Machine$double.eps^0.33
    y_unique <- sort(unique(y))
    if (length(y_unique) > 1L) delta <- pmax(min(diff(y_unique)) / 5e6, delta)
  }

  ## simple special cases: raw response or Pearson residuals
  if (type %in% c("pearson", "response")) {
    res <- y - mean(pd)
    if (type == "pearson") res <- res/sqrt(variance(pd))
    return(res)
  }

  ## simple special cases: continuous quantile or PIT residuals, replicated if nsim > 1 or length(prob) > 1
  if (all(is_continuous(pd))) {
    res <- cdf(pd, y, elementwise = TRUE)
    if (type == "quantile") {
      res <- qnorm(res)
      ## try to catch infinite quantile residuals due to CDF = 0 or = 1
      if (length(ix <- which(!is.finite(res))) > 0L) {
        res[ix] <- qnorm(cdf(pd[ix], y[ix], elementwise = TRUE, log.p = TRUE), log.p = TRUE)
      }
    }
    if (nc > 1L) {
      res <- rep.int(res, nc)
      res <- matrix(res, nrow = n, ncol = nc, dimnames = list(names(y), paste("r", clab, sep = "_")))
    }
    return(res)
  }
  
  ## otherwise: obtain some PIT flavor for distribution with point masses and optionally transform to normal scale

  ## quantile probabilities simulated or replicated
  prob <- if (random) runif(n * nc) else rep.int(prob, rep.int(n, nc))
  prob <- matrix(prob, nrow = n, ncol = nc, dimnames = list(names(y), paste("r", clab, sep = "_")))

  ## probability integral transform: CDF at y and supremum left of y
  pit_y   <- cdf(pd, y, elementwise = TRUE)
  pit_sup <- cdf(pd, y - delta, elementwise = TRUE)
  
  ## no interpolation necessary for y observations with continous CDF
  if (!all(is_discrete(pd)) && any(icont <- abs(pit_y - pit_sup) < .Machine$double.eps^0.33)) prob[icont, ] <- 1

  ## interpolate between PIT values, drop dimension if possible
  if (nc == 1L) prob <- prob[, 1L]
  res <- (1 - prob) * pit_sup + prob * pit_y

  ## for quantile residuals: map to normal scale
  if (type == "quantile") res <- qnorm(res)

  return(res)
}
