#' Create a GAMLSS Distribution
#'
#' A single class and corresponding methods encompassing all distributions
#' from the \pkg{gamlss.dist} package using the workflow from the
#' \pkg{distributions3} package.
#' 
#' The constructor function \code{GAMLSS} sets up a distribution
#' object, representing a distribution from the GAMLSS (generalized additive
#' model of location, scale, and shape) framework by the corresponding parameters
#' plus a \code{family} attribute, e.g., \code{\link[gamlss.dist]{NO}} for the
#' normal distribution or \code{\link[gamlss.dist]{BI}} for the binomial
#' distribution. There can be up to four parameters, called \code{mu} (often some
#' sort of location parameter), \code{sigma} (often some sort of scale parameter),
#' \code{tau} and \code{nu} (other parameters, e.g., capturing shape, etc.).
#'
#' All parameters can also be vectors, so that it is possible to define a vector
#' of GAMLSS distributions from the same family with potentially different parameters.
#' All parameters need to have the same length or must be scalars (i.e.,
#' of length 1) which are then recycled to the length of the other parameters.
#' 
#' Note that not all distributions use all four parameters, i.e., some use just a
#' subset. In that case, the corresponding arguments in \code{GAMLSS} should be
#' unspecified, \code{NULL}, or \code{NA}.
#' 
#' For the \code{GAMLSS} distribution objects there is a wide range
#' of standard methods available to the generics provided in the \pkg{distributions3}
#' package: \code{\link[distributions3]{pdf}} and \code{\link[distributions3]{log_pdf}}
#' for the (log-)density (PDF), \code{\link[distributions3]{cdf}} for the probability
#' from the cumulative distribution function (CDF), \code{quantile} for quantiles,
#' \code{\link[distributions3]{random}} for simulating random variables,
#' and \code{\link[distributions3]{support}} for the support interval
#' (minimum and maximum). Internally, these methods rely on the usual d/p/q/r
#' functions provided in \pkg{gamlss.dist}, see the manual pages of the individual
#' families. The methods \code{\link[distributions3]{is_discrete}} and
#' \code{\link[distributions3]{is_continuous}} can be used to query whether the
#' distributions are discrete on the entire support or continuous on the entire
#' support, respectively.
#'
#' Additionally, for some families there is also a \code{\link[scoringRules]{crps}}
#' method for computing the continuous ranked probability score (CRPS) via the
#' \pkg{scoringRules} package. This is only available for those families which are
#' supported by both packages.
#' 
#' See the examples below for an illustration of the workflow for the class and methods.
#' 
#' @param family character. Name of a GAMLSS family provided by \pkg{gamlss.dist}, e.g.,
#'   \code{\link[gamlss.dist]{NO}} or \code{\link[gamlss.dist]{BI}} for the normal or
#'   binomial distribution, respectively.
#' @param mu numeric. GAMLSS \code{mu} parameter. Can be a scalar or a vector
#'   or missing if not part of the \code{family}.
#' @param sigma numeric. GAMLSS \code{sigma} parameter. Can be a scalar or a vector
#'   or missing if not part of the \code{family}.
#' @param tau numeric. GAMLSS \code{tau} parameter. Can be a scalar or a vector
#'   or missing if not part of the \code{family}.
#' @param nu numeric. GAMLSS \code{nu} parameter. Can be a scalar or a vector
#'   or missing if not part of the \code{family}.
#' 
#' @return A \code{GAMLSS} distribution object.
#' 
#' @seealso \code{\link[gamlss.dist]{gamlss.family}}
#' 
#' @examples
#' \dontshow{ if(!requireNamespace("gamlss.dist")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("not all packages required for the example are installed")
#'   } else q() }
#' }
#' ## package and random seed
#' library("distributions3")
#' set.seed(6020)
#' 
#' ## three Weibull distributions
#' X <- GAMLSS("WEI", mu = c(1, 1, 2), sigma = c(1, 2, 2))
#' X
#' 
#' ## moments
#' mean(X)
#' variance(X)
#' 
#' ## support interval (minimum and maximum)
#' support(X)
#' is_discrete(X)
#' is_continuous(X)
#' 
#' ## simulate random variables
#' random(X, 5)
#' 
#' ## histograms of 1,000 simulated observations
#' x <- random(X, 1000)
#' hist(x[1, ], main = "WEI(1,1)")
#' hist(x[2, ], main = "WEI(1,2)")
#' hist(x[3, ], main = "WEI(2,2)")
#' 
#' ## probability density function (PDF) and log-density (or log-likelihood)
#' x <- c(2, 2, 1)
#' pdf(X, x)
#' pdf(X, x, log = TRUE)
#' log_pdf(X, x)
#' 
#' ## cumulative distribution function (CDF)
#' cdf(X, x)
#' 
#' ## quantiles
#' quantile(X, 0.5)
#' 
#' ## cdf() and quantile() are inverses
#' cdf(X, quantile(X, 0.5))
#' quantile(X, cdf(X, 1))
#' 
#' ## all methods above can either be applied elementwise or for
#' ## all combinations of X and x, if length(X) = length(x),
#' ## also the result can be assured to be a matrix via drop = FALSE
#' p <- c(0.05, 0.5, 0.95)
#' quantile(X, p, elementwise = FALSE)
#' quantile(X, p, elementwise = TRUE)
#' quantile(X, p, elementwise = TRUE, drop = FALSE)
#' 
#' ## compare theoretical and empirical mean from 1,000 simulated observations
#' cbind(
#'   "theoretical" = mean(X),
#'   "empirical" = rowMeans(random(X, 1000))
#' )
#' @export
GAMLSS <- function(family, mu, sigma, tau, nu) {
  stopifnot(requireNamespace("gamlss.dist"))
  ## get family object
  f <- .gamlss_family(family)

  ## get parameters
  par <- names(f$parameters)
  if(!("mu" %in% par)    && !missing(mu)    && (!is.null(mu)    || !is.na(mu)))    warning(sprintf("'mu' is not a parameter of the '%s' family", family))
  if(!("sigma" %in% par) && !missing(sigma) && (!is.null(sigma) || !is.na(sigma))) warning(sprintf("'sigma' is not a parameter of the '%s' family", family))
  if(!("tau" %in% par)   && !missing(tau)   && (!is.null(tau)   || !is.na(tau)))   warning(sprintf("'tau' is not a parameter of the '%s' family", family))
  if(!("nu" %in% par)    && !missing(nu)    && (!is.null(nu)    || !is.na(nu)))    warning(sprintf("'nu' is not a parameter of the '%s' family", family))
  
  ## set up distribution
  d <- eval(parse(text = sprintf("data.frame(%s)", paste(par, "=", par, collapse = ", "))))
  class(d) <- c("GAMLSS", "distribution")
  attr(d, "family") <- f$family
  return(d)
}

## auxiliary functions for getting family and d/p/q/r function
## from gamlss.dist namespace (FIXME: where to get() for user-defined
## family?
.gamlss_family <- function(x) get(if(is.character(x)) x else attr(x, "family")[1L], asNamespace("gamlss.dist"))()
.gamlss_dpqr <- function(x, type) get(paste0(type, attr(x, "family")[1L]), asNamespace("gamlss.dist"))


#' @rdname GAMLSS
#' @method mean GAMLSS
#' @export
#' @usage NULL
#' @importFrom stats setNames
mean.GAMLSS <- function(x, ...) {
  stopifnot(requireNamespace("gamlss.dist"))
  f <- .gamlss_family(x)
  m <- do.call(f$mean, as.list(x))
  setNames(m, names(x))
}

#' @rdname GAMLSS
#' @method variance GAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 variance
#' @importFrom stats setNames
variance.GAMLSS <- function(x, ...) {
  stopifnot(requireNamespace("gamlss.dist"))
  f <- .gamlss_family(x)
  m <- do.call(f$variance, as.list(x))
  setNames(m, names(x))
}

#' @rdname GAMLSS
#' @method skewness GAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 skewness
skewness.GAMLSS <- function(x, ...) {
  stop("not yet implemented")
}

#' @rdname GAMLSS
#' @method kurtosis GAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 kurtosis
kurtosis.GAMLSS <- function(x, ...) {
  stop("not yet implemented")
}

#' @rdname GAMLSS
#' @method random GAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 random apply_dpqr make_positive_integer
random.GAMLSS <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("gamlss.dist"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) do.call(.gamlss_dpqr(d, "r"), c(list(n = at), as.list(d)))
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' @rdname GAMLSS
#' @importFrom distributions3 pdf apply_dpqr
#' @method pdf GAMLSS
#' @export
#' @usage NULL
pdf.GAMLSS <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("gamlss.dist"))
  FUN <- function(at, d) do.call(.gamlss_dpqr(d, "d"), c(list(x = at), as.list(d), ...))
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname GAMLSS
#' @method log_pdf GAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 log_pdf apply_dpqr
log_pdf.GAMLSS <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("gamlss.dist"))
  FUN <- function(at, d) do.call(.gamlss_dpqr(d, "d"), c(list(x = at), as.list(d), list(log = TRUE)))
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' @rdname GAMLSS
#' @method cdf GAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 cdf apply_dpqr
cdf.GAMLSS <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("gamlss.dist"))
  FUN <- function(at, d) do.call(.gamlss_dpqr(d, "p"), c(list(q = at), as.list(d), ...))
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' @rdname GAMLSS
#' @method quantile GAMLSS
#' @export
#' @usage NULL
#' @importFrom stats quantile
#' @importFrom distributions3 apply_dpqr
quantile.GAMLSS <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("gamlss.dist"))
  FUN <- function(at, d) do.call(.gamlss_dpqr(d, "q"), c(list(p = at), as.list(d), ...))
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' @rdname GAMLSS
#' @method support GAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 support make_support
support.GAMLSS <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("gamlss.dist"))
  s <- quantile(d, probs = c(0, 1), elementwise = FALSE)
  distributions3::make_support(s[, 1L], s[, 2L], d, drop = drop)
}

#' @rdname GAMLSS
#' @method is_discrete GAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 is_discrete
#' @importFrom stats setNames
is_discrete.GAMLSS <- function(d, ...) {
  stopifnot(requireNamespace("gamlss.dist"))
  f <- .gamlss_family(d)
  setNames(rep.int(f$type == "Discrete", length(d)), names(d))
}

#' @rdname GAMLSS
#' @method is_continuous GAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 is_continuous
#' @importFrom stats setNames
is_continuous.GAMLSS <- function(d, ...) {
  stopifnot(requireNamespace("gamlss.dist"))
  f <- .gamlss_family(d)
  setNames(rep.int(f$type == "Continuous", length(d)), names(d))
}

#' @rdname GAMLSS
#' @method format GAMLSS
#' @export
#' @usage NULL
format.GAMLSS <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
  class(x) <- c(paste("GAMLSS", attr(x, "family")[1L]), "distribution")
  NextMethod()
}

#' @rdname GAMLSS
#' @method print GAMLSS
#' @export
#' @usage NULL
print.GAMLSS <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
  class(x) <- c(paste("GAMLSS", attr(x, "family")[1L]), "distribution")
  NextMethod()
}

#' @rdname GAMLSS
#' @exportS3Method scoringRules::crps GAMLSS
#' @usage NULL
#' @importFrom distributions3 apply_dpqr
crps.GAMLSS <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  ## manually match gamlss.dist distributions names with scoringRules
  f <- attr(y, "family")[1L]
  FUN <- switch(EXPR = f,
    "NO"  = function(at, d) scoringRules::crps_norm(y = at, mean = d$mu, sd = d$sigma),
    "LO"  = function(at, d) scoringRules::crps_logis(y = at, mean = d$mu, sd = d$sigma),
    "TF"  = function(at, d) scoringRules::crps_t(y = at, location = d$mu, scale = d$sigma, df = d$nu),
    "LNO" = function(at, d) scoringRules::crps_lnorm(y = at, meanlog = d$mu, sdlog = d$sigma),
    "PO"  = function(at, d) scoringRules::crps_pois(y = at, lambda = d$mu),
    "BI"  = function(at, d) scoringRules::crps_binom(y = at, prob = d$mu, size = 1L), ## FIXME: size?
    "NBI" = function(at, d) scoringRules::crps_nbinom(y = at, size = 1/d$sigma, mu = d$mu),
    "NBII"= function(at, d) scoringRules::crps_nbinom(y = at, size = d$mu/d$sigma, mu = d$m),
    "EXP" = function(at, d) scoringRules::crps_exp(y = at, rate = 1/d$mu),
    "GA"  = function(at, d) scoringRules::crps_gamma(y = at, shape = 1/d$sigma^2, scale = d$mu * d$sigma^2),
    "BE"  = function(at, d) scoringRules::crps_beta(y = at, shape1 = d$mu * (1 - d$sigma^2)/(d$sigma^2), shape2 = (1 - d$mu) * (1 - d$sigma^2)/(d$sigma^2)),
    NULL ## FIXME: scoringRules also has gpd, gev, lapl, cnorm, tnorm, ...
  )
  if(is.null(FUN)) {
    stop(sprintf("crps() not supported yet for the '%s' family", f))
  } else {
    distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
  }
}



#' Extracting Fitted or Predicted Probability Distributions from gamlss Models
#' 
#' Methods for \pkg{gamlss} model objects for extracting
#' fitted (in-sample) or predicted (out-of-sample) probability distribution
#' objects.
#' 
#' To facilitate making probabilistic forecasts based on \code{\link[gamlss]{gamlss}}
#' model objects, the \code{\link[distributions3]{prodist}} method extracts fitted or
#' predicted probability \code{distribution} objects. Internally, the
#' \code{\link[gamlss]{predictAll}} method from the \pkg{gamlss} package is
#' used first to obtain the distribution parameters (\code{mu}, \code{sigma}, \code{tau},
#' \code{nu}, or a subset thereof). Subsequently, the corresponding \code{distribution}
#' object is set up using \code{\link{GAMLSS}}, enabling the workflow provided by
#' the \pkg{distributions3} package.
#' 
#' Note that these probability distributions only reflect the random variation in
#' the dependent variable based on the model employed (and its associated
#' distributional assumption for the dependent variable). This does not capture
#' the uncertainty in the parameter estimates.
#' 
#' @param object A model object of class \code{\link[gamlss]{gamlss}}.
#' @param ... Arguments passed on to \code{\link[gamlss]{predictAll}}, 
#' e.g., \code{newdata}.
#' 
#' @return An object inheriting from \code{distribution}.
#' 
#' @seealso \code{\link{GAMLSS}} \code{\link[gamlss]{predictAll}}
#' 
#' @keywords distribution
#' 
#' @examples
#' \dontshow{ if(!requireNamespace("gamlss")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("not all packages required for the example are installed")
#'   } else q() }
#' }
#' ## packages, code, and data
#' library("gamlss")
#' library("distributions3")
#' data("cars", package = "datasets")
#' 
#' ## fit heteroscedastic normal GAMLSS model
#' m <- gamlss(dist ~ pb(speed), ~ pb(speed), data = cars, family = "NO")
#' 
#' ## obtain predicted distributions for three levels of speed
#' d <- prodist(m, newdata = data.frame(speed = c(10, 20, 30)))
#' print(d)
#' 
#' ## obtain quantiles (works the same for any distribution object 'd' !!)
#' quantile(d, 0.5)
#' quantile(d, c(0.05, 0.5, 0.95), elementwise = FALSE)
#' quantile(d, c(0.05, 0.5, 0.95), elementwise = TRUE)
#' 
#' ## visualization
#' plot(dist ~ speed, data = cars)
#' nd <- data.frame(speed = 0:240/4)
#' nd$dist <- prodist(m, newdata = nd)
#' nd$fit <- quantile(nd$dist, c(0.05, 0.5, 0.95))
#' matplot(nd$speed, nd$fit, type = "l", lty = 1, col = "slategray", add = TRUE)
#' 
#' ## moments
#' mean(d)
#' variance(d)
#' 
#' ## simulate random numbers
#' random(d, 5)
#' 
#' ## density and distribution
#' pdf(d, 50 * -2:2)
#' cdf(d, 50 * -2:2)
#' 
#' ## further diagnostics: graphical and scores
#' pithist(m)
#' qqrplot(m)
#' proscore(m, type = c("LogLik", "CRPS", "MAE", "MSE"), aggregate = TRUE)
#' 
#' ## note that proscore can replicate logLik() value
#' proscore(m, aggregate = sum)
#' logLik(m)
#' 
#' ## Poisson example
#' data("FIFA2018", package = "distributions3")
#' m2 <- gamlss(goals ~ pb(difference), data = FIFA2018, family = "PO")
#' d2 <- prodist(m2, newdata = data.frame(difference = 0))
#' print(d2)
#' quantile(d2, c(0.05, 0.5, 0.95))
#' @export
prodist.gamlss <- function(object, ...) {
  stopifnot(requireNamespace("gamlss"))
  d <- gamlss::predictAll(object, ...)
  d$y <- NULL
  class(d) <- c("GAMLSS", "distribution")
  return(d)  
}
