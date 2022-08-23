#' Methods for Evaluating the CRPS of Various Distributions
#'
#' Methods to the \code{\link[scoringRules]{crps}} generic function from
#' the \pkg{scoringRules} package for evaluating the (continuous) ranked probability
#' score (CRPS) of many distribution objects from the \pkg{distributions3} package.
#'
#' Methods for most univariate distributions are provided for which \pkg{scoringRules}
#' provides a CRPS function \emph{and} \pkg{distributions3} provides a distribution
#' object. For example, the \code{crps} method for \code{Normal} objects simply
#' leverages \code{\link[scoringRules]{crps_norm}} for objects of class
#' \code{\link[distributions3]{Normal}}.
#'
#' @param y A distribution object, e.g., as created by
#'   \code{\link[distributions3]{Normal}} or \code{\link[distributions3]{Binomial}}.
#' @param x A vector of elements whose CRPS should be determined given the
#'   distribution \code{y}.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise logical. Should each distribution in \code{y} be evaluated
#'   at all elements of \code{x} (\code{elementwise = FALSE}, yielding a matrix)?
#'   Or, if \code{y} and \code{x} have the same length, should the evaluation be
#'   done element by element (\code{elementwise = TRUE}, yielding a vector)? The
#'   default of \code{NULL} means that \code{elementwise = TRUE} is used if the
#'   lengths match and otherwise \code{elementwise = FALSE} is used.
#' @param ... Currently not used.
#'
#' @return In case of a single distribution object, either a numeric
#'   vector of \code{length(x)} (if \code{drop = TRUE}, default) or a matrix with
#'   \code{length(x)} columns (if \code{drop = FALSE}). In case of a vectorized distribution
#'   object, a matrix with \code{length(x)} columns containing all possible combinations.
#'
#' @examples
#' \dontshow{ if(!requireNamespace("distributions3")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("not all packages required for the example are installed")
#'   } else q() }
#' }
#' set.seed(6020)
#'
#' ## three normal distributions X and observations x
#' library("distributions3")
#' X <- Normal(mu = c(0, 1, 2), sigma = c(2, 1, 1))
#' x <- c(0, 0, 1)
#'
#' ## evaluate crps using infrastructure from scoringRules
#' library("scoringRules")
#' crps(X, x)
#' 
#' ## analogously for Poisson distribution
#' Y <- Poisson(c(0.5, 1, 2))
#' crps(Y, x)
#' @exportS3Method scoringRules::crps Beta
crps.Beta <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_beta(y = at, shape1 = d$alpha, shape2 = d$beta)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps Bernoulli
crps.Bernoulli <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_binom(y = at, prob = d$p, size = 1)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps Binomial
crps.Binomial <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_binom(y = at, prob = d$p, size = d$size)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps Erlang
crps.Erlang <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_gamma(y = at, shape = d$k, rate = d$lambda)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps Exponential
crps.Exponential <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_exp(y = at, rate = d$rate)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps Gamma
crps.Gamma <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_gamma(y = at, shape = d$shape, rate = d$rate)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps GEV
crps.GEV <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_gev(y = at, location = d$mu, scale = d$sigma, shape = d$xi)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps Geometric
crps.Geometric <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_gev(y = at, prob = d$p, size = 1)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps Gumbel
crps.Gumbel <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_gev(y = at, location = d$mu, scale = d$sigma, shape = 0)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps HyperGeometric
crps.HyperGeometric <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_hyper(y = at, m = d$m, n = d$n, k = d$k)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps Logistic
crps.Logistic <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_logis(y = at, location = d$location, scale = d$scale)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps LogNormal
crps.LogNormal <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_lnorm(y = at, meanlog = d$log_mu, sdlog = d$log_sigma)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps NegativeBinomial
crps.NegativeBinomial <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- if ("mu" %in% names(unclass(y))) {
    function(at, d) scoringRules::crps_nbinom(y = at, mu = d$mu, size = d$size)
  } else {
    function(at, d) scoringRules::crps_nbinom(y = at, p = d$p, size = d$size)
  }
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps Normal
crps.Normal <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_norm(y = at, mean = d$mu, sd = d$sigma)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps Poisson
crps.Poisson <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_pois(y = at, lambda = d$lambda)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps StudentsT
crps.StudentsT <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_t(y = at, df = d$df)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

#' @rdname crps.Beta
#' @exportS3Method scoringRules::crps Uniform
crps.Uniform <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("scoringRules"))
  FUN <- function(at, d) scoringRules::crps_unif(y = at, min = d$a, max = d$b)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

