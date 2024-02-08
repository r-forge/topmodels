TruncatedNormal <- function(mu = 0, sigma = 1, left = -Inf, right = Inf) {
  n <- c(length(mu), length(sigma), length(left), length(right))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))
  d <- data.frame(mu = mu, sigma = sigma, left = left, right = right)
  class(d) <- c("TruncatedNormal", "distribution")
  d
}

mean.TruncatedNormal <- function(x, ...) {
  m <- etnorm(mean = x$mu, sd = x$sigma, left = x$left, right = x$right)
  setNames(m, names(x))
}

variance.TruncatedNormal <- function(x, ...) {
  s <- sdtnorm(mean = x$mu, sd = x$sigma, left = x$left, right = x$right)
  setNames(s^2, names(x))
}

skewness.TruncatedNormal <- function(x, ...) {
  stop("not yet implemented")
}

kurtosis.TruncatedNormal <- function(x, ...) {
  stop("not yet implemented")
}

random.TruncatedNormal <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rtnorm(n = at, mean = d$mu, sd = d$sigma, left = d$left, right = d$right)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.TruncatedNormal <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dtnorm(x = at, mean = d$mu, sd = d$sigma, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.TruncatedNormal <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dtnorm(x = at, mean = d$mu, sd = d$sigma, left = d$left, right = d$right, log = TRUE)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.TruncatedNormal <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) ptnorm(q = at, mean = d$mu, sd = d$sigma, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.TruncatedNormal <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qtnorm(at, mean = d$mu, sd = d$sigma, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

crps.TruncatedNormal <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) scoringRules::crps_tnorm(y = at, location = d$mu, scale = d$sigma, lower = d$left, upper = d$right)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

support.TruncatedNormal <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(d$left, d$right, d, drop = drop)
}

is_discrete.TruncatedNormal <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.TruncatedNormal <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}
