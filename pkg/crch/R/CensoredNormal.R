CensoredNormal <- function(mu = 0, sigma = 1, left = -Inf, right = Inf) {
  n <- c(length(mu), length(sigma), length(left), length(right))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))
  d <- data.frame(mu = mu, sigma = sigma, left = left, right = right)
  class(d) <- c("CensoredNormal", "distribution")
  d
}

mean.CensoredNormal <- function(x, ...) {
  m <- ecnorm(mean = x$mu, sd = x$sigma, left = x$left, right = x$right)
  setNames(m, names(x))
}

variance.CensoredNormal <- function(x, ...) {
  s <- sdcnorm(mean = x$mu, sd = x$sigma, left = x$left, right = x$right)
  setNames(s^2, names(x))
}

skewness.CensoredNormal <- function(x, ...) {
  stop("not yet implemented")
}

kurtosis.CensoredNormal <- function(x, ...) {
  stop("not yet implemented")
}

random.CensoredNormal <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rcnorm(n = at, mean = d$mu, sd = d$sigma, left = d$left, right = d$right)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.CensoredNormal <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dcnorm(x = at, mean = d$mu, sd = d$sigma, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.CensoredNormal <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dcnorm(x = at, mean = d$mu, sd = d$sigma, left = d$left, right = d$right, log = TRUE)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.CensoredNormal <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) pcnorm(q = at, mean = d$mu, sd = d$sigma, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.CensoredNormal <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qcnorm(at, mean = d$mu, sd = d$sigma, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

crps.CensoredNormal <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) scoringRules::crps_cnorm(y = at, location = d$mu, scale = d$sigma, lower = d$left, upper = d$right)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

support.CensoredNormal <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(d$left, d$right, d, drop = drop)
}

is_discrete.CensoredNormal <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.CensoredNormal <- function(d, ...) {
  setNames(!is.finite(d$left) & !is.finite(d$right), names(d))
}
