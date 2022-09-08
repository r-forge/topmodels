TruncatedLogistic <- function(location = 0, scale = 1, left = -Inf, right = Inf) {
  n <- c(length(location), length(scale), length(left), length(right))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))
  d <- data.frame(location = location, scale = scale, left = left, right = right)
  class(d) <- c("TruncatedLogistic", "distribution")
  d
}

mean.TruncatedLogistic <- function(x, ...) {
  m <- etlogis(location = x$location, scale = x$scale, left = x$left, right = x$right)
  setNames(m, names(x))
}

variance.TruncatedLogistic <- function(x, ...) {
  stop("not yet implemented")
}

skewness.TruncatedLogistic <- function(x, ...) {
  stop("not yet implemented")
}

kurtosis.TruncatedLogistic <- function(x, ...) {
  stop("not yet implemented")
}

random.TruncatedLogistic <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rtlogis(n = at, location = d$location, scale = d$scale, left = d$left, right = d$right)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.TruncatedLogistic <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dtlogis(x = at, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.TruncatedLogistic <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dtlogis(x = at, location = d$location, scale = d$scale, left = d$left, right = d$right, log = TRUE)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.TruncatedLogistic <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) ptlogis(q = at, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.TruncatedLogistic <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qtlogis(at, location = d$location, scale = d$scale, left = d$left, right = d$right, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

crps.TruncatedLogistic <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) scoringRules::crps_tlogis(y = at, location = d$location, scale = d$scale, lower = d$left, upper = d$right)
  distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
}

support.TruncatedLogistic <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(d$left, d$right, d, drop = drop)
}

is_discrete.TruncatedLogistic <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.TruncatedLogistic <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}
