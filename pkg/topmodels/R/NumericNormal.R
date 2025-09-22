
# ------------------------------- NumericNormal -------------------------------

#' Normal Distribution for Testing Numeric Moment Calculations
#'
#' @param mu numeric vector with means.
#' @param var numeric vector with variances.
#'
#' @return Object of class `c("NumericNormal", "distribution")`.
#'
#' @examples
#' \dontrun{
#' n <- Normal(5, sqrt(4))
#' r <- NumericNormal(5, 4)
#' 
#' ## ---------------- PDF ----------------
#' ## 'NumericNormal' only knows the (exact) CDF of a Gaussian distribution
#' ## and approximates the PDF via numeric derivation ([stats::numericDeriv()]).
#' pdf(n, 3.5)
#' pdf(r, 3.5)
#' 
#' ## Visual example
#' xx <- seq(-5, 15, by = 0.1)
#' plot(xx, pdf(n, xx), main = "Testing numeric CDF->PDF")
#' lines(xx, pdf(r, xx), col = 2, lwd = 3)
#' legend("topleft", legend = c("Normal", "NumericNormal"),
#'        pch = c(1, NA), lty = c(NA, 1), lwd = c(NA, 3), col = c(1, 2))
#' 
#' ## -------------- Quantile -------------
#' ## 'NumericNormal' only knows the CDF from which it calculates the quantiles.
#' quantile(n, 0.3)
#' quantile(r, 0.3)
#' 
#' qq <- c(0.0001, seq(0.01, 0.99, by = 0.005), 0.9999)
#' plot(quantile(n, qq), qq, main = "Testing numeric PDF->Quantile")
#' lines(quantile(r, qq), qq, col = 2, lwd = 3)
#' legend("topleft",
#'        legend = c("Normal", "NumericNormal"),
#'        pch = c(1, NA), lty = c(NA, 1), lwd = c(NA, 3), col = c(1, 2))
#' 
#' ## ----------- Central Moments ---------
#' ## Uses 'gridsize = 500L' by default (number of grid size to approximate the distribution).
#' c(mean     = mean(r),
#'   variance = variance(r),
#'   skewness = skewness(r),
#'   kurtosis = kurtosis(r))
#' 
#' ## Central moments with 'gridsize = 100L' (rougher numeric approximation).
#' c(mean     = mean(r, gridsize = 100),
#'   variance = variance(r, gridsize = 100),
#'   skewness = skewness(r, gridsize = 100),
#'   kurtosis = kurtosis(r, gridsize = 100))
#' }
#'
#' @rdname NumericNormal
#' @export
NumericNormal <- function(mu, var) {
    stopifnot(is.numeric(mu), is.numeric(var),
              length(mu) > 0L, length(var) > 0L,
              all(is.finite(mu)), all(is.finite(var)), all(var > 0))
    structure(data.frame(mu = mu, var = var),
              class = c("NumericNormal", "distribution"))
}


#' @rdname NumericNormal
#' @exportS3Method
support.NumericNormal <- function(d, drop = TRUE, ...) {
    min <- rep(-Inf, length(d))
    max <- rep(+Inf, length(d))
    make_support(min, max, d, drop = drop)
}


#' @rdname NumericNormal
#' @exportS3Method
is_discrete.NumericNormal <- function(d, ...) {
    setNames(rep.int(FALSE, length(d)), names(d))
}


#' @rdname NumericNormal
#' @exportS3Method
is_continuous.NumericNormal <- function(d, ...) {
    setNames(rep.int(TRUE, length(d)), names(d))
}


#' @param d object of class `"NumericNormal"`.
#' @param x numeric, where to evaluate the distribution object.
#'
#' @rdname NumericNormal
#' @exportS3Method
cdf.NumericNormal <- function(d, x, drop = TRUE, elementwise = NULL, ...) {

    stopifnot(is.numeric(x) && all(is.finite(x)))
    drop <- as.logical(drop)[[1L]]
    stopifnot(isTRUE(drop) || isFALSE(drop))
    stopifnot(is.null(drop) || isTRUE(drop) || isFALSE(drop))

    # Check if calculation is performed elementwise or not
    k <- length(x); n <- length(d)
    if (is.null(elementwise))
        elementwise <- (k == 1L & n == 1L) || (k > 1L && k == n)
    if (elementwise && k > 1L && k != n)
        stop(sprintf("lengths of distributions and arguments do not match: %s != %s", n, k))

    if (elementwise) {
        res <- sapply(seq_along(d), function(i) {
                          pnorm(x[i], mean = d[i]$mu, sd = sqrt(d[i]$var))
               })
    } else {
        res <- sapply(seq_along(d), function(i) {
                          pnorm(x, mean = d[i]$mu, sd = sqrt(d[i]$var))
               })
        if (is.matrix(res)) {
            res <- t(res)
            colnames(res) <- paste0("p_", format(x))
        }
        if (drop && any(dim(res) == 1L))
            res <- as.numeric(res)
    }
    return(res)
}

## for testing ## # --- this was only for testing the moments calculation algorithm ---
## for testing ## #' @rdname NumericNormal
## for testing ## #' @exportS3Method
## for testing ## quantile.NumericNormal <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
## for testing ## 
## for testing ##     stopifnot(is.numeric(x) && all(is.finite(x)))
## for testing ##     drop <- as.logical(drop)[[1L]]
## for testing ##     stopifnot(isTRUE(drop) || isFALSE(drop))
## for testing ##     stopifnot(is.null(drop) || isTRUE(drop) || isFALSE(drop))
## for testing ## 
## for testing ##     # Check if calculation is performed elementwise or not
## for testing ##     k <- length(x); n <- length(d)
## for testing ##     if (is.null(elementwise))
## for testing ##         elementwise <- (k == 1L & n == 1L) || (k > 1L && k == n)
## for testing ##     if (elementwise && k > 1L && k != n)
## for testing ##         stop(sprintf("lengths of distributions and arguments do not match: %s != %s", n, k))
## for testing ## 
## for testing ##     if (elementwise) {
## for testing ##         res <- sapply(seq_along(d), function(i) {
## for testing ##                           qnorm(x[i], mean = d[i]$mu, sd = sqrt(d[i]$var))
## for testing ##                })
## for testing ##     } else {
## for testing ##         res <- sapply(seq_along(d), function(i) {
## for testing ##                           qnorm(x, mean = d[i]$mu, sd = sqrt(d[i]$var))
## for testing ##                })
## for testing ##         if (is.matrix(res)) {
## for testing ##             res <- t(res)
## for testing ##             colnames(res) <- paste0("p_", format(x))
## for testing ##         }
## for testing ##         if (drop && any(dim(res) == 1L))
## for testing ##             res <- as.numeric(res)
## for testing ##     }
## for testing ##     return(res)
## for testing ## }
# ----------------------------- End NumericNormal --------------------------------

# ------------------------------- NumericPoisson -------------------------------
#' @param lambda vector of (non-negative) means.
#' @rdname NumericNormal
#' @export
NumericPoisson <- function(lambda) {
    stopifnot(is.numeric(lambda), length(lambda) > 0L,
              all(is.finite(lambda)), all(lambda > 0))
    structure(data.frame(lambda = lambda),
              class = c("NumericPoisson", "distribution"))
}


#' @rdname NumericNormal
#' @exportS3Method
support.NumericPoisson <- function(d, drop = TRUE, ...) {
    min <- rep(0L, length(d))
    max <- rep(+Inf, length(d))
    make_support(min, max, d, drop = drop)
}


#' @rdname NumericNormal
#' @exportS3Method
is_discrete.NumericPoisson <- function(d, ...) {
    setNames(rep.int(TRUE, length(d)), names(d))
}


#' @rdname NumericNormal
#' @exportS3Method
is_continuous.NumericPoisson <- function(d, ...) {
    setNames(rep.int(FALSE, length(d)), names(d))
}


#' @param d object of class `"NumericPoisson"`.
#' @param x numeric, where to evaluate the distribution object.
#'
#' @rdname NumericNormal
#' @exportS3Method
cdf.NumericPoisson <- function(d, x, drop = TRUE, elementwise = NULL, ...) {

    stopifnot(is.numeric(x) && all(is.finite(x)))
    drop <- as.logical(drop)[[1L]]
    stopifnot(isTRUE(drop) || isFALSE(drop))
    stopifnot(is.null(drop) || isTRUE(drop) || isFALSE(drop))

    # Check if calculation is performed elementwise or not
    k <- length(x); n <- length(d)
    if (is.null(elementwise))
        elementwise <- (k == 1L & n == 1L) || (k > 1L && k == n)
    if (elementwise && k > 1L && k != n)
        stop(sprintf("lengths of distributions and arguments do not match: %s != %s", n, k))

    if (elementwise) {
        res <- sapply(seq_along(d), function(i) ppois(x[i], lambda = d[i]$lambda))
    } else {
        res <- sapply(seq_along(d), function(i) ppois(x, lambda = d[i]$lambda))
        if (is.matrix(res)) {
            res <- t(res)
            colnames(res) <- paste0("p_", format(x))
        }
        if (drop && any(dim(res) == 1L))
            res <- as.numeric(res)
    }
    return(res)
}

# ----------------------------- End NumericPoisson -------------------------------

#' @rdname NumericNormal
#' @exportS3Method
pdf.distribution <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
    if (!hasS3method("cdf", class(d)))
        stop("no S3 method 'cdf' found for object of class: ", paste(class(d), collapse = ", "))

    stopifnot(is.numeric(x) && all(is.finite(x)))
    x <- as.numeric(x) # Required for numericDeriv
    drop <- as.logical(drop)[[1L]]
    stopifnot(isTRUE(drop) || isFALSE(drop))
    stopifnot(is.null(drop) || isTRUE(drop) || isFALSE(drop))

    ## Check if calculation is performed elementwise or not
    k <- length(x); n <- length(d)
    if (is.null(elementwise))
        elementwise <- (k == 1L & n == 1L) || (k > 1L && k == n)
    if (elementwise && k > 1L && k != n)
        stop(sprintf("lengths of distributions and arguments do not match: %s != %s", n, k))

    ## Setting up results matrix of dimension 'n x k'
    res <- matrix(NA, nrow = n, ncol = k, dimnames = list(NULL, paste0("d_", format(x))))

    ## Discrete distribution?
    discrete <- all(distributions3::is_discrete(d))

    ## For discrete distributions: Take difference F(floor(x_i)) - F(floor(x_i) - 1) and
    ## and F(x_i) for x_i == 0.
    if (discrete) {

        ## Elementwise: one 'x[i]' per distribution 'x[i]' or
        ## same probability 'x[1L]' for a set of distributions 'x[i]'.
        if (elementwise || k == 1L) {
            x <- rep(floor(x), length(d))
            is_zero <- x == 0
            if (any(is_zero)) {
                res[is_zero, 1L]  <- cdf(d[is_zero], 0, drop = TRUE, elementwise = TRUE)
            }
            if (any(!is_zero)) {
                res[!is_zero, 1L] <- cdf(d[!is_zero], x[!is_zero],     drop = TRUE, elementwise = TRUE) -
                                     cdf(d[!is_zero], x[!is_zero] - 1, drop = TRUE, elementwise = TRUE)
            }
        } else {
            x <- floor(x)
            for (i in seq_len(k)) {
                tmpx <- rep(x[i], length(d))
                if (x[i] == 0) {
                    res[, i] <- cdf(d, tmpx, drop = TRUE, elementwise = TRUE)
                } else {
                    res[, i] <- cdf(d, tmpx,     drop = TRUE, elementwise = TRUE) -
                                cdf(d, tmpx - 1, drop = TRUE, elementwise = TRUE)
                }
            }
        }

    ## For non-discrete (continuous) distributions, calculate numeric derivative
    } else {

        ## Elementwise: one 'x[i]' per distribution 'x[i]' or
        ## same probability 'x[1L]' for a set of distributions 'x[i]'.
        if (elementwise || k == 1L) {
            # Numeric approximation of the CDF
            myenv   <- new.env()
            myenv$d <- d; myenv$x <- if (k == 1L) x[[1L]] else x[i]
            nd <- tryCatch(numericDeriv(quote(cdf(d, x)), "x", myenv),
                           error = function(e) stop("problem calculating numeric derivative (cdf->pdf)"))
            res[, 1L] <- attr(nd, "gradient")[, 1L]
        } else {
            for (i in seq_len(k)) {
                myenv   <- new.env()
                myenv$d <- d; myenv$x <- x[i]
                numericDeriv(quote(cdf(d, x)), "x", myenv)
                nd <- tryCatch(numericDeriv(quote(cdf(d, x)), "x", myenv),
                               error = function(e) stop("problem calculating numeric derivative (cdf->pdf): ", e))
                res[, i] <- attr(nd, "gradient")[, 1L]
            }
        }
    }

    ## Drop dimensions if requested/possible
    if (drop & any(dim(res) == 1L)) dim(res) <- NULL
    return(res)
}


#' @param probs numeric. Vector of probabilities.
#' @param drop logical. Should the result be simplified to a vector if possible?
#' @param elementwise `NULL` (default; auto-detect) or logical. Should the
#'        quantiles defined by `probs` be evaluated at each distribution (`x`)? 
#' @param lower,upper numeric. Lower and upper end points for the interval to
#'        be searched, forwarded to [stats::uniroot()].
#' @param tol numeric. Desired accuracy for [stats::uniroot()].
#' @param ... currently ignored.
#'
#' @rdname NumericNormal
#' @exportS3Method
quantile.distribution <- function(x, probs, drop = TRUE, elementwise = NULL,
                                  lower = -1 / sqrt(.Machine$double.eps),
                                  upper = +1 / sqrt(.Machine$double.eps),
                                  tol   = .Machine$double.eps^0.5, ...) {

    if (!hasS3method("cdf", class(x)))
        stop("no S3 method 'cdf' found for object of class: ", paste(class(x), collapse = ", "))

    stopifnot(is.numeric(lower), length(lower) == 1L)
    stopifnot(is.numeric(upper), length(upper) == 1L)
    stopifnot(lower < upper)
    stopifnot(is.numeric(tol), length(tol) == 1L, tol >= .Machine$double.eps, tol < 0.01)

    # Check if calculation is performed elementwise or not
    k <- length(probs); n <- length(x)
    if (is.null(elementwise))
        elementwise <- (k == 1L & n == 1L) || (k > 1L && k == n)
    if (elementwise && k > 1L && k != n)
        stop(sprintf("lengths of distributions and arguments do not match: %s != %s", n, k))

    ## Numeric approximation of quantiles
    inverse_cdf <- function(d, p, lower, upper, tol, ...) {
        uniroot(function(y, d, ...) cdf(d, y, ...) - p, d = d, lower = lower, upper = upper, tol = tol, ...)$root
    }

    # Evaluate support; define lower/upper bounds for uniroot
    sup <- support(d = x, drop = FALSE)
    sup[sup[, "min"] < lower, "min"] <- lower
    sup[sup[, "max"] > upper, "max"] <- upper

    ## Setting up results matrix of dimension 'n x k' or 'n x 1' (elementwise)
    if (elementwise) {
        res <- matrix(NA, nrow = n, ncol = 1L)
    } else {
        res <- matrix(NA, nrow = n, ncol = k, dimnames = list(NULL, paste0("q_", format(probs))))
    }

    for (i in seq_len(n)) {
        ## Elementwise: one probability 'probs[i]' per distribution 'x[i]' or
        ## same probability 'probs[1L]' for a set of distributions 'x[i]'.
        if (elementwise || k == 1L) {
            res[i, 1L] <- inverse_cdf(d = x[i], if (k == 1) probs[1L] else probs[i],
                                      lower = sup[[i, "min"]], upper = sup[[i, "max"]], tol = tol)
        ## Calculate quantiles for each probability 'probs[j]' for each distribution 'x[i]';
        ## Scoping variable 'tol'.
        } else {
            res[i, ] <- sapply(probs, function(p, d, l, u) inverse_cdf(d = d, p = p, lower = l, upper = u, tol = tol),
                               d = x[i], l = sup[[i, "min"]], u = sup[[i, "max"]])
        }
    }

    ## Discrete distribution? Round result
    if (all(distributions3::is_discrete(x))) res <- round(res)

    ## Drop dimensions if requested/possible
    if (drop & any(dim(res) == 1L)) dim(res) <- NULL
    return(res)
}


#' Method for Numerically Evaluate Central Moments of Probability Distributions
#'
#' Method used to evaluate (approximate) the central moments (mean, variance, skewness, and kurtosis)
#' for probability distributions for which only the cummulative distribution function (CDF) and -
#' potentially - the quantile function is provided.
#'
#' For discrete distributions spanning a range less than `gridsize` the PDF is calculated
#' at $i = \{0, 1, 2, 3, \dots\}$ by differenciating the CDF provided which is then used
#' to calculate the central moments.
#'
#' For continuous distributions as well as discrete distributions spanning a
#' wide range of values (larger than `gridsize`) a grid with `gridsize`
#' intervals is created. Given the distribution provides a quantile function,
#' this grid is specified on a (mostly) uniform grid on the quantile scale. If
#' no quantile function is provided, the $0.01$ and $99.99$ percentile are
#' calculated approximated via the CDF, between which a uniform grid is
#' spanned. For each interval the density is approximated using numeric
#' forward differences
#'
#' * $f(x[j]) = (F(x[i+1]) - F(x_i)) / (x[i+1] - x[i])$
#'
#' at each $x[j] = (x[i+1] + x[i]) * 0.5$ with interval width $x_[i+1] - x[i]$.
#' The densities $f(x[j])$ and interval mids $x_j$ are used to calculate
#' weighted moments using the trapezoidal rule.
#'
#' @param x  Object of class 'distribution'
#' @param what single integer, controls what the C code returns. 1L (mean) 2L
#'        (variance) 3L (skewness) 4L (kurtosis).
#' @param gridsize positive size of the grid used to approximate the CDF for
#'   the numerical calculation of the CRPS.
#' @param batchsize maximum batch size. Used to split the input into batches.
#'   Lower values reduce required memory but may increase computation time.
#' @param applyfun an optional \code{\link[base]{lapply}}-style function with arguments
#'   \code{function(X, FUN, \dots)}. It is used to compute the CRPS for each element
#'   of \code{y}. The default is to use the basic \code{lapply}
#'   function unless the \code{cores} argument is specified (see below).
#' @param cores numeric. If set to an integer the \code{applyfun} is set to····
#'   \code{\link[parallel]{mclapply}} with the desired number of \code{cores},
#'   except on Windows where \code{\link[parallel]{parLapply}} with
#'   \code{makeCluster(cores)} is used.
#' @param method character. Should the grid be set up on the observation scale
#'   and \code{method = "cdf"} be used to compute the corresponding probabilities?
#'   Or should the grid be set up on the probability scale and \code{method = "quantile"}
#'   be used to compute the corresponding observations? By default, \code{"cdf"}
#'   is used for discrete observations whose range is smaller than the \code{gridsize}
#'   and \code{"quantile"} otherwise.
#' @param ... [mean.distribution()], [variance.distribution()], [skewness.distribution()] and
#'   [kurtosis.distribution()] forward the additional arguments (`...`) to
#'   [distribution_calculate_moments()]; all other functions/methods ignore additional arguments.
#'
#' @return A (potentially named) numeric vector of length `length(x)` with the requested central moment.
#'
#' @rdname NumericNormal
distribution_calculate_moments <- function(x, what, gridsize = 500L, batchsize = 1e4L, applyfun = NULL, cores = NULL, method = NULL, ...) {

  ## sanity checks
  stopifnot(inherits(x, "distribution"))
  stopifnot(is.integer(what), length(what) == 1L, what >= 1L && what <= 4L)
  stopifnot(is.numeric(gridsize), length(gridsize) == 1L, gridsize >= 2L)
  gridsize <- as.integer(gridsize)
  stopifnot(is.numeric(batchsize), length(batchsize) == 1L, batchsize >= 1L)
  stopifnot(is.null(cores) || is.numeric(cores))
  stopifnot(is.null(applyfun) || is.function(applyfun))
  if (is.numeric(cores)) {
    cores <- as.integer(cores)
    stopifnot(length(cores) == 1L, cores >= 1L)
  }

  ## basic properties:
  ## n = number of distributions
  n <- length(x)

  ## handle zero-length distribution vector
  if (n == 0L) return(numeric(0L))

  ## Names of 'x' (if any) to be added later
  xnam <- names(x)

  ## define apply functions for parallelization
  if(is.null(applyfun)) {
    applyfun <- if(is.null(cores)) {
      lapply
    } else {
      if(.Platform$OS.type == "windows") {
        cl_cores <- parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl_cores))
        function(X, FUN, ...) parallel::parLapply(cl = cl_cores, X, FUN, ...)
      } else {
        function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = cores)
      }
    }
  } 

  ## calculate batches
  ## batch_n:   number of batches required
  ## batch_id:  integer vector of same length as x
  batch_n  <- (function(cores, N, b) {
                  rval <- max(ifelse(is.null(cores), 0, cores), N %/% b + as.integer(N %% b != 0))
                  if (!is.null(cores) && rval %% cores != 0) rval <- (rval %/% cores + 1) * cores
                  return(rval)
                })(cores, length(x), batchsize)
  batch_id <- if (batch_n == 1) rep(1L, length(x)) else rep(seq_len(batch_n), each = ceiling(length(x) / batch_n))[seq_along(x)]

  ## selecting default method
  xrange <- range(support(x), na.rm = TRUE)
  if (!is.finite(xrange[1L])) xrange[1L] <- min(quantile(x, 0.0001), na.rm = TRUE)
  if (!is.finite(xrange[2L])) xrange[2L] <- max(quantile(x, 0.9999), na.rm = TRUE)
  discrete <- all(distributions3::is_discrete(x))

  ## decide on method to use. If discrete and the range is <= gridsize
  ## we use the discrete CDF.
  ## Else we span up a grid across the range where we can calculate the
  ## density (numerical differences of the CDF). If a quantile function is
  ## available the quantile function is used to span a nice grid, else we
  ## use a univariate grid (method = "cdf") as we expect the approximated
  ## inverse cdf (i.e., the quantile function) will be inefficient/slow.
  if (is.null(method) && discrete) {
    method <- if (diff(xrange) <= gridsize) "cdf" else "quantile"
  } else if (is.null(method)) {
    method <- if (!hasS3method("quantile", class(x)[!class(x) == "distribution"])) "cdf" else "quantile"
  }


  ## Discrete and 'cdf': Calculate true PDF for all counts from 0 to max(xrange)
  if (discrete && method == "cdf") {
    q <- seq(min(xrange), max(max(xrange) * 1.5, gridsize), by = 1.0) # Must be numeric!

    ## Scoping `batch_id`, `x`
    batch_fn <- function(i) {
        ## Index of `x` falling into current batch `i`
        idx  <- which(batch_id == i)

        ## Calculate CDF, convert to PDF
        d <- sapply(idx, function(i) cdf(x[i], q, drop = TRUE, elementwise = FALSE)) # <- transposed
        d <- cbind(d[1, ], t(apply(d, MARGIN = 2, diff)))

        ## Mean
        mean <- sapply(seq_along(idx), function(i) weighted.mean(q, d[idx[i], ]))
        if (what == 1L) return(mean)

        ## Variance
        variance <- sapply(seq_along(idx), function(i) weighted.mean((q - mean[i])^2, d[idx[i], ]))
        if (what == 2L) return(variance)

        ## Skewness or kurtosis
        pow <- if (what == 3) 3.0 else 4.0
        offset <- if (what == 3) 0.0 else 3.0 # To get excess kurtosis
        return(sapply(seq_along(idx), function(i) weighted.mean((q - mean[i])^pow, d[idx[i], ]) / variance[i]^(pow/2)) - offset)
    }
    rval <- do.call(c, applyfun(seq_len(batch_n), batch_fn))
  }

  ## Discrete and 'quantile': Approximation oft he discrete distribution using a grid
  cat("  ---- is discrete? ", discrete, "\n")
  cat("       method selected to be \"", method, "\"\n", sep = "")
  if (discrete && method == "quantile") {
      stop("This must be implemented (discrete + method = 'quantile')")
  }

  ## cdf method: set up observations and compute probabilities via cdf()
  if (!discrete && method == "cdf") {
    ## Drawing one set of quantiles; calculate probabilities for all distributions
    q <- if (discrete && (diff(xrange) <= gridsize)) {
      seq(xrange[1L], xrange[2L] + 1, by = 1.0)
    } else {
      seq(xrange[1L], xrange[2L], length.out = gridsize)
    }
    qdelta <- diff(q)         # Interval width
    qmid   <- q[-1] - diff(q) # Interval mid

    ## Scoping `batch_id`, `x`
    batch_fn <- function(i) {
        ## Index of `x` falling into current batch `i`
        idx  <- which(batch_id == i)

        ## Getting xlimits
        xlims <- quantile(x[idx], c(0.0001, 0.9999), drop = FALSE, elementwise = FALSE)

        ## Calculate (uniform) grid; note that the matrix is transposed to
        ## avoid flipping it back and force unnecessarily.
        grd       <- apply(xlims, MARGIN = 1, function(z) seq(z[1L], z[2L], length.out = gridsize + 1))
        ## Grid widths and mid points
        grd_width <- grd[2L, ] - grd[1L, ]
        grd_mids  <- t(apply(grd, MARGIN = 2, function(z) (z[-1] + z[-length(z)]) / 2))

        ## Approximate density
        d <- t(sapply(idx, function(i) diff(cdf(x[i], grd[, i], drop = TRUE, elementwise = FALSE)) / grd_width[i]))

        ## - q: grid points (quantiles)
        ## - d: density at these grid points
        ## Based on these two matrices the C code calculates a weighted mean
        .Call("c_d2moments_numeric",
              at       = as.numeric(grd_mids),
              pat      = as.numeric(d),
              dim      = dim(grd_mids),       # dimension of matrices 'at/pat'.
              discrete = FALSE,
              what     = what,                # mean (1L), variance (2L), skewness (3L), kurtosis (4L)
              PACKAGE  = "topmodels")
    }
    rval <- do.call(c, applyfun(seq_len(batch_n), batch_fn))
  }

  ## quantile method: set up probabilities and compute observations via quantile()
  if (!discrete && method == "quantile") {
    ## Drawing one set of probabilities; calculate quantiles for all distributions
    p <- c(0.001, 0.01, 0.1, 1L:(gridsize - 1L), gridsize - c(0.1, 0.01, 0.001)) / gridsize

    ## Scoping `batch_id`, `x`, `p`
    batch_fn <- function(i) {
        ## Index of `x`/`y` falling into current batch `i`
        idx  <- which(batch_id == i)


        ## Calculate quantile-based grid; note that the matrix is transposed to
        ## avoid flipping it back and force unnecessarily.
        grd <- sapply(idx, function(i) quantile(x[i], p, drop = FALSE, elementwise = FALSE))
        ## Grid widths and mid points
        grd_width <- apply(grd, MARGIN = 2, function(z) diff(z))
        grd_mids  <- t(apply(grd, MARGIN = 2, function(z) (z[-1] + z[-length(z)]) / 2))

        ## Approximate density
        d <- t(sapply(idx, function(i) diff(cdf(x[i], grd[, i], drop = TRUE, elementwise = FALSE)) / grd_width[, i]))

        ## - q: grid points (quantiles)
        ## - d: density at these grid points
        ## Based on these two matrices the C code calculates a weighted mean
        .Call("c_d2moments_numeric",
              at       = as.numeric(grd_mids),
              pat      = as.numeric(d),
              dim      = dim(grd_mids),       # dimension of matrices 'at/pat'.
              discrete = FALSE,
              what     = what,                # mean (1L), variance (2L), skewness (3L), kurtosis (4L)
              PACKAGE  = "topmodels")
    }
    rval <- do.call(c, applyfun(seq_len(batch_n), batch_fn))
  }

  ## handle dimensions
  return(setNames(rval, xnam))
}

#' @param x object of class `c("NumericNormal", "distribution")`.
#' @param gridsize integer, number of grid points used for approximation. Defaults to `500L`.
#'
#' @rdname NumericNormal
#' @exportS3Method
mean.distribution <- function(x, ...) {
    distribution_calculate_moments(x = x, what = 1L, ...)
}


#' @rdname NumericNormal
#' @exportS3Method
variance.distribution <- function(x, ...) {
    distribution_calculate_moments(x = x, what = 2L, ...)
}


#' @rdname NumericNormal
#' @exportS3Method
skewness.distribution <- function(x, ...) {
    distribution_calculate_moments(x = x, what = 3L, ...)
}


#' @rdname NumericNormal
#' @exportS3Method
kurtosis.distribution <- function(x, ...) {
    distribution_calculate_moments(x = x, what = 4L, ...)
}


