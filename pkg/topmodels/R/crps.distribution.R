

#' @param m positive size of the grid used to approximate the CDF for
#'        the numerical calculation of the CRPS.
#' @param batchsize maximum batch size. Used to split the input into batches.
#'        Lower values reduce required memory but may increase computation time.
#' @param applyfun \code{NULL} or a suitable apply function. If \code{NULL},
#'        \code{\link[base]{lapply}} is used. If set, parallelization is not possible anymore.
#' @param cores \code{NULL} or positive integer, number of cores for parallelization.
#'
#' @details For univariate distributions where \pkg{distributions3} provides
#' a distribution object but \pkg{scoringRules} does not offer a CRPS function,
#' a numeric approximation is used. \code{\link[topmodels]{crps.distribution}};
#' see section 'Numeric approximation' for more information.
#'
#' The argument \code{cores} can be used to enable parallelization. When used, \code{y}
#' is split into \code{B} equidistant batches; at least \code{B = cores} batches or
#' a multiple of \code{cores} with a maximum size of \code{batchsize}. On systems running
#' MS Windows \code{parlapply} is used, else \code{mclapply} if \code{applyfun = NULL}.
#' Note that parallelization is ignored if a custom \code{applyfun} has been specified.
#'
#' @section Numeric approximation:
#' The numerical approximation/numerical integration to calculate the CRPS
#' supports continuous, discrete, and mixed distributions.
#'
#' For continuous and mixed distributions, an equidistant grid of \code{m + 5}
#' probabilities is drawn for which the corresponding quantiles for each
#' distribution \code{y} is calculated (includes the observation itself). The
#' calculation of the CRPS is done using trapezodial approximation for the
#' numeric integration.  For discrete distributions \code{m} quantiles are
#' drawn and the corresponding probabilities calculated for each distribution
#' \code{y} (includes the observation itself) and the CRPS calculated using
#' numeric integration.  If a grid of size \code{m} is not sufficient to cover
#' the required range, the method falls back to the procedure used for
#' continuous distributions to approximate the CRPS.
#'
#' The numeric approximation requires to set up a matrix of dimension
#' \code{length(y) * (m + 5)} (or \code{length(y) * (m + 1)}) which is very
#' memory intensive. Thus, the data is split into equidistant batches with a
#' maximal size of \code{batchsize}, limiting the size of the matrix to a
#' maximum of \code{batchsize * (m + 5)}. A smaller value of \code{batchsize}
#' will reduce memory footprint but will slightly increase computation time.
#'
#' The error (deviation between numerical approximation and analytic solution)
#' has been shown to be in the order of \code{1e-2} for a series of distributions
#' tested. Accuracy can be increased by increasing \code{m} and will be lower
#' for a smaller \code{m}.
#'
#' @useDynLib topmodels, .registration = TRUE
#' @exportS3Method scoringRules::crps distributon
#' @rdname crps.Beta
crps.distribution <- function(y, x, drop = TRUE, elementwise = NULL, m = 500, batchsize = 1e4, applyfun = NULL, cores = NULL, ...) {
  ## essentially follow apply_dpqr() but try to exploit specific structure of CRPS

  ## Checking sanity
  stopifnot(inherits(y, "distribution"), is.numeric(x))
  stopifnot(is.null(drop) || isTRUE(drop) || isFALSE(drop))
  stopifnot(is.null(elementwise) || isTRUE(elementwise) || isFALSE(elementwise))
  stopifnot(is.numeric(m), length(m) == 1, m >= 2); m <- as.integer(m)
  stopifnot(is.numeric(batchsize), length(batchsize) == 1, batchsize >= 1)

  ## Sanity check for optional arguments for parallelization
  stopifnot(is.null(cores) || is.numeric(cores))
  stopifnot(is.null(applyfun) || is.function(applyfun))
  if (is.numeric(cores)) { cores <- as.integer(cores); stopifnot(length(cores) == 1, cores >= 1) }

  ## basic properties:
  ## rows n = number of distributions
  ## columns k = number of arguments x
  rnam <- names(y)
  n <- length(y)
  k <- length(x)

  ## determine the dimension of the return value:
  ## * elementwise = FALSE: n x k matrix,
  ##   corresponding to all combinations of 'd' and 'x'
  ## * elementwise = TRUE: n vector,
  ##   corresponding to combinations of each element in 'd' with only the corresponding element in 'x'
  ##   only possible if n = k
  ## * elementwise = NULL: guess the type (default),
  ##   only use TRUE if n = k > 1, and FALSE otherwise
  if(is.null(elementwise)) elementwise <- k > 1L && k == n && is.null(dim(x))
  if(elementwise && k > 1L && k != n) stop(
    sprintf("lengths of distributions and arguments do not match: %s != %s", n, k))

  ## "x" names (if not dropped)
  anam <- if ((k == 1L || n == 1L) && drop) {
    NULL
  } else {
    distributions3:::make_suffix(x, digits = pmax(3L, getOption("digits") - 3L))
  }

  ## handle different types of "x"
  if (k == 0L) {
    return(matrix(numeric(0L), nrow = n, ncol = 0L, dimnames = list(rnam, NULL)))
  } else if (k == 1L) {
    x <- rep.int(as.vector(x), n)
  } else if (elementwise) {
    k <- 1L
  } else {
    x <- as.vector(x)
    k <- length(x)
  }

  ## columns names (if not dropped)
  cnam <- if ((k == 1L || n == 1L) && drop) {
    NULL
  } else if (length(anam) > k) {
    "crps"
  } else {
    paste("crps", anam, sep = "_")
  }

  ## handle zero-length distribution vector
  if (n == 0L) return(matrix(numeric(0L), nrow = 0L, ncol = k, dimnames = list(NULL, cnam)))

  ## Define apply functions for parallelization; code by Achim Zeileis
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

  ## Calculate batches
  ## batch_n:   number of batches required
  ## batch_id:  integer vector of same length as y
  batch_n  <- (function(cores, N, b) {
                  rval <- max(ifelse(is.null(cores), 0, cores), N %/% b + as.integer(N %% b != 0))
                  if (!is.null(cores) && rval %% cores != 0) rval <- (rval %/% cores + 1) * cores
                  return(rval)
                })(cores, length(y), batchsize)
  batch_id <- if (batch_n == 1) rep(1L, length(y)) else rep(seq_len(batch_n), each = ceiling(length(y) / batch_n))[seq_along(y)]

  ## Calling C function for numeric CRPS calculation
  discrete <- all(distributions3::is_discrete(y))
  ## If all distributions are discrete: Calculate required grid size.
  ## If this exceeds `m` continuous approximation will be used.
  ## TODO(R): Split data; use approximation only where grid exceeds m but use
  ##          the discrete (more accurate) approximation for the others?
  if (discrete) {
    xrange <- pmin(range(quantile(y, c(0.001, 0.999)), na.rm = TRUE), range(support(y)))
    ## If grid size larger m; set `discrete = FALSE` and continue; fallback to discrete case
    if (diff(xrange) > m) {
        discrete <- FALSE ## Falling back to continuous mode
        warning("grid size too large; falling back to continuous CRPS approximation")
    } else {
      ## Drawing one set of quantiles; calculate probabilities for all distributions
      q <- seq(xrange[1], xrange[2] + 1, by = 1.0)
      ## Not elementwise: Same observation(s) `x` for all distributions `y`
      if (isFALSE(elementwise)) {
          # Scoping `y`, `q`, `batch_n`, `applyfun`
          # Inside; calling the batch function used for possible parallelization/batching
          fn <- function(x) {
              # Scoping `y`, `q`, `batch_id`
              batch_fn <- function(i, z) {
                  idx  <- which(batch_id == i) ## Index of `x`/`y` falling into current batch `i`
                  p    <- cdf(y[idx], q)       ## Calculating quantiles at `q` for `y[idx]`
                  px   <- cdf(y[idx], z)       ## Probability at `x[idx]`
                  .Call("c_CRPS_numeric", as.numeric(z), px, p, q, FALSE, PACKAGE = "topmodels")
              }
              return(do.call(c, applyfun(seq_len(batch_n), batch_fn, z = x)))
          }
          # Using `lapply()` to iterate over observations `x` by calling `fn(x)`
          rval <- do.call(c, lapply(x, fn))
      ## Multiple distributions `y`, one (or multiple) observations `x` evaluated for each `y`
      } else {
          ## Scoping `batch_id`, `y`, `x`, `q`
          batch_fn <- function(i) {
              idx  <- which(batch_id == i) ## Index of `x`/`y` falling into current batch `i`
              p    <- cdf(y[idx], q)       ## Calculating quantiles at `p` for `y[idx]`
              px   <- cdf(y[idx], x[idx])  ## Probability at `x[idx]`
              .Call("c_CRPS_numeric", as.numeric(x[idx]), px, p, q, FALSE, PACKAGE = "topmodels")
          }
          rval <- do.call(c, applyfun(seq_len(batch_n), batch_fn))
      }
    }
  }
  ## Continuous mode
  if (!discrete) {
    ## Drawing one set of probabilities; calculate quantiles for all distributions
    p <- c(0.001, 0.01, 0.1, 1:(m-1), m - c(0.1, 0.01, 0.001)) / m
    ## Not elementwise: Same observation(s) `x` for all distributions `y`
    if (isFALSE(elementwise)) {
        # Scoping `y`, `p`, `batch_n`, `applyfun`
        # Inside; calling the batch function used for possible parallelization/batching
        fn <- function(x) {
            # Scoping `y`, `p`, `batch_id`
            batch_fn <- function(i, z) {
                idx  <- which(batch_id == i) ## Index of `x`/`y` falling into current batch `i`
                q    <- quantile(y[idx], p)  ## Calculating quantiles at `p` for `y[idx]`
                px   <- cdf(y[idx], z)  ## Probability at `x[idx]`
                .Call("c_CRPS_numeric", as.numeric(z), px, p, q, TRUE, PACKAGE = "topmodels")
            }
            return(do.call(c, applyfun(seq_len(batch_n), batch_fn, z = x)))
        }
        # Using `lapply()` to iterate over observations `x` by calling `fn(x)`
        rval <- do.call(c, lapply(x, fn))
    ## Multiple distributions `y`, one (or multiple) observations `x` evaluated for each `y`
    } else {
        ## Scoping `batch_id`, `y`, `x`, `p`
        batch_fn <- function(i) {
            idx  <- which(batch_id == i) ## Index of `x`/`y` falling into current batch `i`
            q    <- quantile(y[idx], p)  ## Calculating quantiles at `p` for `y[idx]`
            px   <- cdf(y[idx], x[idx])  ## Probability at `x[idx]`
            .Call("c_CRPS_numeric", as.numeric(x[idx]), px, p, q, TRUE, PACKAGE = "topmodels")
        }
        rval <- do.call(c, applyfun(seq_len(batch_n), batch_fn))
    }
  }

  ## handle dimensions
  if (k == 1L && drop) {
    rval <- as.vector(rval)
    names(rval) <- rnam
  } else if (n == 1L && drop) {
    rval <- as.vector(rval)
  } else {
    dim(rval) <- c(n, k)
    dimnames(rval) <- list(rnam, cnam)
  }

  return(rval)
}
