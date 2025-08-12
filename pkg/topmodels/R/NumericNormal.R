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
#' ## ---------------- CDF ----------------
#' ## 'NumericNormal' only knows the (exact) pdf of a Gaussian distribution
#' ## and approximates the CDF via numeric integration (default grid size n = 101L).
#' cdf(n, 3.5)
#' cdf(r, 3.5)
#'
#' ## Visual example
#' xx <- seq(-5, 15, by = 0.1)
#' plot(xx, cdf(n, xx), main = "Testing numeric PDF->CDF")
#' lines(xx, cdf(r, xx), col = 2, lwd = 3)
#' lines(xx, cdf(r, xx, n = 5L), col = 4, lwd = 3)
#' legend("topleft", legend = c("Normal", "NumericNormal n = 101", "NumericNormal n = 5"),
#'        pch = c(1, NA, NA), lty = c(NA, 1, 1), lwd = c(NA, 3, 3), col = c(1, 2, 4))
#'
#' ## -------------- Quantile -------------
#' ## 'NumericNormal' only knows the PDF from which it calculates the CDF
#' ## to retrieve the quantiles.
#' quantile(n, 0.3)
#' quantile(r, 0.3)
#'
#' qq <- c(0.0001, seq(0.01, 0.99, by = 0.005), 0.9999)
#' plot(quantile(n, qq), qq, main = "Testing numeric PDF->Quantile")
#' lines(quantile(r, qq), qq, col = 2, lwd = 3)
#' lines(quantile(r, qq, n = 5L), qq, col = 4, lwd = 3)
#' legend("topleft",
#'        legend = c("Normal", "NumericNormal n = 101", "NumericNormal n = 5"),
#'        pch = c(1, NA, NA), lty = c(NA, 1, 1), lwd = c(NA, 3, 3), col = c(1, 2, 4))
#'
#' ## Central moments
#' c(mean = mean(r), variance = variance(r),
#'   skewness = skewness(r), kurtosis = kurtosis(r))
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

#' @param d object of class `"NumericNormal"`.
#' @param x numeric, where to evaluate the distribution object.
#'
#' @rdname NumericNormal
#' @exportS3Method
pdf.NumericNormal <- function(d, x, ...) {
    # Just for testing-simplifying the elementwise case
    stopifnot(is.numeric(x) && all(is.finite(x)))

    res <- sapply(seq_along(d), function(i) {
                      dnorm(x, mean = d[i]$mu, sd = sqrt(d[i]$var))
           })
    if (is.matrix(res)) {
        res <- t(res)
        colnames(res) <- paste0("d_", format(x))
    }
    return(res)
}

# TODO(R): This is only for testing; sets up a grid for each
#          distribution alongside a matrix with the corresponding PDF.
#          Used for implementing the basic C functions allowing us to
#          retrieve CDF, quantiles, and moments from a distribution
#          of which we only know the PDF.
NumericNormal_get_grid <- function(d, n) {
    stopifnot(is.numeric(n))
    n <- as.integer(n)[1L]
    stopifnot(is.finite(n) && n >= 5L) # For testing

    # Must guess useful range to calculate the density
    # which is used for numerical integration to calculate the cdf.
    # These are the grid points.
    get_grid <- function(i, n) {
        matrix(seq(d[i]$mu - 4 * sqrt(d[i]$var),
                   d[i]$mu + 4 * sqrt(d[i]$var),
                   length.out = n), nrow = 1L)
    }
    get_dens <- function(i) {
        matrix(pdf(d[i], at[i, ]))
    }
    at  <- t(sapply(seq_along(d), get_grid, n = n))
    pat <- t(sapply(seq_along(d), get_dens))
    return(list(at = at, pat = pat))
}

#' @param n integer, number of grid points used for numerical integration.
#'
#' @rdname NumericNormal
#' @exportS3Method
cdf.NumericNormal <- function(d, x, ..., n = 101L) {

    # TODO(R): This is solely for testing purposes
    dev <- NumericNormal_get_grid(d, n)
    at <- dev$at; pat <- dev$pat; rm(dev)

    # Important, 'probs' must be sorted for the C code to work properly.
    # To be able to restore the original order we keep the order
    # on probsorder
    xorder <- order(x)
    x      <- as.numeric(sort(x))

    res <- .Call("c_d2pq_numeric",
                 at       = as.numeric(at),     # grid points
                 pat      = as.numeric(pat),    # matrix of dimension length(d) x length(at)
                 dim      = dim(at),            # dimension of matrices at/pat
                 x        = x,                  # point(s) to calculate the cdf for
                 discrete = FALSE,
                 what     = "cdf",
                 PACKAGE  = "topmodels")

    if (length(d) == length(x) || length(d) == 1L) {
        return(res[xorder])
    } else {
        return(matrix(res, nrow = length(d))[, xorder])
    }

}

#' @param x object of class `c("NumericNormal", "distribution")`.
#' @param probs numeric vector, probabilities.
#' @param n integer, number of grid points used for numerical integration.
#'
#' @rdname NumericNormal
#' @exportS3Method
quantile.NumericNormal <- function(x, probs, ..., n = 101L) {

    # TODO(R): This is solely for testing purposes
    dev <- NumericNormal_get_grid(x, n)
    at <- dev$at; pat <- dev$pat; rm(dev)

    # Important, 'probs' must be sorted for the C code to work properly.
    # To be able to restore the original order we keep the order
    # on probsorder
    probsorder <- order(probs)
    probs      <- as.numeric(sort(probs))

    res <- .Call("c_d2pq_numeric",
                 at       = as.numeric(at),     # grid points
                 pat      = as.numeric(pat),    # matrix of dimension length(d) x length(at)
                 dim      = dim(at),            # dimension of matrices at/pat
                 x        = probs,              # point(s) to calculate the cdf for
                 discrete = FALSE,
                 what     = "quantile",
                 PACKAGE  = "topmodels")

    if (length(x) == length(probs) || length(x) == 1L) {
        return(res[probsorder])
    } else {
        return(matrix(res, nrow = length(x))[, probsorder])
    }

}



#' @param x object of class `c("NumericNormal", "distribution")`.
#' @param \dots ignored.
#'
#' @rdname NumericNormal
#' @exportS3Method
mean.NumericNormal <- function(x, ..., n = 101L) {

    # TODO(R): This is solely for testing purposes
    dev <- NumericNormal_get_grid(x, n)
    at <- dev$at; pat <- dev$pat; rm(dev)

    res <- .Call("c_d2moments_numeric",
                 at       = as.numeric(at),     # grid points
                 pat      = as.numeric(pat),    # matrix of dimension length(d) x length(at)
                 dim      = dim(at),            # dimension of matrices at/pat
                 discrete = FALSE,
                 what     = 1L,                 # 1L = mean
                 PACKAGE  = "topmodels")

    return(res)
}


#' @param x object of class `c("NumericNormal", "distribution")`.
#' @param \dots ignored.
#'
#' @rdname NumericNormal
#' @exportS3Method
variance.NumericNormal <- function(x, ..., n = 101L) {

    # TODO(R): This is solely for testing purposes
    dev <- NumericNormal_get_grid(x, n)
    at <- dev$at; pat <- dev$pat; rm(dev)

    res <- .Call("c_d2moments_numeric",
                 at       = as.numeric(at),     # grid points
                 pat      = as.numeric(pat),    # matrix of dimension length(d) x length(at)
                 dim      = dim(at),            # dimension of matrices at/pat
                 discrete = FALSE,
                 what     = 2L,                 # 2L = variance
                 PACKAGE  = "topmodels")

    return(res)
}


#' @param x object of class `c("NumericNormal", "distribution")`.
#' @param \dots ignored.
#'
#' @rdname NumericNormal
#' @exportS3Method
skewness.NumericNormal <- function(x, ..., n = 101L) {

    # TODO(R): This is solely for testing purposes
    dev <- NumericNormal_get_grid(x, n)
    at <- dev$at; pat <- dev$pat; rm(dev)

    res <- .Call("c_d2moments_numeric",
                 at       = as.numeric(at),     # grid points
                 pat      = as.numeric(pat),    # matrix of dimension length(d) x length(at)
                 dim      = dim(at),            # dimension of matrices at/pat
                 discrete = FALSE,
                 what     = 3L,                 # 3L = skewness
                 PACKAGE  = "topmodels")

    return(res)
}


#' @param x object of class `c("NumericNormal", "distribution")`.
#' @param \dots ignored.
#'
#' @rdname NumericNormal
#' @exportS3Method
kurtosis.NumericNormal <- function(x, ..., n = 101L) {

    # TODO(R): This is solely for testing purposes
    dev <- NumericNormal_get_grid(x, n)
    at <- dev$at; pat <- dev$pat; rm(dev)

    res <- .Call("c_d2moments_numeric",
                 at       = as.numeric(at),     # grid points
                 pat      = as.numeric(pat),    # matrix of dimension length(d) x length(at)
                 dim      = dim(at),            # dimension of matrices at/pat
                 discrete = FALSE,
                 what     = 4L,                 # 4L = kurtosis
                 PACKAGE  = "topmodels")

    return(res)
}
