# -------------------------------------------------------
# Testing numeric approximations for the Poisson distribution
# given we only have cdf.cdfPoisson, is_discrete.cdfPoisson, and
# support.cdfPoisson.
# -------------------------------------------------------

if (interactive()) { library("tinytest"); library("topmodels") }

suppressPackageStartupMessages(library("distributions3"))

## ------------- custom Poisson distribution (cdfPoisson) ----------------

## Constructor function for new 'cdfPoisson' distribution
cdfPoisson <- function(lambda) {
    d <- data.frame(lambda = lambda)
    class(d) <- c("cdfPoisson", "distribution")
    return(d)
}

## Additional S3 methods required (currently unregistered, thus invisible to tinytest)
cdf.cdfPoisson         <- getS3method("cdf", class = "Poisson")
is_discrete.cdfPoisson <- getS3method("is_discrete", class = "Poisson")
support.cdfPoisson     <- getS3method("support", class = "Poisson")

## Setting up single-distribution object, testing constructor function
expect_silent(p1 <- cdfPoisson(lambda = 5),
    info = "Creator function must be silent")
expect_identical(p1,
    structure(data.frame(lambda = 5), class = c("cdfPoisson", "distribution")),
    info = "Object returned by cdfPoisson not as expected.")


## --------------- Checking behavior when S3 methods are missing -----------------
## Currently cdf.cdfPoisson is not defined which we need (in these tests)
## to approximate pdf/quantile. In this case pdf() and quantile() must throw errors.
expect_error(pdf(p1, x = 0),
             pattern = "^no S3 method 'cdf' found for object of class: cdfPoisson$",
             info = "No 'cdf.cdfPoisson' method available must throw an error.")
expect_error(quantile(p1, probs = 0.5),
             pattern = "^no S3 method 'cdf' or 'quantile' found for object of class: cdfPoisson$",
             info = "No method 'cdf.cdfPoisson' or 'quantile.cdfPoisson' available must throw an error.")
# ... registering pdf (required for tinytest)
registerS3method("cdf", "cdfPoisson", cdf.cdfPoisson, envir = asNamespace("topmodels"))

## Now the cdf method is registered which is used to approximate cdf/quantile.
## However, is_discrete (required) still missing.
expect_error(pdf(p1, x = 0),
             pattern = "^S3 method 'is_discrete' missing for object of class: cdfPoisson$",
             info = "Method is_discrete is required, else an error must be thrown.")
## ... registering is_discrete
registerS3method("is_discrete", "cdfPoisson", is_discrete.cdfPoisson, envir = asNamespace("topmodels"))

## Same (error) should happen as support is missing
expect_error(pdf(p1, x = 0),
             pattern = "^S3 method 'support' missing for object of class: cdfPoisson$",
             info = "Method support is required, else an error must be thrown.")
## ... registering support
registerS3method("support", "cdfPoisson", support.cdfPoisson, envir = asNamespace("topmodels"))
## ------------- End of checking behavior when S3 methods are missing -------------


## Testing required S3 methods needed
expect_identical(is_discrete(p1), TRUE,
    info = "is_discrete() must return FALSE.")
expect_identical(support(p1), c(min = 0, max = Inf),
    info = "support() returns incorrect object.")
expect_identical(support(p1, drop = FALSE),
    matrix(c(0, Inf), nrow = 1, dimnames = list(NULL, c("min", "max"))),
    info = "support(..., drop = FALSE) returns incorrect object.")

## Setting up named object with three distributions
expect_silent(p3 <- setNames(cdfPoisson(lambda = c(2, 3.5, 5)), LETTERS[1:3]),
    info = "Creator function must be silent")
expect_identical(p3,
    structure(data.frame(lambda = c(2, 3.5, 5)),
              class = c("cdfPoisson", "distribution"),
              row.names = LETTERS[1:3]),
    info = "Object returned by cdfPoisson not as expected.")
expect_identical(is_discrete(p3), setNames(rep(TRUE, 3L), LETTERS[1:3]),
    info = "is_discrete() must named vector containing TRUE.")
expect_identical(support(p3),
    matrix(c(0, Inf), nrow = 3, ncol = 2, byrow = TRUE,
           dimnames = list(LETTERS[1:3], c("min", "max"))),
    info = "support(..., drop = FALSE) returns incorrect object.")

## Setting up 'analytic' Poisson distributions used for comparison
P1 <- Poisson(lambda = 5)
P3 <- setNames(Poisson(lambda = c(2, 3.5, 5)), LETTERS[1:3])


# --------------- cdf ---------------
expect_identical(cdf(p1, x = 2), cdf(P1, x = 2),
        info = "cdf of cdfPoisson and Poisson must be identical.")
expect_identical(cdf(p3, x = 2), cdf(P3, x = 2),
        info = "cdf of cdfPoisson and Poisson must be identical.")
expect_identical(cdf(p3, x = 2, drop = FALSE), cdf(P3, x = 2, drop = FALSE),
        info = "cdf of cdfPoisson and Poisson must be identical.")
expect_identical(cdf(p3, x = 1:5), cdf(P3, x = 1:5),
        info = "cdf of cdfPoisson and Poisson must be identical.")
expect_identical(cdf(p3, x = 1:5, lower.tail = FALSE), cdf(P3, x = 1:5, lower.tail = FALSE),
        info = "cdf of cdfPoisson and Poisson w/ lower.tail = FALSE must be identical.")

# --------------- pdf ---------------
# The first call uses the default (numDeriv::grad if installed, else stats::numericDeriv),
# we then test both methods individually using the undocumented deriv.method dev argument.
expect_equal(pdf(p1, x = 2), pdf(P1, x = 2),
        info = "pdf of cdfPoisson and Poisson should be nearly identical.")
expect_equal(pdf(p3, x = 1:3), pdf(P3, x = 1:3),
        info = "pdf of cdfPoisson using stats::numericDeriv and Poisson should be nearly identical.")
expect_identical(pdf(p3, x = 2), pdf(p3, x = rep(2, 3L)),
        info = "Argument 'x' not properly replicated.")

# Testing return objects when using w/ and w/o drop, elementwise, ...
expect_equal(pdf(p1, x = 2, drop = FALSE), pdf(P1, x = 2, drop = FALSE))
expect_equal(pdf(p3, x = 1:3, drop = FALSE), pdf(P3, x = 1:3, drop = FALSE))
expect_equal(pdf(p3, x = 1:3), pdf(P3, x = 1:3))
expect_equal(pdf(p3, x = 1:3, drop = FALSE, elementwise = TRUE), pdf(P3, x = 1:3, drop = FALSE, elementwise = TRUE))

# log(pdf)
expect_equal(pdf(p3, x = 1:5, log = TRUE), pdf(P3, x = 1:5, log = TRUE),
        info = "pdf of cdfPoisson using and Poisson should be nearly identical.")

# ------------- quantile ------------
p <- c(0, 0.00001, 0.001, seq(0.01, 0.99, by = 0.01), 0.999, 0.99999, 1)
expect_equal(quantile(p1, probs = p), quantile(P1, probs = p),
        info = "quantile of cdfPoisson and Poisson should be nearly identical.")

expect_equal(quantile(p3, probs = c(0.25, 0.5, 0.75)),
             quantile(P3, probs = c(0.25, 0.5, 0.75)),
        info = "quantile of cdfPoisson and Poisson should be nearly identical.")
expect_equal(quantile(p3, probs = c(0.25, 0.5, 0.75), drop = FALSE),
             quantile(P3, probs = c(0.25, 0.5, 0.75), drop = FALSE),
        info = "quantile of cdfPoisson and Poisson should be nearly identical.")
expect_equal(quantile(p3, probs = c(0.25, 0.5, 0.75), elementwise = FALSE),
             quantile(P3, probs = c(0.25, 0.5, 0.75), elementwise = FALSE),
        info = "quantile of cdfPoisson and Poisson should be nearly identical.")


# ------- random number gen ---------

## Given this is random by nature I am only checking the return objects

expect_silent(r1 <- random(p1),                   info = "random(p1) should run silent.")
expect_inherits(r1, "numeric",                    info = "random(p1) should return numeric.")
expect_identical(length(r1), 1L,                  info = "random(p1) should return single numeric.")
expect_null(names(r1),                            info = "random(p1) should return unnamed vector.")
rm(r1)

expect_silent(r1 <- random(p1, drop = FALSE),     info = "random(p1, drop = FALSE) should run silent.")
expect_true(is.numeric(r1),                       info = "random(p1, drop = FALSE) should return numeric.")
expect_inherits(r1, "matrix",                     info = "random(p1, drop = FALSE) should return matrix.")
expect_identical(dim(r1), c(1L, 1L),              info = "random(p1, drop = FALSE) should return matrix of dimension 1x1.")
expect_identical(dimnames(r1), list(NULL, "r_1"), info = "dimension names of random(p1, drop = FALSE) incorrect.")
rm(r1)

expect_silent(r1 <- random(p1, 50),               info = "random(p1, 50) should run silent.")
expect_inherits(r1, "numeric",                    info = "random(p1) should return numeric.")
expect_identical(length(r1), 50L,                 info = "random(p1) should return numeric of length 50L.")
expect_null(names(r1),                            info = "random(p1) should return unnamed vector.")
rm(r1)

expect_silent(r3 <- random(p3),                   info = "random(p3) should run silent.")
expect_inherits(r3, "numeric",                    info = "random(p3) should return numeric.")
expect_identical(length(r3), 3L,                  info = "random(p3) should return numeric of length 3.")
expect_identical(names(r3), names(p3),            info = "random(p3) should return named vector.")
rm(r3)

expect_silent(r3 <- random(p3, 5, drop = FALSE),  info = "random(p3, 5, drop = FALSE) should run silent.")
expect_true(is.numeric(r3),                       info = "random(p3, 5, drop = FALSE) should return numeric.")
expect_inherits(r3, "matrix",                     info = "random(p3, 5, drop = FALSE) should return matrix.")
expect_identical(dim(r3), c(3L, 5L),              info = "random(p3, 5, drop = FALSE) should return matrix of dimension 3x5.")
expect_identical(dimnames(r3), list(names(p3), paste("r", 1:5, sep = "_")), info = "dimension names of random(p3, 5, drop = FALSE) incorrect.")
rm(r3)

# --------- central moments (default; cdf) ---------
# These tests are written to be using the default methods which, here, is 
# method = 'cdf' as cdf.cdfPoisson is available. We later enforce/test
# method = 'quantile'.

# Using a farily small gridsize to make these tests decently fast, on cost of
# the precision. Thus the tolerance is sometimes fairly big. By default,
# gridsize = 500 is used.

# Testing mean (default cdf method)
expect_equal(mean(p3, gridsize = 50), mean(P3),
        tolerance = 1e-7,
        info = "mean of cdfPoisson and Poisson should be nearly identical")
expect_equal(mean(p3, gridsize = 50), mean(p3, gridsize = 50, method = "cdf"),
        info = "Expected mean to use method='cdf' by default, but got different results.")

# Testing variance (default cdf method)
expect_equal(variance(p3), variance(P3),
        tolerance = 1e-5,
        info = "variance of cdfPoisson and Poisson should be nearly identical")
expect_equal(variance(p3, gridsize = 50), variance(p3, gridsize = 50, method = "cdf"),
        info = "Expected variance to use method='cdf' by default, but got different results.")

# Testing skewness (default cdf method)
expect_equal(skewness(p3), skewness(P3),
        tolerance = 1e-5,
        info = "skewness of cdfPoisson and Poisson should be nearly identical")
expect_equal(skewness(p3, gridsize = 50), skewness(p3, gridsize = 50, method = "cdf"),
        info = "Expected skewness to use method='cdf' by default, but got different results.")

# Testing kurtosis (default cdf method)
expect_equal(kurtosis(p3), kurtosis(P3),
        tolerance = 1e-4,
        info = "kurtosis of cdfPoisson and Poisson should be nearly identical")
expect_equal(kurtosis(p3, gridsize = 50), kurtosis(p3, gridsize = 50, method = "cdf"),
        info = "Expected kurtosis to use method='cdf' by default, but got different results.")


# --------- central moments (enforcing quantile method) -------
# Next we enforce 'method = "quantile"' for distributions with
# a fairly small parameter (lambda) where the method should fall 
# back to the 'cdf' method as the 'quantile' method may result in
# inaccurate results (mainly for third and forth order moments).

# Mean (enforcing quantile, falling back to 'cdf')
expect_warning(m <- mean(p3, method = "quantile"),
               pattern = "switching to method = 'cdf'",
               info = "Expected warning switching back to cdf method.")
expect_equal(m, mean(P3),
             tolerance = 1e-7,
             info = "Approximated and analytic mean expected to be equal.")
rm(m)

# Variance (enforcing quantile, falling back to 'cdf')
expect_warning(v <- variance(p3, method = "quantile"),
               pattern = "switching to method = 'cdf'",
               info = "Expected warning switching back to cdf method.")
expect_equal(v, variance(P3),
             tolerance = 1e-5,
             info = "Approximated and analytic variance expected to be equal.")
rm(v)

# Skewness (enforcing quantile, falling back to 'cdf')
expect_warning(s <- mean(p3, method = "quantile"),
               pattern = "switching to method = 'cdf'",
               info = "Expected warning switching back to cdf method.")
expect_equal(s, mean(P3),
             tolerance = 1e-5,
             info = "Approximated and analytic mean expected to be equal.")
rm(s)

# Kurtosis (enforcing quantile, falling back to 'cdf')
expect_warning(k <- kurtosis(p3, method = "quantile"),
               pattern = "switching to method = 'cdf'",
               info = "Expected warning switching back to cdf method.")
expect_equal(k, kurtosis(P3),
             tolerance = 1e-4,
             info = "Approximated and analytic kurtosis expected to be equal.")
rm(k)


# --------- central moments (making use of quantiles) ---------
# Enforcing 'method = "quantile"' where it is actually used (the
# approximated quantiles (cdf -> quantile).

# New set of distributions with large lambda to test 'quantile' approximation,
# where the discrete distribution is handled as a continuous one. Only valid/useful
# if the range to be evaluated exceeds the gridsize.
p3x <- cdfPoisson(lambda = 1:3 * 100)
P3x <- Poisson(lambda = 1:3 * 100)

# Mean (once w/ default method, once enforcing quantile)
expect_silent(m <- mean(p3x, gridsize = 100L),
               info = "Expected mean(..., method = \"quantile\") to be executed silently.")
expect_equal(m, mean(P3x),
             tolerance = 1e-4,
             info = "Approximated and analytic mean expected to be equal.")
rm(m)
expect_silent(m <- mean(p3x, method = "quantile", gridsize = 100L),
               info = "Expected mean(..., method = \"quantile\") to be executed silently.")
expect_equal(m, mean(P3x),
             tolerance = 1e-2,
             info = "Approximated and analytic mean expected to be equal.")
rm(m)

# Variance (once w/ default method, once enforcing quantile)
expect_silent(v <- variance(p3x, gridsize = 100L),
               info = "Expected variance(..., method = \"quantile\") to be executed silently.")
expect_equal(v, variance(P3x),
             tolerance = 1e-2,
             info = "Approximated and analytic variance expected to be equal.")
rm(v)
expect_silent(v <- variance(p3x, method = "quantile", gridsize = 100L),
               info = "Expected variance(..., method = \"quantile\") to be executed silently.")
expect_equal(v, variance(P3x),
             tolerance = 1e-1,
             info = "Approximated and analytic variance expected to be equal.")
rm(v)

# Skewness (once w/ default method, once enforcing quantile)
expect_silent(s <- skewness(p3x, gridsize = 100L),
               info = "Expected skewness(..., method = \"quantile\") to be executed silently.")
expect_equal(s, skewness(P3x),
             tolerance = 0.25,
             info = "Approximated and analytic skewness expected to be equal.")
rm(s)
expect_silent(s <- skewness(p3x, method = "quantile", gridsize = 100L),
               info = "Expected skewness(..., method = \"quantile\") to be executed silently.")
expect_equal(s, skewness(P3x),
             tolerance = 0.25,
             info = "Approximated and analytic skewness expected to be equal.")
rm(s)


# Kurtosis (once w/ default method, once enforcing quantile)
expect_silent(k <- kurtosis(p3x, gridsize = 100L),
               info = "Expected kurtosis(..., method = \"quantile\") to be executed silently.")
expect_equal(k, kurtosis(P3x),
             tolerance = 1e-1,
             info = "Approximated and analytic kurtosis expected to be equal.")
rm(k)
expect_silent(k <- kurtosis(p3x, method = "quantile", gridsize = 100L),
               info = "Expected kurtosis(..., method = \"quantile\") to be executed silently.")
expect_equal(k, kurtosis(P3x),
             tolerance = 0.25,
             info = "Approximated and analytic kurtosis expected to be equal.")
rm(k)
