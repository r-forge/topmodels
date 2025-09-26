# -------------------------------------------------------
# Testing numeric approximations for the Normal distribution
# given we only have cdf.cdfNormal, is_discrete.cdfNormal, and
# support.cdfNormal.
# -------------------------------------------------------

if (interactive()) { library("tinytest"); library("topmodels") }

suppressPackageStartupMessages(library("distributions3"))

## ------------- custom Normal distribution (cdfNormal) ----------------

## Constructor function for new 'cdfNormal' distribution
cdfNormal <- function(mu, sigma) {
    d <- data.frame(mu = mu, sigma = sigma)
    class(d) <- c("cdfNormal", "distribution")
    return(d)
}

## Additional S3 methods required (currently unregistered, thus invisible to tinytest)
cdf.cdfNormal         <- getS3method("cdf", class = "Normal")
is_discrete.cdfNormal <- getS3method("is_discrete", class = "Normal")
support.cdfNormal     <- getS3method("support", class = "Normal")

## Setting up single-distribution object, testing constructor function
expect_silent(n1 <- cdfNormal(mu = 3, sigma = 2),
    info = "Creator function must be silent")
expect_identical(n1,
    structure(data.frame(mu = 3, sigma = 2), class = c("cdfNormal", "distribution")),
    info = "Object returned by cdfNormal not as expected.")

## --------------- Checking behavior when S3 methods are missing -----------------
## Currently cdf.cdfNormal is not defined which we need (in these tests)
## to approximate pdf/quantile. In this case pdf() and quantile() must throw errors.
expect_error(pdf(n1, x = 0),
             pattern = "^no S3 method 'cdf' found for object of class: cdfNormal$",
             info = "No 'cdf.cdfNormal' method available must throw an error.")
expect_error(quantile(n1, probs = 0.5),
             pattern = "^no S3 method 'cdf' or 'quantile' found for object of class: cdfNormal$",
             info = "No method 'cdf.cdfNormal' or 'quantile.cdfNormal' available must throw an error.")
# ... registering pdf (required for tinytest)
registerS3method("cdf", "cdfNormal", cdf.cdfNormal, envir = asNamespace("topmodels"))

## Now the cdf method is registered which is used to approximate cdf/quantile.
## However, is_discrete (required) still missing.
expect_error(pdf(n1, x = 0),
             pattern = "^S3 method 'is_discrete' missing for object of class: cdfNormal$",
             info = "Method is_discrete is required, else an error must be thrown.")
## ... registering is_discrete
registerS3method("is_discrete", "cdfNormal", is_discrete.cdfNormal, envir = asNamespace("topmodels"))

## Same (error) should happen as support is missing
expect_error(pdf(n1, x = 0),
             pattern = "^S3 method 'support' missing for object of class: cdfNormal$",
             info = "Method support is required, else an error must be thrown.")
## ... registering support
registerS3method("support", "cdfNormal", support.cdfNormal, envir = asNamespace("topmodels"))
## ------------- End of checking behavior when S3 methods are missing -------------


## Testing required S3 methods needed
expect_identical(is_discrete(n1), FALSE,
    info = "is_discrete() must return FALSE.")
expect_identical(support(n1), c(min = -Inf, max = Inf),
    info = "support() returns incorrect object.")
expect_identical(support(n1, drop = FALSE),
    matrix(c(-Inf, Inf), nrow = 1, dimnames = list(NULL, c("min", "max"))),
    info = "support(..., drop = FALSE) returns incorrect object.")

## Setting up named object with three distributions
expect_silent(n3 <- setNames(cdfNormal(mu = 1:3, sigma = c(1, 2.5, 5)), LETTERS[1:3]),
    info = "Creator function must be silent")
expect_identical(n3,
    structure(data.frame(mu = 1:3, sigma = c(1, 2.5, 5)),
              class = c("cdfNormal", "distribution"),
              row.names = LETTERS[1:3]),
    info = "Object returned by cdfNormal not as expected.")
expect_identical(is_discrete(n3), setNames(rep(FALSE, 3L), LETTERS[1:3]),
    info = "is_discrete() must named vector containing FALSE.")
expect_identical(support(n3),
    matrix(c(-Inf, Inf), nrow = 3, ncol = 2, byrow = TRUE,
           dimnames = list(LETTERS[1:3], c("min", "max"))),
    info = "support(..., drop = FALSE) returns incorrect object.")

## Setting up 'analytic' Normal distributions used for comparison
N1 <- Normal(mu = 3, sigma = 2)
N3 <- setNames(Normal(mu = 1:3, sigma = c(1, 2.5, 5)), LETTERS[1:3])


# --------------- cdf ---------------
expect_identical(cdf(n1, x = 2), cdf(N1, x = 2),
        info = "cdf of cdfNormal and Normal must be identical.")
expect_identical(cdf(n3, x = 2), cdf(N3, x = 2),
        info = "cdf of cdfNormal and Normal must be identical.")
expect_identical(cdf(n3, x = 2, drop = FALSE), cdf(N3, x = 2, drop = FALSE),
        info = "cdf of cdfNormal and Normal must be identical.")
expect_identical(cdf(n3, x = 1:5), cdf(N3, x = 1:5),
        info = "cdf of cdfNormal and Normal must be identical.")
expect_identical(cdf(n3, x = 1:5, lower.tail = FALSE), cdf(N3, x = 1:5, lower.tail = FALSE),
        info = "cdf of cdfNormal and Normal w/ lower.tail = FALSE must be identical.")

# --------------- pdf ---------------
# The first call uses the default (numDeriv::grad if installed, else stats::numericDeriv),
# we then test both methods individually using the undocumented deriv.method dev argument.
expect_equal(pdf(n1, x = 2), pdf(N1, x = 2),
        info = "pdf of cdfNormal and Normal should be nearly identical.")

expect_equal(pdf(n1, x = 2, deriv.method = "numericDeriv"), pdf(N1, x = 2),
        info = "pdf of cdfNormal using stats::numericDeriv and Normal should be nearly identical.")
expect_equal(pdf(n1, x = 2, deriv.method = "grad"), pdf(N1, x = 2),
        info = "pdf of cdfNormal using numDeriv::grad and Normal should be nearly identical.")

expect_equal(pdf(n3, x = 1:3, deriv.method = "numericDeriv"), pdf(N3, x = 1:3),
        info = "pdf of cdfNormal using stats::numericDeriv and Normal should be nearly identical.")
expect_equal(pdf(n3, x = 1:3, deriv.method = "numericDeriv"), pdf(N3, x = 1:3),
        info = "pdf of cdfNormal using stats::numericDeriv and Normal should be nearly identical.")

expect_identical(pdf(n3, x = 2), pdf(n3, x = rep(2, 3L)),
        info = "Argument 'x' not properly replicated.")

# Testing return objects when using w/ and w/o drop, elementwise, ...
expect_equal(pdf(n1, x = 2, drop = FALSE), pdf(N1, x = 2, drop = FALSE))
expect_equal(pdf(n3, x = 1:3, drop = FALSE), pdf(N3, x = 1:3, drop = FALSE))
expect_equal(pdf(n3, x = 1:3), pdf(N3, x = 1:3))
expect_equal(pdf(n3, x = 1:3, drop = FALSE, elementwise = TRUE), pdf(N3, x = 1:3, drop = FALSE, elementwise = TRUE))

# log(pdf)
expect_equal(pdf(n3, x = 1:5, log = TRUE, deriv.method = "numericDeriv"), pdf(N3, x = 1:5, log = TRUE),
        tolerance = 1e-7, # slightly less precise than the numDeriv::grad method
        info = "pdf of cdfNormal using stats::numericDeriv and Normal should be nearly identical.")
expect_equal(pdf(n3, x = 1:5, log = TRUE, deriv.method = "grad"), pdf(N3, x = 1:5, log = TRUE),
        info = "pdf of cdfNormal using numDeriv::grad and Normal should be nearly identical.")
expect_equal(pdf(n3, 1:3, log = TRUE), pdf(N3, 1:3, log = TRUE),
        info = "Return should be a named vector.")

# ------------- quantile ------------
p <- c(0, 0.00001, 0.001, seq(0.01, 0.99, by = 0.01), 0.999, 0.99999, 1)
expect_equal(quantile(n1, probs = p), quantile(N1, probs = p),
        info = "quantile of cdfNormal and Normal should be nearly identical.")

expect_equal(quantile(n3, probs = c(0.25, 0.5, 0.75)),
             quantile(N3, probs = c(0.25, 0.5, 0.75)),
        info = "quantile of cdfNormal and Normal should be nearly identical.")
expect_equal(quantile(n3, probs = c(0.25, 0.5, 0.75), drop = FALSE),
             quantile(N3, probs = c(0.25, 0.5, 0.75), drop = FALSE),
        info = "quantile of cdfNormal and Normal should be nearly identical.")
expect_equal(quantile(n3, probs = c(0.25, 0.5, 0.75), elementwise = FALSE),
             quantile(N3, probs = c(0.25, 0.5, 0.75), elementwise = FALSE),
        info = "quantile of cdfNormal and Normal should be nearly identical.")


# ------- random number gen ---------

## Given this is random by nature I am only checking the return objects

expect_silent(r1 <- random(n1),                   info = "random(n1) should run silent.")
expect_inherits(r1, "numeric",                    info = "random(n1) should return numeric.")
expect_identical(length(r1), 1L,                  info = "random(n1) should return single numeric.")
expect_null(names(r1),                            info = "random(n1) should return unnamed vector.")
rm(r1)

expect_silent(r1 <- random(n1, drop = FALSE),     info = "random(n1, drop = FALSE) should run silent.")
expect_true(is.numeric(r1),                       info = "random(n1, drop = FALSE) should return numeric.")
expect_inherits(r1, "matrix",                     info = "random(n1, drop = FALSE) should return matrix.")
expect_identical(dim(r1), c(1L, 1L),              info = "random(n1, drop = FALSE) should return matrix of dimension 1x1.")
expect_identical(dimnames(r1), list(NULL, "r_1"), info = "dimension names of random(n1, drop = FALSE) incorrect.")
rm(r1)

expect_silent(r1 <- random(n1, 50),               info = "random(n1, 50) should run silent.")
expect_inherits(r1, "numeric",                    info = "random(n1) should return numeric.")
expect_identical(length(r1), 50L,                 info = "random(n1) should return numeric of length 50L.")
expect_null(names(r1),                            info = "random(n1) should return unnamed vector.")
rm(r1)

expect_silent(r3 <- random(n3),                   info = "random(n3) should run silent.")
expect_inherits(r3, "numeric",                    info = "random(n3) should return numeric.")
expect_identical(length(r3), 3L,                  info = "random(n3) should return numeric of length 3.")
expect_identical(names(r3), names(n3),            info = "random(n3) should return named vector.")
rm(r3)

expect_silent(r3 <- random(n3, 5, drop = FALSE),  info = "random(n3, 5, drop = FALSE) should run silent.")
expect_true(is.numeric(r3),                       info = "random(n3, 5, drop = FALSE) should return numeric.")
expect_inherits(r3, "matrix",                     info = "random(n3, 5, drop = FALSE) should return matrix.")
expect_identical(dim(r3), c(3L, 5L),              info = "random(n3, 5, drop = FALSE) should return matrix of dimension 3x5.")
expect_identical(dimnames(r3), list(names(n3), paste("r", 1:5, sep = "_")), info = "dimension names of random(n3, 5, drop = FALSE) incorrect.")
rm(r3)

# --------- central moments ---------

# Using a farily small gridsize to make these tests decently fast, on cost of
# the precision. Thus the tolerance is sometimes fairly big. By default, gridsize = 500 is used.

# Testing mean

# default method falls to 'cdf' as we have no quantile function for MyMean
expect_equal(mean(n3), mean(N3),
        info = "mean of cdfNormal and Normal should be nearly identical")
expect_equal(mean(n3, gridsize = 50), mean(n3, gridsize = 50, method = "cdf"),
        tolerance = 1e-2,
        info = "Expected mean to use method='cdf' by default, but got different results.")
# Testing if method = "quantile" gives similar results
expect_equal(mean(n3, gridsize = 50, method = "quantile"), mean(N3),
        tolerance = 1e-2,
        info = "Results using mean(..., method = \"quantile\") are not correct.")

# Testing variance

# default method falls to 'cdf' as we have no quantile function for MyMean
expect_equal(variance(n3), variance(N3),
        tolerance = 1e-2,
        info = "variance of cdfNormal and Normal should be nearly identical")
expect_equal(variance(n3, gridsize = 50), variance(n3, gridsize = 50, method = "cdf"),
        tolerance = 1e-2,
        info = "Expected variance to use method='cdf' by default, but got different results.")
# Testing if method = "quantile" gives similar results
expect_equal(variance(n3, gridsize = 50, method = "quantile"), variance(N3),
        tolerance = 1e-1,
        info = "Results using variance(..., method = \"quantile\") are not correct.")

# Testing skewness

# default method falls to 'cdf' as we have no quantile function for MyMean
expect_equal(skewness(n3), skewness(N3),
        info = "skewness of cdfNormal and Normal should be nearly identical")
expect_equal(skewness(n3, gridsize = 50), skewness(n3, gridsize = 50, method = "cdf"),
        info = "Expected skewness to use method='cdf' by default, but got different results.")
# Testing if method = "quantile" gives similar results
expect_equal(skewness(n3, gridsize = 50, method = "quantile"), skewness(N3),
        info = "Results using skewness(..., method = \"quantile\") are not correct.")

# Testing kurtosis

# default method falls to 'cdf' as we have no quantile function for MyMean
expect_equal(kurtosis(n3), kurtosis(N3),
        tolerance = 1e-1,
        info = "kurtosis of cdfNormal and Normal should be nearly identical")
expect_equal(kurtosis(n3, gridsize = 50), kurtosis(n3, gridsize = 50, method = "cdf"),
        info = "Expected kurtosis to use method='cdf' by default, but got different results.")
# Testing if method = "quantile" gives similar results
expect_equal(kurtosis(n3, gridsize = 50, method = "quantile"), kurtosis(N3),
        tolerance = 1e-1,
        info = "Results using kurtosis(..., method = \"quantile\") are not correct.")

