# -------------------------------------------------------
# Testing numeric approximations for the Poisson distribution
# given we only have cdf.pdfPoisson, is_discrete.pdfPoisson, and
# support.pdfPoisson.
# -------------------------------------------------------

if (interactive()) { library("tinytest"); library("topmodels") }

suppressPackageStartupMessages(library("distributions3"))

## ------------- custom Poisson distribution (pdfPoisson) ----------------

## Constructor function for new 'pdfPoisson' distribution
pdfPoisson <- function(lambda) {
    d <- data.frame(lambda = lambda)
    class(d) <- c("pdfPoisson", "distribution")
    return(d)
}

## Additional S3 methods required
pdf.pdfPoisson         <- getS3method("pdf", class = "Poisson")
is_discrete.pdfPoisson <- getS3method("is_discrete", class = "Poisson")
support.pdfPoisson     <- getS3method("support", class = "Poisson")

## Registering S3 methods (required for testing)
registerS3method("pdf",         "pdfPoisson", pdf.pdfPoisson,         envir = asNamespace("topmodels"))
registerS3method("is_discrete", "pdfPoisson", is_discrete.pdfPoisson, envir = asNamespace("topmodels"))
registerS3method("support",     "pdfPoisson", support.pdfPoisson,     envir = asNamespace("topmodels"))


## Setting up single-distribution object
expect_silent(p1 <- pdfPoisson(lambda = 5),
    info = "Creator function must be silent")
expect_identical(p1,
    structure(data.frame(lambda = 5), class = c("pdfPoisson", "distribution")),
    info = "Object returned by pdfPoisson not as expected.")
expect_identical(is_discrete(p1), TRUE,
    info = "is_discrete() must return FALSE.")
expect_identical(support(p1), c(min = 0, max = Inf),
    info = "support() returns incorrect object.")
expect_identical(support(p1, drop = FALSE),
    matrix(c(0, Inf), nrow = 1, dimnames = list(NULL, c("min", "max"))),
    info = "support(..., drop = FALSE) returns incorrect object.")

## Setting up named object with three distributions
expect_silent(p3 <- setNames(pdfPoisson(lambda = c(2, 3.5, 5)), LETTERS[1:3]),
    info = "Creator function must be silent")
expect_identical(p3,
    structure(data.frame(lambda = c(2, 3.5, 5)),
              class = c("pdfPoisson", "distribution"),
              row.names = LETTERS[1:3]),
    info = "Object returned by pdfPoisson not as expected.")
expect_identical(is_discrete(p3), setNames(rep(TRUE, 3L), LETTERS[1:3]),
    info = "is_discrete() must named vector containing TRUE.")
expect_identical(support(p3),
    matrix(c(0, Inf), nrow = 3, ncol = 2, byrow = TRUE,
           dimnames = list(LETTERS[1:3], c("min", "max"))),
    info = "support(..., drop = FALSE) returns incorrect object.")

## Setting up 'analytic' Poisson distributions used for comparison
P1 <- Poisson(lambda = 5)
P3 <- setNames(Poisson(lambda = c(2, 3.5, 5)), LETTERS[1:3])


# --------------- pdf ---------------
expect_identical(pdf(p1, x = 2), pdf(P1, x = 2),
        info = "pdf of pdfPoisson and Poisson must be identical.")
expect_identical(pdf(p3, x = 2), pdf(P3, x = 2),
        info = "pdf of pdfPoisson and Poisson must be identical.")
expect_identical(pdf(p3, x = 2, drop = FALSE), pdf(P3, x = 2, drop = FALSE),
        info = "pdf of pdfPoisson and Poisson must be identical.")
expect_identical(pdf(p3, x = 1:5), pdf(P3, x = 1:5),
        info = "pdf of pdfPoisson and Poisson must be identical.")
expect_identical(pdf(p3, x = 1:5, log = TRUE), pdf(P3, x = 1:5, log = TRUE),
        info = "pdf of pdfPoisson and Poisson w/ log = TRUE must be identical.")

# --------------- cdf ---------------
# Uses the pdf() function and iteratively sums up the pdf (discrete count data)
# to evaluate the CDF.
expect_equal(cdf(p1, x = 2), cdf(P1, x = 2),
        info = "pdf of pdfPoisson and Poisson should be nearly identical.")
expect_equal(cdf(p3, x = 1:3), cdf(P3, x = 1:3),
        info = "pdf of pdfPoisson using stats::numericDeriv and Poisson should be nearly identical.")
expect_identical(cdf(p3, x = 2), cdf(p3, x = rep(2, 3L)),
        info = "Argument 'x' not properly replicated.")

# Testing return objects when using w/ and w/o drop, elementwise, ...
expect_equal(cdf(p1, x = 2, drop = FALSE), cdf(P1, x = 2, drop = FALSE))
expect_equal(cdf(p3, x = 1:3, drop = FALSE), cdf(P3, x = 1:3, drop = FALSE))
expect_equal(cdf(p3, x = 1:3), cdf(P3, x = 1:3))
expect_equal(cdf(p3, x = 1:3, drop = FALSE, elementwise = TRUE), cdf(P3, x = 1:3, drop = FALSE, elementwise = TRUE))

# log(pdf)
cdf(p3, x = 1:5, lower.tail = FALSE)
cdf(p3, x = 1:5, lower.tail = TRUE)
cdf(P3, x = 1:5, lower.tail = FALSE)
cdf(P3, x = 1:5, lower.tail = TRUE)
expect_equal(cdf(p3, x = 1:5, lower.tail = FALSE), cdf(P3, x = 1:5, lower.tail = FALSE),
        info = "cdf of pdfPoisson and Poisson using lower.tail = FALSE should be nearly identical.")

# ------------- quantile ------------
p <- c(0, 0.00001, 0.001, seq(0.01, 0.99, by = 0.01), 0.999, 0.99999, 1)
expect_equal(quantile(p1, probs = p), quantile(P1, probs = p),
        info = "quantile of pdfPoisson and Poisson should be nearly identical.")

expect_equal(quantile(p3, probs = c(0.25, 0.5, 0.75)),
             quantile(P3, probs = c(0.25, 0.5, 0.75)),
        info = "quantile of pdfPoisson and Poisson should be nearly identical.")
expect_equal(quantile(p3, probs = c(0.25, 0.5, 0.75), drop = FALSE),
             quantile(P3, probs = c(0.25, 0.5, 0.75), drop = FALSE),
        info = "quantile of pdfPoisson and Poisson should be nearly identical.")
expect_equal(quantile(p3, probs = c(0.25, 0.5, 0.75), elementwise = FALSE),
             quantile(P3, probs = c(0.25, 0.5, 0.75), elementwise = FALSE),
        info = "quantile of pdfPoisson and Poisson should be nearly identical.")


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

# --------- central moments ---------



# Using a farily small gridsize to make these tests decently fast, on cost of
# the precision. Thus the tolerance is sometimes fairly big. By default, gridsize = 500 is used.

# Testing mean

# default method falls to 'cdf' as we have no quantile function for MyMean
expect_equal(mean(p3, gridsize = 50), mean(P3),
        tolerance = 1e-7,
        info = "mean of pdfPoisson and Poisson should be nearly identical")
expect_equal(mean(p3, gridsize = 50), mean(p3, gridsize = 50, method = "cdf"),
        info = "Expected mean to use method='cdf' by default, but got different results.")
# Testing if method = "quantile" gives similar results; typically we would like to
# use method = 'cdf' but if the range of 'x' to cover the distribution(s) is larger than
# the grid size we fall back to 'quantile' and handle the discrete distribution as a
# continuous distribution; this is less precise though.
expect_equal(mean(p3, gridsize = 50, method = "quantile"), mean(P3),
        tolerance = 1e-1,
        info = "Results using mean(..., method = \"quantile\") are not correct.")

# Testing variance

# default method falls to 'cdf' as we have no quantile function for MyMean
expect_equal(variance(p3, gridsize = 50), variance(p3, gridsize = 50, method = "cdf"),
        info = "Expected variance to use method='cdf' by default, but got different results.")
# TODO(R): Incorrect # # Testing if method = "quantile" gives similar results
# TODO(R): Incorrect # expect_equal(variance(p3, gridsize = 50, method = "quantile"), variance(P3),
# TODO(R): Incorrect #         info = "Results using variance(..., method = \"quantile\") are not correct.")

# Testing skewness

# default method falls to 'cdf' as we have no quantile function for MyMean
expect_equal(skewness(p3), skewness(P3),
        tolerance = 1e-5,
        info = "skewness of pdfPoisson and Poisson should be nearly identical")
expect_equal(skewness(p3, gridsize = 50), skewness(p3, gridsize = 50, method = "cdf"),
        info = "Expected skewness to use method='cdf' by default, but got different results.")
# TODO(R): Incorrect # # Testing if method = "quantile" gives similar results
# TODO(R): Incorrect # expect_equal(skewness(p3, gridsize = 50, method = "quantile"), skewness(P3),
# TODO(R): Incorrect #         info = "Results using skewness(..., method = \"quantile\") are not correct.")

# Testing kurtosis

# default method falls to 'cdf' as we have no quantile function for MyMean
expect_equal(kurtosis(p3), kurtosis(P3),
        tolerance = 1e-4,
        info = "kurtosis of pdfPoisson and Poisson should be nearly identical")
expect_equal(kurtosis(p3, gridsize = 50), kurtosis(p3, gridsize = 50, method = "cdf"),
        info = "Expected kurtosis to use method='cdf' by default, but got different results.")
# TODO(R): Incorrect # # Testing if method = "quantile" gives similar results
# TODO(R): Incorrect # expect_equal(kurtosis(p3, gridsize = 50, method = "quantile"), kurtosis(P3),
# TODO(R): Incorrect #         info = "Results using kurtosis(..., method = \"quantile\") are not correct.")

# TODO(R): moments for quantile.pdfPoisson(..., method = 'q') must be fixed; the code
#          below is just for testing. Remove once that issue is resolved.
#f <- function(...) { devtools::document("../../"); devtools:: load_all("../../") }
#f()
#
#rbind(analytic     = mean(P3),
#      via_cdf      = mean(p3, m = "c"),
#      via_quantile = mean(p3, m = "q", gridsize = 101))
#
#rbind(analytic     = variance(P3),
#      via_cdf      = variance(p3, m = "c"),
#      via_quantile = variance(p3, m = "q", gridsize = 101))
#
#rbind(analytic     = skewness(P3),
#      via_cdf      = skewness(p3, m = "c"),
#      via_quantile = skewness(p3, m = "q", gridsize = 101))
#
#rbind(analytic     = kurtosis(P3),
#      via_cdf      = kurtosis(p3, m = "c"),
#      via_quantile = kurtosis(p3, m = "q", gridsize = 101))

