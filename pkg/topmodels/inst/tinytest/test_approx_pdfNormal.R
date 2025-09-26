# -------------------------------------------------------
# Testing numeric approximations for the Normal distribution
# given we only have pdf.pdfNormal, is_discrete.pdfNormal, and
# support.pdfNormal.
# -------------------------------------------------------

if (interactive()) { library("tinytest"); library("topmodels") }

suppressPackageStartupMessages(library("distributions3"))

## ------------- custom Normal distribution (pdfNormal) ----------------

## Constructor function for new 'pdfNormal' distribution
pdfNormal <- function(mu, sigma) {
    d <- data.frame(mu = mu, sigma = sigma)
    class(d) <- c("pdfNormal", "distribution")
    return(d)
}

## Additional S3 methods required (currently unregistered, thus invisible to tinytest)
pdf.pdfNormal <- getS3method("pdf", class = "Normal")
is_discrete.pdfNormal   <- getS3method("is_discrete", class = "Normal")
support.pdfNormal       <- getS3method("support", class = "Normal")

## Setting up single-distribution object, testing constructor function
expect_silent(n1 <- pdfNormal(mu = 3, sigma = 2),
    info = "Creator function must be silent")
expect_identical(n1,
    structure(data.frame(mu = 3, sigma = 2), class = c("pdfNormal", "distribution")),
    info = "Object returned by pdfNormal not as expected.")


## --------------- Checking behavior when S3 methods are missing -----------------
## Currently pdf.pdfNormal is not defined which we need (in these tests)
## to approximate cdf/quantile. In this case cdf() and quantile() must throw errors.
expect_error(cdf(n1, x = 0),
             pattern = "^no S3 method 'pdf' found for object of class: pdfNormal$",
             info = "No 'pdf.pdfNormal' method available must throw an error.")
expect_error(quantile(n1, probs = 0.5),
             pattern = "^no S3 method 'cdf' or 'quantile' found for object of class: pdfNormal$",
             info = "No method 'cdf.pdfNormal' or 'quantile.pdfNormal' available must throw an error.")
## ... registering pdf (required for tinytest)
registerS3method("pdf", "pdfNormal", pdf.pdfNormal, envir = asNamespace("topmodels"))

## Now the pdf method is registered which is used to approximate cdf/quantile.
## However, is_discrete (required) still missing.
expect_error(cdf(n1, x = 0),
             pattern = "^S3 method 'is_discrete' missing for object of class: pdfNormal$",
             info = "Method is_discrete is required, else an error must be thrown.")
## ... registering is_discrete (required for tinytest)
registerS3method("is_discrete", "pdfNormal", is_discrete.pdfNormal, envir = asNamespace("topmodels"))

## Same (error) should happen as support is missing
expect_error(cdf(n1, x = 0),
             pattern = "^S3 method 'support' missing for object of class: pdfNormal$",
             info = "Method support is required, else an error must be thrown.")
## ... registering support (required for tinytest)
registerS3method("support", "pdfNormal", support.pdfNormal, envir = asNamespace("topmodels"))
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
expect_silent(n3 <- setNames(pdfNormal(mu = 1:3, sigma = c(1, 2.5, 5)), LETTERS[1:3]),
    info = "Creator function must be silent")
expect_identical(n3,
    structure(data.frame(mu = 1:3, sigma = c(1, 2.5, 5)),
              class = c("pdfNormal", "distribution"),
              row.names = LETTERS[1:3]),
    info = "Object returned by pdfNormal not as expected.")
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
expect_identical(pdf(n1, x = 2), pdf(N1, x = 2),
        info = "pdf of pdfNormal and Normal must be identical.")
expect_identical(pdf(n3, x = 2), pdf(N3, x = 2),
        info = "pdf of pdfNormal and Normal must be identical.")
expect_identical(pdf(n3, x = 2, drop = FALSE), pdf(N3, x = 2, drop = FALSE),
        info = "pdf of pdfNormal and Normal must be identical.")
expect_identical(pdf(n3, x = 1:5), pdf(N3, x = 1:5),
        info = "pdf of pdfNormal and Normal must be identical.")
expect_identical(pdf(n3, x = 1:5, log = TRUE), pdf(N3, x = 1:5, log = TRUE),
        info = "pdf of pdfNormal and Normal w/ log = TRUE must be identical.")

# --------------- cdf ---------------
# The cdf is based on 'stats::integrate' and should (numerically) give us the
# same results.
expect_equal(cdf(n1, x = 2), cdf(N1, x = 2),
        info = "pdf of pdfNormal and Normal should be nearly identical.")
expect_equal(cdf(n1, x = 1:5), cdf(N1, x = 1:5),
        info = "pdf of pdfNormal and Normal should be nearly identical.")

expect_equal(cdf(n3, x = 3), cdf(N3, x = 3),
        info = "pdf of pdfNormal and Normal should be nearly identical.")
expect_equal(cdf(n3, x = c(1, 5, 10)), cdf(N3, x = c(1, 5, 10)),
        info = "pdf of pdfNormal and Normal should be nearly identical.")

expect_identical(cdf(n3, x = 2), cdf(n3, x = rep(2, 3L)),
        info = "Argument 'x' not properly replicated.")

# Testing return objects when using w/ and w/o drop, elementwise, ...
expect_equal(cdf(n1, x = 2, drop = FALSE), cdf(N1, x = 2, drop = FALSE))
expect_equal(cdf(n3, x = 1:3, drop = FALSE), cdf(N3, x = 1:3, drop = FALSE))
expect_equal(cdf(n3, x = 1:3), cdf(N3, x = 1:3))
expect_equal(cdf(n3, x = 1:3, drop = FALSE, elementwise = TRUE), cdf(N3, x = 1:3, drop = FALSE, elementwise = TRUE))

# lower.tail
expect_equal(cdf(n3, x = 1:5, lower.tail = FALSE), cdf(N3, x = 1:5, lower.tail = FALSE),
        info = "cdf of pdfNormal using lower.tail = FALSE should be equal.")

# ------------- quantile ------------
# Approximating the quantile function based on the PDF for continuous distributions
# is currently not implemented; expecting error message.
expect_error(quantile(n3, probs = 0.5),
             pattern = "^approximation for quantile function via pdf for continuous distributions currently not possible$",
             info = "Expected an error message as numeric pdf->quantile for continuous distributions not implemented.")

### p <- c(0, 0.00001, 0.001, seq(0.01, 0.99, by = 0.01), 0.999, 0.99999, 1)
### expect_equal(quantile(n1, probs = p), quantile(N1, probs = p),
###         info = "quantile of pdfNormal and Normal should be nearly identical.")
### 
### expect_equal(quantile(n3, probs = c(0.25, 0.5, 0.75)),
###              quantile(N3, probs = c(0.25, 0.5, 0.75)),
###         info = "quantile of pdfNormal and Normal should be nearly identical.")
### expect_equal(quantile(n3, probs = c(0.25, 0.5, 0.75), drop = FALSE),
###              quantile(N3, probs = c(0.25, 0.5, 0.75), drop = FALSE),
###         info = "quantile of pdfNormal and Normal should be nearly identical.")
### expect_equal(quantile(n3, probs = c(0.25, 0.5, 0.75), elementwise = FALSE),
###              quantile(N3, probs = c(0.25, 0.5, 0.75), elementwise = FALSE),
###         info = "quantile of pdfNormal and Normal should be nearly identical.")


# ------- random number gen ---------
# As the approximation for the quantile function (based on the pdf) is not yet
# implemented for continuous distributions, random() will not work (is based
# on the quantile function); expecting the same error message as for quantile(...) for now.
expect_error(random(n3, 3L),
             pattern = "^approximation for quantile function via pdf for continuous distributions currently not possible$",
             info = "Expected an error message as numeric pdf->quantile currently missing for continuous distributions (which random is based on).")

### expect_silent(r1 <- random(n1),                   info = "random(n1) should run silent.")
### expect_inherits(r1, "numeric",                    info = "random(n1) should return numeric.")
### expect_identical(length(r1), 1L,                  info = "random(n1) should return single numeric.")
### expect_null(names(r1),                            info = "random(n1) should return unnamed vector.")
### rm(r1)
### 
### expect_silent(r1 <- random(n1, drop = FALSE),     info = "random(n1, drop = FALSE) should run silent.")
### expect_true(is.numeric(r1),                       info = "random(n1, drop = FALSE) should return numeric.")
### expect_inherits(r1, "matrix",                     info = "random(n1, drop = FALSE) should return matrix.")
### expect_identical(dim(r1), c(1L, 1L),              info = "random(n1, drop = FALSE) should return matrix of dimension 1x1.")
### expect_identical(dimnames(r1), list(NULL, "r_1"), info = "dimension names of random(n1, drop = FALSE) incorrect.")
### rm(r1)
### 
### expect_silent(r1 <- random(n1, 50),               info = "random(n1, 50) should run silent.")
### expect_inherits(r1, "numeric",                    info = "random(n1) should return numeric.")
### expect_identical(length(r1), 50L,                 info = "random(n1) should return numeric of length 50L.")
### expect_null(names(r1),                            info = "random(n1) should return unnamed vector.")
### rm(r1)
### 
### expect_silent(r3 <- random(n3),                   info = "random(n3) should run silent.")
### expect_inherits(r3, "numeric",                    info = "random(n3) should return numeric.")
### expect_identical(length(r3), 3L,                  info = "random(n3) should return numeric of length 3.")
### expect_identical(names(r3), names(n3),            info = "random(n3) should return named vector.")
### rm(r3)
### 
### expect_silent(r3 <- random(n3, 5, drop = FALSE),  info = "random(n3, 5, drop = FALSE) should run silent.")
### expect_true(is.numeric(r3),                       info = "random(n3, 5, drop = FALSE) should return numeric.")
### expect_inherits(r3, "matrix",                     info = "random(n3, 5, drop = FALSE) should return matrix.")
### expect_identical(dim(r3), c(3L, 5L),              info = "random(n3, 5, drop = FALSE) should return matrix of dimension 3x5.")
### expect_identical(dimnames(r3), list(names(n3), paste("r", 1:5, sep = "_")), info = "dimension names of random(n3, 5, drop = FALSE) incorrect.")
### rm(r3)

# --------- central moments ---------
# As the approximation for the quantile function (based on the pdf) is not yet
# implemented for continuous distributions, numeric approximation of central
# moments not yet possible; expecting errors.

expect_error(mean(n1),
             pattern = "^approximation for quantile function via pdf for continuous distributions currently not possible$",
             info = "Currently missing pdf->quantile, thus mean.distribution should throw error.")
expect_error(variance(n1),
             pattern = "^approximation for quantile function via pdf for continuous distributions currently not possible$",
             info = "Currently missing pdf->quantile, thus variance.distribution should throw error.")
expect_error(skewness(n1),
             pattern = "^approximation for quantile function via pdf for continuous distributions currently not possible$",
             info = "Currently missing pdf->quantile, thus skewness.distribution should throw error.")
expect_error(kurtosis(n1),
             pattern = "^approximation for quantile function via pdf for continuous distributions currently not possible$",
             info = "Currently missing pdf->quantile, thus kurtosis.distribution should throw error.")

### # Testing mean
### 
### # default method falls to 'cdf' as we have no quantile function for MyMean
### expect_equal(mean(n3), mean(N3),
###         info = "mean of pdfNormal and Normal should be nearly identical")
### expect_equal(mean(n3, gridsize = 50), mean(n3, gridsize = 50, method = "cdf"),
###         tolerance = 1e-2,
###         info = "Expected mean to use method='cdf' by default, but got different results.")
### # Testing if method = "quantile" gives similar results
### expect_equal(mean(n3, gridsize = 50, method = "quantile"), mean(N3),
###         tolerance = 1e-2,
###         info = "Results using mean(..., method = \"quantile\") are not correct.")
### 
### # Testing variance
### 
### # default method falls to 'cdf' as we have no quantile function for MyMean
### expect_equal(variance(n3), variance(N3),
###         tolerance = 1e-2,
###         info = "variance of pdfNormal and Normal should be nearly identical")
### expect_equal(variance(n3, gridsize = 50), variance(n3, gridsize = 50, method = "cdf"),
###         tolerance = 1e-2,
###         info = "Expected variance to use method='cdf' by default, but got different results.")
### # Testing if method = "quantile" gives similar results
### expect_equal(variance(n3, gridsize = 50, method = "quantile"), variance(N3),
###         tolerance = 1e-1,
###         info = "Results using variance(..., method = \"quantile\") are not correct.")
### 
### # Testing skewness
### 
### # default method falls to 'cdf' as we have no quantile function for MyMean
### expect_equal(skewness(n3), skewness(N3),
###         info = "skewness of pdfNormal and Normal should be nearly identical")
### expect_equal(skewness(n3, gridsize = 50), skewness(n3, gridsize = 50, method = "cdf"),
###         info = "Expected skewness to use method='cdf' by default, but got different results.")
### # Testing if method = "quantile" gives similar results
### expect_equal(skewness(n3, gridsize = 50, method = "quantile"), skewness(N3),
###         info = "Results using skewness(..., method = \"quantile\") are not correct.")
### 
### # Testing kurtosis
### 
### # default method falls to 'cdf' as we have no quantile function for MyMean
### expect_equal(kurtosis(n3), kurtosis(N3),
###         tolerance = 1e-1,
###         info = "kurtosis of pdfNormal and Normal should be nearly identical")
### expect_equal(kurtosis(n3, gridsize = 50), kurtosis(n3, gridsize = 50, method = "cdf"),
###         info = "Expected kurtosis to use method='cdf' by default, but got different results.")
### # Testing if method = "quantile" gives similar results
### expect_equal(kurtosis(n3, gridsize = 50, method = "quantile"), kurtosis(N3),
###         tolerance = 1e-1,
###         info = "Results using kurtosis(..., method = \"quantile\") are not correct.")

