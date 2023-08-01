# --------------------------------------------------------------------
# Testing bamlss family distrubiton3 functionality
# --------------------------------------------------------------------

# TODO(R): Possible improvements
# - Currently families accept invalid parameters (e.g.,
#   the gaussian_bamlss accepts negative sigmas)
if (interactive()) { library("tinytest") }
source("fam_basetest.R")

suppressPackageStartupMessages(library("topmodels"))
suppressPackageStartupMessages(library("distributions3"))
suppressPackageStartupMessages(library("bamlss"))


# --------------------------------------------------------------------
# Set of parameters for testing
# --------------------------------------------------------------------
par <- list(mu = c(4, 5, 6, 2.82134), sigma = c(1, 1.5, 10, 0.01))
expect_identical(names(par), c("mu", "sigma"))
expect_true(length(par) == 2L && all(sapply(par, is.numeric)))

# --------------------------------------------------------------------
# gaussian_bamlss: Normal distribution
# --------------------------------------------------------------------

expect_error(bamlss::BAMLSS("gaussian"),
             info = "No parameters")
expect_error(bamlss::BAMLSS("gaussian", unlist(sapply(par, function(x) x[[1]]))),
             info = "Parameters not list")

# Testing basic functionality.
expect_silent(bamlss::BAMLSS("gaussian", par))
expect_silent(bamlss::BAMLSS("gaussian_bamlss", par))
expect_silent(bamlss::BAMLSS(bamlss.family("gaussian"), par))
expect_silent(d <- bamlss::BAMLSS(bamlss.family("gaussian_bamlss"), par))

# Testing return 'd'
expect_identical(class(d), c("BAMLSS", "distribution"))
expect_identical(dim(d), NULL)
expect_identical(length(d), length(par[[1]]))

# Leaving out one parameter
for (n in names(par)) {
    tmp <- par; tmp[[n]] <- NULL
    expect_error(bamlss::BAMLSS("gaussian", tmp), info = "Parameter missing")
}
 
# Arithmetic mean
expect_identical(mean(d), par$mu)

# Variance
expect_equal(variance(d), par$sigma^2)

# Skewness
expect_error(skewness(d))

# Kurtosis
expect_error(kurtosis(d))

# Random number generator
expect_silent(r <- random(d))
expect_true(is.numeric(r) && length(r) == length(d) && is.null(dim(r)))

# Probability density function
expect_equal(pdf(d, 3),             dnorm(3, par$mu, par$sigma))
expect_equal(pdf(d, 3, log = TRUE), dnorm(3, par$mu, par$sigma, log = TRUE))
expect_equal(log_pdf(d, 3),         dnorm(3, par$mu, par$sigma, log = TRUE))

# Cumulative distribution function
expect_equal(cdf(d, 3),             pnorm(3, par$mu, par$sigma))

# Quantiles
expect_equal(quantile(d, 0.01), qnorm(0.01, par$mu, par$sigma))
expect_equal(quantile(d, 0.50), qnorm(0.50, par$mu, par$sigma))
expect_equal(quantile(d, 0.99), qnorm(0.99, par$mu, par$sigma))

# Support
names(d) <- LETTERS[seq_len(length(d))] # Used later
expect_silent(supp <- support(d))
expect_true(is.matrix(supp) && identical(dim(supp), c(length(d), length(par))))
expect_identical(colnames(supp), c("min", "max"))
expect_identical(rownames(supp), names(d))
expect_true(all(supp[, "min"] == -Inf))
expect_true(all(supp[, "max"] == Inf))

# is_*
expect_warning(id <- is_discrete(d))
expect_warning(ic <- is_continuous(d))
expect_identical(id, setNames(rep(FALSE, length(d)), names(d)))
expect_identical(ic, setNames(rep(FALSE, length(d)), names(d)))

get_regex <- function(name, par) {
    reg_float <- "[\\+\\-\\.0-9]+"
    reg_par   <- sprintf("%s\\s=\\s{1,}%s", names(par), reg_float)
    sprintf("^BAMLSS %s distribution \\(%s\\)$", name, paste(reg_par, collapse = ", "))
}
reg <- get_regex("gaussian", par)
expect_true(all(grepl(reg, format(d))))
# Not testing print 


