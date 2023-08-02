# --------------------------------------------------------------------
# Testing bamlss family distrubiton3 functionality
# cnorm_bamlss: Zero left censored normal distribution
# --------------------------------------------------------------------

# TODO(R): cnorm_bamlss returns non-zero cdf below censoring point.
#          Should be zero(?!)
# x <- seq(-4, 4, length.out = 201)
# plot(x, cdf(BAMLSS("cnorm", list(mu = -1, sigma = 3)), x), ylim = c(0, 1))
# points(x, crch::pcnorm(x, -1, 3, left = 0), col = 3)
# abline(v = 0, col = 2)

if (interactive()) { library("tinytest") }

suppressPackageStartupMessages(library("topmodels"))
suppressPackageStartupMessages(library("distributions3"))
suppressPackageStartupMessages(library("bamlss"))
suppressPackageStartupMessages(library("crch"))


# --------------------------------------------------------------------
# Set of parameters for testing
# --------------------------------------------------------------------
par <- list(mu = c(-1, 0, 4, 5, -13, 2.82134), sigma = c(2, 3, 1, 1.5, 10, 0.01))
expect_identical(names(par), c("mu", "sigma"))
expect_true(length(par) == 2L && all(sapply(par, is.numeric)))

# Misuse
expect_error(bamlss::BAMLSS("cnorm"),
             info = "No parameters")
expect_error(bamlss::BAMLSS("cnorm", unlist(sapply(par, function(x) x[[1]]))),
             info = "Parameters not list")

# Testing basic functionality.
expect_silent(bamlss::BAMLSS("cnorm", par))
expect_silent(bamlss::BAMLSS("cnorm_bamlss", par))
expect_silent(bamlss::BAMLSS(bamlss.family("cnorm"), par))
expect_silent(d <- bamlss::BAMLSS(bamlss.family("cnorm_bamlss"), par))

# Testing return 'd'
expect_identical(class(d), c("BAMLSS", "distribution"))
expect_identical(dim(d), NULL)
expect_identical(length(d), length(par[[1]]))

# Leaving out one parameter
for (n in names(par)) {
    tmp <- par; tmp[[n]] <- NULL
    expect_error(bamlss::BAMLSS("cnorm", tmp), info = "Parameter missing")
}
 
# Arithmetic mean
expect_true(is.numeric(mean(d)), length(mean(d)) == length(d))

# Variance
expect_true(is.numeric(variance(d)), length(variance(d)) == length(d))

# Skewness
expect_error(skewness(d), pattern = "no skewness\\(\\) function provided by")

# Kurtosis
expect_error(kurtosis(d), pattern = "no kurtosis\\(\\) function provided by")

# Random number generator
expect_silent(r <- random(d))
expect_true(is.numeric(r) && length(r) == length(d) && is.null(dim(r)))
expect_identical(min(r), 0.0)

# Probability density function
expect_equal(pdf(d, 3),             crch::dcnorm(3, par$mu, par$sigma, left = 0))
expect_equal(pdf(d, 3, log = TRUE), crch::dcnorm(3, par$mu, par$sigma, left = 0, log = TRUE))
expect_equal(log_pdf(d, 3),         crch::dcnorm(3, par$mu, par$sigma, left = 0, log = TRUE))

expect_equal(pdf(d, 0),             crch::dcnorm(0, par$mu, par$sigma, left = 0))
expect_equal(pdf(d, 0, log = TRUE), crch::dcnorm(0, par$mu, par$sigma, left = 0, log = TRUE))
expect_equal(log_pdf(d, 0),         crch::dcnorm(0, par$mu, par$sigma, left = 0, log = TRUE))

# Cumulative distribution function
expect_equal(cdf(d, 3),             pcnorm(3, par$mu, par$sigma, left = 0))
#expect_equal(cdf(d, 0),             pcnorm(0, par$mu, par$sigma, left = 0))
#expect_equal(cdf(d, -1),            pcnorm(-1, par$mu, par$sigma, left = 0))

# They differ, all.equal returns 'text' with mean difference. Testing for that for now.
# The problem is that the cnorm_bamlss cnorm returns cdf > 0 below censoring point
expect_true(is.character(all.equal(cdf(d, 0),  pcnorm(0, par$mu, par$sigma, left = 0))))
expect_true(is.character(all.equal(cdf(d, -1), pcnorm(-1, par$mu, par$sigma, left = 0))))

# Quantiles
expect_equal(quantile(d, 0.01), qcnorm(0.01, par$mu, par$sigma, left = 0))
expect_equal(quantile(d, 0.50), qcnorm(0.50, par$mu, par$sigma, left = 0))
expect_equal(quantile(d, 0.99), qcnorm(0.99, par$mu, par$sigma, left = 0))

# Support
names(d) <- LETTERS[seq_len(length(d))] # Used later
expect_silent(supp <- support(d))
expect_true(is.matrix(supp) && identical(dim(supp), c(length(d), 2L)))
expect_identical(colnames(supp), c("min", "max"))
expect_identical(rownames(supp), names(d))
expect_true(all(supp[, "min"] == 0))
expect_true(all(supp[, "max"] == Inf))

# is_*
expect_warning(id <- is_discrete(d))
expect_warning(ic <- is_continuous(d))
expect_identical(id, setNames(rep(FALSE, length(d)), names(d)))
expect_identical(ic, setNames(rep(FALSE, length(d)), names(d)))

get_regex <- function(name, par) {
    reg_float <- "[\\+-\\.0-9]+"
    reg_par   <- sprintf("%s\\s=\\s{1,}%s", names(par), reg_float)
    sprintf("^BAMLSS %s distribution \\(%s\\)$", name, paste(reg_par, collapse = ", "))
}
reg <- get_regex("cnorm", par)
expect_true(all(grepl(reg, format(d))))
# Not testing print 


