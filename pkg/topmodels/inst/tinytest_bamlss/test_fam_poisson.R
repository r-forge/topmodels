# --------------------------------------------------------------------
# Testing bamlss family distrubiton3 functionality
# poisson_bamlss: Poisson count data distribution
# --------------------------------------------------------------------

if (interactive()) { library("tinytest") }

suppressPackageStartupMessages(library("topmodels"))
suppressPackageStartupMessages(library("distributions3"))
suppressPackageStartupMessages(library("bamlss"))


# --------------------------------------------------------------------
# Set of parameters for testing
# --------------------------------------------------------------------
par <- list(lambda = c(0, 0.1, 2.82134, 4, 6, 20))
expect_identical(names(par), c("lambda"))
expect_true(length(par) == 1L && all(sapply(par, is.numeric)))

# Misuse
expect_error(bamlss::BAMLSS("poisson"),
             info = "No parameters")
expect_error(bamlss::BAMLSS("poisson", unlist(sapply(par, function(x) x[[1]]))),
             info = "Parameters not list")

# Testing basic functionality.
expect_silent(bamlss::BAMLSS("poisson", par))
expect_silent(bamlss::BAMLSS("poisson_bamlss", par))
expect_silent(bamlss::BAMLSS(bamlss.family("poisson"), par))
expect_silent(d <- bamlss::BAMLSS(bamlss.family("poisson_bamlss"), par))

# Testing return 'd'
expect_identical(class(d), c("BAMLSS", "distribution"))
expect_identical(dim(d), NULL)
expect_identical(length(d), length(par[[1]]))

# Leaving out one parameter
for (n in names(par)) {
    tmp <- par; tmp[[n]] <- NULL
    expect_error(bamlss::BAMLSS("poisson", tmp), info = "Parameter missing")
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
expect_true(is.integer(r) && length(r) == length(d) && is.null(dim(r)))
expect_identical(min(r), 0L)

# Probability density function
expect_equal(pdf(d, 3),             dpois(3, par$lambda))
expect_equal(pdf(d, 3, log = TRUE), dpois(3, par$lambda, log = TRUE))
expect_equal(log_pdf(d, 3),         dpois(3, par$lambda, log = TRUE))

expect_equal(pdf(d, 0),             dpois(0, par$lambda))
expect_equal(pdf(d, 0, log = TRUE), dpois(0, par$lambda, log = TRUE))
expect_equal(log_pdf(d, 0),         dpois(0, par$lambda, log = TRUE))

# Cumulative distribution function
expect_equal(cdf(d, 3),             ppois(3,  par$lambda))
expect_equal(cdf(d, 0),             ppois(0,  par$lambda))
expect_equal(cdf(d, -1),            ppois(-1, par$lambda))

# Quantiles
expect_equal(quantile(d, 0.01), qpois(0.01, par$lambda))
expect_equal(quantile(d, 0.50), qpois(0.50, par$lambda))
expect_equal(quantile(d, 0.99), qpois(0.99, par$lambda))

# Support
names(d) <- LETTERS[seq_len(length(d))] # Used later
expect_silent(supp <- support(d))
expect_true(is.matrix(supp) && identical(dim(supp), c(length(d), 2L)))
expect_identical(colnames(supp), c("min", "max"))
expect_identical(rownames(supp), names(d))
expect_true(all(supp[, "min"] == 0))
expect_identical(unname(supp[, "max"]), ifelse(par$lambda == 0, 0L, Inf))

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
reg <- get_regex("poisson", par)
expect_true(all(grepl(reg, format(d))))
# Not testing print 


