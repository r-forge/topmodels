# -------------------------------------------------------
# Checking CRPS.distribution for Poisson distribution
# against the analytical solution from scoringRules
# -------------------------------------------------------

if (interactive()) library("tinytest")

library("scoringRules")
library("distributions3")


# -------------------------------------------------------
# Single distrubiton, evaluate at many observations 'x'
# -------------------------------------------------------
set.seed(1234)
x <- rpois(500, lambda = 22.5)                             # Sequnece where to evaluate the CRPS

expect_silent(d <- Poisson(22.5),                          info = "Setting up Poisson(22.5)")
expect_inherits(d, "distribution",                         info = "Checking return type of Poisson(22.5)(")
crps_analyt <- crps(d, x)                                  # Ground truth

expect_silent(crps <- topmodels:::crps.distribution(d, x),
              info = "Silent execution of crps.distribution(d, x)")
expect_true(is.vector(crps),                               info = "Return is vector")
expect_true(is.double(crps),                               info = "Return is double")
expect_identical(length(crps), length(x),                  info = "Matrix of length(x)")
expect_equal(crps_analyt, crps, tol = 0.005,               info = "Result within tolerance")

# Drop = FALSE
expect_silent(crpsM <- topmodels:::crps.distribution(d, x, drop = F),
              info = "Silent execution of crps(d, x, drop = FALSE)")
expect_true(is.matrix(crpsM),                              info = "Return is matrix")
expect_true(is.double(crpsM),                              info = "Return is double")
expect_identical(dim(crpsM), c(1L, 500L),                  info = "Matrix of dimension (1 x 500)")
expect_equal(crps_analyt, as.numeric(crpsM), tol = 0.005,  info = "Matrix of dimension (1 x 500)")
expect_equal(crps, as.numeric(crpsM),                      info = "Same result with drop = TRUE and drop = FALSE")

# -------------------------------------------------------
# Multiple (50) distributions, evaluate at one specific observation each
# -------------------------------------------------------
set.seed(1234)

lambda <- runif(50, 0.1, 7)^2
x      <- sample(0:100, 50, replace = TRUE)

expect_silent(d <- Poisson(lambda),                        info = "Setting up Poisson(lambda)")
expect_inherits(d, "distribution",                         info = "Checking return type of Poisson(lambda)")
crps_analyt <- crps(d, x)                                  # Ground truth

expect_silent(crps <- topmodels:::crps.distribution(d, x),
              info = "Silent execution of crps.distribution(d, x)")
expect_true(is.vector(crps),                               info = "Return is vector")
expect_true(is.double(crps),                               info = "Return is double")
expect_identical(length(crps), length(d),                  info = "Vector of length(d)")
expect_equal(crps_analyt, crps, tol = 0.005,               info = "Results within tolerance")


# -------------------------------------------------------
# Testing fallback to discrete approximation when
# grid is too large (0:max larger than m)
# -------------------------------------------------------
set.seed(1234)

expect_silent(d <- Poisson(c(1:2, 1000)),                  info = "Setting up Poisson(c(1:2, 1000))")
expect_inherits(d, "distribution",                         info = "Checking return type of Poisson(c(1:2, 10000))")
crps_analyt <- crps(d, 1:2)                   # Ground truth

expect_warning(crps <- topmodels:::crps.distribution(d, 1:2),
               pattern = "^grid size too large; falling back to continuous CRPS approximation$",
               info = "Show fallback warning")
expect_inherits(crps, "matrix",                            info = "Return is matrix")
expect_true(is.double(crps),                               info = "Return is double")
expect_identical(dim(crps), dim(crps_analyt),              info = "Matrix of dimension (2 x 3)")
expect_true(all(abs(crps_analyt - crps) < 0.005),          info = "Results within tolerance")
# TODO(R): Currently names don't match thus not equal
#expect_equal(crps_analyt, crps, tol = 0.005,               info = "Results within tolerance")
