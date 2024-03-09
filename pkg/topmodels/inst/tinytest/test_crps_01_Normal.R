# -------------------------------------------------------
# Checking CRPS.distribution for Normal distribution
# against the analytical solution from scoringRules
# -------------------------------------------------------

if (interactive()) library("tinytest")

library("scoringRules")
library("distributions3")


# -------------------------------------------------------
# Single distrubiton, evaluate at many observations 'x'
# -------------------------------------------------------
set.seed(1234)

x <- runif(500, 4 - 50, 4 + 50)               # Sequnece where to evaluate the CRPS

expect_silent(d <- Normal(4, 2),                     info = "Setting up Normal(4, 2)")
expect_true(inherits(d, "distribution"),             info = "Checking return type of Normal(4, 2)")
crps_analyt <- crps(d, x)                            # Ground truth

expect_silent(crps <- crps.distribution(d, x),
              info = "Silent execution of crps.distribution(d, x)")
expect_true(is.vector(crps),                         info = "Return is vector")
expect_true(is.double(crps),                         info = "Return is double")
expect_identical(length(crps), length(x),            info = "Vector of length(x)")
expect_true(all(abs(crps_analyt - crps) < 0.005),    info = "Error tolerance < 0.005")
expect_equal(crps_analyt, crps, tol = 0.005,         info = "Result within tolerance")

# Drop = FALSE
expect_silent(crpsM <- crps.distribution(d, x, drop = FALSE),
              info = "Silent execution of crps.distribution(d, x, drop = FALSE)")
expect_true(is.matrix(crpsM),                        info = "Return is matrix")
expect_true(is.double(crpsM),                        info = "Return is double")
expect_identical(dim(crpsM), c(1L, 500L),            info = "Matrix of dimension (1 x 500)")
expect_equal(crps_analyt, as.numeric(crpsM), tol = 0.005, info = "Result within tolerance")

# -------------------------------------------------------
# Multiple (50) distributions, evaluate at one specific observation each
# -------------------------------------------------------
set.seed(1234)
mu <- rnorm(50, -10, 10)
sd <- exp(runif(50, -5, 5))
x  <- mu + rnorm(50, sd = 10)

expect_silent(d <- Normal(mu, sd),                   info = "Setting up Normal(mu, sd)")
expect_true(inherits(d, "distribution"),             info = "Checking return type of Normal(mu, sd)")
crps_analyt <- crps(d, x)                            # Ground truth

expect_silent(crps <- crps.distribution(d, x),
              info = "Silent execution of crps.distribution(d, x)")
expect_true(is.vector(crps),                         info = "Return is vector")
expect_true(is.double(crps),                         info = "Return is double")
expect_identical(length(crps), length(d),            info = "Vector of length(d)")
expect_equal(crps_analyt, crps, tol = 0.005,         info = "Results within tolerance")


# -------------------------------------------------------
# 20 Distributons; 20 observations (elementwise = TRUE))
# 20 Distributions; 3 observations (elementwise = FALSE)
# Testing batching and parallelization
# -------------------------------------------------------
set.seed(1234)
mu <- rnorm(20, -10, 10)
sd <- exp(runif(20, -5, 5))
x  <- mu + rnorm(20, sd = 10)

expect_silent(d <- Normal(mu, sd),                   info = "Setting up Normal(mu, sd)")
expect_true(inherits(d, "distribution"),             info = "Checking return type of Normal(mu, sd)")
crps_analyt <- crps(d, x)                            # Ground truth

# With elementwise = TRUE mode (auto)
expect_silent(crps <- crps.distribution(d, x, batchsize = 5),
              info = "crps.distribution(d, x), batchsize = 5 (must crate 4 batches)")
expect_true(is.vector(crps),                         info = "Return is vector")
expect_true(is.double(crps),                         info = "Return is double")
expect_identical(length(crps), length(d),            info = "Vector of length(d)")
expect_equal(crps_analyt, crps, tol = 0.005,         info = "Results within tolerance")

if (require("parallel")) {
    expect_silent(crps <- crps.distribution(d, x, cores = 2),
                  info = "crps.distribution(d, x), cores = 2 (must parallelize + batch)")
    expect_true(is.vector(crps),                       info = "Return is vector")
    expect_true(is.double(crps),                       info = "Return is double")
    expect_identical(length(crps), length(d),          info = "Vector of length(d)")
    expect_equal(crps_analyt, crps, tol = 0.005,       info = "Results within tolerance")
}

# With elementwise = FALSE mode (auto)
x <- c(-10, 0, 10)
crps_analyt <- crps(d, x)                     # Ground truth
expect_silent(crps <- crps.distribution(d, x, batchsize = 5),
              info = "crps.distribution(d, x), batchsize = 5 (must crate 4 batches)")
expect_true(is.matrix(crps),                          info = "Return is matrix")
expect_true(is.double(crps),                          info = "Return is double")
expect_identical(dim(crps_analyt), dim(crps),         info = "Return dimension")
expect_true(all(abs(crps_analyt - crps) < 0.005),     info = "Results within tolerance")
# TODO(R): Currently names don't match thus not equal
#expect_equal(crps_analyt, crps, tol = 0.005,          info = "Results within tolerance")

if (require("parallel")) {
    expect_silent(crps <- crps.distribution(d, x, cores = 2),
                  info = "crps.distribution(d, x), cores = 2 (must parallelize + batch)")
    expect_true(is.matrix(crps),                          info = "Return is matrix")
    expect_true(is.double(crps),                          info = "Return is double")
    expect_identical(dim(crps_analyt), dim(crps),         info = "Return dimension")
    expect_true(all(abs(crps_analyt - crps) < 0.005),     info = "Results within tolerance")
    # TODO(R): Currently names don't match thus not equal
    #expect_equal(crps_analyt, crps, tol = 0.005,          info = "Results within tolerance")
}
