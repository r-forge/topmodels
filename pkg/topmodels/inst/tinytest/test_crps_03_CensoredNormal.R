# -------------------------------------------------------
# Checking CRPS.distribution for CensoredNormal distribution
# against the analytical solution from scoringRules.
# Requires 'crch' (for CensoredNormal distribution).
# -------------------------------------------------------

if (interactive()) library("tinytest")

library("scoringRules")
library("distributions3")

if (require("crch")) {

    # -------------------------------------------------------
    # Single distrubiton, evaluate at many observations 'x'
    # -------------------------------------------------------
    set.seed(1234)
    expect_silent(d <- CensoredNormal(4, 2, left = 0),
                  info = "Setting up CensoredNormal(4, 2, left = 0)")

    x <- runif(10, 4 - 50, 4 + 50)                   # Sequnece where to evaluate the CRPS
    crps_analyt <- crps(d, x)                        # Ground truth

    expect_true(inherits(d, "distribution"),
                info = "Checking return type of CensoredNormal(4, 2, left = 0)")
    expect_silent(crps <- topmodels:::crps.distribution(d, x),
                info = "Silent execution of crps.distribution(d, x)")
    expect_true(is.vector(crps),                              info = "Return is vector")
    expect_true(is.double(crps),                              info = "Return is double")
    expect_identical(length(crps), length(x),                 info = "Vector of length(x)")
    expect_equal(crps_analyt, crps, tol = 0.005,              info = "Result within tolerance")

    # Drop = FALSE
    expect_silent(crpsM <- topmodels:::crps.distribution(d, x, drop = F),
                  info = "Silent execution of crps(d, x, drop = FALSE)")
    expect_true(is.matrix(crpsM),                             info = "Return is matrix")
    expect_true(is.double(crpsM),                             info = "Return is double")
    expect_identical(dim(crpsM), c(1L, 10L),                  info = "Matrix of dimension (1 x 10)")
    expect_equal(crps_analyt, as.numeric(crpsM), tol = 0.005, info = "Result within tolerance")


    # -------------------------------------------------------
    # Multiple (50) distributions, evaluate at one specific observation each
    # -------------------------------------------------------
    set.seed(1234)
    mu <- rnorm(100, -10, 10)
    sd <- exp(runif(100, -5, 5))
    x  <- mu + rnorm(100, sd = 10)
    left  <- sample(c(-Inf, -Inf, -5, 0, 5), 100, replace = TRUE)
    right <- sample(c(Inf, Inf, 10, 15, 20), 100, replace = TRUE)

    ##table(left, right)
    expect_silent(d <- CensoredNormal(mu, sd, left = left, right = right),
                  info = "Setting up CensoredNormal()")
    expect_true(inherits(d, "distribution"),                  info = "Checking return type of Normal(mu, sd)")
    crps_analyt <- crps(d, x)                                 # Ground truth

    expect_silent(crps <- topmodels:::crps.distribution(d, x),
                  info = "Silent execution of crps.distribution(d, x)")
    expect_true(is.vector(crps),                              info = "Return is vector")
    expect_true(is.double(crps),                              info = "Return is double")
    expect_identical(length(crps), length(d),                 info = "Vector of length(d)")
    expect_equal(crps_analyt, crps, tol = 0.005,              info = "Result within tolerance")
}
