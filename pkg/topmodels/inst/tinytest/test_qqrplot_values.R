# --------------------------------------------------------------------
# Testing usage of the `qqrplot()` return values (calculation)
# --------------------------------------------------------------------

if (interactive()) { rm(list = objects()); library("devtools"); library("tinytest"); library("topmodels") }

library("digest")
library("crch")
library("ggplot2")
library("tibble")

# --------------------------------------------------------------------
# Setting up the data sets/models used to test the function
# --------------------------------------------------------------------
data("CrabSatellites", package = "countreg")
CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]

# Different regression models (lm, censored lm, poisson count data model)
expect_silent(m1 <- lm(dist ~ speed, data = cars))
expect_silent(m2 <- crch(dist ~ speed | speed, left = 3, data = cars))
expect_silent(m3 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson))

# test_qqrplot_usage.R compares objects returned when using class = "data.frame"
# and class = "tibble", thus it is enough to test only one; sticking to data.frame here.
seed <- 6020

set.seed(seed)
expect_silent(q1u <- qqrplot(m1, class = "data.frame", plot = FALSE, scale = "uniform"))
set.seed(seed)
expect_silent(q2u <- qqrplot(m2, class = "data.frame", plot = FALSE, scale = "uniform"))
set.seed(seed)
expect_silent(q3u <- qqrplot(m3, class = "data.frame", plot = FALSE, scale = "uniform"))

set.seed(seed)
expect_silent(q1n <- qqrplot(m1, class = "data.frame", plot = FALSE, scale = "normal"))
set.seed(seed)
expect_silent(q2n <- qqrplot(m2, class = "data.frame", plot = FALSE, scale = "normal"))
set.seed(seed)
expect_silent(q3n <- qqrplot(m3, class = "data.frame", plot = FALSE, scale = "normal"))


# --------------------------------------------------------------------
# All objects must have identical names (identical order)
# --------------------------------------------------------------------
expect_identical(names(q1u), names(q2u), info = "testing names of return")
expect_identical(names(q1u), names(q3u), info = "testing names of return")
expect_identical(names(q1u), names(q1n), info = "testing names of return")
expect_identical(names(q1u), names(q2n), info = "testing names of return")
expect_identical(names(q1u), names(q3n), info = "testing names of return")


# --------------------------------------------------------------------
# Dimension of the objects; must match the original data set
# --------------------------------------------------------------------
expect_identical(NROW(cars), NROW(q1u),           info = "testing dimension 1 (number of rows/observations)")
expect_identical(NROW(cars), NROW(q1n),           info = "testing dimension 1 (number of rows/observations)")
expect_identical(NROW(cars), NROW(q2u),           info = "testing dimension 1 (number of rows/observations)")
expect_identical(NROW(cars), NROW(q2n),           info = "testing dimension 1 (number of rows/observations)")
expect_identical(NROW(CrabSatellites), NROW(q3u), info = "testing dimension 1 (number of rows/observations)")
expect_identical(NROW(CrabSatellites), NROW(q3n), info = "testing dimension 1 (number of rows/observations)")


# --------------------------------------------------------------------
# As we are using a fixed seed; check md5 sum of the returned objects.
# In case a test fails it is not precise but shows that there are differences.
# TODO: (RS2ML) Not sure if this works across different OSs!!!
# --------------------------------------------------------------------
expect_identical(digest(q1u), "641a9d725b724d97757c7a1b83ca6883")
expect_identical(digest(q1n), "950d8cf9f6064b7368698fd44c69fbb6")
expect_identical(digest(q2u), "9429373875f9dbb8286a4d9308d4ef0a")
expect_identical(digest(q2n), "4eb11af5b60908c41a0f7b00c6de6df9")
expect_identical(digest(q3u), "89a51471444630b56275d60dc9db5d5a")
expect_identical(digest(q3n), "69ee4fb089738de8a03899c47922f063")

# Thus, altenratively, gong for colSums(abs()) for all finite values
# which is a much weaker test.
fn <- function(x) as.numeric(sapply(x, function(y) ifelse(all(is.na(y)), NA, sum(abs(y[is.finite(y)])))))
expect_equal(fn(q1u), c(23.8823086240537, 25, NA, NA, NA))
expect_equal(fn(q1n), c(38.4240155605576, 39.6761730673991, NA, NA, NA))
expect_equal(fn(q2u), c(23.6139646074072, 25, NA, NA, NA))
expect_equal(fn(q2n), c(38.5067890113182, 39.6761730673991, NA, NA, NA))
expect_equal(fn(q3u), c(74.4564233315605, 86.5, 72.0636565265156, 78.8877997788859, 86.5))
expect_equal(fn(q3n), c(241.192716338567, 137.831905004867, 251.110835542692, 239.737769513666, 137.831905004867))
rm(fn)


# --------------------------------------------------------------------
# Per definition, uniform quantiles (observed and expected) must
# be within [0, 1].
# --------------------------------------------------------------------
expect_true(all(sapply(list(q1u, q2u, q3u), function(x) all(x$observed >= 0 & x$expected <= 1))),
            info = "testing for observed and expected quantiles for uniform distribution to lie in [0,1]")


# LINEAR MODEL ----------------------------
# Calculate residuals with ML corrected standard deviation (scale)
m1_eps   <- residuals(m1)
m1_scale <- sqrt(1 / length(m1_eps) * sum((m1_eps - mean(m1_eps))^2))
m1_prob  <- pnorm(model.response(model.frame(m1)),
                  mean = fitted(m1),
                  sd   = m1_scale)
expect_equivalent(q1n$observed, qnorm(m1_prob))
expect_equivalent(q1u$observed, m1_prob)

# NOTE: (RS2ML) Comparing ordered qq-(r)-residuals does not work as
#       observations are not unique, thus the order itself is not unique.
##expect_equivalent(q1u$expected, ppoints(length(eps))[order(order(eps))])
##plot(q1u$expected, ppoints(length(eps))[order(order(eps))])
expect_equivalent(sort(q1n$expected), qnorm(ppoints(length(m1_eps))))
expect_equivalent(sort(q1u$expected), ppoints(length(m1_eps)))

rm(m1_eps, m1_scale, m1_prob)


# CRCH MODEL ------------------------------
# Same for the censored model
m2_eps  <- residuals(m2, type = "standardized")
m2_prob <- pcnorm(model.response(model.frame(m2)),
                mean = fitted(m2, "location"),
                sd   = fitted(m2, "scale"),
                left = m2$cens$left, right = m2$cens$right)
expect_equal(q2n$observed, qnorm(m2_prob))
expect_equal(sort(q2n$expected), qnorm(ppoints(length(m2_eps))))

expect_equal(q2u$observed, m2_prob)
expect_equal(sort(q2u$expected), qunif(ppoints(length(m2_eps))))

rm(m2_eps, m2_prob)


# POISSON MODEL ---------------------------
m3_eps   <- residuals(m3)
m3_delta <- min(diff(sort(unique(model.response(model.frame(m3)))))) / 5e6
m3_prob1 <- ppois(model.response(model.frame(m3)) - m3_delta, lambda = fitted(m3))
m3_prob2 <- ppois(model.response(model.frame(m3)),            lambda = fitted(m3))

# Oberved and expected are based on randomized quantile reiduals.
# Observed but must lie in ppois(y - delta | lambda) and ppois(y | lambda).
expect_true(all(q3n$observed >= qnorm(m3_prob1) & q3n$observed <= qnorm(m3_prob2)))
expect_true(all(q3u$observed >= m3_prob1        & q3u$observed <= m3_prob2))

# TODO: (RS2ML) How to test the expected randomized quantile residuals?
#       How to test the simulated confidenc intervals?
####### m3_exp <- ppoints(length(m3_eps))[order(order(m3_eps))]
####### expect_true(all(abs(q3u$expected - m3_exp) < q3u$simint_expected))




