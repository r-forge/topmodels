# --------------------------------------------------------------------
# Testing usage of the `qqrplot()` return values (calculation)
# --------------------------------------------------------------------

if (interactive()) { rm(list = objects()); library("devtools"); library("tinytest"); library("topmodels") }

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
# Per definition, uniform quantiles (observed and expected) must
# be within [0, 1].
# --------------------------------------------------------------------
expect_true(all(sapply(list(q1u, q2u, q3u), function(x) all(x$observed >= 0 & x$expected <= 1))),
            info = "testing for observed and expected quantiles for uniform distribution to lie in [0,1]")


# LINEAR MODEL ----------------------------
# Calculate residuals with ML corrected standard deviation (scale)
eps   <- residuals(m1)
##delta <- min(diff(sort(unique(model.response(model.frame(m1)))))) / 5e6
scale <- sqrt(1 / length(eps) * sum((eps - mean(eps))^2))
expect_equivalent(q1n$observed, residuals(m1) / scale)
expect_equivalent(q1u$observed, qunif(pnorm(residuals(m1) / scale)))

# TODO: (RS2ML) This fails as order(order(eps)) is not identical as two
#       data points are identical and the order is, thus, no longer unique.
##expect_equivalent(q1u$expected, ppoints(length(eps))[order(order(eps))])
##plot(q1u$expected, ppoints(length(eps))[order(order(eps))])

expect_equivalent(sort(q1n$expected), qnorm(ppoints(length(eps))))
expect_equivalent(sort(q1u$expected), ppoints(length(eps)))


## CRCH MODEL ------------------------------
## Same for the censored model
#eps   <- residuals(m2, type = "standardized")
#idx   <- which(y > m2$cens$left)
#expect_equivalent(q2n$observed[idx], eps[idx])                      # observed [non-censored]
#expect_equivalent(sort(q2n$expected), qnorm(ppoints(length(eps))))  # expected
## TODO: (RS2ML) This fails as order(order(eps)) is not identical as two
##       data points are identical and the order is, thus, no longer unique.
###expect_equivalent(q2n$expected, qnorm(ppoints(length(eps)))[order(order(eps))])  # expected
#plot(q2u$expected, ppoints(length(eps))[order(order(eps))])
#
## ------> da ist was mit meiner standardized order falsch oder? Ist das nicht einfach der CENS?
#
#
#expect_equivalent(q2n$expected, qnorm(ppoints(length(eps)))[order(order(eps))])  # expected
#
#plot(q2n$expected, qcnorm(ppoints(length(eps)), left = m2$cens$left)[order(order(eps))])  # expected
#k <- data.frame(q = q2n$expected, order = order(order(eps)), reto = qnorm(ppoints(length(eps)))[order(order(eps))])
#subset(k, order == 1)
#which(k$q < -2)
#k[1,] 
#plot(reto ~ q, data = k, col = ifelse(order == 1, 'red', 'black'))
#
#head(q2n)
#
#y     <- model.response(model.frame(m2))
#
#expect_equivalent(qunif(pnorm(eps)), q2u$observed) # TODO: (RS) Hm ..
#
#head(q2n)
#hist(q2n$expected)
#f <- sort(q2n$expected)









