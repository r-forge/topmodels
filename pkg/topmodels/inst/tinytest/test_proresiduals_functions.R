# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `proresiduals.R`
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test `proresiduals.default()` with continuous model
# --------------------------------------------------------------------
m <- lm(dist ~ speed, data = cars)

expect_equal(
  proresiduals(m),
  qnorm(pnorm(setNames(cars$dist, rownames(cars)),
    mean = fitted(m), sd = sqrt(mean(residuals(m)^2))))
)


# --------------------------------------------------------------------
# Test `proresiduals.default()` with discrete model
# --------------------------------------------------------------------
m2 <- glm(dist ~ speed, data = cars, family = poisson)

expect_equal(
  proresiduals(m2, prob = 1),
  qnorm(ppois(setNames(cars$dist, rownames(cars)), lambda = fitted(m2)))
)

expect_equal(
  proresiduals(m2, prob = 0),
  qnorm(ppois(setNames(cars$dist, rownames(cars)) - 1, fitted(m2)))
)

# --------------------------------------------------------------------
# Test `proresiduals.default()` with censored model
# --------------------------------------------------------------------
suppressPackageStartupMessages(require("crch"))
cars$dist30 <- pmax(cars$dist, 30)
m3 <- crch(dist30 ~ speed | speed, left = 30, data = cars)

expect_equal(
  proresiduals(m3, prob = 1),
  qnorm(predict(m3, type = "probability", at = cars$dist30))
)

expect_equal(
  proresiduals(m3, prob = 0),
  qnorm(predict(m3, type = "probability", at = cars$dist30 - sqrt(.Machine$double.eps)))
)
