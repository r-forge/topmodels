# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `qresiduals.R`
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test `qresiduals.default()` with default arguments
# --------------------------------------------------------------------
m <- lm(dist ~ speed, data = cars)

expect_equal(
  qresiduals(m),
  qnorm(pnorm(setNames(cars$dist, rownames(cars)),
    mean = fitted(m), sd = sqrt(mean(residuals(m)^2))))
)


# --------------------------------------------------------------------
# Test `qresiduals.default()` with censor point
# --------------------------------------------------------------------
suppressPackageStartupMessages(require("crch"))
m2 <- crch(dist ~ speed | speed, left = 30, data = cars)

expect_equal(
  qresiduals(m2),
  qnorm(predict(m2, type = "probability", at = cars$dist))
)

# TODO: Improve/extend tests for censoring/truncation
