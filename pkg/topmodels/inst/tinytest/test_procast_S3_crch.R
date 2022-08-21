# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `procast.R` FOR S3-CLASS `crch`
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test procast.crch argument = `type`
# --------------------------------------------------------------------
suppressPackageStartupMessages(require("crch"))
m1 <- crch(dist ~ speed | speed, data = cars)

expect_equal(
  procast(m1, type = "quantile", at = 0.5, drop = TRUE), 
  predict(m1, type = "quantile", at = 0.5),
  check.attributes = FALSE
)

expect_equal(
  procast(m1, type = "mean", drop = TRUE), 
  predict(m1, type = "response"),
  check.attributes = FALSE
)

expect_equal(
  procast(m1, type = "density", at = 1, drop = TRUE),
  predict(m1, type = "density", at = 1),
  check.attributes = FALSE
)

expect_equal(
  procast(m1, type = "probability", at = 1, drop = TRUE),
  predict(m1, type = "probability", at = 1),
  check.attributes = FALSE
)


# --------------------------------------------------------------------
# Test procast.crch argument = `drop`
# --------------------------------------------------------------------
expect_equal(
  procast(m1, type = "quantile", at = 0.5, drop = FALSE),
  data.frame(
    quantile = predict(m1, type = "quantile", at = 0.5)
  ),
  check.attributes = FALSE
)


# --------------------------------------------------------------------
# Test procast.crch argument = `at`
# --------------------------------------------------------------------
m2 <- crch(dist ~ speed | speed, data = cars)
nd <- data.frame(speed = c(10, 15, 20))

expect_equal(
  procast(m2, at = c(0.25, 0.5, 0.75), elementwise = FALSE),
  data.frame(
    q_0.25 = predict(m2, type = "quantile", at = 0.25),
    q_0.5  = predict(m2, type = "quantile", at = 0.50),
    q_0.75 = predict(m2, type = "quantile", at = 0.75)
  ), 
  check.attributes = FALSE
)

expect_equal(
  procast(m2, type = "quantile", at = c(0.25, 0.5, 0.75), elementwise = FALSE, drop = TRUE),
  predict(m2, type = "quantile", at = c(0.25, 0.5, 0.75)),
  check.attributes = FALSE
)

