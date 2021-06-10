# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `procast.R` FOR S3-CLASS `crch`
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test procast.crch argument = `type`
# --------------------------------------------------------------------
require("crch")
m1 <- crch(dist ~ speed | speed, data = cars)

expect_equal(
  procast(m1, drop = TRUE), 
  predict(m1, type = "quantile"),
  check.attributes = FALSE
)

expect_equal(
  procast(m1, type = "location", drop = TRUE), 
  predict(m1, type = "location"),
  check.attributes = FALSE
)

expect_equal(
  procast(m1, type = "scale", drop = TRUE),
  predict(m1, type = "scale"),
  check.attributes = FALSE
)

expect_equal(
  procast(m1, type = "parameter"),
  predict(m1, type = "parameter"),
  check.attributes = FALSE
)

expect_equal(
  procast(m1, type = "density", drop = TRUE),
  predict(m1, type = "density"),
  check.attributes = FALSE
)

expect_equal(
  procast(m1, type = "probability", drop = TRUE),
  predict(m1, type = "probability"),
  check.attributes = FALSE
)

expect_error(
  procast(m1, type = "score", drop = TRUE)
)


# --------------------------------------------------------------------
# Test procast.crch argument = `drop`
# --------------------------------------------------------------------
expect_equal(
  procast(m1, type = "quantile", drop = FALSE),
  data.frame(
    quantile = predict(m1, type = "quantile")
  ),
  check.attributes = FALSE
)

## TODO: (ML) Implement test for condition (!is.null(dim(rval)) && NCOL(rval) == 1L)


# --------------------------------------------------------------------
# Test procast.crch argument = `at`
# --------------------------------------------------------------------
m2 <- crch(dist ~ speed | speed, data = cars)
nd <- data.frame(speed = c(10, 15, 20))

expect_equal(
  procast(m2, at = c(0.25, 0.5, 0.75)),
  data.frame(
    q_0.25 = predict(m2, type = "quantile", at = c(0.25)),
    q_0.5  = predict(m2, type = "quantile", at = c(0.5)),
    q_0.75 = predict(m2, type = "quantile", at = c(0.75))
  ), 
  check.attributes = FALSE
)

expect_equal(
  procast(m2, at = rbind(c(0.25, 0.5, 0.75)), drop = TRUE),
  predict(m2, type = "quantile", at = c(0.25, 0.5, 0.75)),
  check.attributes = FALSE
)

qnt1 <- procast(m2, type = "quantile", newdata = nd, at = "function")
expect_equal(
  as.numeric(qnt1(c(0.25, 0.5, 0.75))),  # TODO: (ML) Should it really be a named vector?
  qnorm(
    c(0.25, 0.5, 0.75),
    predict(m2, newdata = nd, type = "location"),
    predict(m2, newdata = nd, type = "scale")
  )
)

expect_error(procast(m2, newdata = nd, at = "list"))

