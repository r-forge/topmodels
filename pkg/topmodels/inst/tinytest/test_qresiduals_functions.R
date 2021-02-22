# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `qresiduals.R`
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test `qresiduals.default()` with default arguments
# --------------------------------------------------------------------
m <- lm(dist ~ speed, data = cars)

expect_equal(
  as.numeric(qresiduals(m)),  # TODO: (ML) Should it really be a named vector?
  qnorm(pnorm(
    cars$dist, 
    m$fitted.values, 
    summary(m)$sigma * sqrt(df.residual(m) / nobs(m))
  ))
)


# --------------------------------------------------------------------
# Test `qresiduals.default()` with censor point
# --------------------------------------------------------------------
library("crch")
m2 <- crch(dist ~ speed | speed, left = 30, data = cars)
idx <- which(cars$dist <= 30)

expect_equal(
  qresiduals(m2)[-idx],  # TODO: (ML) Should it really be a named vector?
  qnorm(predict(m2, type = "probability", at = cars$dist))[-idx]
)

# TODO: Improve/extent tests for censoring/truncation


