# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `qresiduals.R'
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test `qresiduals.default()' with default arguments
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


