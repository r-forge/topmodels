# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `procast.R` FOR S3-CLASS `lm`
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test procast.lm argument = `type`
# --------------------------------------------------------------------
m <- lm(dist ~ speed, data = cars)

expect_equal(
  procast(m, drop = TRUE), 
  qnorm(0.5, m$fitted.values, sd = summary(m)$sigma * sqrt(df.residual(m) / nobs(m)))
)

expect_equal(
  procast(m, type = "location", drop = TRUE), 
  m$fitted.values
)

expect_equal(
  unique(procast(m, type = "scale", drop = TRUE)),
  summary(m)$sigma * sqrt(df.residual(m) / nobs(m))
)

expect_equal(
  procast(m, type = "parameter"),
  data.frame(
    mu = m$fitted.values, 
    sigma = summary(m)$sigma * sqrt(df.residual(m) / nobs(m))
  )
)

expect_equal(
  procast(m, type = "density", drop = TRUE),
  dnorm(
    0.5,
    m$fitted.values, 
    summary(m)$sigma * sqrt(df.residual(m) / nobs(m))
  )
)

expect_equal(
  procast(m, type = "probability", drop = TRUE),
  pnorm(
    0.5,
    m$fitted.values, 
    summary(m)$sigma * sqrt(df.residual(m) / nobs(m))
  )
)

## TODO: (ML) Implement type = `score`


# --------------------------------------------------------------------
# Test procast.lm argument = `drop`
# --------------------------------------------------------------------
expect_equal(
  procast(m, drop = FALSE),
  data.frame(
    quantile = qnorm(0.5, m$fitted.values, sd = summary(m)$sigma * sqrt(df.residual(m) / nobs(m)))
  )
)

## TODO: (ML) Implement test for condition (!is.null(dim(rval)) && NCOL(rval) == 1L)


# --------------------------------------------------------------------
# Test procast.lm argument = `newdata`
# --------------------------------------------------------------------
m <- lm(dist ~ speed, data = cars)
nd <- data.frame(speed = c(10, 15, 20))

expect_equal(
  procast(m, type = "parameter", newdata = nd),
  data.frame(
    mu = predict(m, newdata = nd), 
    sigma = summary(m)$sigma * sqrt(df.residual(m) / nobs(m))
  )
)


# --------------------------------------------------------------------
# Test procast.lm argument = `at`
# --------------------------------------------------------------------
m <- lm(dist ~ speed, data = cars)
nd <- data.frame(speed = c(10, 15, 20))

expect_equal(
  procast(m, at = c(0.25, 0.5, 0.75)),
  data.frame(
    q_0.25 = qnorm(0.25, m$fitted.values, sd = summary(m)$sigma * sqrt(df.residual(m) / nobs(m))),
    q_0.5  = qnorm(0.5, m$fitted.values, sd = summary(m)$sigma * sqrt(df.residual(m) / nobs(m))),
    q_0.75 = qnorm(0.75, m$fitted.values, sd = summary(m)$sigma * sqrt(df.residual(m) / nobs(m)))
  )
)

expect_equal(
  procast(m, at = rbind(c(0.25, 0.5, 0.75))),
  data.frame(
    q_0.25 = qnorm(0.25, m$fitted.values, sd = summary(m)$sigma * sqrt(df.residual(m) / nobs(m))),
    q_0.5  = qnorm(0.5, m$fitted.values, sd = summary(m)$sigma * sqrt(df.residual(m) / nobs(m))),
    q_0.75 = qnorm(0.75, m$fitted.values, sd = summary(m)$sigma * sqrt(df.residual(m) / nobs(m)))
  )
)

qnt1 <- procast(m, newdata = nd, at = "function")
expect_equal(
  as.numeric(qnt1(c(0.25, 0.5, 0.75))),  # TODO: (ML) Should it really be a named vector?
  qnorm(
    c(0.25, 0.5, 0.75), 
    predict(m, newdata = nd),
    summary(m)$sigma * sqrt(df.residual(m) / nobs(m))
  )
)

qnt2 <- procast(m, newdata = nd, at = "list")
expect_equal(
  qnt2[[2]](c(0.25, 0.5, 0.75)),
  qnorm(
    c(0.25, 0.5, 0.75), 
    predict(m, newdata = nd)[2],
    summary(m)$sigma * sqrt(df.residual(m) / nobs(m))
  )
)


# --------------------------------------------------------------------
# Test procast.lm argument = `na.action`
# --------------------------------------------------------------------
m <- lm(dist ~ speed, data = cars)
nd <- data.frame(speed = c(10, NA, 20))

expect_equal(
  procast(m, type = "parameter", newdata = nd),
  data.frame(
    mu = predict(m, newdata = nd), 
    sigma = summary(m)$sigma * sqrt(df.residual(m) / nobs(m))
  )
)

expect_equal(
  procast(m, type = "parameter", newdata = nd, na.action = na.omit),
  data.frame(
    mu = predict(m, newdata = nd, na.action = na.omit), 
    sigma = summary(m)$sigma * sqrt(df.residual(m) / nobs(m))
  )
)