# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `newresponse.R`
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test for `newresponse.default()`
# --------------------------------------------------------------------
## Test for numeric response w/o NAs
m1 <- lm(dist ~ speed, data = cars)

expect_equal(
  as.numeric(newresponse(m1)$dist),
  cars$dist
)

expect_equal(
  as.numeric(newresponse(m1, newdata = cars[1:3, ])$dist),
  cars$dist[1:3]
)


## Test for numeric response w/ NAs in covariates
nd1 <- cars
nd1[c(3, 5, 7), "speed"] <- NA
expect_equal(
  as.numeric(newresponse(m1, newdata = nd1)$dist),
  cars$dist
)

## Test for numeric response w/ NAs in response
nd2 <- cars
nd2[c(3, 5, 7), "dist"] <- NA
expect_equal(
  as.numeric(newresponse(m1, newdata = nd2)$dist),
  nd2$dist
)


# --------------------------------------------------------------------
# Test for `newresponse.glm()`
# --------------------------------------------------------------------
## Test for numeric response w/o NAs
m2 <- glm(dist ~ speed, data = cars)

expect_equal(
  as.numeric(newresponse(m2)$dist),
  cars$dist
)

expect_equal(
  as.numeric(newresponse(m2, newdata = cars[1:3, ])$dist),
  cars$dist[1:3]
)


## Test for numeric response w/ NAs in covariates
nd1 <- cars
nd1[c(3, 5, 7), "speed"] <- NA
expect_equal(
  as.numeric(newresponse(m2, newdata = nd1)$dist),
  cars$dist
)

## Test for numeric response w/ NAs in response
nd2 <- cars
nd2[c(3, 5, 7), "dist"] <- NA
expect_equal(
  as.numeric(newresponse(m2, newdata = nd2)$dist),
  nd2$dist
)
