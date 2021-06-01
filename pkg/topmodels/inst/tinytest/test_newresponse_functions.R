# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `newresponse.R`
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test for `newresponse.default()`
# --------------------------------------------------------------------
## Test for numeric response w/o NAs
m1 <- lm(dist ~ speed, data = cars)

expect_equal(
  as.numeric(newresponse(m1)),  # TODO: (ML) Do we really want a named vector?
  cars$dist
)

expect_equal(
  as.numeric(newresponse(m1, newdata = cars[1:3, ])),  # TODO: (ML) Do we really want a named vector?
  cars$dist[1:3]
)


## Test for numeric response w/ NAs in covariates
nd1 <- cars
nd1[c(3, 5, 7), "speed"] <- NA
expect_equal(
  as.numeric(newresponse(m1, newdata = nd1, na.action = na.pass)),  # TODO: (ML) Do we really want a named vector?
  cars$dist
)

expect_equal(
  as.numeric(attr(newresponse(m1, newdata = nd1, na.action = na.pass), "weights")),
  ifelse(is.na(nd1$dist), NA, 1)
)


## Test for numeric response w/ NAs in response
nd2 <- cars
nd2[c(3, 5, 7), "dist"] <- NA
expect_equal(
  as.numeric(newresponse(m1, newdata = nd2, na.action = na.pass)),  # TODO: (ML) Do we really want a named vector?
  nd2$dist
)

expect_equal(
  as.numeric(attr(newresponse(m1, newdata = nd2, na.action = na.pass), "weights")),
  ifelse(is.na(nd2$dist), NA, 1)
)



# --------------------------------------------------------------------
# Test for `newresponse.glm()`
# --------------------------------------------------------------------
## Test for numeric response w/o NAs
m2 <- glm(dist ~ speed, data = cars)

expect_equal(
  as.numeric(newresponse(m2)),  # TODO: (ML) Do we really want a named vector?
  cars$dist
)

expect_equal(
  as.numeric(newresponse(m2, newdata = cars[1:3, ])),  # TODO: (ML) Do we really want a named vector?
  cars$dist[1:3]
)


## Test for numeric response w/ NAs in covariates
nd1 <- cars
nd1[c(3, 5, 7), "speed"] <- NA
expect_equal(
  as.numeric(newresponse(m2, newdata = nd1, na.action = na.pass)),  # TODO: (ML) Do we really want a named vector?
  cars$dist
)

expect_equal(
  as.numeric(attr(newresponse(m2, newdata = nd1, na.action = na.pass), "weights")),
  ifelse(is.na(nd1$dist), NA, 1)
)


## Test for numeric response w/ NAs in response
nd2 <- cars
nd2[c(3, 5, 7), "dist"] <- NA
expect_equal(
  as.numeric(newresponse(m2, newdata = nd2, na.action = na.pass)),  # TODO: (ML) Do we really want a named vector?
  nd2$dist
)

expect_equal(
  as.numeric(attr(newresponse(m2, newdata = nd2, na.action = na.pass), "weights")),
  ifelse(is.na(nd2$dist), NA, 1)
)
