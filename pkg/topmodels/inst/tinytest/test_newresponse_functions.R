# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `newresponse.R'
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test for newresponse.default()
# --------------------------------------------------------------------
m <- lm(dist ~ speed, data = cars)
expect_equal(
  as.numeric(newresponse(m)),  # TODO: (ML) Do we really want a named vector?
  cars$dist
)

expect_equal(
  as.numeric(newresponse(m, newdata = cars[1:3, ])),  # TODO: (ML) Do we really want a named vector?
  cars$dist[1:3]
)

