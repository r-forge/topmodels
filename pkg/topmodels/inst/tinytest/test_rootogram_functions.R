# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `rootogram.R`
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# # Test NA handling in rootogram()
# --------------------------------------------------------------------
m1 <- lm(dist ~ speed, data = cars)

## Test w/ NAs in covariates
nd1 <- cars
nd1[c(3, 5, 7), "speed"] <- NA
expect_silent(
  rootogram(m1, newdata = nd1)
)

## Test w/ NAs in response
nd2 <- cars
nd2[c(3, 5, 7), "dist"] <- NA
expect_silent(
  rootogram(m1, newdata = nd2)
)

## Test w/ NAs in covaries and response
nd3 <- cars
nd3[c(10, 11), "speed"] <- NA
nd3[c(3, 5, 7), "dist"] <- NA
expect_silent(
  rootogram(m1, newdata = nd3)
)
