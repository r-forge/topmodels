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

# --------------------------------------------------------------------
# Run some tests FIXME: (ML) Implement more systematic testing
# --------------------------------------------------------------------
require("crch")
require("ggplot2")
data("CrabSatellites", package = "countreg")
CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]

m1 <- lm(dist ~ speed, data = cars)
m2 <- crch(dist ~ speed | speed, left = 3, data = cars)
m3 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)

nd_cars <- cars[1:10, ]

expect_silent( 
  rootogram(m1, 
    newdata = nd_cars, class = "tibble", scale = "raw", style = "hanging")
)

expect_silent( 
  rootogram(m1, 
    newdata = nd_cars, class = "data.frame", style = "standing",
    xlab = "test", ylab = "test", main = "test", plot = "ggplot2"
  )
)

r2 <- rootogram(m2, style = "suspended", plot = FALSE)
r3 <- rootogram(m3, plot = FALSE)

expect_message(
  autoplot(c(r2, r3), pch = 3, legend = TRUE)
)

#expect_silent( #FIXME: (ML)
#  plot(c(r2, r3), pch = 3)
#)
