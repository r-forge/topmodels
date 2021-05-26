# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `reliagramt.R`
# --------------------------------------------------------------------

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
  reliagram(m1, 
    newdata = nd_cars, class = "tibble", breaks = seq(0, 1, by = 0.25), quantiles = 0.2, 
    confint = "black", confint_level = 0.9, confint_nboot = 100,
    confint_seed = 3, xlab = "test", ylab = "test", main = "test", plot = "base"
  )
)

expect_silent( 
  reliagram(m1, 
    newdata = nd_cars, class = "tibble", breaks = seq(0, 1, by = 0.25), quantiles = 0.2, 
    confint = "black", confint_level = 0.9, confint_nboot = 100,
    confint_seed = 3, xlab = "test", ylab = "test", main = "test", plot = "ggplot2"
  )
)

r2 <- reliagram(m2, quantiles = c(0.2, 0.6), plot = FALSE)
r3 <- reliagram(m3, quantiles = c(0.2, 0.6), plot = FALSE, minimum = c(3, 6))

expect_silent(
  autoplot(r2, single_graph = TRUE, confint = FALSE, pch = 3, legend = TRUE, minimum = 8)
)

expect_silent(
  plot(r3, single_graph = TRUE, confint = FALSE)
)

expect_silent(
  lines(r2)
)

expect_silent(
  lines(r2, confint = TRUE, ref = TRUE, minimum = 8)
)




