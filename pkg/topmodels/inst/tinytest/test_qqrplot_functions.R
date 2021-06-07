# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `qqrplot.R`
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
  qqrplot(m1, 
    newdata = nd_cars, class = "tibble", trafo = qnorm, nsim = 1L, delta = NULL, 
    confint = "black", confint_level = 0.9, confint_nsim = 20, confint_seed = 3, 
    xlab = "test", ylab = "test", main = "test", plot = "base"
  )
)

expect_silent( 
  qqrplot(m2, 
    newdata = nd_cars, class = "tibble", trafo = qnorm, nsim = 1L, delta = 1E-10,
    confint = "black", confint_level = 0.9, confint_nsim = 20, confint_seed = 3, 
    xlab = "test", ylab = "test", main = "test", plot = "base"
  )
)

q2 <- qqrplot(m2, nsim = 50L, plot = FALSE)
q3 <- qqrplot(m3, nsim = 50L, plot = TRUE)

expect_silent(
  autoplot(c(q2, q3), single_graph = TRUE, confint = FALSE, pch = 3, legend = TRUE)
)

expect_silent(
  plot(c(q2, q3), single_graph = TRUE, confint = FALSE)
)

expect_silent(
  points(q3, col = 2)
)

expect_silent(
  points(q3, confint = TRUE, ref = TRUE)
)




