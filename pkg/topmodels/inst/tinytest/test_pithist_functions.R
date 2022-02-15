# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `pithist.R`
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
  pithist(m1, 
    newdata = nd_cars, class = "tibble", breaks = seq(0, 1, by = 0.25), freq = FALSE, 
    confint_col = "black", confint_level = 0.9, confint_type = "exact",
    xlab = "test", ylab = "test", main = "test", plot = "base"
  )
)

expect_silent( 
  pithist(m1, 
    newdata = nd_cars, class = "tibble", breaks = seq(0, 1, by = 0.25), freq = TRUE, 
    confint_col = "black", confint_level = 0.9, confint_type = "exact",
    xlab = "test", ylab = "test", main = "test", plot = "ggplot2"
  )
)

p2 <- pithist(m2, confint_type = "approximation", plot = FALSE)
p3 <- pithist(m3, nsim = 50L, plot = TRUE)

expect_silent(
  autoplot(c(p2, p3), single_graph = TRUE, confint = FALSE, pch = 3, legend = TRUE)
)

expect_silent(
  plot(c(p2, p3), single_graph = TRUE, confint = FALSE)
)

expect_silent(
  lines(p3, col = 2)
)

expect_silent(
  lines(p3, confint = TRUE, ref = TRUE)
)




