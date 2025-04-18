##### --------------------------------------------------------------------
##### TESTS FOR FUNCTIONS WITHIN `qqrplot.R`
##### --------------------------------------------------------------------
####
##### --------------------------------------------------------------------
##### Run some tests FIXME: (ML) Implement more systematic testing
##### --------------------------------------------------------------------
####suppressPackageStartupMessages(library("crch"))
####suppressPackageStartupMessages(library("ggplot2"))
####
####data("CrabSatellites", package = "countreg")
####
####m1 <- lm(dist ~ speed, data = cars)
####m2 <- crch(dist ~ speed | speed, left = 3, data = cars)
####m3 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
####
####nd_cars <- cars[1:10, ]
####
####expect_silent( 
####  qqrplot(m1, 
####    newdata = nd_cars, class = "tibble", scale = "normal", nsim = 1L, delta = NULL, 
####    simint_col = "black", simint_level = 0.9, simint_nrep = 20, 
####    xlab = "test", ylab = "test", main = "test", plot = "base"
####  )
####)
####
####expect_silent( 
####  qqrplot(m2, 
####    newdata = nd_cars, class = "tibble", scale = "normal", nsim = 1L, delta = 1E-10,
####    simint_col = "black", simint_level = 0.9, simint_nrep = 20, 
####    xlab = "test", ylab = "test", main = "test", plot = "base"
####  )
####)
####
####q2 <- qqrplot(m2, nsim = 50L, plot = FALSE)
####q3 <- qqrplot(m3, nsim = 50L, plot = TRUE)
####
####expect_message(
####  autoplot(c(q2, q3), single_graph = TRUE, pch = 3, legend = TRUE)
####)
####
####expect_message(
####  autoplot(c(q2, q3), single_graph = FALSE, pch = 3, legend = TRUE)
####)
####
####dev.new()
####expect_silent({
####  plot(c(q2, q3), single_graph = TRUE, simint = FALSE)
####  points(q3, col = 2)
####  points(q3, simint = TRUE)
####})
####dev.off()
