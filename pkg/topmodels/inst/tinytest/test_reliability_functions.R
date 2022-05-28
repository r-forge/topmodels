# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `reliagramt.R`
# --------------------------------------------------------------------
## FIXME: (ML) Implement more systematic testing

# --------------------------------------------------------------------
# Run some tests for computation 
# --------------------------------------------------------------------
suppressPackageStartupMessages(require("crch"))
suppressPackageStartupMessages(require("ggplot2"))

data("CrabSatellites", package = "countreg")
CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]

m1 <- lm(dist ~ speed, data = cars)
m2 <- crch(dist ~ speed | speed, left = 3, data = cars)
m3 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)

nd_cars <- cars[1:10, ]

## test computation and base plot
expect_silent( 
  reliagram(m1, 
    newdata = nd_cars, class = "tibble", breaks = seq(0, 1, by = 0.25), quantiles = 0.2, 
    confint = "black", confint_level = 0.9, confint_nboot = 100,
    confint_seed = 3, xlab = "test", ylab = "test", main = "test", plot = "base"
  )
)

## test computation and ggplot2 plot
expect_silent( 
  reliagram(m1, 
    newdata = nd_cars, class = "tibble", breaks = seq(0, 1, by = 0.25), quantiles = 0.2, 
    confint = "black", confint_level = 0.9, confint_nboot = 100,
    confint_seed = 3, xlab = "test", ylab = "test", main = "test", plot = "ggplot2"
  )
)

## test computation
r2 <- reliagram(m2, plot = FALSE)
r3 <- reliagram(m3, quantiles = c(0.2, 0.6), plot = FALSE, minimum = c(3, 6))

# --------------------------------------------------------------------
# Run some tests for `plot()`
# --------------------------------------------------------------------
## test autoplot for single reliagram
expect_silent(
  plot(r2, single_graph = TRUE, confint = FALSE, pch = 3, minimum = 8)
)

## test autoplot for single reliagram w/ two quantiles
expect_silent(
  plot(r3, single_graph = TRUE, confint = FALSE)
)

## test autoplot for two reliagrams
expect_silent(
  plot(c(r2, r3), single_graph = TRUE, confint = FALSE, pch = 3, minimum = 8)
)

## test `lines()` w/o confint and ref
expect_silent(
  lines(r2)
)

## test `lines()` w/ confint and ref
expect_silent(
  lines(r2, confint = TRUE, ref = TRUE, minimum = 8)
)

# --------------------------------------------------------------------
# Run some tests for `autoplot()`
# --------------------------------------------------------------------
## test autoplot for single reliagram
expect_silent(
  autoplot(r2, single_graph = TRUE, confint = FALSE, pch = 3, legend = TRUE, minimum = 8)
)

## test autoplot for single reliagram w/ two quantiles
expect_silent(
  autoplot(r3, single_graph = TRUE, confint = FALSE)
)

## test autoplot for two reliagrams
expect_silent(
  autoplot(c(r2, r3), single_graph = TRUE, confint = FALSE, pch = 3, legend = TRUE, minimum = 8)
)


