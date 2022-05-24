# --------------------------------------------------------------------
# Testing the calculation of the pithist method.
# Partial tests only.
# --------------------------------------------------------------------

if (interactive()) { library("devtools"); library("tinytest"); library("topmodels") }

library("crch")
library("ggplot2")
library("tibble")

# --------------------------------------------------------------------
# Setting up the data sets/models used to test the function
# --------------------------------------------------------------------
data("CrabSatellites", package = "countreg")
CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]

# Different regression models (lm, censored lm, poisson count data model)
expect_silent(m1 <- lm(dist ~ speed, data = cars))
expect_silent(m2 <- crch(dist ~ speed | speed, left = 3, data = cars))
expect_silent(m3 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson))

# Newdata for predictions
expect_silent(nd_cars <- cars[1:10, ])

# Create pithist objects needed (forcing class = 'data.frame')
bk4 <- seq(0, 1, by = 0.25)
expect_silent(p1 <- pithist(m1, class = "data.frame", breaks = bk4, plot = FALSE))
expect_silent(p2 <- pithist(m2, class = "data.frame", breaks = bk4, plot = FALSE))
expect_silent(p3 <- pithist(m3, class = "data.frame", breaks = bk4, plot = FALSE))

bk10 <- seq(0, 1, by = 0.1)
expect_silent(p1_10 <- pithist(m1, class = "data.frame", breaks = bk10, plot = FALSE))


# --------------------------------------------------------------------
# Manually calculate the PIT histogram
# --------------------------------------------------------------------
manual_pit_lm <- function(x, breaks, eps = 0) { ######sqrt(.Machine$double.eps)) {
    location <- predict(x)
    scale    <- rep(sd(residuals(x)), length(location))
    obs      <- model.response(model.frame(x))
    pitresid <- pnorm(obs, mean = location, sd = scale)
    bin      <- cut(pitresid + eps, breaks = breaks, include.lowest = TRUE, left = TRUE)
    print(hist(pitresid, breaks = breaks, plot = FALSE))
    data.frame(observed  = as.numeric(table(bin)) / (length(obs) * diff(breaks)),
               expected  = qbinom(p = .5, size = length(obs), prob = diff(breaks)) / (length(obs) * diff(breaks)),
               mid       = head(breaks, -1) + diff(breaks) / 2,
               width     = diff(breaks))
}

# Calculate PIT result manually with 4 bins
expect_silent(p1_manual <- manual_pit_lm(m1, bk4))
for (n in names(p1_manual)) expect_equal(p1[, n], p1_manual[, n])

# TODO(R): Problem, cannot exactely reproduce the observed frequencies.
#          I am fairly sure it has someting to do with eps/left/right.
# Same but with 10 bins
expect_silent(p1_10_manual <- manual_pit_lm(m1, bk10))
#for (n in names(p1_10_manual)) expect_equal(p1_10[, n], p1_10_manual[, n])
#cbind(pithist = p1_10$observed, manual = p1_10_manual$observed, diff = round(p1_10$observed - p1_10_manual$observed, 4))









