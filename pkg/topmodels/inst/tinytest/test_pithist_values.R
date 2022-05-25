# --------------------------------------------------------------------
# Testing the calculation of the pithist method.
# Partial tests only.
# --------------------------------------------------------------------

if (interactive()) { rm(list = objects()); library("devtools"); library("tinytest"); library("topmodels") }

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
expect_silent(p2_10 <- pithist(m2, class = "data.frame", breaks = bk10, plot = FALSE))
expect_silent(p3_10 <- pithist(m3, class = "data.frame", breaks = bk10, plot = FALSE))


# --------------------------------------------------------------------
# Manually calculate the PIT histogram
# --------------------------------------------------------------------

# LINEAR MODEL -------------------------
manual_pit_lm <- function(x, breaks) {
    location <- predict(x)
    #scale    <- summary(x)$sigma * sqrt(df.residual(x) / nobs(x)) # ML CONS
    eps      <- residuals(x)
    scale    <- sqrt(1 / length(eps) * sum((eps - mean(eps))^2))
    obs      <- model.response(model.frame(x))
    pitresid <- pnorm(obs, mean = location, sd = scale)
    bin      <- cut(pitresid, breaks = breaks, include.lowest = TRUE, left = TRUE)
    data.frame(observed  = as.numeric(table(bin)) / (length(obs) * diff(breaks)),
               expected  = qbinom(p = .5, size = length(obs), prob = diff(breaks)) / (length(obs) * diff(breaks)),
               mid       = head(breaks, -1) + diff(breaks) / 2,
               width     = diff(breaks))
}

# Calculate PIT result manually with 4 bins
expect_silent(p1_manual <- manual_pit_lm(m1, bk4))
for (n in names(p1_manual)) expect_equal(p1[, n], p1_manual[, n])
# Same but with 10 bins
expect_silent(p1_10_manual <- manual_pit_lm(m1, bk10))
for (n in names(p1_10_manual)) expect_equal(p1_10[, n], p1_10_manual[, n])


# CENSORED GAUSSIAN MODEL --------------
manual_pit_crch <- function(x, breaks) {
    obs      <- model.response(model.frame(x))
    pitresid <- predict(x, type = "probability", at = model.response(model.frame(m2)))
    bin      <- cut(pitresid, breaks = breaks, include.lowest = TRUE, left = FALSE)
    #bin      <- cut(pitresid, breaks = breaks, include.lowest = TRUE, left = TRUE)
    data.frame(observed  = as.numeric(table(bin)) / (length(obs) * diff(breaks)),
               expected  = qbinom(p = .5, size = length(obs), prob = diff(breaks)) / (length(obs) * diff(breaks)),
               mid       = head(breaks, -1) + diff(breaks) / 2,
               width     = diff(breaks))
}

# Calculate PIT result manually with 4 bins
expect_silent(p2_manual <- manual_pit_crch(m2, bk4))
for (n in names(p2_manual)) expect_equal(p2[, n], p2_manual[, n])

# Same but with 10 bins
expect_silent(p2_10_manual <- manual_pit_crch(m2, bk10))
for (n in names(p2_10_manual)) expect_equal(p2_10[, n], p2_10_manual[, n])


# POISSON MODEL ------------------------
manual_pit_poisson <- function(x, breaks) {
    obs      <- model.response(model.frame(x))
    delta <- min(diff(sort(unique(obs)))) / 5e6 # Taken from package
    fn <- function(b, obs, lambda, delta) {
        p2 <- ppois(obs,         lambda)
        p1 <- ppois(obs - delta, lambda)
        mean(pmax(0, pmin(1, (b - p1) / (p2 - p1))))
    }
    res <- sapply(breaks, fn, obs = obs, lambda = fitted(x), delta = delta)
    data.frame(observed  = (diff(res) + c(res[0], rep(0, length(res) - 1))) / diff(breaks),
               expected  = qbinom(p = .5, size = length(obs), prob = diff(breaks)) / (length(obs) * diff(breaks)),
               mid       = head(breaks, -1) + diff(breaks) / 2,
               width     = diff(breaks))
}

# Calculate PIT result manually with 4 bins
expect_silent(p3_manual <- manual_pit_poisson(m3, bk4))
for (n in names(p3_manual)) expect_equal(p3[, n], p3_manual[, n])

# Same but with 10 bins
expect_silent(p3_10_manual <- manual_pit_poisson(m3, bk10))
for (n in names(p3_10_manual)) expect_equal(p3_10[, n], p3_10_manual[, n])






