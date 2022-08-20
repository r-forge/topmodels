# --------------------------------------------------------------------
# Testing the calculation of the pithist method.
# Partial tests only.
# --------------------------------------------------------------------

if (interactive()) { rm(list = objects()); library("devtools"); library("tinytest"); library("topmodels") }

suppressPackageStartupMessages(library("digest"))
suppressPackageStartupMessages(library("crch"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tibble"))

# --------------------------------------------------------------------
# Setting up the data sets/models used to test the function
# --------------------------------------------------------------------
data("CrabSatellites", package = "countreg")

# Different regression models (lm, censored lm, poisson count data model)
expect_silent(m1 <- lm(dist ~ speed, data = cars))
expect_silent(m2 <- crch(dist ~ speed | speed, left = 3, data = cars))
expect_silent(m3 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson))

# Newdata for predictions
expect_silent(nd_cars <- cars[1:10, ])

# Create pithist objects needed (forcing class = 'data.frame')
# All on scale = "uniform"
bk4 <- seq(0, 1, by = 0.25)
expect_silent(p1_uniform <- pithist(m1, scale = "uniform", class = "data.frame", breaks = bk4, plot = FALSE))
expect_silent(p2_uniform <- pithist(m2, scale = "uniform", class = "data.frame", breaks = bk4, plot = FALSE))
expect_silent(p3_uniform <- pithist(m3, scale = "uniform", class = "data.frame", breaks = bk4, plot = FALSE))

bk10 <- seq(0, 1, by = 0.1)
expect_silent(p1_uniform10 <- pithist(m1, scale = "uniform", class = "data.frame", breaks = bk10, plot = FALSE))
expect_silent(p2_uniform10 <- pithist(m2, scale = "uniform", class = "data.frame", breaks = bk10, plot = FALSE))
expect_silent(p3_uniform10 <- pithist(m3, scale = "uniform", class = "data.frame", breaks = bk10, plot = FALSE))


# --------------------------------------------------------------------
# Testing the entire objects uisng sha1 hasing. Not precise in case it
# fails but covers the entire object.
# --------------------------------------------------------------------
expect_identical(digest::sha1(p1_uniform, digits = 10),    "bd8e42cb9ee24db09adf1372adefb3cc73848168")
expect_identical(digest::sha1(p2_uniform, digits = 10),    "93f9b7c74f1324ae34d7ba003dd37632ae842ccd")
expect_identical(digest::sha1(p3_uniform, digits = 10),    "c7aacccc28c18857aaf2fe9e7606c7def127130c")
expect_identical(digest::sha1(p1_uniform10, digits = 10), "a7df92257fac687ae019d1dea8242078e6f68ec7")
expect_identical(digest::sha1(p2_uniform10, digits = 10), "bf74794d30af290ec04fb914801eee2d334b1ebf")
expect_identical(digest::sha1(p3_uniform10, digits = 10), "0e26567a0a5205fb6568687a582ea4e27929663f")

# --------------------------------------------------------------------
# Manually calculate the PIT histogram
# --------------------------------------------------------------------

# LINEAR MODEL -------------------------
manual_uniform_pit_lm <- function(x, breaks) {
    location <- predict(x)
    #scale    <- summary(x)$sigma # ML CONS
    eps      <- residuals(x)
    scale    <- sqrt(1 / df.residual(x) * sum((eps - mean(eps))^2))
    obs      <- model.response(model.frame(x))
    pitresid <- pnorm(obs, mean = location, sd = scale)
    bin      <- cut(pitresid, breaks = breaks, include.lowest = TRUE, left = TRUE)
    data.frame(observed  = as.numeric(table(bin)) / (length(obs) * diff(breaks)),
               expected  = qbinom(p = .5, size = length(obs), prob = diff(breaks)) / (length(obs) * diff(breaks)),
               mid       = head(breaks, -1) + diff(breaks) / 2,
               width     = diff(breaks))
}

# Calculate PIT result manually with 4 bins
expect_silent(p1_uniform_manual <- manual_uniform_pit_lm(m1, bk4))
for (n in names(p1_uniform_manual)) expect_equal(p1_uniform[, n], p1_uniform_manual[, n])
# Same but with 10 bins
expect_silent(p1_uniform10_manual <- manual_uniform_pit_lm(m1, bk10))
for (n in names(p1_uniform10_manual)) expect_equal(p1_uniform10[, n], p1_uniform10_manual[, n])

rm(p1_uniform, p1_uniform10, p1_uniform_manual, p1_uniform10_manual, manual_uniform_pit_lm)



# CENSORED GAUSSIAN MODEL --------------
manual_uniform_pit_crch <- function(x, breaks) {
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
expect_silent(p2_uniform_manual <- manual_uniform_pit_crch(m2, bk4))
for (n in names(p2_uniform_manual)) expect_equal(p2_uniform[, n], p2_uniform_manual[, n])

# Same but with 10 bins
expect_silent(p2_uniform10_manual <- manual_uniform_pit_crch(m2, bk10))
for (n in names(p2_uniform10_manual)) expect_equal(p2_uniform10[, n], p2_uniform10_manual[, n])

rm(p2_uniform, p2_uniform10, p2_uniform_manual, p2_uniform10_manual, manual_uniform_pit_crch)


# POISSON MODEL ------------------------
manual_uniform_pit_poisson <- function(x, breaks) {
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
expect_silent(p3_uniform_manual <- manual_uniform_pit_poisson(m3, bk4))
for (n in names(p3_uniform_manual)) expect_equal(p3_uniform[, n], p3_uniform_manual[, n])

# Same but with 10 bins
expect_silent(p3_uniform10_manual <- manual_uniform_pit_poisson(m3, bk10))
for (n in names(p3_uniform10_manual)) expect_equal(p3_uniform10[, n], p3_uniform10_manual[, n])

rm(p3_uniform, p3_uniform10, p3_uniform_manual, p3_uniform10_manual, manual_uniform_pit_poisson)

# --------------------------------------------------------
# Repeating the same (but only with 10 breaks)
# for scale = "normal".
# --------------------------------------------------------

bk11 <- seq(-5.5, 5.5, by = 1)
expect_silent(p1_normal <- pithist(m1, scale = "normal", class = "data.frame", breaks = bk11, plot = FALSE))

# LINEAR MODEL -------------------------
manual_normal_pit_lm <- function(x, breaks) {
    location <- predict(x)
    #scale    <- summary(x)$sigma # ML CONS
    eps      <- residuals(x)
    scale    <- sqrt(1 / df.residual(x) * sum((eps - mean(eps))^2))
    obs      <- model.response(model.frame(x))
    pitresid <- qnorm(pnorm(obs, mean = location, sd = scale))
    return(pitresid)
    bin      <- cut(pitresid, breaks = breaks, include.lowest = TRUE, left = TRUE)
    width    <- diff(pnorm(breaks))
    data.frame(observed  = as.numeric(table(bin)) / (length(obs) * diff(breaks)),
               expected  = qbinom(p = .5, size = length(obs), prob = width) / (length(obs) * width),
               mid       = head(breaks, -1) + diff(breaks) / 2,
               width     = diff(breaks))
}
#k <- manual_normal_pit_lm(m1, breaks = bk11)
#k
#p1_normal
#plot(p1_normal10)
###???
###???# Calculate PIT result manually with 4 bins
###???expect_silent(p1_uniform_manual <- manual_uniform_pit_lm(m1, bk4))
###???for (n in names(p1_uniform_manual)) expect_equal(p1_uniform[, n], p1_uniform_manual[, n])
###???# Same but with 10 bins
###???expect_silent(p1_uniform10_manual <- manual_uniform_pit_lm(m1, bk10))
###???for (n in names(p1_uniform10_manual)) expect_equal(p1_uniform10[, n], p1_uniform10_manual[, n]) 
###???rm(p1_uniform, p1_uniform10, p1_uniform_manual, p1_uniform10_manual, manual_uniform_pit_lm)





