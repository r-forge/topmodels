# --------------------------------------------------------------------
# Testing usage of the `pithist()` methods
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
expect_silent(p1 <- pithist(m1, class = "data.frame", plot = FALSE))
expect_silent(p2 <- pithist(m2, class = "data.frame", plot = FALSE))
expect_silent(p3 <- pithist(m3, class = "data.frame", plot = FALSE))


# --------------------------------------------------------------------
# Pithist summary method
# * Testing functionality of `summary.pithist`
# --------------------------------------------------------------------
tmp_additional_variables <- c("confint_lwr", "confint_upr")

expect_silent(p1_summary <- summary(p1)); expect_identical(class(p1_summary), class(p1))
expect_silent(p2_summary <- summary(p2)); expect_identical(class(p2_summary), class(p2))
expect_silent(p3_summary <- summary(p3)); expect_identical(class(p3_summary), class(p3))
expect_identical(names(p1_summary), c(names(p1), tmp_additional_variables))
expect_identical(names(p2_summary), c(names(p2), tmp_additional_variables))
expect_identical(names(p3_summary), c(names(p3), tmp_additional_variables))

rm(tmp_additional_variables)


# --------------------------------------------------------------------
# Combining multiple pithist objects
# * Testing functionality of `c.pithist` and `rbind.pithist`
# * Ensure nrow(ABC) = nrow(A) + nrow(B) + nrow(C)
# * Ensure $group contains the number of groups used
# * Ensure object ABC contains the same names as A/B/C + group + confint_lwr + confint_upr
#   (latter two from summary)
# --------------------------------------------------------------------
tmp_additional_variables <- c("confint_lwr", "confint_upr", "group")

# Combining two pit histograms
expect_silent(p1p2 <- c(p1, p2))
expect_identical(nrow(p1p2), nrow(p1) + nrow(p2))
expect_identical(ncol(p1p2), ncol(p2) + length(tmp_additional_variables))
expect_identical(names(p1p2), c(names(p1), tmp_additional_variables))
expect_identical(sort(unique(p1p2$group)), 1:2)
expect_identical(c(p1, p2), rbind(p1, p2))
rm(p1p2)

# Combining three pit histograms
expect_silent(p1p2p3 <- c(p1, p2, p3))
expect_identical(nrow(p1p2p3), nrow(p1) + nrow(p2) + nrow(p3))
expect_identical(ncol(p1p2p3), ncol(p2) + length(tmp_additional_variables))
expect_identical(names(p1p2p3), c(names(p1), tmp_additional_variables))
expect_identical(sort(unique(p1p2p3$group)), 1:3)
expect_identical(c(p1, p2, p3), rbind(p1, p2, p3))

# NOTE: We can also row-bind the objects returned by summary;
# `c.pithist()` ensures the order of the variables is the same, 
# thus we should get identical results.
expect_silent(s1s2s3 <- c(p1_summary, p2_summary, p3_summary))
identical(p1p2p3, s1s2s3)
expect_silent(s1s2s3 <- rbind(p1_summary, p2_summary, p3_summary))
identical(p1p2p3, s1s2s3)

rm(p1p2p3, s1s2s3)
rm(tmp_additional_variables)


# TODO: Plotting not yet covered by tests






