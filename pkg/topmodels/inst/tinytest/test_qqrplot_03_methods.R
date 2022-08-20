# --------------------------------------------------------------------
# Testing usage of the `qqrplot()` methods
# --------------------------------------------------------------------

if (interactive()) { rm(list = objects()); library("devtools"); library("tinytest"); library("topmodels") }

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

# Create qqrplot objects needed (forcing class = 'data.frame')
expect_silent(q1 <- qqrplot(m1, class = "data.frame", plot = FALSE))
expect_silent(q2 <- qqrplot(m2, class = "data.frame", plot = FALSE))
expect_silent(q3 <- qqrplot(m3, class = "data.frame", plot = FALSE))


# --------------------------------------------------------------------
# Pithist summary method
# * Testing functionality of `summary.qqrplot`
# --------------------------------------------------------------------
expect_silent(q1_summary <- summary(q1))
expect_silent(q2_summary <- summary(q2))
expect_silent(q3_summary <- summary(q3))
expect_identical(q1, q1_summary)
expect_identical(q2, q2_summary)
expect_identical(q3, q3_summary)


# --------------------------------------------------------------------
# Combining multiple qqrplot objects
# * Testing functionality of `c.qqrpot` and `rbind.qqrpot`
# * Ensure nrow(ABC) = nrow(A) + nrow(B) + nrow(C)
# * Ensure $group contains the number of groups used
# * Ensure object ABC contains the same names as A/B/C + group + confint_lwr + confint_upr
#   (latter two from summary)
# --------------------------------------------------------------------
tmp_additional_variables <- c("group")

# Combining two pit histograms
expect_silent(q1q2 <- c(q1, q2))
expect_identical(class(q1q2), class(q1))
expect_identical(nrow(q1q2), nrow(q1) + nrow(q2))
expect_identical(ncol(q1q2), ncol(q2) + length(tmp_additional_variables))
expect_identical(names(q1q2), c(names(q1), tmp_additional_variables))
expect_identical(sort(unique(q1q2$group)), 1:2)
expect_identical(c(q1, q2), rbind(q1, q2))
rm(q1q2)

# Combining three Q-Q residual plots
expect_silent(q1q2q3 <- c(q1, q2, q3))
expect_identical(class(q1q2q3), class(q1))
expect_identical(nrow(q1q2q3), nrow(q1) + nrow(q2) + nrow(q3))
expect_identical(ncol(q1q2q3), ncol(q2) + length(tmp_additional_variables))
expect_identical(names(q1q2q3), c(names(q1), tmp_additional_variables))
expect_identical(sort(unique(q1q2q3$group)), 1:3)
expect_identical(c(q1, q2, q3), rbind(q1, q2, q3))

# NOTE: We can also row-bind the objects returned by summary;
# `c.qqrplot()` ensures the order of the variables is the same, 
# thus we should get identical results.
expect_silent(s1s2s3 <- c(q1_summary, q2_summary, q3_summary))
identical(q1q2q3, s1s2s3)
expect_silent(s1s2s3 <- rbind(q1_summary, q2_summary, q3_summary))
identical(q1q2q3, s1s2s3)

# NOTE: We can also row-bind the objects returned by summary;
# `c.qqrplot()` ensures the order of the variables is the same, 
# thus we should get identical results.
expect_silent(s1s2s3 <- c(q1_summary, q2_summary, q3_summary))
identical(q1q2q3, s1s2s3)
expect_silent(s1s2s3 <- rbind(q1_summary, q2_summary, q3_summary))
identical(q1q2q3, s1s2s3)

rm(q1q2q3, s1s2s3)

# Testing the same method for c("qqrplot", "tbl_df", ...) objects; only for the
# tribble-combo (combining tbl_p1, tbl_p2 and tbl_p3).
expect_silent(tbl_q1 <- qqrplot(m1, class = "tibble", plot = FALSE))
expect_silent(tbl_q2 <- qqrplot(m2, class = "tibble", plot = FALSE))
expect_silent(tbl_q3 <- qqrplot(m3, class = "tibble", plot = FALSE))
expect_silent(tbl_q1q2q3 <- c(tbl_q1, tbl_q2, tbl_q3))
expect_identical(class(tbl_q1), class(tbl_q1q2q3))
expect_identical(nrow(tbl_q1q2q3), nrow(tbl_q1) + nrow(tbl_q2) + nrow(tbl_q3))
expect_identical(ncol(tbl_q1q2q3), ncol(tbl_q2) + length(tmp_additional_variables))
expect_identical(names(tbl_q1q2q3), c(names(tbl_q1), tmp_additional_variables))
expect_identical(sort(unique(tbl_q1q2q3$group)), 1:3)
expect_identical(c(tbl_q1, tbl_q2, tbl_q3), rbind(tbl_q1, tbl_q2, tbl_q3))

rm(tbl_q1, tbl_q2, tbl_q3, tbl_q1q2q3)

rm(tmp_additional_variables)


# TODO: Plotting not yet covered by tests






