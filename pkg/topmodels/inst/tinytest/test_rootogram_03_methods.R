# --------------------------------------------------------------------
# Testing usage of the `rootogram()` methods
# --------------------------------------------------------------------

if (interactive()) { library("devtools"); library("tinytest"); library("topmodels") }

suppressPackageStartupMessages(library("crch"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tibble"))

# --------------------------------------------------------------------
# Setting up the data sets/models used to test the function
# --------------------------------------------------------------------
data("CrabSatellites", package = "countreg")
CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]

# Different regression models (lm, censored lm, poisson count data model)
expect_silent(m1 <- lm(dist ~ speed, data = cars))
expect_silent(m2 <- crch(dist ~ speed | speed, left = 3, data = cars))
expect_silent(m3 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson))

# Create rootogram objects needed (forcing class = 'data.frame')
expect_silent(r1 <- rootogram(m1, class = "data.frame", plot = FALSE))
expect_silent(r2 <- rootogram(m2, class = "data.frame", plot = FALSE))
expect_silent(r3 <- rootogram(m3, class = "data.frame", plot = FALSE))


# --------------------------------------------------------------------
# Pithist summary method
# * Testing functionality of `summary.rootogram`
# --------------------------------------------------------------------
tmp_additional_variables <- c("ymin", "ymax", "confint_lwr", "confint_upr")
expect_silent(r1_summary <- summary(r1)); expect_identical(class(r1_summary), class(r1))
expect_silent(r2_summary <- summary(r2)); expect_identical(class(r2_summary), class(r2))
expect_silent(r3_summary <- summary(r3)); expect_identical(class(r3_summary), class(r3))
expect_identical(names(r1_summary), c(names(r1), tmp_additional_variables))
expect_identical(names(r2_summary), c(names(r2), tmp_additional_variables))
expect_identical(names(r3_summary), c(names(r3), tmp_additional_variables))

rm(tmp_additional_variables)


# --------------------------------------------------------------------
# Combining multiple rootogram objects
# * Testing functionality of `c.rootogram` and `rbind.rootogram`
# * Ensure nrow(ABC) = nrow(A) + nrow(B) + nrow(C)
# * Ensure $group contains the number of groups used
# * Ensure object ABC contains the same names as A/B/C + group + confint_lwr + confint_upr
#   (latter two from summary)
# --------------------------------------------------------------------
tmp_additional_variables <- c("group")

# Combining two rootograms
expect_silent(r1r2 <- c(r1, r2))
expect_identical(class(r1r2), class(r1))
expect_identical(nrow(r1r2), nrow(r1) + nrow(r2))
expect_identical(ncol(r1r2), ncol(r2) + length(tmp_additional_variables))
expect_identical(names(r1r2), c(names(r1), tmp_additional_variables))
expect_identical(sort(unique(r1r2$group)), 1:2)
expect_identical(c(r1, r2), rbind(r1, r2))
rm(r1r2)

# Combining three rootograms
expect_silent(r1r2r3 <- c(r1, r2, r3))
expect_identical(class(r1r2r3), class(r1))
expect_identical(nrow(r1r2r3), nrow(r1) + nrow(r2) + nrow(r3))
expect_identical(ncol(r1r2r3), ncol(r2) + length(tmp_additional_variables))
expect_identical(names(r1r2r3), c(names(r1), tmp_additional_variables))
expect_identical(sort(unique(r1r2r3$group)), 1:3)
expect_identical(c(r1, r2, r3), rbind(r1, r2, r3))

# NOTE: Combining rootogramd does not work.
# TODO: (RS2ML) Is this intended?
expect_error(c(r1_summary, r2_summary, r3_summary), info = "combining rootograms not possible")
##### expect_silent(s1s2s3 <- c(r1_summary, r2_summary, r3_summary))
##### identical(r1r2r3, s1s2s3)
##### expect_silent(s1s2s3 <- rbind(r1_summary, r2_summary, r3_summary))
##### identical(r1r2r3, s1s2s3)

# Testing the same method for c("rootogram", "tbl_df", ...) objects; only for the
# tribble-combo (combining tbl_r1, tbl_r2 and tbl_r3).
expect_silent(tbl_r1 <- rootogram(m1, class = "tibble", plot = FALSE))
expect_silent(tbl_r2 <- rootogram(m2, class = "tibble", plot = FALSE))
expect_silent(tbl_r3 <- rootogram(m3, class = "tibble", plot = FALSE))
expect_silent(tbl_r1r2r3 <- c(tbl_r1, tbl_r2, tbl_r3))
expect_identical(class(tbl_r1), class(tbl_r1r2r3))
expect_identical(nrow(tbl_r1r2r3), nrow(tbl_r1) + nrow(tbl_r2) + nrow(tbl_r3))
expect_identical(ncol(tbl_r1r2r3), ncol(tbl_r2) + length(tmp_additional_variables))
expect_identical(names(tbl_r1r2r3), c(names(tbl_r1), tmp_additional_variables))
expect_identical(sort(unique(tbl_r1r2r3$group)), 1:3)
expect_identical(c(tbl_r1, tbl_r2, tbl_r3), rbind(tbl_r1, tbl_r2, tbl_r3))

rm(tbl_r1, tbl_r2, tbl_r3, tbl_r1r2r3)

rm(tmp_additional_variables)


# TODO: Plotting not yet covered by tests






