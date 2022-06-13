# --------------------------------------------------------------------
# Testing usage of the `pithist()` methods
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
expect_identical(class(p1p2), class(p1))
expect_identical(nrow(p1p2), nrow(p1) + nrow(p2))
expect_identical(ncol(p1p2), ncol(p2) + length(tmp_additional_variables))
expect_identical(names(p1p2), c(names(p1), tmp_additional_variables))
expect_identical(sort(unique(p1p2$group)), 1:2)
expect_identical(c(p1, p2), rbind(p1, p2))
rm(p1p2)

# Combining three pit histograms
expect_silent(p1p2p3 <- c(p1, p2, p3))
expect_identical(class(p1p2p3), class(p1))
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

# NOTE: We can also row-bind the objects returned by summary;
# `c.pithist()` ensures the order of the variables is the same, 
# thus we should get identical results.
expect_silent(s1s2s3 <- c(p1_summary, p2_summary, p3_summary))
identical(p1p2p3, s1s2s3)
expect_silent(s1s2s3 <- rbind(p1_summary, p2_summary, p3_summary))
identical(p1p2p3, s1s2s3)

rm(p1p2p3, s1s2s3)

# Testing the same method for c("pithist", "tbl_df", ...) objects; only for the
# tribble-combo (combining tbl_p1, tbl_p2 and tbl_p3).
expect_silent(tbl_p1 <- pithist(m1, class = "tibble", plot = FALSE))
expect_silent(tbl_p2 <- pithist(m2, class = "tibble", plot = FALSE))
expect_silent(tbl_p3 <- pithist(m3, class = "tibble", plot = FALSE))
expect_silent(tbl_p1p2p3 <- c(tbl_p1, tbl_p2, tbl_p3))
expect_identical(class(tbl_p1), class(tbl_p1p2p3))
expect_identical(nrow(tbl_p1p2p3), nrow(tbl_p1) + nrow(tbl_p2) + nrow(tbl_p3))
expect_identical(ncol(tbl_p1p2p3), ncol(tbl_p2) + length(tmp_additional_variables))
expect_identical(names(tbl_p1p2p3), c(names(tbl_p1), tmp_additional_variables))
expect_identical(sort(unique(tbl_p1p2p3$group)), 1:3)
expect_identical(c(tbl_p1, tbl_p2, tbl_p3), rbind(tbl_p1, tbl_p2, tbl_p3))

rm(tbl_p1, tbl_p2, tbl_p3, tbl_p1p2p3)

rm(tmp_additional_variables)


# Testing combine/rbind with two objects of different scales
expect_silent(p1_unif <- pithist(m1, class = "data.frame", scale = "normal", plot = FALSE))
expect_silent(p1_norm <- pithist(m1, class = "data.frame", scale = "uniform", plot = FALSE))
expect_error(c(p1_unif, p1_norm), pattern = "^Can't combine pit histograms which are on different scales.$")

rm(p1_unif, p1_norm)


# Testing combine/rbind with two objects with different titles
expect_silent(p1_A <- pithist(m1, class = "data.frame", plot = FALSE, main = "Model A", freq = TRUE))
expect_silent(p1_B <- pithist(m1, class = "data.frame", plot = FALSE, main = "Model B", freq = TRUE))
expect_silent(p1_C <- pithist(m1, class = "data.frame", plot = FALSE, main = "Model C", freq = FALSE))
expect_identical(attr(p1_A, "main"), "Model A")
expect_identical(attr(p1_B, "main"), "Model B")
expect_identical(attr(p1_C, "main"), "Model C")

# Combining A + B, both with freq = TRUE
expect_silent(p1_AB <- c("model 1" = p1_A, "model 2" = p1_B))
expect_identical(attr(p1_AB, "main"), c("model 1", "model 2"))
expect_true(attr(p1_AB, "freq"))
expect_identical(attr(p1_AB, "ylab"), c("model 1" = "Frequency", "model 2" = "Frequency"),
                 info = "combined pit histogram with freq = TRUE and freq = TRUE should have ylab = \"Frequency\"")

# Combining B + C with freq = TRUE and freq = FALSE; should take first one
expect_message(p1_BC <- c("model 2" = p1_B, "model 3" = p1_C),
               pattern = "^\\s\\* as arg `freq`'s definition is not unique, using solely the first \\(\"TRUE\"\\)\\s$")
expect_true(attr(p1_BC, "freq"), info = "combined pit histogram with freq = TRUE and freq = FALSE should result in freq = TRUE")
expect_identical(attr(p1_BC, "ylab"), c("model 2" = "Frequency", "model 3" = "Frequency"),
                 info = "combined pit histogram with freq = TRUE and freq = FALSE should have ylab = \"Frequency\"")

expect_message(p1_CB <- c("model 3" = p1_C, "model 2" = p1_B),
               pattern = "^\\s\\* as arg `freq`'s definition is not unique, using solely the first \\(\"FALSE\"\\)\\s$")
expect_false(attr(p1_CB, "freq"), info = "combined pit histogram with freq = FALSE and freq = TRUE should result in freq = FALSE")
expect_identical(attr(p1_CB, "ylab"), c("model 3" = "Density", "model 2" = "Density"),
                 info = "combined pit histogram with freq = FALSE and freq = TRUE should have ylab = \"Density\"")



rm(p1_A, p1_B)

# TODO: Plotting not yet covered by tests






