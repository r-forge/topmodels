# --------------------------------------------------------------------
# Testing usage of the `qqrplot()` function.
# --------------------------------------------------------------------

if (interactive()) { rm(list = objects()); library("devtools"); library("tinytest"); library("topmodels") }

suppressPackageStartupMessages(library("crch"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tibble"))

# Global seed for the tests in this file
seed <- 123

# --------------------------------------------------------------------
# Setting up the data sets/models used to test the function
# --------------------------------------------------------------------
data("CrabSatellites", package = "countreg")
CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]

# Different regression models (lm, censored lm, poisson count data model)
expect_silent(m1 <- lm(dist ~ speed, data = cars))
expect_silent(m2 <- crch(dist ~ speed | speed, left = 3, data = cars))
expect_silent(m3 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson))

# --------------------------------------------------------------------
# Sanity checks and allowed parameters
# --------------------------------------------------------------------
# Main argument must be a model object
expect_error(qqrplot(1),            info = "Main object is not a model object")
expect_error(qqrplot(NA),           info = "Main object is not a model object")

# newdata
expect_error(qqrplot(m1, newdata = 1),     info = "newdata of wrong class")
expect_error(qqrplot(m1, newdata = NA),    info = "newdata of wrong class")

# plot
expect_error(qqrplot(m1, plot = 1),                 info = "numeric not allowed")
expect_error(qqrplot(m1, plot = logical(0)),        info = "zero-length logical not allowed")
expect_error(qqrplot(m1, plot = rep(TRUE, 2)),      info = "logical length > 1 not allowed")
expect_error(qqrplot(m1, plot = "foo"),             info = "option not allowed")
expect_error(qqrplot(m1, plot = character(0)),      info = "character of length != 1 not allowed")
expect_error(qqrplot(m1, plot = c("base", "ggplot2")), info = "character of length != 1 not allowed")

# class, scale, type, and style
expect_error(qqrplot(m1, class = 1),                info = "class must be character")
expect_error(qqrplot(m1, class = character(0)),     info = "zero-length class not allowed")
expect_error(qqrplot(m1, class = "foo"),            info = "invalid argument for class")

# detrend
expect_error(qqrplot(m1, detrend = "TRUE"),            info = "detrend must be logical")
expect_error(qqrplot(m1, detrend = 1),                 info = "detrend must be logical")
expect_error(qqrplot(m1, detrend = c(TRUE, FALSE)),    info = "detrend must be of length 1")
expect_error(qqrplot(m1, detrend = logical(0)),        info = "detrend must be of length 1")

## TODO(RS): Currently testing that "type == "proportional" fails,
##           must be removed in the future when this option has been added.
#expect_error(qqrplot(m1, type = "proportional"),    info = "test 'not implemented' option")

## scale
expect_error(qqrplot(m1, scale = 1),                info = "scale must be character")
expect_error(qqrplot(m1, scale = character(0)),     info = "zero-length scale not allowed")
expect_error(qqrplot(m1, scale = "foo"),            info = "invalid argument for scale")

# simint
expect_error(qqrplot(m1, nsim = "3"),               info = "simint must be NULL, TRUE, or FALSE")
expect_error(qqrplot(m1, nsim = c(TRUE, FALSE)),    info = "simint must be NULL, TRUE, or FALSE")
expect_error(qqrplot(m1, nsim = 2:3),               info = "simint must be of length 1")
expect_error(qqrplot(m1, nsim = 0),                 info = "simint must be positive")

# delta: forwarded to pitresiduals()
expect_error(qqrplot(m1, delta = TRUE),              info = "delta must be numeric")
expect_error(qqrplot(m1, delta = "1"),               info = "delta must be numeric")
expect_error(qqrplot(m1, delta = 1:3),               info = "delta must be of length 1")
expect_error(qqrplot(m1, delta = 0),                 info = "delta must be > 0.0")

# simint_level
expect_error(qqrplot(m1, simint_level = "0.5"),     info = "simint_level must be numeric")
expect_error(qqrplot(m1, simint_level = -1e10),     info = "simint_level must be >= 0")
expect_error(qqrplot(m1, simint_level = 1+1e10),    info = "simint_level must be <= 1")
expect_error(qqrplot(m1, simint_level = 1:2 / 10),  info = "simint_level must be length 1")

# simint_nrep: forwarded to pitresiduals()
expect_error(qqrplot(m1, simint_nrep = "3"),        info = "simint_nrep must be numeric")
expect_error(qqrplot(m1, simint_nrep = 10:20),      info = "simint_nrep must be length 1")
expect_error(qqrplot(m1, simint_nrep = 0),          info = "simint_nrep must be >= 1")

# confint
expect_error(qqrplot(m1, confint = "TRUE"),          info = "confint must be logical")
expect_error(qqrplot(m1, confint = 1),               info = "confint must be logical")
expect_error(qqrplot(m1, confint = c(TRUE, FALSE)),  info = "confint must be of length 1")
expect_error(qqrplot(m1, confint = logical(0)),      info = "confint must be of length 1")

# xlab/ylab/main
expect_error(qqrplot(m1, xlab = 3),                 info = "xlab must be character")
expect_error(qqrplot(m1, xlab = character(0)),      info = "xlab must be length 1")
expect_error(qqrplot(m1, xlab = LETTERS[1:2]),      info = "xlab must be length 1")
expect_error(qqrplot(m1, ylab = 3),                 info = "xlab must be character")
expect_error(qqrplot(m1, ylab = character(0)),      info = "xlab must be length 1")
expect_error(qqrplot(m1, ylab = LETTERS[1:2]),      info = "xlab must be length 1")
expect_error(qqrplot(m1, main = 3),                 info = "main must be character or NULL")
expect_error(qqrplot(m1, main = character(0)),      info = "main must be length 1 if character")
expect_error(qqrplot(m1, main = LETTERS[1:2]),      info = "main must be length 1 if character")



# --------------------------------------------------------------------
# Basic usage; testing return objects
# --------------------------------------------------------------------

# Return class data.frame
set.seed(seed)
expect_silent(q1 <- qqrplot(m1, class = "data.frame"))
# ---> throws a warning when plotted; should be silent when not.
set.seed(seed)
expect_silent(q2 <- qqrplot(m2, class = "data.frame", plot = FALSE))
set.seed(seed)
expect_warning(q2 <- qqrplot(m2, class = "data.frame"),
               "Removed 1 rows containing non-finite values \\((stat_qqrplot_ref|stat_qqrplot_confint)\\)\\.$")
set.seed(seed)
expect_silent(q3 <- qqrplot(m3, class = "data.frame"))

# Return class tibble
set.seed(seed)
expect_silent(tbl_q1 <- qqrplot(m1, class = "tibble"))
# ---> throws a warning when plotted; should be silent when not.
set.seed(seed)
expect_silent(tbl_q2 <- qqrplot(m2, class = "tibble", plot = FALSE))
set.seed(seed)
expect_warning(tbl_q2 <- qqrplot(m2, class = "tibble"),
               "Removed 1 rows containing non-finite values \\((stat_qqrplot_ref|stat_qqrplot_confint)\\)\\.$")
set.seed(seed)
expect_silent(tbl_q3 <- qqrplot(m3, class = "tibble"))

set.seed(Sys.time()) # Resetting seed

# Check if we get a data.frame in return
expect_true(all(sapply(list(q1, q2, q3), function(x) inherits(x, "data.frame"))),
            info = "qqrplot(..., class = 'data.frame') did not return a data.frame")

# Check if we get a tbl_df in return
expect_true(all(sapply(list(tbl_q1, tbl_q2, tbl_q3), function(x) inherits(x, "tbl_df"))),
            info = "qqrplot(..., class = 'tibble') did not return a tbl_df")

## Check that all objects have class "qqrplot" as first main class
tmp <- list(q1, q2, q3, tbl_q1, tbl_q2, tbl_q3)
expect_true(all(sapply(tmp, function(x) class(x)[1] == "qqrplot")),
            info = "Missing \"qqrplot\" as main class for at least one object returned by qqrplot()")
expected_names <- c("observed", "expected", "simint_observed_lwr", "simint_observed_upr", "simint_expected")
expect_true(all(sapply(tmp, function(x) all(sort(expected_names) == sort(names(x))))),
            info = "Unexpected/missing variables in object returned by qqrplot()")

# Simint should be all NA for the two models with discrete distributions
tmp2 <- list(q1, q2, tbl_q1, tbl_q2)
expect_true(all(sapply(tmp2, function(x) all(sapply(x[grepl("^simint_", names(x))], function(k) all(is.na(k)))))),
            info = "simint_observed_lwr, simint_observed_upr should, and simint_expected should all be NA!")
rm(tmp2)

# Test that the attributes contain the same information.
# tbl_df and data.frame return names/row.names in different order; thus we sort them first
# and remove the class attribute as it must differ between the two objects.
tmp_get_attr <- function(x, drop = "class") {
    tmp <- attributes(x)[!names(attributes(x)) %in% drop]
    tmp[sort(names(tmp))]
}
expect_identical(tmp_get_attr(q1), tmp_get_attr(tbl_q1))
expect_identical(tmp_get_attr(q2), tmp_get_attr(tbl_q2))
expect_identical(tmp_get_attr(q3), tmp_get_attr(tbl_q3))
# Same is true for q1/q2 except the 'main' title differs (name of the original object)
expect_identical(tmp_get_attr(q1, drop = c("class", "main")), tmp_get_attr(q2, drop = c("class", "main")))
# q3 contains a different data set; 'row.names' must differ, rest should be identical except
# 'simint' which is auto-set to TRUE for models with discrete distribution (Poisson glm)
expect_identical(tmp_get_attr(q1, drop = c("class", "main", "simint", "row.names")),
                 tmp_get_attr(q3, drop = c("class", "main", "simint", "row.names")))
rm(tmp_get_attr)


# --------------------------------------------------------------------
# Testing default values to make sure they do not change.
# We need to force three things:
# * class: to not be dependent on whether tibble has been loaded already
# * ylab: no lazy evaluation, thus setting to `"Density"` by default
# * main: overwrite main to get the identical object
# --------------------------------------------------------------------
tmp_qqrplot_with_defaults <- function(object, main, ...) {
    qqrplot(object,
            newdata = NULL,
            plot = TRUE,
            class = "data.frame",
            detrend = FALSE,
            scale = c("normal", "uniform"),
            nsim = 1L,
            delta = NULL,
            simint = TRUE,
            simint_level = 0.95,
            simint_nrep = 250,

            ## plotting arguments
            confint = TRUE,
            ref = TRUE,
            xlab = "Theoretical quantiles",
            ylab = "Quantile residuals", # <- forced
            main = main, ...)
}

# Note: we need seeding
set.seed(seed); expect_silent(q1_default <- tmp_qqrplot_with_defaults(m1, main = "m1"))
set.seed(seed); expect_identical(q1, q1_default)
set.seed(seed); expect_warning(q2_default <- tmp_qqrplot_with_defaults(m2, main = "m2"),
                               "Removed 1 rows containing non-finite values \\((stat_qqrplot_ref|stat_qqrplot_confint)\\)\\.$")
set.seed(seed); expect_equivalent(q2, q2_default)
set.seed(seed); expect_silent(q3_default <- tmp_qqrplot_with_defaults(m3, main = "m3"))
set.seed(seed); expect_identical(q3, q3_default)
rm(tmp_qqrplot_with_defaults, q1_default, q2_default, q3_default)

set.seed(Sys.time()) # Resetting seed

# Given we know that p1, p2, p3 follow the defaults: Test/check
# a series of attributes.
tmp <- list(q1, q2, q3)
expect_true(all(sapply(tmp, function(x) isFALSE(attr(x, "detrend")))))
expect_true(all(sapply(tmp, function(x) identical(attr(x, "scale"),  "normal"))))
expect_true(all(sapply(tmp, function(x) isTRUE(attr(x, "confint")))))
expect_true(all(sapply(tmp, function(x) isTRUE(attr(x, "ref")))))
expect_true(all(sapply(tmp, function(x) identical(attr(x, "xlab"), "Theoretical quantiles"))))
expect_true(all(sapply(tmp, function(x) identical(attr(x, "ylab"), "Quantile residuals"))))
expect_true(isFALSE(attr(q1, "simint")))
expect_true(isFALSE(attr(q2, "simint")))
expect_true(isTRUE(attr(q3,  "simint")))


