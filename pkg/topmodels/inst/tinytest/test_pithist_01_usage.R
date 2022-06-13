# --------------------------------------------------------------------
# Testing usage of the `pithist()` function.
# --------------------------------------------------------------------

if (interactive()) { rm(list = objects()); library("devtools"); library("tinytest"); library("topmodels") }

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

# --------------------------------------------------------------------
# Sanity checks and allowed parameters
# --------------------------------------------------------------------
# Main argument must be a model object
expect_error(pithist(1),            info = "Main object is not a model object")
expect_error(pithist(NA),           info = "Main object is not a model object")

# newdata: forwarded to `pitresiduals()`
expect_error(pithist(m1, newdata = 3),              info = "newdata of wrong class")
expect_error(pithist(m1, newdata = NA),             info = "newdata of wrong class")

# plot
expect_error(pithist(m1, plot = 1),                 info = "numeric not allowed")
expect_error(pithist(m1, plot = logical(0)),        info = "zero-length logical not allowed")
expect_error(pithist(m1, plot = rep(TRUE, 2)),      info = "logical length > 1 not allowed")
expect_error(pithist(m1, plot = "foo"),             info = "option not allowed")
expect_error(pithist(m1, plot = character(0)),      info = "character of length != 1 not allowed")
expect_error(pithist(m1, plot = c("base", "ggplot2")), info = "character of length != 1 not allowed")

# class, scale, type, and style
expect_error(pithist(m1, class = 1),                info = "class must be character")
expect_error(pithist(m1, class = character(0)),     info = "zero-length class not allowed")
expect_error(pithist(m1, class = "foo"),            info = "invalid argument for class")
expect_error(pithist(m1, scale = 1),                info = "scale must be character")
expect_error(pithist(m1, scale = character(0)),     info = "zero-length scale not allowed")
expect_error(pithist(m1, scale = "foo"),            info = "invalid argument for scale")
expect_error(pithist(m1, type = 1),                 info = "type must be character")
expect_error(pithist(m1, type = character(0)),      info = "zero-length type not allowed")
expect_error(pithist(m1, type = "foo"),             info = "invalid argument for type")
expect_error(pithist(m1, style = 1),                info = "style must be character")
expect_error(pithist(m1, style = character(0)),     info = "zero-length style not allowed")
expect_error(pithist(m1, style = "foo"),            info = "invalid argument for style")

# TODO(RS): Currently testing that "type == "proportional" fails,
#           must be removed in the future when this option has been added.
expect_error(pithist(m1, type = "proportional"),    info = "test 'not implemented' option")

# Testing 'breaks'
expect_error(pithist(m1, breaks = "foo"),           info = "breaks non-numeric")
expect_error(pithist(m1, breaks = TRUE),            info = "breaks non-numeric")
expect_error(pithist(m1, breaks = matrix(1, nrow = 3, ncol = 1)), info = "breaks not dimension-less numeric vector")
expect_error(pithist(m1, breaks = numeric(0)),      info = "breaks invalid length (must be 2 or more)")
expect_error(pithist(m1, breaks = .99),             info = "breaks invalid; if single numeric it must be >= 1")

# simint; not tested, tested by pitresiduals() if used

# simint_level
expect_error(pithist(m1, simint_level = "0.5"),     info = "simint_level must be numeric")
expect_error(pithist(m1, simint_level = -1e10),     info = "simint_level must be >= 0")
expect_error(pithist(m1, simint_level = 1+1e10),    info = "simint_level must be <= 1")
expect_error(pithist(m1, simint_level = 1:2 / 10),  info = "simint_level must be length 1")

# simint_nrep
expect_error(pithist(m1, simint_nrep = "3"),        info = "simint_nrep must be numeric")
expect_error(pithist(m1, simint_nrep = 10:20),      info = "simint_nrep must be length 1")
expect_error(pithist(m1, simint_nrep = 0),          info = "simint_nrep must be >= 1")

# delta: forwarded to pitresiduals()
expect_error(pithist(m1, delta = TRUE),              info = "delta must be numeric")
expect_error(pithist(m1, delta = "1"),               info = "delta must be numeric")
expect_error(pithist(m1, delta = 1:3),               info = "delta must be of length 1")
expect_error(pithist(m1, delta = 0),                 info = "delta must be > 0.0")

# freq
expect_error(pithist(m1, freq = "TRUE"),            info = "freq must be logical")
expect_error(pithist(m1, freq = 1),                 info = "freq must be logical")
expect_error(pithist(m1, freq = c(TRUE, FALSE)),    info = "freq must be of length 1")
expect_error(pithist(m1, freq = logical(0)),        info = "freq must be of length 1")

# expected, confint
expect_error(pithist(m1, expected = "TRUE"),         info = "expected must be logical")
expect_error(pithist(m1, expected = 1),              info = "expected must be logical")
expect_error(pithist(m1, expected = c(TRUE, FALSE)), info = "expected must be of length 1")
expect_error(pithist(m1, expected = logical(0)),     info = "expected must be of length 1")
expect_error(pithist(m1, confint = "TRUE"),          info = "confint must be logical")
expect_error(pithist(m1, confint = 1),               info = "confint must be logical")
expect_error(pithist(m1, confint = c(TRUE, FALSE)),  info = "confint must be of length 1")
expect_error(pithist(m1, confint = logical(0)),      info = "confint must be of length 1")

# xlab/ylab/main
expect_error(pithist(m1, xlab = 3),                 info = "xlab must be character")
expect_error(pithist(m1, xlab = character(0)),      info = "xlab must be length 1")
expect_error(pithist(m1, xlab = LETTERS[1:2]),      info = "xlab must be length 1")
expect_error(pithist(m1, ylab = 3),                 info = "xlab must be character")
expect_error(pithist(m1, ylab = character(0)),      info = "xlab must be length 1")
expect_error(pithist(m1, ylab = LETTERS[1:2]),      info = "xlab must be length 1")
expect_error(pithist(m1, main = 3),                 info = "main must be character or NULL")
expect_error(pithist(m1, main = character(0)),      info = "main must be length 1 if character")
expect_error(pithist(m1, main = LETTERS[1:2]),      info = "main must be length 1 if character")

# There is a warning when `freq = TRUE` is used in combination with
# non-equidistant breaks.
expect_warning(pithist(m1, freq = TRUE, breaks = c(0, 0.1, 0.5, 1)),
               "^For non-equidistant breaks `freq = FALSE` must be used and has been set accordingly.$",
               info = "freq = TRUE and non-equidistant breaks trigger warning")
expect_silent(pithist(m1, freq = FALSE, breaks = c(0, 0.1, 0.5, 1)),
               info = "freq = FALSE and non-equidistant breaks trigger warning is silent")



# --------------------------------------------------------------------
# Basic usage; testing return objects
# --------------------------------------------------------------------
expect_silent(p1 <- pithist(m1, class = "data.frame"))
expect_silent(p2 <- pithist(m2, class = "data.frame"))
expect_silent(p3 <- pithist(m3, class = "data.frame"))
expect_silent(tbl_p1 <- pithist(m1, class = "tibble"))
expect_silent(tbl_p2 <- pithist(m2, class = "tibble"))
expect_silent(tbl_p3 <- pithist(m3, class = "tibble"))

# Check if we get a data.frame in return
expect_true(all(sapply(list(p1, p2, p3), function(x) inherits(x, "data.frame"))),
            info = "pithist(..., class = 'data.frame') did not return a data.frame")

# Check if we get a tbl_df in return
expect_true(all(sapply(list(tbl_p1, tbl_p2, tbl_p3), function(x) inherits(x, "tbl_df"))),
            info = "pithist(..., class = 'tibble') did not return a tbl_df")

# Check that all objects have class "pithist" as first main class
tmp <- list(p1, p2, p3, tbl_p1, tbl_p2, tbl_p3)
expect_true(all(sapply(tmp, function(x) class(x)[1] == "pithist")),
            info = "Missing \"pithist\" as main class for at least one object returned by pithist()")
expected_names <- c("observed", "expected", "mid", "width", "simint_lwr", "simint_upr")
expect_true(all(sapply(tmp, function(x) all(sort(expected_names) == sort(names(x))))),
            info = "Unexpected/missing variables in object returned by pithist()")
# As we are not asking for confidence intervals, simint must be empty (NA)
expect_true(all(sapply(tmp, function(x) all(is.na(x$simint_lwr)) && all(is.na(x$simint_upr)))),
            info = "simint_lwr and simint_upr should all be NA!")
expect_true(all(sapply(tmp, function(x) all(diff(x$mid) > 0))),
            info = "PIT midpoints should be unique and increasing!")
expect_true(all(sapply(tmp, function(x) all.equal(x$width, rep(min(x$width), length(x$width))))),
            info = "PIT width should be equal for all bins!")
expect_true(all(sapply(tmp, function(x) all.equal(x$expected, rep(min(x$expected), length(x$expected))))),
            info = "PIT expected should be equal for all bins!")
rm(tmp)

# Given the data sets used we know what the attr(x, "count") value must be.
# Only test on p1/p2/p3; testing them against tbl_p1, tbl_p2, tbl_p3 below.
expect_identical(attr(p1, "count"), nrow(cars))
expect_identical(attr(p2, "count"), nrow(cars))
expect_identical(attr(p3, "count"), nrow(CrabSatellites))

# Test that the attributes contain the same information.
# tbl_df and data.frame return names/row.names in different order; thus we sort them first
# and remove the class attribute as it must differ between the two objects.
tmp_get_attr <- function(x, drop = "class") {
    tmp <- attributes(x)[!names(attributes(x)) %in% drop]
    tmp[sort(names(tmp))]
}
expect_identical(tmp_get_attr(p1), tmp_get_attr(tbl_p1))
expect_identical(tmp_get_attr(p2), tmp_get_attr(tbl_p2))
expect_identical(tmp_get_attr(p3), tmp_get_attr(tbl_p3))
# Same is true for p1/p2 except the 'main' title differs (name of the original object)
expect_identical(tmp_get_attr(p1, drop = c("class", "main")), tmp_get_attr(p2, drop = c("class", "main")))
# p3 contains a different data set; thus 'counts' and 'row.names' must differ, but the rest is the same
expect_identical(tmp_get_attr(p1, drop = c("class", "main", "counts", "row.names")),
                 tmp_get_attr(p3, drop = c("class", "main", "counts", "row.names")))
rm(tmp_get_attr)


# --------------------------------------------------------------------
# Testing default values to make sure they do not change.
# We need to force three things:
# * class: to not be dependent on whether tibble has been loaded already
# * ylab: no lazy evaluation, thus setting to `"Density"` by default
# * main: overwrite main to get the identical object
# --------------------------------------------------------------------
tmp_pithist_with_defaults <- function(object, main, ...) {
    pithist(object,
            newdata = NULL,
            plot = TRUE,
            class = "data.frame", # <--- forced
            scale = "uniform",
            breaks = NULL,
            type = "expected",
            nsim = 1L,
            delta = NULL,
            simint = NULL,
            simint_level = 0.95,
            simint_nrep = 250,
            style = "bar",
            freq = FALSE,
            expected = TRUE,
            confint = TRUE,
            xlab = "PIT",
            ylab = "Density",    # <---- forced
            main = main,
       ...
     )
}
expect_silent(p1_default <- tmp_pithist_with_defaults(m1, main = "m1"))
expect_identical(p1, p1_default)
expect_silent(p2_default <- tmp_pithist_with_defaults(m2, main = "m2"))
expect_identical(p2, p2_default)
expect_silent(p3_default <- tmp_pithist_with_defaults(m3, main = "m3"))
expect_identical(p3, p3_default)
rm(tmp_pithist_with_defaults, p1_default, p2_default, p3_default)

# Given we know that p1, p2, p3 follow the defaults: Test/check
# a series of attributes.
tmp <- list(p1, p2, p3)
expect_true(all(sapply(tmp, function(x) identical(attr(x, "scale"), "uniform"))))
expect_true(all(sapply(tmp, function(x) identical(attr(x, "type"),  "expected"))))
expect_true(all(sapply(tmp, function(x) identical(attr(x, "style"), "bar"))))
expect_true(all(sapply(tmp, function(x) identical(attr(x, "xlab"), "PIT"))))
expect_true(all(sapply(tmp, function(x) identical(attr(x, "ylab"), "Density"))))
expect_true(all(sapply(tmp, function(x) is.character(attr(x, "main")) && length(attr(x, "main")) == 1L)))

expect_true(all(sapply(tmp, function(x) is.numeric(attr(x, "counts")) && length(attr(x, "counts")) == 1L)))

expect_true(all(sapply(tmp, function(x) isFALSE(attr(x, "simint")))))
expect_true(all(sapply(tmp, function(x) isFALSE(attr(x, "freq")))))
expect_true(all(sapply(tmp, function(x) isTRUE(attr(x, "expected")))))
expect_true(all(sapply(tmp, function(x) isTRUE(attr(x, "confint")))))


# Frequency or density?
expect_silent(p_density   <- pithist(m1, plot = FALSE, class = "data.frame", freq = FALSE))
expect_silent(p_frequency <- pithist(m1, plot = FALSE, class = "data.frame", freq = TRUE))
expect_identical(attr(p_density,   "ylab"), "Density")
expect_identical(attr(p_frequency, "ylab"), "Frequency")
expect_equal(p_density$observed * p_density$width * nrow(model.frame(m1)),
             p_frequency$observed, info = "checking observed frequency versus density")
expect_equal(p_density$expected * p_density$width * nrow(model.frame(m1)),
             p_frequency$expected, info = "checking expected frequency versus density")



