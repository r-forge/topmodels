# --------------------------------------------------------------------
# Testing usage of the `pithist()` function.
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Setting up the data sets/models used to test the function
# --------------------------------------------------------------------
library("crch")
library("ggplot2")
library("tibble")

data("CrabSatellites", package = "countreg")
CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]

# Different regression models (lm, censored lm, poisson count data model)
expect_silent(m1 <- lm(dist ~ speed, data = cars))
expect_silent(m2 <- crch(dist ~ speed | speed, left = 3, data = cars))
expect_silent(m3 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson))

# Newdata for predictions
expect_silent(nd_cars <- cars[1:10, ])


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
            info = "pithist(..., class = 'data.frame') did not return a tbl_df")

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


# --------------------------------------------------------------------
# --------------------------------------------------------------------





