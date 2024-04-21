# --------------------------------------------------------------------
# Testing usage of the `rootogram()` function.
# --------------------------------------------------------------------

if (interactive()) { library("devtools"); library("tinytest"); library("topmodels") }

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

# --------------------------------------------------------------------
# Sanity checks and allowed parameters
# --------------------------------------------------------------------
# Main argument must be a model object
expect_error(rootogram(1),          info = "Main object is not a model object")
expect_error(rootogram(NA),         info = "Main object is not a model object")


# newdata
expect_error(rootogram(m1, newdata = 3),              info = "newdata of wrong class")
expect_error(rootogram(m1, newdata = NA),             info = "newdata of wrong class")

# plot
expect_error(rootogram(m1, plot = 1),                 info = "numeric not allowed")
expect_error(rootogram(m1, plot = logical(0)),        info = "zero-length logical not allowed")
expect_error(rootogram(m1, plot = rep(TRUE, 2)),      info = "logical length > 1 not allowed")
expect_error(rootogram(m1, plot = "foo"),             info = "option not allowed")
expect_error(rootogram(m1, plot = character(0)),      info = "character of length != 1 not allowed")
expect_error(rootogram(m1, plot = c("base", "ggplot2")), info = "character of length != 1 not allowed")

# class, scale, response_type
expect_error(rootogram(m1, class = 1),                info = "class must be character")
expect_error(rootogram(m1, class = character(0)),     info = "zero-length class not allowed")
expect_error(rootogram(m1, class = "foo"),            info = "invalid argument for class")
expect_error(rootogram(m1, scale = 1),                info = "scale must be character")
expect_error(rootogram(m1, scale = character(0)),     info = "zero-length scale not allowed")
expect_error(rootogram(m1, scale = "foo"),            info = "invalid argument for scale")

# breaks and width
expect_error(rootogram(m1, breaks = c(TRUE, FALSE)),  info = "breaks must be NULL or numeric")
expect_error(rootogram(m1, breaks = TRUE),            info = "breaks must be NULL or numeric")
expect_error(rootogram(m1, breaks = "foo"),           info = "breaks must be NULL or numeric")
expect_error(rootogram(m1, breaks = .99),             info = "breaks invalid; if single numeric it must be >= 1")
expect_error(rootogram(me, breaks = matrix(1:5, nrow = 1)), info = "breaks must have no dimension")
expect_error(rootogram(m1, width = "foo"),            info = "width must be numeric")
expect_error(rootogram(m1, width = TRUE),             info = "width must be numeric")
expect_error(rootogram(m1, width = 1:2),              info = "width must be single numeric")
expect_error(rootogram(m1, width = -0.2),             info = "width must be positive")

# style and scale
expect_error(rootogram(m1, style = 1),                info = "style must be character")
expect_error(rootogram(m1, style = character(0)),     info = "zero-length style not allowed")
expect_error(rootogram(m1, style = "foo"),            info = "invalid argument for style")
expect_error(rootogram(m1, style = letters),          info = "style must be character length 1")
expect_error(rootogram(m1, scale = 1),                info = "scale must be character")
expect_error(rootogram(m1, scale = character(0)),     info = "zero-length scale not allowed")
expect_error(rootogram(m1, scale = "foo"),            info = "invalid argument for scale")
expect_error(rootogram(m1, scale = letters),          info = "scale must be character length 1")

# expected, confint, ref
expect_error(rootogram(m1, expected = 1),              info = "expected must be logical or character")
expect_error(rootogram(m1, expected = c(TRUE, FALSE)), info = "expected must be of length 1")
expect_error(rootogram(m1, expected = logical(0)),     info = "expected must be of length 1")
expect_error(rootogram(m1, expected = character(0)),   info = "expected must be of length 1")
expect_error(rootogram(m1, expected = letters),        info = "expected must be of length 1")
expect_error(rootogram(m1, expected = "foo"),          info = "invalid argument for expected")
expect_error(rootogram(m1, confint = "TRUE"),          info = "confint must be logical")
expect_error(rootogram(m1, confint = 1),               info = "confint must be logical")
expect_error(rootogram(m1, confint = c(TRUE, FALSE)),  info = "confint must be of length 1")
expect_error(rootogram(m1, confint = logical(0)),      info = "confint must be of length 1")
expect_error(rootogram(m1, ref = "TRUE"),              info = "ref must be logical")
expect_error(rootogram(m1, ref = 1),                   info = "ref must be logical")
expect_error(rootogram(m1, ref = c(TRUE, FALSE)),      info = "ref must be of length 1")
expect_error(rootogram(m1, ref = logical(0)),          info = "ref must be of length 1")

# xlab/ylab/main
expect_error(rootogram(m1, xlab = 3),                 info = "xlab must be character")
expect_error(rootogram(m1, xlab = character(0)),      info = "xlab must be length 1")
expect_error(rootogram(m1, xlab = LETTERS[1:2]),      info = "xlab must be length 1")
expect_error(rootogram(m1, ylab = 3),                 info = "xlab must be character")
expect_error(rootogram(m1, ylab = character(0)),      info = "xlab must be length 1")
expect_error(rootogram(m1, ylab = LETTERS[1:2]),      info = "xlab must be length 1")
expect_error(rootogram(m1, main = 3),                 info = "main must be character or NULL")
expect_error(rootogram(m1, main = character(0)),      info = "main must be length 1 if character")
expect_error(rootogram(m1, main = LETTERS[1:2]),      info = "main must be length 1 if character")


# --------------------------------------------------------------------
# Basic usage; testing return objects
# --------------------------------------------------------------------
expect_silent(r1 <- rootogram(m1, class = "data.frame"))
expect_silent(r2 <- rootogram(m2, class = "data.frame"))
expect_silent(r3 <- rootogram(m3, class = "data.frame"))
expect_silent(tbl_r1 <- rootogram(m1, class = "tibble"))
expect_silent(tbl_r2 <- rootogram(m2, class = "tibble"))
expect_silent(tbl_r3 <- rootogram(m3, class = "tibble"))

# Check if we get a data.frame in return
expect_true(all(sapply(list(r1, r2, r3), function(x) inherits(x, "data.frame"))),
            info = "rootogram(..., class = 'data.frame') did not return a data.frame")

# Check if we get a tbl_df in return
expect_true(all(sapply(list(tbl_r1, tbl_r2, tbl_r3), function(x) inherits(x, "tbl_df"))),
            info = "rootogram(..., class = 'tibble') did not return a tbl_df")

# Check that all objects have class "rootogram" as first main class
tmp <- list(r1, r2, r3, tbl_r1, tbl_r2, tbl_r3)
expect_true(all(sapply(tmp, function(x) class(x)[1] == "rootogram")),
            info = "Missing \"rootogram\" as main class for at least one object returned by rootogram()")
expected_names <- c("observed", "expected", "mid", "width")
expect_true(all(sapply(tmp, function(x) all(sort(expected_names) == sort(names(x))))),
            info = "Unexpected/missing variables in object returned by rootogram()")
expect_true(all(sapply(tmp, function(x) all(diff(x$mid) > 0))),
            info = "rootogram midpoints should be unique and increasing!")
expect_true(all(sapply(tmp, function(x) all(x$width > 0))),
            info = "rootogram width should be greater than 0")
expect_true(all(sapply(tmp, function(x) all(is.finite(x$observed) & x$observed >= 0))),
            info = "rootogram observed must be finite greater greater or equal than 0")
expect_true(all(sapply(tmp, function(x) all(is.finite(x$expected) & x$expected >= 0))),
            info = "rootogram expected must be finite greater greater or equal than 0")
rm(tmp)

# Same test but only for non-censored models (which allow for having
# non-uniform bin sizes)
tmp <- list(r1, r3, tbl_r1, tbl_r3)
expect_true(all(sapply(tmp, function(x) all.equal(x$width, rep(min(x$width), length(x$width))))),
            info = "rootogram width should be equal for all bins!")
rm(tmp)

# Test that the attributes contain the same information.
# tbl_df and data.frame return names/row.names in different order; thus we sort them first
# and remove the class attribute as it must differ between the two objects.
tmp_get_attr <- function(x, drop = "class") {
    tmp <- attributes(x)[!names(attributes(x)) %in% drop]
    tmp[sort(names(tmp))]
}
expect_identical(tmp_get_attr(r1), tmp_get_attr(tbl_r1))
expect_identical(tmp_get_attr(r2), tmp_get_attr(tbl_r2))
expect_identical(tmp_get_attr(r3), tmp_get_attr(tbl_r3))

# As these models are based on different data sets (and distributions) some
# attributes must differ, others must be the same. Testing for all ! in "drop"
expect_identical(tmp_get_attr(r1, drop = c("class", "main", "row.names", "xlab")),
                 tmp_get_attr(r2, drop = c("class", "main", "row.names", "xlab")))
expect_identical(tmp_get_attr(r1, drop = c("class", "main", "row.names", "xlab")),
                 tmp_get_attr(r3, drop = c("class", "main", "row.names", "xlab")))
rm(tmp_get_attr)


# --------------------------------------------------------------------
# Testing default values to make sure they do not change.
# We need to force three things:
# * class: to not be dependent on whether tibble has been loaded already
# * ylab: no lazy evaluation, thus setting to `"Density"` by default
# * main: overwrite main to get the identical object
# --------------------------------------------------------------------
tmp_rootogram_with_defaults <- function(object, main, ...) {
    rootogram(object,
              newdata = NULL,
              plot = TRUE,
              class = "data.frame",    # <--- forced
              response_type = NULL,
              breaks = NULL,
              width = NULL,

              ## plotting arguments
              style = c("hanging", "standing", "suspended"),
              scale = c("sqrt", "raw"),
              expected = TRUE,
              confint = TRUE,
              ref = TRUE,
              xlab = NULL,
              ylab = NULL,
              main = main,
              ...
    )
}
expect_silent(r1_default <- tmp_rootogram_with_defaults(m1, main = "m1"))
expect_identical(r1, r1_default)
expect_silent(r2_default <- tmp_rootogram_with_defaults(m2, main = "m2"))
expect_identical(r2, r2_default)
expect_silent(r3_default <- tmp_rootogram_with_defaults(m3, main = "m3"))
expect_identical(r3, r3_default)
rm(tmp_rootogram_with_defaults, r1_default, r2_default, r3_default)


## Given we know that r1, r2, r3 follow the defaults: Test/check
## a series of attributes.
tmp <- list(r1, r2, r3)
expect_true(all(sapply(tmp, function(x) identical(attr(x, "scale"), "sqrt"))))
expect_true(all(sapply(tmp, function(x) identical(attr(x, "style"), "hanging"))))
expect_true(all(sapply(tmp, function(x) identical(attr(x, "ylab"), "sqrt(Frequency)"))))
expect_true(all(sapply(tmp, function(x) is.character(attr(x, "main")) && length(attr(x, "main")) == 1L)))
expect_identical(attr(r1, "xlab"), "dist")
expect_identical(attr(r2, "xlab"), "dist")
expect_identical(attr(r3, "xlab"), "satellites")

expect_true(all(sapply(tmp, function(x) isTRUE(attr(x, "expected")))))
expect_true(all(sapply(tmp, function(x) isTRUE(attr(x, "confint")))))
expect_true(all(sapply(tmp, function(x) isTRUE(attr(x, "ref")))))


# TODO(R): Additional tests for `breaks` if set by user.
# error if too short
# check handling of NA/Inf/-Inf
# error when breaks do not cover any observations
