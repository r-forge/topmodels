# --------------------------------------------------------------------
# Testing bamlss family distrubiton3 functionality
# gaussian_bamlss: Normal distribution
# --------------------------------------------------------------------

if (interactive()) { library("tinytest") }

suppressPackageStartupMessages(library("topmodels"))
suppressPackageStartupMessages(library("distributions3"))
suppressPackageStartupMessages(library("bamlss"))
suppressPackageStartupMessages(library("countreg"))

data("CrabSatellites")

# --------------------------------------------------------------------
# Gaussian model
# --------------------------------------------------------------------

# First model
expect_silent(mod <- bamlss(speed ~ dist, data = cars, family = "gaussian"))

expect_silent(d <- prodist(mod))
expect_identical(class(d), c("BAMLSS", "distribution"))
expect_identical(length(d), nrow(cars))
expect_true(all(grepl("^BAMLSS gaussian distribution", format(d))))

expect_warning(topmodels(mod, ask = FALSE),
               pattern = "no 'type' information provided")
rm(mod, d)

# Second model
expect_silent(mod <- bamlss(speed ~ dist | s(dist), data = cars, family = "gaussian"))

expect_silent(d <- prodist(mod))
expect_identical(class(d), c("BAMLSS", "distribution"))
expect_identical(length(d), nrow(cars))
expect_true(all(grepl("^BAMLSS gaussian distribution", format(d))))

expect_warning(topmodels(mod, ask = FALSE),
               pattern = "no 'type' information provided",
               info = "Check warning")
rm(mod, d)


# --------------------------------------------------------------------
# Zero left-censored
# --------------------------------------------------------------------
set.seed(1)
data <- data.frame(x = pmax(0, rnorm(200, 2, 3)))
devtools::load_all("~/Software/bayesr/pkg/bamlss")
expect_silent(mod <- bamlss(x ~ 1 | 1, data = data, family = "cnorm"))

expect_silent(d <- prodist(mod))
expect_identical(class(d), c("BAMLSS", "distribution"))
expect_identical(length(d), nrow(data))
expect_true(all(grepl("^BAMLSS cnorm distribution", format(d))))

# Currently runs into an error related to the calculation of the
# probs (object2).
expect_error(topmodels(mod, ask = FALSE),
             pattern = "object and object2 are not equal")
#expect_warning(topmodels(mod, ask = FALSE),
#               pattern = "no 'type' information provided",
#               info = "Check warning")
rm(mod, d)


# --------------------------------------------------------------------
# Poisson model
# --------------------------------------------------------------------
expect_silent(mod <- bamlss(satellites ~ weight * width, data = CrabSatellites, family = "poisson"))

expect_silent(d <- prodist(mod))
expect_identical(class(d), c("BAMLSS", "distribution"))
expect_identical(length(d), nrow(CrabSatellites))
expect_true(all(grepl("^BAMLSS poisson distribution", format(d))))

expect_warning(topmodels(mod, ask = FALSE),
               pattern = "no 'type' information provided")
rm(mod, d)



