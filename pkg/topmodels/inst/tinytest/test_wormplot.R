# --------------------------------------------------------------------
# Tests for function 'wormplot'. Note that wormplot is just a
# wrapper around `qqrplot` (detrended qqrplot); thus, all we have
# to test here is that the defaults work correctly.
# --------------------------------------------------------------------

if (interactive()) { library("devtools"); library("tinytest"); library("topmodels") }

library("crch")
library("ggplot2")
library("tibble")

m <- lm(dist ~ speed, data = cars)

# Calling wormplot with function defaults, and once with a named
# list according to the defaults defined in the function to check
# that they have not changed. Seed must be set due to randomization
# used in the calculation.
set.seed(123)
expect_silent(w1 <- wormplot(m, plot = FALSE))

set.seed(123)
expect_silent(w2 <- wormplot(m))

defaults <- list(m,
                 newdata = NULL,
                 plot = TRUE,
                 class = NULL,
                 detrend = TRUE,
                 scale = c("normal", "uniform"),
                 nsim = 1L,
                 delta = NULL,
                 confint = TRUE,
                 simint = TRUE,
                 simint_level = 0.95,
                 simint_nrep = 250,
                 single_graph = FALSE,
                 xlab = "Theoretical quantiles",
                 ylab = "Deviation",
                 main = NULL)

set.seed(123)
expect_silent(w3 <- do.call(wormplot, defaults))


# Performing tests on the return objects
expect_true(attr(w1, "detrend"))
expect_identical(w1, w2)
expect_identical(w1, w3)


# `wormplot()` returns an object of class `"qqrplot"` with 
# the attribute `detrend = TRUE`; else it should return the
# very same as `qqrplot(..., detrend = TRUE)` when we take
# care of 'main' (title).
set.seed(123)
expect_silent(q1 <- qqrplot(m, detrend = TRUE))
w1_tmp <- w1; attr(w1_tmp, "main") <- "test"
q1_tmp <- q1; attr(q1_tmp, "main") <- "test"
expect_identical(w1_tmp, q1_tmp)
rm(w1_tmp, q1_tmp)


# Breaking static seed 
set.seed(Sys.time())

