# --------------------------------------------------------------------
# Testing usage of the `pithist()` with random PIT residuals
# --------------------------------------------------------------------

if (interactive()) { rm(list = objects()); library("devtools"); library("tinytest"); library("topmodels") }

suppressPackageStartupMessages(library("crch"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("digest"))

# --------------------------------------------------------------------
# Setting up the data sets/models used to test the function
# --------------------------------------------------------------------
data("CrabSatellites", package = "countreg")

# Different regression models (lm, censored lm, poisson count data model)
expect_silent(m1 <- lm(dist ~ speed, data = cars))
expect_silent(m2 <- crch(dist ~ speed | speed, left = 3, data = cars))
expect_silent(m3 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson))




# --------------------------------------------------------------------
# Main return value; uniform scale
# --------------------------------------------------------------------
seed <- 23456

# simulated intervals (uniform scale)
set.seed(seed); expect_silent(p1r_uniform <- pithist(m1, plot = FALSE, type = "random", scale = "uniform", class = "data.frame"))
set.seed(seed); expect_silent(p2r_uniform <- pithist(m2, plot = FALSE, type = "random", scale = "uniform", class = "data.frame"))
set.seed(seed); expect_silent(p3r_uniform <- pithist(m3, plot = FALSE, type = "random", scale = "uniform", class = "data.frame"))
set.seed(Sys.time()) # Resetting seed

p_tmp <- list("p1r_uniform" = p1r_uniform, "p2r_uniform" = p2r_uniform, "p3r_uniform" = p3r_uniform)
expect_true(all(sapply(p_tmp, function(x) identical(class(x), c("pithist", "data.frame")))),
            info = "return class of pithist with randomized PIT residuals")
expect_true(all(sapply(p_tmp, function(x) sum(is.na(x)) == 0)),
            info = "no missing values in object returned by pithist with randomized PIT residuals")
expect_true(all(sapply(p_tmp, function(x) all(with(x, simint_upr - simint_lwr) >= 0))),
            info = "simint_upr must be larger or equal to simint_lwr")
expect_true(all(sapply(p_tmp, function(x) all(with(x, simint_upr - simint_lwr) >= 0))),
            info = "simint_upr must be larger or equal to simint_lwr")
expect_equal(p1r_uniform$simint_lwr, p1r_uniform$simint_upr,
             info = "discrete distribution, lower/upper bound of simulated PIT identical")
expect_equal(p1r_uniform$simint_lwr, p1r_uniform$observed,
             info = "discrete distribution, observed and bounds of simulated PIT identical")
expect_equal(p2r_uniform$simint_lwr, p2r_uniform$simint_upr,
             info = "discrete distribution, lower/upper bound of simulated PIT identical")
expect_equal(p2r_uniform$simint_lwr, p2r_uniform$observed,
             info = "discrete distribution, observed and bounds of simulated PIT identical")

# No longer need the two based on the continuous distribution
rm(p1r_uniform, p2r_uniform, p_tmp)

# discrete distribution (Poisson): random interval must contain the observed value
expect_true(all(p3r_uniform$observed > p3r_uniform$simint_lwr),
            info = "lower bound of simulated random PIT must be lower than observed")
expect_true(all(p3r_uniform$observed < p3r_uniform$simint_upr),
            info = "upper bound of sumulated random PIT must be lower than observed")

# given we used a fixed seed the sha1 hash must be ...
expect_identical(sha1(p3r_uniform), "8fa1a2ef5ef61c338d724ab86b22224869267acb",
                 info = "compare result of (pseudo-)random pithist with sha1")

# Testing nsim = 1 (default) against nsim = 10 and nsim = 20; expecting sharper intervals
set.seed(seed); expect_silent(p3r_uniform10 <- pithist(m3, plot = FALSE, type = "random", scale = "uniform", nsim = 10))
set.seed(seed); expect_silent(p3r_uniform20 <- pithist(m3, plot = FALSE, type = "random", scale = "uniform", nsim = 20))
set.seed(Sys.time())

fn <- function(x) with(x, x$simint_upr - x$simint_lwr)
expect_true(all(fn(p3r_uniform) > fn(p3r_uniform10)),
            info = "width of simulated interval with nsim1 > nsim10")
expect_true(all(fn(p3r_uniform10) > fn(p3r_uniform20)),
            info = "width of simulated interval with nsim10 > nsim20")


# ---------------------------------------------------------
# TODO: (RS2ML): Can the simulated interval be lower than 0???

if (any(p3r_uniform[, c("simint_lwr", "simint_upr")] < 0)) warning("NEGATIVE SIMINT IN PITHIST!")
if (any(p3r_uniform10[, c("simint_lwr", "simint_upr")] < 0)) warning("NEGATIVE SIMINT IN PITHIST!")
if (any(p3r_uniform20[, c("simint_lwr", "simint_upr")] < 0)) warning("NEGATIVE SIMINT IN PITHIST!")

rm(p3r_uniform10, p3r_uniform20)

# same seed, same nsim, different simint_levels (confidence level)
set.seed(seed); expect_silent(p3r_uniform94 <- pithist(m3, plot = FALSE, type = "random", scale = "uniform", simint_level = 0.94))
set.seed(seed); expect_silent(p3r_uniform96 <- pithist(m3, plot = FALSE, type = "random", scale = "uniform", simint_level = 0.95))
expect_true(all(p3r_uniform$simint_lwr <= p3r_uniform94$simint_lwr),
            info = "lower bound of simulated interval for simint_level 0.95 <= 0.94")
expect_true(all(p3r_uniform$simint_lwr >= p3r_uniform96$simint_lwr),
            info = "lower bound of simulated interval for simint_level 0.95 >= 0.96")
expect_true(all(p3r_uniform$simint_upr >= p3r_uniform94$simint_upr),
            info = "upper bound of simulated interval for simint_level 0.95 >= 0.94")
expect_true(all(p3r_uniform$simint_upr <= p3r_uniform96$simint_upr),
            info = "upper bound of simulated interval for simint_level 0.95 <= 0.96")

rm(p3r_uniform94, p3r_uniform96, fn)

# same as p3r_uniform but we set simint = FALSE at the same time to not plot simulated confidence
# interval. In this case simint_lwr and simint_upr should be empty while observed should come
# from the simulation (thus identical to p3r_uniform)
set.seed(seed); expect_silent(p3r_uniform_nosimint <- pithist(m3, plot = FALSE, type = "random", scale = "uniform", class = "data.frame", simint = FALSE))
expect_true(all(is.na(p3r_uniform_nosimint$simint_lwr)) && all(is.na(p3r_uniform_nosimint$simint_upr)),
            info = "lower/upper bound must be NA if simint = FALSE")
expect_identical(p3r_uniform[, 1:4], p3r_uniform_nosimint[, 1:4],
            info = "observed, expected, mid, and width must be identical no matter what simint is set to")

rm(p3r_uniform_nosimint)

# TODO: (RS2ML) Any better way to test the impact of simint_nrep?
set.seed(seed); expect_silent(tmp <- pithist(m3, plot = FALSE, nsim = 1, type = "random", scale = "uniform", class = "data.frame", simint_nrep = 100))
expect_true(!identical(p3r_uniform, tmp),
            info = "different number of repetitions for drawing quantiles - results not identical")
expect_identical(sha1(tmp), "aa988c3967e99fe5ce79951621388e56f1c00803")

rm(p3r_uniform, tmp)



# --------------------------------------------------------------------
# Main return value; same tests as above for normal scale
# --------------------------------------------------------------------
seed <- 23456

# simulated intervals (normal scale)
set.seed(seed); expect_silent(p1r_normal <- pithist(m1, plot = FALSE, type = "random", scale = "normal", class = "data.frame"))
set.seed(seed); expect_silent(p2r_normal <- pithist(m2, plot = FALSE, type = "random", scale = "normal", class = "data.frame"))
set.seed(seed); expect_silent(p3r_normal <- pithist(m3, plot = FALSE, type = "random", scale = "normal", class = "data.frame"))
set.seed(Sys.time()) # Resetting seed

p_tmp <- list("p1r_normal" = p1r_normal, "p2r_normal" = p2r_normal, "p3r_normal" = p3r_normal)
expect_true(all(sapply(p_tmp, function(x) identical(class(x), c("pithist", "data.frame")))),
            info = "return class of pithist with randomized PIT residuals")
expect_true(all(sapply(p_tmp, function(x) sum(is.na(x)) == 0)),
            info = "no missing values in object returned by pithist with randomized PIT residuals")
expect_true(all(sapply(p_tmp, function(x) all(with(x, simint_upr - simint_lwr) >= 0))),
            info = "simint_upr must be larger or equal to simint_lwr")
expect_true(all(sapply(p_tmp, function(x) all(with(x, simint_upr - simint_lwr) >= 0))),
            info = "simint_upr must be larger or equal to simint_lwr")
expect_equal(p1r_normal$simint_lwr, p1r_normal$simint_upr,
             info = "discrete distribution, lower/upper bound of simulated PIT identical")
expect_equal(p1r_normal$simint_lwr, p1r_normal$observed,
             info = "discrete distribution, observed and bounds of simulated PIT identical")
expect_equal(p2r_normal$simint_lwr, p2r_normal$simint_upr,
             info = "discrete distribution, lower/upper bound of simulated PIT identical")
expect_equal(p2r_normal$simint_lwr, p2r_normal$observed,
             info = "discrete distribution, observed and bounds of simulated PIT identical")

# No longer need the two based on the continuous distribution
rm(p1r_normal, p2r_normal, p_tmp)

# discrete distribution (Poisson): random interval must contain the observed value
expect_true(all(p3r_normal$observed > p3r_normal$simint_lwr),
            info = "lower bound of simulated random PIT must be lower than observed")
expect_true(all(p3r_normal$observed < p3r_normal$simint_upr),
            info = "upper bound of sumulated random PIT must be lower than observed")


## FIXME: fails on R-Forge
## # given we used a fixed seed the sha1 hash must be ...
## expect_identical(sha1(p3r_normal), "80a66ac3be1f028cb8e37d452b232fe1dfb216e2",
##                  info = "compare result of (pseudo-)random pithist with sha1")

# Testing nsim = 1 (default) against nsim = 10 and nsim = 20; expecting sharper intervals
set.seed(seed); expect_silent(p3r_normal10 <- pithist(m3, plot = FALSE, type = "random", scale = "normal", nsim = 10))
set.seed(seed); expect_silent(p3r_normal20 <- pithist(m3, plot = FALSE, type = "random", scale = "normal", nsim = 20))
set.seed(Sys.time())

fn <- function(x) with(x, x$simint_upr - x$simint_lwr)
expect_true(all(fn(p3r_normal) > fn(p3r_normal10)),
            info = "width of simulated interval with nsim1 > nsim10")
expect_true(all(fn(p3r_normal10) > fn(p3r_normal20)),
            info = "width of simulated interval with nsim10 > nsim20")

# ---------------------------------------------------------
# TODO: (RS2ML): Can the simulated interval be lower than 0???
if (any(p3r_normal[, c("simint_lwr", "simint_upr")] < 0)) warning("NEGATIVE SIMINT IN PITHIST!")
if (any(p3r_normal10[, c("simint_lwr", "simint_upr")] < 0)) warning("NEGATIVE SIMINT IN PITHIST!")
if (any(p3r_normal20[, c("simint_lwr", "simint_upr")] < 0)) warning("NEGATIVE SIMINT IN PITHIST!")

rm(p3r_normal10, p3r_normal20)

# same seed, same nsim, different simint_levels (confidence level)
set.seed(seed); expect_silent(p3r_normal94 <- pithist(m3, plot = FALSE, type = "random", scale = "normal", simint_level = 0.94))
set.seed(seed); expect_silent(p3r_normal96 <- pithist(m3, plot = FALSE, type = "random", scale = "normal", simint_level = 0.95))
expect_true(all(p3r_normal$simint_lwr <= p3r_normal94$simint_lwr),
            info = "lower bound of simulated interval for simint_level 0.95 <= 0.94")
expect_true(all(p3r_normal$simint_lwr >= p3r_normal96$simint_lwr),
            info = "lower bound of simulated interval for simint_level 0.95 >= 0.96")
expect_true(all(p3r_normal$simint_upr >= p3r_normal94$simint_upr),
            info = "upper bound of simulated interval for simint_level 0.95 >= 0.94")
expect_true(all(p3r_normal$simint_upr <= p3r_normal96$simint_upr),
            info = "upper bound of simulated interval for simint_level 0.95 <= 0.96")

rm(p3r_normal94, p3r_normal96, fn)

# same as p3r_normal but we set simint = FALSE at the same time to not plot simulated confidence
# interval. In this case simint_lwr and simint_upr should be empty while observed should come
# from the simulation (thus identical to p3r_normal)
set.seed(seed); expect_silent(p3r_normal_nosimint <- pithist(m3, plot = FALSE, type = "random", scale = "normal", class = "data.frame", simint = FALSE))
expect_true(all(is.na(p3r_normal_nosimint$simint_lwr)) && all(is.na(p3r_normal_nosimint$simint_upr)),
            info = "lower/upper bound must be NA if simint = FALSE")
expect_identical(p3r_normal[, 1:4], p3r_normal_nosimint[, 1:4],
            info = "observed, expected, mid, and width must be identical no matter what simint is set to")

rm(p3r_normal_nosimint)

## FIXME: fails on R-Forge
## set.seed(seed); expect_silent(tmp <- pithist(m3, plot = FALSE, nsim = 1, type = "random", scale = "normal", class = "data.frame", simint_nrep = 100))
## expect_true(!identical(p3r_normal, tmp),
##             info = "different number of repetitions for drawing quantiles - results not identical")
## expect_identical(sha1(tmp), "7883be6c78488cbee41b9add5206bf40697c5597")
## 
## rm(p3r_normal, tmp)
