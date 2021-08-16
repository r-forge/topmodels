# -------------------------------------------------------------------
# - NAME:   demonstrate_misspecifications.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2021-08-03
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2021-08-16 on thinkmoritz
# -------------------------------------------------------------------
library("topmodels")
library("crch")
library("rmutil")  # for laplace dist
library("sn") # for skewed N dist

# -------------------------------------------------------------------
# (1) MISSPECIFIED RESPONSE: POISSON vs NEGBIN
# -------------------------------------------------------------------
## artificial data from negative binomial (mu = 3, theta = 2)
## and Poisson (mu = 3) distribution
set.seed(1090)
y <- rnbinom(100, mu = 3, size = 2)
x <- rpois(100, lambda = 3)

plot(dpois(seq(1, 10), lambda = 3), type = "l")
lines(dnbinom(seq(1, 10), mu = 3, size = 2), type = "l", col = 2)
legend("topright", c("y ~ nbinom", "x ~ rpois"), lty = 1, col = c(1, 2))

## NEW
d <- data.frame(x = runif(500))
d$mu <- 0.5 + 1 * d$x
d$yp <- rpois(500, lambda = d$mu)
d$ynb <- rnbinom(500, mu = d$mu, size = 1)
m1 <- glm(yp ~ x, data = d, family = poisson)
m2 <- glm(ynb ~ x, data = d, family = poisson)
m3 <- glm.nb(ynb ~ x, data = d)
##

## correctly specified Poisson model fit (mu = 3.34)
m1_corr <- glm(x ~ 1, family = poisson)

## misspecified Poisson model fit (mu = 3.32):
## wrong assumption of underlying distribution -> underdispersive model fit (residuals)
## compare Kleiber and Zeileis: 
## “wave-like” pattern around the horizontal reference line: the data exhibit too
## many small counts, notably zeros, as well as too many large counts for a
## Poisson model to provide an adequate fit. In summary, the patterns encountered
## in the bottom row of Figure 1 reflect a substantial amount of overdispersion
## that is not captured by the fitted Poisson distribution.
## [U-shaped PIT, S-shaped QQ-plot]
m1_missp <- glm(y ~ 1, family = poisson)
     
par(mfrow = c(2, 5))
topmodels(m1_corr, single_page = TRUE, spar = FALSE)
topmodels(m1_missp, single_page = TRUE, spar = FALSE)
par(mfrow = c(1, 1))

par(mfrow = c(3, 5))
topmodels(m1, single_page = TRUE, spar = FALSE)
topmodels(m2, single_page = TRUE, spar = FALSE)
topmodels(m3, single_page = TRUE, spar = FALSE)
par(mfrow = c(1, 1))


# -------------------------------------------------------------------
# (2) MISSPECIFIED RESPONSE: HEAVY TAILS (N vs t_4)
# -------------------------------------------------------------------
## artificial data from a t_4 distribution
set.seed(1090)
y <- rt(1000, 4)

curve(dt(x, df = 4), from = -6, to = 6, ylim = c(0, 0.5))
curve(dnorm, from = -6, to = 6, add = TRUE, col = 2)
legend("topleft", c("y ~ t_4", "norm"), lty = 1, col = c(1, 2))

## correct t fit
m2_corr <- crch(y ~ 1, dist = "student")

## incorrect normal fit
## wrong assumpion of underlying distribution: tails too light -> overdispersive model fit (residuals)
## [inverse U-shaped PIT, inverse S-shaped QQ-plot]
m2_missp <- crch(y ~ 1, dist = "gaussian")

root2_corr <- rootogram(m2_corr,  breaks = 40, plot = FALSE)
root2_missp <- rootogram(m2_missp,  breaks = 40, plot = FALSE)
plot(c(root2_corr, root2_missp), xlim = c(-6, 6), ylim = c(-2, 14))
par(mfrow = c(2, 5))
topmodels(m2_corr, single_page = TRUE, spar = FALSE)
topmodels(m2_missp, single_page = TRUE, spar = FALSE)
par(mfrow = c(1, 1))


# -------------------------------------------------------------------
# (3) MISSPECIFIED RESPONSE: OVER-/UNDERDISPERSION
# -------------------------------------------------------------------
## artificial data from a Gaussian, Uniform
## and Laplace (mu = 0, s = 0.4) distribution
set.seed(1090)
y1 <- rnorm(1000)
y2 <- runif(1000, min = -6, max = 6)
y3 <- rlaplace(1000, m = 0, s = 0.4)

curve(dnorm, from = -6, to = 6, ylim = c(0, 1))
curve(dunif(x, min = -6, max = 6), from = -6, to = 6, add = TRUE, col = 2)
curve(dlaplace, from = -6, to = 6, add = TRUE, col = 3)
legend("topleft", c("y1 ~ norm", "y2 ~ unif", "y3 ~ laplace"), lty = 1, col = c(1, 2, 3))

## correct gaussian fit
m3_corr <- crch(y1 ~ 1, dist = "gaussian")

## incorrect fit (underdispersive model fit)
## [U-shaped PIT, S-shaped QQ-plot, "thin tails" worm plot]
m3_missp_u <- crch(y2 ~ 1, dist = "gaussian")

## incorrect fit (overdispersive model fit)
## [inverse U-shaped PIT, inverse S-shaped QQ-plot, "fat tails" worm plot]
m3_missp_lp <- crch(y3 ~ 1, dist = "gaussian")

par(mfrow = c(3, 5))
topmodels(m3_corr, single_page = TRUE, spar = FALSE)
topmodels(m3_missp_u, single_page = TRUE, spar = FALSE)
topmodels(m3_missp_lp, single_page = TRUE, spar = FALSE)
par(mfrow = c(1, 1))


# -------------------------------------------------------------------
# (4) MISSPECIFIED RESPONSE: SKEWED DISTRIBUTION
# -------------------------------------------------------------------
## artificial data from a Gaussian, Exponential and Negative Exponential
## distribution
set.seed(1090)
#y1 <- rnorm(1000, mean = 3)
#y2 <- rweibull(1000, shape = 1.5)
#y3 <- rweibull(1000, shape = 12)
d4 <- data.frame(yn = rnorm(1000, mean = 0))
d4$yn_rs <- rsn(1000, dp = cp2dp(c(0, 1, 0.9), family = "SN"))
d4$yn_ls <- rsn(1000, dp = cp2dp(c(0, 1, -0.9), family = "SN"))

#curve(dnorm(x, mean = 2.5), from = -2, to = 8, ylim = c(0, 1))
#curve(dweibull(x, shape = 1.5), from = -2, to = 8, col = 2, add = TRUE)
#curve(dweibull(x, shape = 12, scale = 6), from = -2, to = 8, col = 3, add = TRUE)
#legend("topleft", c("y1 ~ norm", "y2 ~ weibull_1.5", "y3 ~ weibull_12"), lty = 1, col = c(1, 2, 3))

curve(dnorm(x, mean = 0), from = -8, to = 8, ylim = c(0, 1), col = 1)
curve(dsn(x, dp = cp2dp(c(0, 1, 0.9), family = "SN")), from = -8, to = 8, add = TRUE, col = 2)
curve(dsn(x, dp = cp2dp(c(0, 1, -0.9), family = "SN")), from = -8, to = 8, add = TRUE, col = 3)
legend("topleft", c("y1 ~ norm", "y2 ~ right skewed", "y3 ~ left skewed"), lty = 1, col = c(1, 2, 3))

## correct gaussian fit
m4_corr <- crch(yn ~ 1, data = d4, dist = "gaussian")

## incorrect fit: right-skewed residuals
## [curved (positive skewed) QQ-Plot, U-shape wormplot]
m4_missp1 <- crch(yn_rs ~ 1, data = d4, dist = "gaussian")

## incorrect fit: left-skewed residuals
## [curved (negative skewed) QQ-Plot, inverse U-shape wormplot]
m4_missp2 <- crch(yn_ls ~ 1, data = d4, dist = "gaussian")

par(mfrow = c(3, 4))
topmodels(m4_corr, which = c(1, 2, 4, 5), single_page = TRUE, spar = FALSE)
topmodels(m4_missp1, which = c(1, 2, 4, 5), single_page = TRUE, spar = FALSE)
topmodels(m4_missp2, which = c(1, 2, 4, 5), single_page = TRUE, spar = FALSE)
par(mfrow = c(1, 1))

# -------------------------------------------------------------------
# (4b) MISSPECIFIED RESPONSE: SKEWED DISTRIBUTION
# -------------------------------------------------------------------
library("topmodels")
library("crch")
library("sn")
set.seed(1090)
d4 <- data.frame(yn = rnorm(1000, mean = 0))
d4$yn_rs <- rsn(1000, dp = cp2dp(c(0, 1, 0.9), family = "SN"))
d4$yn_ls <- rsn(1000, dp = cp2dp(c(0, 1, -0.9), family = "SN"))

curve(dnorm(x, mean = 0), from = -8, to = 8, ylim = c(0, 1), col = 1)
curve(dsn(x, dp = cp2dp(c(0, 1, 0.9), family = "SN")), from = -8, to = 8, add = TRUE, col = 2)
curve(dsn(x, dp = cp2dp(c(0, 1, -0.9), family = "SN")), from = -8, to = 8, add = TRUE, col = 3)
legend("topleft", c("y1 ~ norm", "y2 ~ right skewed", "y3 ~ left skewed"), lty = 1, col = c(1, 2, 3))

m4_corr <- crch(yn ~ 1, data = d4, dist = "gaussian")
m4_missp1 <- crch(yn_rs ~ 1, data = d4, dist = "gaussian")
m4_missp2 <- crch(yn_ls ~ 1, data = d4, dist = "gaussian")

par(mfrow = c(3, 4))
topmodels(m4_corr, which = c(1, 2, 4, 5), single_page = TRUE, spar = FALSE)
topmodels(m4_missp1, which = c(1, 2, 4, 5), single_page = TRUE, spar = FALSE)
topmodels(m4_missp2, which = c(1, 2, 4, 5), single_page = TRUE, spar = FALSE)
par(mfrow = c(1, 1))

# -------------------------------------------------------------------
# (4c) MISSPECIFIED RESPONSE: SKEWED DISTRIBUTION
# -------------------------------------------------------------------
library("topmodels")
library("crch")
library("sn")

## artificial data from Normal distribution: no skewness, right skewed, and left skewed
set.seed(1090)
d4 <- data.frame(
  x = runif(1000),
  z <- runif(1000)
)
d4$mu <- 0.1 + 0.85 * d4$x
d4$sigma <- exp(0.7 + 0.1 * d4$z)

d4$yn <- rsn(n = 1000, xi = d4$mu, omega = d4$sigma, alpha = 0, tau = 0)
d4$yn_rs <- rsn(n = 1000, xi = d4$mu, omega = d4$sigma, alpha = 3, tau = 0)
d4$yn_ls <- rsn(n = 1000, xi = d4$mu, omega = d4$sigma, alpha = -3, tau = 0)

curve(dsn(x, xi = 0, omega = 1, alpha = 0, tau = 0), from = -8, to = 8, col = 1, ylim = c(0, 1))
curve(dsn(x, xi = 0, omega = 1.6, alpha = 3, tau = 0), from = -8, to = 8, add = TRUE, col = 2)
curve(dsn(x, xi = 0, omega = 1.6, alpha = -3, tau = 0), from = -8, to = 8, add = TRUE, col = 3)
legend("topleft", c("y1 ~ norm", "y2 ~ right skewed", "y3 ~ left skewed"), lty = 1, col = c(1, 2, 3))

## correct gaussian fit
m4_corr <- crch(yn ~ x | z, data = d4, dist = "gaussian")

## incorrect fit: right-skewed residuals
## [curved (positive skewed) QQ-Plot, U-shape wormplot]
m4_missp1 <- crch(yn_rs ~ x | z, data = d4, dist = "gaussian")

## incorrect fit: left-skewed residuals
## [curved (negative skewed) QQ-Plot, inverse U-shape wormplot]
m4_missp2 <- crch(yn_ls ~ x | z, data = d4, dist = "gaussian")

par(mfrow = c(3, 4))
topmodels(m4_corr, which = c(1, 2, 4, 5), single_page = TRUE, spar = FALSE)
topmodels(m4_missp1, which = c(1, 2, 4, 5), single_page = TRUE, spar = FALSE)
topmodels(m4_missp2, which = c(1, 2, 4, 5), single_page = TRUE, spar = FALSE)
par(mfrow = c(1, 1))

# -------------------------------------------------------------------
# (5) MANUS'S EXAMPLE
# -------------------------------------------------------------------
ens_mean <- rnorm(5000, 0.35, 6.91)
ens_logsd <- rnorm(5000, -0.56, 0.43)

mean <- 6.5 + 1 * ens_mean
logsd <- 0.9 + 1.3 * ens_logsd
y <- rlogis(5000, mean, exp(logsd))

d <- data.frame(y = y, x1 = ens_mean, z1 = ens_logsd)

curve(dnorm, from = -5, to = 5)
curve(dlogis, from = -5, to = 5, col = 2, add = TRUE)
curve(dt(x, df = 2), from = -5, to = 5, col = 3, add = TRUE)
legend("topleft", c("dnorm", "dlogis", "dt_2"), lty = 1, col = c(1, 2, 3))

## wrong assumption of underlying distribution: tails too light -> overdispersive model fit (residuals)
m5_gauss <- crch(y ~ x1 | z1, dist = "gaussian", data = d)

## correct model fit
m5_logis <- crch(y ~ x1 | z1, dist = "logistic", data = d)

## wrong assumption of underlying distribution, but almost not visible
m5_student <- crch(y ~ x1 | z1, dist = "student", data = d)

par(mfrow = c(3, 5))
topmodels(m5_gauss, single_page = TRUE, spar = FALSE)
topmodels(m5_logis, single_page = TRUE, spar = FALSE)
topmodels(m5_student, single_page = TRUE, spar = FALSE)
par(mfrow = c(1, 1))

# -------------------------------------------------------------------
# (6) LEAVING OUT RELEVANT REGRESSORS: MISSPECIFIED MEAN AND VARIANCE
# -------------------------------------------------------------------
ens_mean <- rnorm(5000, 10, 6.91)
ens_logsd <- rnorm(5000, -0.56, 0.43)

mean <- 6.5 + 4 * ens_mean
logsd <- 0.9 + 1.3 * ens_logsd
y <- rnorm(5000, mean, exp(logsd))

d <- data.frame(y = y, x1 = ens_mean, z1 = ens_logsd)

## correct model fit
m6_corr <- crch(y ~ x1 | z1, dist = "gaussian", data = d)

## misspecified mean
m6_missp1 <- crch(y ~ 1 | z1, dist = "gaussian", data = d)

## misspecified variance
m6_missp2 <- crch(y ~ x1 | 1, dist = "gaussian", data = d)

par(mfrow = c(3, 5))
topmodels(m6_corr, single_page = TRUE, spar = FALSE)
topmodels(m6_missp1, single_page = TRUE, spar = FALSE)
topmodels(m6_missp2, single_page = TRUE, spar = FALSE)
par(mfrow = c(1, 1))

