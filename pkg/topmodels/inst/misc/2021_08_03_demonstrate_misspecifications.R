# -------------------------------------------------------------------
# - NAME:   demonstrate_misspecifications.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2021-08-03
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2021-08-03 on thinkmoritz
# -------------------------------------------------------------------
library("topmodels")
library("crch")
library("rmutil")

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
## [U-shaped PIT, S-shaped QQ-plot]
m3_missp_u <- crch(y2 ~ 1, dist = "gaussian")

## incorrect fit (overdispersive model fit)
## [inverse U-shaped PIT, inverse S-shaped QQ-plot]
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
y1 <- rnorm(1000, mean = 3)
y2 <- rweibull(1000, shape = 1.5)
y3 <- rweibull(1000, shape = 12)

curve(dnorm(x, mean = 2.5), from = -2, to = 8, ylim = c(0, 1))
curve(dweibull(x, shape = 1.5), from = -2, to = 8, col = 2, add = TRUE)
curve(dweibull(x, shape = 12, scale = 6), from = -2, to = 8, col = 3, add = TRUE)
legend("topleft", c("y1 ~ norm", "y2 ~ weibull_1.5", "y3 ~ weibull_12"), lty = 1, col = c(1, 2, 3))

## correct gaussian fit
m4_corr <- crch(y1 ~ 1, dist = "gaussian")

## incorrect fit: right-skewed residuals
## [curved (positive skewed) QQ-Plot]
m4_missp1 <- crch(y2 ~ 1, dist = "gaussian")

## incorrect fit: left-skewed residuals
## [curved (negative skewed) QQ-Plot]
m4_missp2 <- crch(y3 ~ 1, dist = "gaussian")

par(mfrow = c(3, 5))
topmodels(m4_corr, single_page = TRUE, spar = FALSE)
topmodels(m4_missp1, single_page = TRUE, spar = FALSE)
topmodels(m4_missp2, single_page = TRUE, spar = FALSE)
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
ens_mean <- rnorm(5000, 0.35, 6.91)
ens_logsd <- rnorm(5000, -0.56, 0.43)

mean <- 6.5 + 1 * ens_mean
logsd <- 0.9 + 1.3 * ens_logsd
y <- rlogis(5000, mean, exp(logsd))

d <- data.frame(y = y, x1 = ens_mean, z1 = ens_logsd)

## correct model fit
m6_corr <- crch(y ~ x1 | z1, dist = "logistic", data = d)

## misspecified mean
m6_missp1 <- crch(y ~ 1 | z1, dist = "logistic", data = d)

## misspecified variance
m6_missp2 <- crch(y ~ x1 | 1, dist = "logistic", data = d)

par(mfrow = c(3, 5))
topmodels(m6_corr, single_page = TRUE, spar = FALSE)
topmodels(m6_missp1, single_page = TRUE, spar = FALSE)
topmodels(m6_missp2, single_page = TRUE, spar = FALSE)
par(mfrow = c(1, 1))

