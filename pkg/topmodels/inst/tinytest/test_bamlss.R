# --------------------------------------------------------------------
# Testing bamlss family distrubiton3 functionality
# --------------------------------------------------------------------

# TODO(R): Possible improvements
# - Currently families accept invalid parameters (e.g.,
#   the gaussian_bamlss accepts negative sigmas)


if (interactive()) { library("devtools"); library("tinytest"); library("topmodels") }

suppressPackageStartupMessages(library("crch"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("bamlss"))

#     ALD_bamlss(..., tau = 0.5, eps = 0.01)
#     beta_bamlss(...)
#     binomial_bamlss(link = "logit", ...)
#     cnorm_bamlss(...)
#     cox_bamlss(...)
#     dw_bamlss(...)
#     DGP_bamlss(...)
#     dirichlet_bamlss(...)
#     ELF_bamlss(..., tau = 0.5)
#     gaussian_bamlss(...)
#     gaussian2_bamlss(...)
#     Gaussian_bamlss(...)
#     gamma_bamlss(...)
#     logNN_bamlss(...)
#     multinomial_bamlss(...)
#     mvnorm_bamlss(k = 2, ...)
#     mvnormAR1_bamlss(k = 2, ...)
#     poisson_bamlss(...)
#     gpareto_bamlss(...)
#     glogis_bamlss(...)
#     AR1_bamlss(...)
#     beta1_bamlss(ar.start, ...)
#     nbinom_bamlss(...)
#     ztnbinom_bamlss(...)
#     lognormal_bamlss(...)
#     weibull_bamlss(...)
#     Sichel_bamlss(...)
#     GEV_bamlss(...)
#     gumbel_bamlss(...)
#     mix_bamlss(f1, f2, ...)
#     ZANBI_bamlss(f1, f2, ...)


# --------------------------------------------------------------------
# gaussian_bamlss: Normal distribution
# --------------------------------------------------------------------

basic_check <- function(name, par) {
    for (n in c(name, gsub("_bamlss$", "", name))) {
        expect_silent(bamlss::BAMLSS(n, par))
        expect_silent(bamlss::BAMLSS(bamlss::bamlss.family(n), par))
    }
}
f <- function(x) {
    expect_identical(x, 3L)
}
f(3L)
f(3)
expect_identical(3, 3)

# Set of parameters for testing
par <- list(mu = c(4, 5, 6), sigma = c(1, 1.5, 10))
expect_identical(names(par), c("mu", "sigma"))
expect_true(length(par) == 2L && all(sapply(par, is.numeric)))

# Creating distribution object
basic_check("gaussian_bamlss", par)

## Expecting errors
#expect_error(bamlss::BAMLSS(f),
#             info = "Parameters missing")
#expect_error(bamlss::BAMLSS(f, list()),
#             info = "Parameters missing")
#expect_error(bamlss::BAMLSS(f, list(sigma = 3)),
#             info = "Parameter mu missing")
#expect_error(bamlss::BAMLSS(f, list(mu = 3)),
#             info = "Parameter sigma missing")
#expect_error(bamlss::BAMLSS(f, c(mu = 3, sigma = 1)),
#             info = "Parameters not as list")
#
## Proper use
#expect_silent(bamlss::BAMLSS("gaussian", par),
#              info = "Creating bamlss distribution object")
#expect_silent(bamlss::BAMLSS("gaussian_bamlss", par),
#              info = "Creating bamlss distribution object")
#expect_silent(d <- bamlss::BAMLSS(f, par),
#              info = "Creating bamlss distribution object") # Store distribution object
#
## Testing return
#expect_identical(class(d), c("BAMLSS", "distribution"))
#expect_identical(dim(d), NULL)
#expect_identical(length(d), length(par[[1]]))
#
## Testing expectation, variance, dprq
#expect_silent(rval <- list(mean = mean(d), median = median(d),
#                           variance = variance(d)))
#expect_identical(rval$mean, par$mu)
#expect_identical(rval$mean, par$mu)
#expect_equal(rval$variance, par$sigma^2)
#
#
#
#
#
#
