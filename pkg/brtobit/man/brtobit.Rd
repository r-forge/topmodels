\name{brtobit}
\alias{brtobit}
\alias{brtobit_fit}
\alias{brtobit_control}

\alias{coef.brtobit}
\alias{vcov.brtobit}
\alias{fitted.brtobit}
\alias{predict.brtobit}

\alias{print.brtobit}
\alias{summary.brtobit}
\alias{print.summary.brtobit}
\alias{logLik.brtobit}
\alias{model.matrix.brtobit}
\alias{model.frame.brtobit}

\alias{bread.brtobit}
\alias{getSummary.brtobit}

\encoding{UTF-8}

\title{Bias-Reduced Tobit Regression}

\description{
  Fitting tobit regression models with bias-reduced estimation (rather than
  plain maximum likelihood). 
}
\usage{
brtobit(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE,
  control = brtobit_control(\dots), \dots)

brtobit_fit(x, y, control = brtobit_control())

brtobit_control(fsmaxit = 100, start = NULL, epsilon = 1e-08, type = "BR", ...)
}
\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x1 + x2} where
    \code{y} is the response and \code{x1} and \code{x2} are regressor variables
    for the location of the latent Gaussian distribution.}
  \item{data}{an optional data frame containing the variables occurring in the
    formulas.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y}{for \code{brtobit}: logical. If \code{TRUE} the model matrix and
    response vector used for fitting are returned as components of the returned value.
    For \code{brtobit_fit}: \code{x} is a design matrix with regressors for the
    location and \code{y} is a vector of observations.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
  \item{control, fsmaxit, start, epsilon}{a list of control parameters passed for the Fisher scoring
    optimization.}
  \item{type}{character. Should bias-reduced (BR) or plain maximum likelihood (ML)
    estimation be used?}
}

\details{
  \code{brtobit} fits tobit regression models with bias-reduced (BR) estimation
  as introduced by Köll et al. (2021). The model assumes an underlying latent Gaussian variable:

  \deqn{y_i^* \sim \mathcal{N}(\mu_i, \sigma^2)}{y* ~ N(mu, sigma^2)}
  
  which is only observed if positive and zero otherwise: \eqn{y_i = \max(0, y_i^*)}{y = max(0, y*)}.
  The latent mean \eqn{\mu_i}{mu} is linked to a linear predictor

  \deqn{\mu_i = x_i^\top \beta}{mu = x'b}

  and the latent variance \eqn{\sigma^2}{sigma^2} is assumed to be constant.

  \code{brtobit_fit} is the lower level function where the actual fitting takes place.

  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"brtobit"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{predict}}, 
  \code{\link[stats]{model.frame}}, \code{\link[stats]{model.matrix}},
  \code{\link[sandwich]{bread}} (from the \pkg{sandwich} package), and
  \code{\link[memisc]{getSummary}} (from the \pkg{memisc} package, enabling \code{\link[memisc]{mtable}}).
  
  In the future we intend to extend the implementation to heteroscedastic tobit
  models in the \code{\link[crch]{crch}} package (Messner, Mayr, Zeileis 2016).
}

\value{
  \code{brtobit} returns an object of class \code{"brtobit"}, i.e., a list with components as follows.
  \code{brtobit_fit} returns an unclassed list with components up to \code{converged}.
  \item{coefficients}{vector of estimated regression coefficients (plus the variance),}
  \item{bias}{bias estimate,}
  \item{vcov}{covariance matrix of all parameters in the model,}
  \item{loglik}{the log-likelihood of the fitted model,}
  \item{df}{number of estimated parameters,}
  \item{nobs}{number of observations,}
  \item{grad}{gradient vector,}
  \item{control}{list of control parameters,}
  \item{iterations}{number of iterations,}
  \item{converged}{logical indicating whether the Fisher scoring optimization converged,}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}  
  \item{terms}{terms objects for the model,}
  \item{levels}{levels of the categorical regressors,}
  \item{contrasts}{contrasts corresponding to \code{levels} from the
    respective models,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the numeric response vector (if \code{y = TRUE}),}
  \item{x}{model matrix (if \code{x = TRUE}).}
}

\references{
  Köll S, Kosmidis I, Kleiber C, Zeileis A (2021).
    \dQuote{Bias Reduction as a Remedy to the Consequences of Infinite Estimates in Poisson and Tobit Regression.}
    arXiv:2101.12345, arXiv.org E-Print Archive. \url{https://arxiv.org/abs/2101.12345}

  Messner JW, Mayr GJ, Zeileis A (2016). Heteroscedastic Censored and
    Truncated Regression with crch.
    \emph{The R Journal}, \bold{8}(1), 173--181.
    \url{https://journal.R-project.org/archive/2016-1/messner-mayr-zeileis.pdf}.
}

\seealso{\code{\link[crch]{crch}}}

\examples{
## artificial data generating process from Koell et al. (2021)
dgp <- function(n = 100, coef = c(1, 1, -10, 2), prob = 0.25) {
  x2 <- runif(n, -1, 1)
  x3 <- rbinom(n, size = 1, prob = ifelse(x2 > 0, prob, 1 - prob))
  y <- rnorm(n, mean = coef[1] + coef[2] * x2 + coef[3] * x3, sd = sqrt(coef[4]))
  y[y <= 0] <- 0
  data.frame(y, x2, x3)
}

set.seed(2020-10-29)
d <- dgp()

## models
m22_ml <- brtobit(y ~ x2 + x3, data = d, type = "ML", fsmaxit = 28)
m22_br <- brtobit(y ~ x2 + x3, data = d, type = "BR")
m2_all <- brtobit(y ~ x2, data = d, type = "ML")
m2_sub <- update(m2_all, subset = x3 == 0)

if(require("memisc")) {

## Table 2
mtable("ML" = m22_ml, "BR" = m22_br, "ML/sub" = m2_sub, "ML/SST" = m2_all,
  summary.stats = c("Log-likelihood", "N"))

}
}

\keyword{regression}
