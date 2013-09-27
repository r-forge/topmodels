\name{hxlr}
\encoding{latin1}

\alias{hxlr}
\alias{print.hxlr}
\alias{summary.hxlr}
\alias{print.summary.hxlr}
\alias{logLik.hxlr}

\title{Heteroskedastic Extended Logistic Regression}

\description{
  This is a wrapper function for \code{\link[ordinal]{clm}} (from package \pkg{ordinal}) to fit
  (heteroskedastic) extended logistic regression (HXLR) models (Messner et al. 2013).
}
\usage{
  hxlr(formula, data, subset, na.action, weights, thresholds, control, ...)
}
\arguments{
  \item{formula}{A formula expression of the form y ~ x | z where y is the response, x are predictor variables for the location and z are predictor variables for the scale of the latent distribution. Response can either be a continuous variable or a factor.}
  \item{data}{An optional data frame containing the variables occurring in the formulas.}
  \item{subset}{An optional vector specifying a subset of observations to be used in the fitting process.}
  \item{na.action}{A function which indicates what should happen when the data contain ‘NA’s. Default is na.omit}
  \item{weights}{Optional case weights in fitting.}
  \item{thresholds}{Vector of (transformed) thresholds that are used to cut the continuous response into categories. Data frames or matrices with multiple columns are allowed as well. Then each column is used as predictor variable for the intercept model.}
  \item{control}{A list of control parameters passed to \code{\link{optim}}. Default is \code{\link{control.hxlr}}}
  \item{\dots}{Additional parameters passed to \code{\link[ordinal]{clm}}. }
}

\details{
Extended logistic regression (Wilks 2009) extends binary logistic regression to multi-category
responses by including the thresholds, that are used to cut a continuous variable into categories,
in the regression equation. Heteroskedastic extended logistic regression (Messner et al. 2013)
extends this model further and allows to add additional predictor variables that are used to
predict the scale of the latent logistic distribution.
}

\value{
  An object of class \code{"hxlr"}, i.e., a list with the following elements.
  \item{coefficients}{List of converted coefficients for location and scale.}
  \item{coefficients.CLM}{List of CLM coefficients for intercept model, predictors, and scale model.}
  \item{fitted.values}{List of fitted location and scale parameters.}
  \item{optim}{Output from optimization.}
  \item{method}{Optimization method used for \code{\link{optim}}}
  \item{method}{list of control parameters used for \code{\link{optim}}}
  \item{start}{Starting values of coefficients used in the optimization.}
  \item{weights}{case weights used for fitting.}
  \item{nobs}{Number of observations in data set.}
  \item{loglik}{Log-likelihood.}
  \item{vcov}{Covariance matrix.}
  \item{converged}{Logical variable whether optimization has converged or not.}
  \item{iterations}{Number of iterations in optimization.}
  \item{call}{Function call.}
  \item{scale}{Formula expression for scale model.}
  \item{terms}{The terms objects used.}
  \item{levels}{list of levels of the factors used in fitting for location and scale respectively.}
  \item{thresholds}{Thresholds.}
}

\references{
  Messner JW, Mayr GJ, Zeileis A, Wilks DS (2013).
    Extending Extended Logistic Regression to Effectively Utilize the Ensemble Spread.
    Working Paper 2013-21. Working Papers in Economics and Statistics,
    Research Platform Empirical and Experimental Economics, Universität Innsbruck.
    \url{http://EconPapers.RePEc.org/RePEc:inn:wpaper:2013-21}

  Wilks DS (2009).
    Extending Logistic Regression to Provide Full-Probability-Distribution MOS Forecasts.
    \emph{Meteorological Applications}, \bold{368}, 361--368.
}

\seealso{\code{\link{predict.hxlr}}, \code{\link[ordinal]{clm}}}

\examples{
data("rainIbk")
q <- unique(quantile(rainIbk$rain, seq(0.1, 0.9, 0.1)))

## fit ordinary extended logistic regression with ensemble mean as 
## predictor variable
XLR <- hxlr(sqrt(rain) ~ sqrtensmean, data = rainIbk, thresholds = sqrt(q))
## print
XLR
## summary
summary(XLR)


## fit ordinary extended logistic regression with ensemble mean 
## and standard deviation as predictor variables
XLRS <- hxlr(sqrt(rain) ~ sqrtensmean + sqrtenssd, data = rainIbk, 
thresholds = sqrt(q))
## fit heteroskedastic extended logistic regression with ensemble 
## standard deviation as predictor for the scale
HXLR <- hxlr(sqrt(rain) ~ sqrtensmean | sqrtenssd, 
data = rainIbk, thresholds = sqrt(q))

## compare AIC of different models
AIC(XLR, XLRS, HXLR)

## XLRS and HXLR are nested in XLR -> likelihood ratiotests
library("lmtest")
lrtest(XLR, XLRS)
lrtest(XLR, HXLR)

\dontrun{
##############################################################################
## Cross-validation and bootstrapping RPS for different models 
## (like in Messner 2013)
N <- NROW(rainIbk)
## function that returns model fits
fits <- function(data, weights = rep(1, N)) {
  list(
    "XLR"    = hxlr(sqrt(rain) ~ sqrtensmean,
data = data, weights = weights, thresholds = sqrt(q)),
    "XLR:S"  = hxlr(sqrt(rain) ~ sqrtensmean + sqrtenssd,
data = data, weights = weights, thresholds = sqrt(q)),
    "XLR:SM" = hxlr(sqrt(rain) ~ sqrtensmean + I(sqrtensmean*sqrtenssd), 
data = data, weights = weights, thresholds = sqrt(q)),
    "HXLR"   = hxlr(sqrt(rain) ~ sqrtensmean | sqrtenssd,
data = data, weights = weights, thresholds = sqrt(q)),
    "HXLR:S" = hxlr(sqrt(rain) ~ sqrtensmean + sqrtenssd | sqrtenssd, 
data = data, weights = weights, thresholds = sqrt(q))
  )
}


## cross validation
id <- sample(1:10, N, replace = TRUE)
obs <- NULL
pred <- list(NULL)
for(i in 1:10) {
  ## splitting into test and training data set
  trainIndex <- which(id != i)     
  testIndex <- which(id == i)                        
  ## weights that are used for fitting the models
  weights <- as.numeric(table(factor(trainIndex, levels = c(1:N))))
  ## testdata
  testdata <- rainIbk[testIndex,]
  ## observations    
  obs <- c(obs, rainIbk$rain[testIndex])
  ## estimation
  modelfits <- fits(rainIbk, weights)
  ## Prediction
  pred2 <- lapply(modelfits, predict, newdata = testdata, type = "cumprob")
  pred <- mapply(rbind, pred, pred2, SIMPLIFY = FALSE)
}
names(pred) <- c(names(modelfits))

## function to compute RPS
rps <- function(pred, obs) {
  OBS <- NULL
  for(i in 1:N) OBS <- rbind(OBS, rep(0:1, c(obs[i] - 1, length(q) - obs[i] + 1))) 
  apply((OBS-pred)^2, 1, sum)
}
## compute rps
RPS <- lapply(pred, rps, obs = as.numeric(cut(obs, c(-Inf, q, Inf))))

## bootstrapping mean rps 
rpsall <- NULL
for(i in 1:250) {
  index <- sample(length(obs), replace = TRUE)
  rpsall <- rbind(rpsall, sapply(RPS, function(x) mean(x[index])))
}
  
rpssall <- 1 - rpsall/rpsall[,1]
boxplot(rpssall[,-1], ylab = "RPSS", main = "RPSS relative to XLR")
abline(h = 0, lty = 2)
}
}

\keyword{regression}
