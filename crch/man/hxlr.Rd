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
  hxlr(formula, scale = ~ NULL, data, weights = NULL, thresholds, ...)
}

\arguments{
  \item{formula}{A formula expression of the form response ~ predictors.
    Response has to be an (ordered) factor, usually derived through cutting a continuous
    variable into categories.}
  \item{scale}{Optional formula expression of the form ~ predictors. These predictors are
    used to predict the scale (dispersion) of the latent logistic distribution.}
  \item{data}{An optional data frame containing the variables occurring in the formulas.}
  \item{weights}{Optional case weights in fitting.}
  \item{thresholds}{Vector of (transformed) thresholds that were used to cut the continuous
    variable into categories. Data frames or matrices with multiple columns are allowed as well.
    Then each column is used as predictor variable for the intercept model.}
  \item{\dots}{Additional parameters passed to \code{\link[ordinal]{clm}}. Not all tested!}
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
  \item{coefficients}{List of coefficients for intercept model, predictors, and scale model.}
  \item{coefficients2}{List of converted coefficients for location and scale.}
  \item{optim}{Output from optimization.}
  \item{loglik}{Log-likelihood.}
  \item{start}{Starting values of coefficients used in the optimization.}
  \item{thresholds}{Thresholds.}
  \item{data}{Data frame containing the variables occurring in the formulas.}
  \item{formula}{Formula}
  \item{scale}{Formula expression for scale model.}
  \item{nobs}{Number of observations in data set.}
  \item{call}{Function call.}
  \item{converged}{Logical variable whether optimization has converged or not.}
  \item{vcov}{Covariance matrix.}
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
data("rrr_ibk")
rrr_ibk <- na.omit(rrr_ibk)
q <- unique(quantile(rrr_ibk$rrr, seq(0.1, 0.9, 0.1)))
rrr_ibk$rrr.cat <- cut(rrr_ibk$rrr, c(-Inf, q, Inf))

## fit ordinary extended logistic regression with ensemble mean as 
## predictor variable
XLR <- hxlr(rrr.cat ~ sqrtensmean, data = rrr_ibk, thresholds = sqrt(q))
## print
XLR
## summary
summary(XLR)


## fit ordinary extended logistic regression with ensemble mean 
## and standard deviation as predictor variables
XLRS <- hxlr(rrr.cat ~ sqrtensmean + sqrtenssd, data = rrr_ibk, 
thresholds = sqrt(q))
## fit heteroskedastic extended logistic regression with ensemble 
## standard deviation as predictor for the scale
HXLR <- hxlr(rrr.cat ~ sqrtensmean, scale = ~sqrtenssd, 
data = rrr_ibk, thresholds = sqrt(q))

## compare AIC of different models
AIC(XLR, XLRS, HXLR)

## XLRS and HXLR are nested in XLR -> likelihood ratiotests
library("lmtest")
lrtest(XLR, XLRS)
lrtest(XLR, HXLR)


##############################################################################
## Cross-validation and bootstrapping RPS for different models 
## (like in Messner 2013)
N <- NROW(rrr_ibk)
## function that returns model fits
fits <- function(data, weights = rep(1, N)) {
  list(
    "XLR"    = hxlr(rrr.cat ~ sqrtensmean,
data = data, weights = weights, thresholds = sqrt(q)),
    "XLR:S"  = hxlr(rrr.cat ~ sqrtensmean + sqrtenssd,
data = data, weights = weights, thresholds = sqrt(q)),
    "XLR:SM" = hxlr(rrr.cat ~ sqrtensmean + I(sqrtensmean*sqrtenssd), 
data = data, weights = weights, thresholds = sqrt(q)),
    "HXLR"   = hxlr(rrr.cat ~ sqrtensmean, scale = ~ sqrtenssd,
data = data, weights = weights, thresholds = sqrt(q)),
    "HXLR:S" = hxlr(rrr.cat ~ sqrtensmean + sqrtenssd, scale = ~ sqrtenssd, 
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
  testdata <- rrr_ibk[testIndex,]
  ## observations    
  obs <- c(obs, rrr_ibk$rrr.cat[testIndex])
  ## estimation
  modelfits <- fits(rrr_ibk, weights)
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
RPS <- lapply(pred, rps, obs = obs)

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

\keyword{regression}
