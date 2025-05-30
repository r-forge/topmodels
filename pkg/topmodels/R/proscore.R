#' Scoring Probabilistic Forecasts
#' 
#' Generic function and default method for computing various kinds of scores
#' for fitted or predicted probability distributions from (regression) models.
#' 
#' The function \code{proscore} provides a unified framework for scoring
#' probabilistic forecasts (in-sample or out-of-sample). The following scores
#' are currently available, using the following notation: \eqn{Y} is the predicted
#' random variable with cumulative distribution function \eqn{F(\cdot)}{F(.)} and probability
#' density function \eqn{f(\cdot)}{f(.)}. The corresponding expectation and variance are denoted
#' by \eqn{E(Y)} and \eqn{V(Y)}. The actual observation is \eqn{y}.
#'
#' \bold{Log-score:} Also known as logarithmic score. This is the negative log-likelihood
#' where the negative sign has the effect that smaller values indicate a better fit.
#'
#' \deqn{
#'   -\log f(y)
#' }{
#'   -log(f(y))
#' }
#'
#' \bold{Log-likelihood:} Also known as log-density. Clearly, this is equivalent to the
#' log-score above but using the conventional sign where bigger values indicate a better fit.
#'
#' \deqn{
#'   \log f(y)
#' }{
#'   log(f(y))
#' }
#'
#' \bold{Continuous ranked probability score (CRPS):}
#'
#' \deqn{
#'   \int_{-\infty}^{\infty} \left( F(x) - 1(x \ge y) \right)^2 \mbox{d} x
#' }{
#'   int { F(x) - 1(x >= y) }^2 dx
#' }
#'
#' where \eqn{1(\cdot)}{1(.)} denotes the indicator function.
#'
#' In case of a discrete rather than a continuous distribution, the ranked probability
#' score (RPS) is defined analogously using the sum rather than the integral. In other
#' words it is then the sum of the squared deviations between the predicted cumulative
#' probabilities \eqn{F(x)} and the ideal step function for the actual observation \eqn{y}.
#'
#' \bold{Mean absolute error (MAE):}
#'
#' \deqn{
#'   \left| y - E(Y) \right|
#' }{
#'   |y - E(Y)|
#' }
#'
#' \bold{Mean squared error (MSE):}
#'
#' \deqn{
#'   \left( y - E(Y) \right)^2
#' }{
#'   (y - E(Y))^2
#' }
#'
#' \bold{Dawid-Sebastiani score (DSS):}
#'
#' \deqn{
#'   \frac{\left( y - E(Y) \right)^2}{V(Y)} + \log(V(Y))
#' }{
#'   (y - E(Y))^2/V(Y) + log(V(Y))
#' }
#'
#' Internally, the default \code{proscore} method first computes the fitted/predicted
#' probability distribution object using \code{\link[distributions3]{prodist}}
#' (corresponding to \eqn{Y} above) and then obtains the corresponding observation
#' \eqn{y} using \code{\link{newresponse}}. Subsequently, the scores are evaluated
#' using either the \code{\link[distributions3]{log_pdf}} method,
#' \code{\link[scoringRules]{crps}} method, or simply the \code{mean}. Finally,
#' the resulting individual scores per observation can be returned as a full
#' data frame, or aggregated (e.g., by using \code{mean}, \code{sum}, or \code{summary}, etc.).
#'
#' @aliases proscore
#' @param object a fitted model object. For the \code{default} method this
#' needs to have a \code{\link[distributions3]{prodist}} and a
#' \code{\link{newresponse}} method.
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict and from which to obtain the response variable.
#' If omitted, the original observations are used.
#' @param na.action function determining what should be done with missing
#' values in \code{newdata}.  The default is to employ \code{NA}.
#' @param type character specifying the type of score to compute. Avaible types:
#' \code{"logs"} (or equivalently \code{"log-score"}),
#' \code{"loglikelihood"} (or equivalently \code{"log_pdf"}),
#' \code{"CRPS"} (or equivalently \code{"RPS"}), \code{"MAE"}, \code{"MSE"},
#' \code{"DSS"} (or equivalently \code{"Dawid-Sebastiani"}).
#' Upper or lower case spellings can be used interchangably, hyphens or underscores
#' can be included or omitted. Setting \code{type = NULL} yields all available scores.
#' @param aggregate logical or function to be used for aggregating scores across
#' observations. Setting \code{aggregate = TRUE} (the default) corresponds to using \code{mean}.
#' @param drop logical. Should scores be returned in a data frame (default)
#' or (if possible) dropped to a vector?
#' @param \dots further parameters passed to the \code{aggregate} function (if any).
#'
#' @return Either a \code{data.frame} of scores (if \code{drop = FALSE}, default) or
#' a named numeric vector (if \code{drop = TRUE} and the scores are not a matrix).
#' The names are the \code{type} specified by the user (i.e., are not canonicalized
#' by partial matching etc.).
#'
#' @keywords regression
#' @examples
#' ## Poisson regression model for FIFA 2018 data:
#' ## number of goals scored by each team in each game, explained by
#' ## predicted ability difference of the competing teams
#' data("FIFA2018", package = "distributions3")
#' m <- glm(goals ~ difference, data = FIFA2018, family = poisson)
#' 
#' ## default: in-sample mean log-score and CRPS
#' proscore(m)
#'
#' ## element-wise score using a new data set (final of the tournament)
#' final <- tail(FIFA2018, 2)
#' proscore(m, newdata = final, aggregate = FALSE)
#'
#' ## replicate in-sample log-likelihood
#' proscore(m, type = "loglik", aggregate = sum)
#' logLik(m)
#'
#' ## compute mean of all available scores
#' proscore(m, type = NULL)
#'
#' ## upper vs. lower case spelling is matched internally but preserved in output
#' proscore(m, type = c("logs", "crps"))
#' proscore(m, type = c("Log-score", "CRPS"))
#'
#' ## least-squares regression for speed and breaking distance of cars
#' data("cars", package = "datasets")
#' m <- lm(dist ~ speed, data = cars)
#'
#' ## replicate in-sample log-likelihood and residual sum of squares
#' ## (aka deviance) by taking the sum (rather than the mean) of the
#' ## log-density and squared errors, respectively
#' proscore(m, type = c("loglik", "MSE"), aggregate = sum)
#' logLik(m)
#' deviance(m)

#' @export 
proscore <- function(object, newdata = NULL, ...) {
  UseMethod("proscore")
}

#' @rdname proscore
#' @importFrom distributions3 prodist log_pdf
#' @export 
proscore.default <- function(object, newdata = NULL, na.action = na.pass, type = c("logs", "crps"), aggregate = TRUE, drop = FALSE, ...)
{
  ## match type
  if (is.null(type)) type <- c("logs", "loglikelihood", "crps", "mae", "mse", "dss")
  otype <- type
  type <- gsub("-|_", "", tolower(type))
  type <- sapply(type, match.arg, c(
    "logs", "logscore",
    "loglikelihood", "logpdf",
    "crps", "rps",
    "mae",
    "mse",
    "dss", "dawidsebastiani"))
  if(any(type %in% c("logscore"))) type[type %in% c("logscore")] <- "logs"
  if(any(type %in% c("logpdf"))) type[type %in% c("logpdf")] <- "loglikelihood"
  if(any(type %in% c("rps"))) type[type %in% c("rps")] <- "crps"
  if(any(type %in% c("dawidsebastiani"))) type[type %in% c("dawidsebastiani")] <- "dss"
  if(any(dup <- duplicated(type))) {
    otype <- otype[!dup]
    type <- type[!dup]
  }

  ## FIXME: how to handle 'size' in binomial family?

  ## extract probability distribution object
  pd <- if(is.null(newdata)) {
    prodist(object)
  } else {
    prodist(object, newdata = newdata, na.action = na.action)
  }
  
  ## extract newresponse
  y <- newresponse(object, newdata = newdata, na.action = na.action)
  y[["(weights)"]] <- NULL
  if (ncol(y) > 1L) stop("multivariate responses not supported yet")
  y <- y[[1L]]
  
  ## evaluate type of proscore
  ps <- list()
  if("logs" %in% type) ps$logs <- -log_pdf(pd, y, drop = TRUE)
  if("loglikelihood" %in% type) ps$loglikelihood <- log_pdf(pd, y, drop = TRUE)
  if("crps" %in% type) {
    crps_fun <- if (requireNamespace("scoringRules")) scoringRules::crps else crps.distribution
    ps$crps <- drop(crps_fun(pd, y))
  }
  if("mae" %in% type) ps$mae <- drop(abs(mean(pd) - y))
  if("mse" %in% type) ps$mse <- drop((mean(pd) - y)^2)
  if("dss" %in% type) ps$dss <- drop((mean(pd) - y)^2/variance(pd) + log(variance(pd)))

  ## FIXME: negative inteval-censored log-likelihood (NILL)  
  ## delta <- 1
  ## if("nill" %in% type) ps$nill <- log(delta) - log(cdf(pd, y + delta/2) - cdf(pd, y - delta/2))

  ps <- ps[type]

  ## aggregate if desired
  if(isTRUE(aggregate)) aggregate <- mean
  if(is.function(aggregate)) ps <- lapply(ps, aggregate, ...)

  ## convert to data frame if drop = FALSE
  ps <- do.call("cbind", ps)
  colnames(ps) <- otype
  if(drop) {
    if(!is.null(dim(ps)) && any(dim(ps) == 1L)) ps <- drop(ps)
  } else {
    if(!inherits(ps, "data.frame")) ps <- as.data.frame(ps)
  }
  
  return(ps)
}
