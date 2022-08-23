#' Scoring Probabilistic Forecasts
#' 
#' Generic function and default method for computing various kinds of scores
#' for fitted or predicted probability distributions from (regression) models.
#' 
#' The function \code{proscore} provides a unified framework for scoring
#' probabilistic forecasts (in-sample or out-of-sample). The following scores
#' are currently available, using the following notation: \eqn{Y} is the predicted
#' random variable with cumulative distribution function \eqn{F(y)} and probability
#' density function \eqn{f(y)}. The actual observation is \eqn{y}.
#'
#' **Log-likelihood:** Also known as log-score, logarithmic score, or log-density.
#'
#' \deqn{
#'   \log f_Y(y)
#' }{
#'   log(f(y))
#' }
#'
#' **Continuous ranked probability score (CRPS):**
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
#' **Mean absolute error (MAE):**
#'
#' \deqn{
#'   \left( E(Y) - y \right)^2
#' }{
#'   (E(Y) - y)^2
#' }
#'
#' where \eqn{E(Y)} is the expectation of the predicted random variable \eqn{Y}.
#'
#' **Mean squared error (MSE):**
#'
#' \deqn{
#'   \left( E(Y) - y \right)^2
#' }{
#'   (E(Y) - y)^2
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
#' \code{"loglikelihood"} (or equivalently \code{"logs"} or \code{"log_pdf"}),
#' \code{"CRPS"} (or equivalently \code{"RPS"}), \code{"MAE"}, and \code{"MSE"}.
#' Upper or lower case spellings can be used interchangably.
#' @param drop logical. Should scores be returned in a data frame (default)
#' or (if possible) dropped to a vector.
#' @param aggregate logical or function to be used for aggregating scores across
#' observations. Setting \code{aggregate = TRUE} corresponds to using \code{mean}.
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
#' ## in-sample log-likelihood
#' proscore(m, type = "loglik", aggregate = sum)
#' logLik(m)
#'
#' ## compute mean of all available scores
#' proscore(m, type = c("LogS", "CRPS", "MAE", "MSE"), aggregate = TRUE)
#' ## note that "LogS" and "loglik" above are both matched to using
#' ## the log-likelihood but the user-supplied spelling is preserved
#'
#' ## prediction using a new data set (final of the tournament)
#' final <- tail(FIFA2018, 2)
#' proscore(m, newdata = final, type = c("LogLik", "CRPS", "MAE", "MSE"))
#' @export 
proscore <- function(object, newdata = NULL, na.action = na.pass, type = "loglikelihood", drop = FALSE, aggregate = FALSE, ...) {
  UseMethod("proscore")
}

#' @rdname proscore
#' @importFrom distributions3 prodist log_pdf
#' @export 
proscore.default <- function(object, newdata = NULL, na.action = na.pass, type = "loglikelihood", drop = FALSE, aggregate = FALSE, ...)
{
  ## match type
  otype <- type
  type <- sapply(tolower(type), match.arg, c(
    "loglikelihood", "logs", "log_pdf",
    "crps", "rps",
    "mae",
    "mse"))
  if(any(type %in% c("logs", "log_pdf"))) type[type %in% c("logs", "log_pdf")] <- "loglikelihood"
  if(any(type %in% c("rps"))) type[type %in% c("rps")] <- "crps"
  if(any(dup <- duplicated(type))) {
    otype <- otype[!dup]
    type <- type[!dup]
  }

  if("crps" %in% type) stopifnot(requireNamespace("scoringRules"))

  ## FIXME: how to handle 'size' in binomial family?
  ## extract probability distribution object
  pd <- if(is.null(newdata)) {
    distributions3::prodist(object)
  } else {
    distributions3::prodist(object, newdata = newdata, na.action = na.action)
  }
  
  ## extract newresponse
  y <- newresponse(object, newdata = newdata, na.action = na.action)
  
  ## evaluate type of proscore
  ps <- list()
  if("loglikelihood" %in% type) ps$loglikelihood <- distributions3::log_pdf(pd, y, drop = TRUE)
  if("crps" %in% type) ps$crps <- scoringRules::crps(pd, y, drop = TRUE)
  if("mae" %in% type) ps$mae <- abs(mean(pd, drop = TRUE) - y)
  if("mse" %in% type) ps$mse <- (mean(pd, drop = TRUE) - y)^2
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
