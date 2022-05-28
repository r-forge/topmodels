#' Probability Integral Transform (PIT) Residuals
#' 
#' Generic function and methods for computing PIT residuals.
#' 
#' For regression models with a continuous response distribution, PIT residuals
#' are simply the predictive cumulative distribution (CDF) evaluated at the
#' observations (Dawid, 1984). For discrete distributions, a random value is drawn
#' from the range of probabilities corresponding to each observation similar to
#' the approach of Dunn and Smyth (1996) to gain (randomized) quantile residuals
#' (\code{\link{qresiduals}}). PIT values have been used under various names, but
#' to emphasize their similar properties to residuals we follow Warton (2017) and
#' refer to them as PIT residuals.
#' 
#' With the method \code{pitresiduals}, the PIT residuals can also be transformed
#' from the probability scale to another distribution scale. Supported quantile
#' scales are uniformly (\code{uniform}) and normally (\code{normal}) distributed.
#' Probabilites can be either a vector or a 2-column matrix of probabilities. The
#' latter offers to either draw \code{"random"} samples from the distribution or
#' compute corresponding \code{"quantile"}s such as the median etc. For
#' (randomized) quantile residuals (on the normal scale), as suggested by Dunn and Smyth (1996), 
#' the \code{scale} must be set to \code{"normal"} or 
#' \code{\link{qresiduals}} can be called directly.
#' 
#' @aliases pitresiduals pitresiduals.default
#' @param object an object. For the \code{default} method this needs to be
#' either a specification of probabilities (vector or 2-dimensional matrix of
#' probabilities) or an object from which the these can be obtained with
#' \code{\link{procast}}.
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param scale On which scale should the PIT residuals be shown; on the probability scale 
#' (\code{"uniform"}) or on the normal scale (\code{"normal"}).
#' @param type character specifying whether - in the case of discrete response
#' distributions - randomized quantile residuals or their corresponding
#' quantiles should be computed.
#' @param nsim numeric. The number of simulated randomized quantile residuals
#' per observation (for \code{type = "numeric"}).
#' @param delta numeric. The minimal difference to compute the range of
#' proabilities corresponding to each observation according to get (randomized)
#' quantile residuals.  For \code{NULL}, the minimal observed difference in the
#' resonse divided by \code{5e-6} is used.
#' @param prob numeric. The probabilities at which quantile residuals should be
#' computed (for \code{type = "quantile"}), defaulting to the median
#' (\code{0.5}).
#' @param \dots further parameters passed to methods.
#' @return A vector or matrix of PIT residuals.
#' @seealso \code{\link[stats]{qnorm}}, \code{\link{qqrplot}}
#' @references Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile
#' Residuals.} \emph{Journal of Computational and Graphical Statistics},
#' \bold{5}(3), 236--244. \doi{10.2307/1390802}
#' 
#' Dawid AP (1984). \dQuote{Present Position and Potential Developments: Some
#' Personal Views: Statistical Theory: The Prequential Approach.} \emph{Journal of
#' the Royal Statistical Society: Series A (General)}, \bold{147}(2), 278--92.
#' \doi{10.2307/2981683}.
#' 
#' Warton DI, Thibaut L, Wang YA (2017) \dQuote{The Pit-Trap--a `Model-Free'
#' Bootstrap Procedure for Inference About Regression Models with Discrete,
#' Multivariate Responses}. \emph{PLOS ONE}, \bold{12}(7), 1--18.
#' \doi{10.1371/journal.pone.0181790}.
#' 
#' @keywords regression
#' @examples
#' 
#' ## linear regression models (homoscedastic Gaussian response)
#' m <- lm(dist ~ speed, data = cars)
#' pitresiduals(m)
#' 
#' @export
pitresiduals <- function(object, ...) {
  UseMethod("pitresiduals")
}

#' @rdname pitresiduals
#' @method pitresiduals default
#' @export
pitresiduals.default <- function(object, 
                                 newdata = NULL, 
                                 scale   = c("uniform", "normal"),
                                 type    = c("random", "quantile"), 
                                 nsim    = 1L, 
                                 delta   = NULL,
                                 prob    = NULL,
                                 ...) {

  # Sanity checks
  stopifnot(is.numeric(nsim), length(nsim) == 1, nsim > sqrt(.Machine$double.eps))
  stopifnot(is.null(delta) || (is.numeric(delta) && length(delta) == 1 && delta > 0.0))

  ## match arguments
  scale <- match.arg(scale)
  type  <- match.arg(type)

  ## check and set prob to default
  if (type == "random" & !is.null(prob)) {
    warning("argument `prob` will be ignored for `type = 'random'`")
    prob <- NULL
  } else if (type == "quantile" & is.null(prob)) {
    prob <- 0.5
  } else if (type == "quantile" & !is.null(prob)) {
    stopifnot(is.numeric(prob) & is.vector(prob))
  } 

  if (!is.object(object) || all(class(object) == "data.frame") || is.numeric(object)) {
    # TODO: (ML) is this case already caught by `is.object()`?
    # NOTE: (RS) A vector with a class is also an object, but a vector
    #            with class is no longer a (plain) vector.
    #            x <- 1:4
    #            y <- 1:4; class(y) <- "foo"
    #            | object  | is.object() | is.vector() | is.numeric() |
    #            |   x     |    FALSE    |    TRUE     |    TRUE      |
    #            |   y     |    TRUE     |    FALSE    |    TRUE      |

    # TODO: (RS2ML) FYI; stop takes up a series of objects which can be combined to
    #       a character, thus stop(pate0("a\n", "b\n")) is identical to stop("a\n", "b\n").
    #       You basically save a paste0() call here.
    stop("`object` must be a (model) 00 object:\n" ,
      "  * you have supplied either a base object,\n",
      "  * an object of the single class `data.frame` or\n",
      "  * a numeric.")

  } else {

    ## FIXME: (ML) Implement correct NA handling: 
    ## * In newresponse() for y
    ## * Which is used for at w/i procast()
    y <- newresponse(object, newdata = newdata)
  
    ## to keep the attributes
    ## FIXME: (ML) in countreg its `at = cbind(y - 1L, y)`, why the difference?
    ## FIXME: (ML) Increased difference, otherwise did not work for binom and pois
    ##at <- cbind(y - .Machine$double.eps^0.4, y)
    ##at <- cbind(y - 1L, y)
    if (is.null(delta)) delta <- min(diff(sort(unique(y)))) / 5e6
    at <- cbind(y - delta, y)

    attr(at, "nobs")    <- attr(y, "nobs")
    attr(at, "n")       <- attr(y, "n")
    attr(at, "weights") <- attr(y, "weights")
    object <- procast(object, newdata = newdata, at = at, type = "probability")

    # TODO: (ML) There is no `try()` environment, which errors can be caught
    # TODO: (RS2ML) Not yet checked whats going on, we could use
    #       a tryCatch() here if needed (tbd).
    #       Someting along x <- tryCatch(<command>,
    #                                    warning = function(w) w,
    #                                    error   = function(e) stop("something went wrong"))
    if (inherits(object, "try-error")) {
      stop("could not obtain probability integral transform from 'object'")
    }
  }

  ## preprocess supplied probabilities
  nc <- NCOL(object)
  nr <- NROW(object)
  if (nc > 2L) stop("quantiles must either be 1- or 2-dimensional")
  if (nc == 2L) {
    if (type == "random") {
      object <- matrix(
        runif(nr * nsim, min = rep(object[, 1L], nsim), max = rep(object[, 2L], nsim)),
        nrow = nr, ncol = nsim, dimnames = list(rownames(object), paste("r", 1L:nsim, sep = "_"))
      )
    } else {
      ## FIXME: (Z) probably needs to be done on the transformed rather than the uniform scale...
      ## TODO: (ML) 
      ## * Otherwise akward features for heavily skewed distributions: observational vs. probability scale
      ## * Compare also "middle-point quantile residuals (MQRs) in Feng et al. (2020)
      nam <- rownames(object)

      ## FIXME: (ML) Alternative computation, which is test below.
      object2 <- sapply(prob, function(x) qunif(x, min = object[, 1L], max = object[, 2L]))

      object <- object[, 1L]  %*% t(1 - prob) + object[, 2L] %*% t(prob)

      stopifnot(all.equal(object, object2))  # FIXME: (ML) Test alternative computation

      dimnames(object) <- list(nam, paste("q", prob, sep = "_"))
    }
    nc <- NCOL(object)
  }
  if (!is.null(dim(object)) & nc == 1L) object <- drop(as.matrix(object))
  # FIXME: (ML) object can be a data.frame, so make sure drop works by converting to matrix

  ## compute quantile residuals  
  if (scale == "normal") object <- qnorm(object)
  return(object)
}
