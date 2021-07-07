#' (Randomized) Quantile Residuals
#' 
#' Generic function and methods for computing (randomized) quantile residuals.
#' 
#' (Randomized) quantile residuals have been suggested by Dunn and Smyth
#' (1996). For regression models with a continuous response distribution this
#' simply computes theoretical standard normal quantiles corresponding to the
#' probability integral transform of the fitted distribution. For discrete
#' distributions, a random theoretical normal quantile is drawn from the range
#' of probabilities corresponding to each observation. Hence, in \code{qqrplot}
#' the default is to use \code{trafo = qnorm} but other transformations can
#' also be used, specifically using the uniform probability scale (via
#' \code{trafo = NULL} or equivalently \code{qunif} or \code{identity}).
#' 
#' The default \code{qresiduals} method can compute randomized quantile
#' residuals from a vector (which essentially just calls
#' \code{\link[stats]{qnorm}}) or a 2-column matrix of probabilities. The
#' latter offers to either draw \code{"random"} samples from the distribution
#' or compute corresponding \code{"quantile"}s such as the median etc.
#' 
#' @aliases qresiduals qresiduals.default
#' @param object an object. For the \code{default} method this needs to be
#' either a specification of probabilities (vector or 2-dimensional matrix of
#' probabilities) or an object from which the these can be obtained with
#' \code{\link{procast}}.
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param trafo function for tranforming residuals from probability scale to a
#' different distribution scale (default: Gaussian).
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
#' @return A vector or matrix of quantile residuals.
#' @note Note that there is also a \code{\link[statmod]{qresiduals}} function
#' in the \pkg{statmod} package that is not generic and always returns a single
#' random quantile residual.
#' @seealso \code{\link[stats]{qnorm}}, \code{\link{qqrplot}}
#' @references Dunn KP, Smyth GK (1996). \dQuote{Randomized Quantile
#' Residuals.} \emph{Journal of Computational and Graphical Statistics},
#' \bold{5}, 1--10. \doi{10.2307/1390802}
#' @keywords regression
#' @examples
#' 
#' ## linear regression models (homoscedastic Gaussian response)
#' m <- lm(dist ~ speed, data = cars)
#' qresiduals(m)
#' 
#' @export
qresiduals <- function(object, ...) {
  UseMethod("qresiduals")
}

#' @rdname qresiduals
#' @method qresiduals default
#' @export
qresiduals.default <- function(object, 
                               newdata = NULL, 
                               trafo = qnorm, 
                               type = c("random", "quantile"), 
                               nsim = 1L, 
                               delta = NULL,
                               prob = NULL,
                               ...) {

  ## type of residual for discrete distribution (if any)
  type <- match.arg(type)

  ## check and set prob to default
  if (type == "random" & !is.null(prob)) {
    warning("argument `prob` will be ignored for `type = 'random'`")
    prob <- NULL
  } else if (type == "quantile" & is.null(prob)) {
    prob <- 0.5
  } else if (type == "quantile" & !is.null(prob)) {
    stopifnot(is.numeric(prob) & is.vector(prob))
  } 

  if(!is.object(object) || 
    all(class(object) == "data.frame") || 
    is.numeric(object)) { # TODO: (ML) is this case already caught by `is.object()`?

    stop(paste0("`object` must be a (model) 00 object:\n" ,
      "  * you have supplied either a base object,\n",
      "  * an object of the single class `data.frame` or\n",
      "  * a numeric."
    ))

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
    if(is.null(delta)) {
      delta <- min(diff(sort(unique(y))))
      delta <- delta / 5e6
    } else {
      stopifnot(is.numeric(delta) && length(delta) == 1)
    }

    at <- cbind(y - delta, y)

    attr(at, "nobs") <-     attr(y, "nobs")
    attr(at, "n") <-        attr(y, "n")
    attr(at, "weights") <-  attr(y, "weights")
    object <- procast(object, newdata = newdata, at = at, type = "probability")

    # TODO: (ML) There is no `try()` environment, which errors can be caught
    if (inherits(object, "try-error")) {
      stop("could not obtain probability integral transform from 'object'")
    }
  }

  ## preprocess supplied probabilities
  nc <- NCOL(object)
  nr <- NROW(object)
  if(nc > 2L) stop("quantiles must either be 1- or 2-dimensional")
  if(nc == 2L) {
    if(type == "random") {
      object <- matrix(
        runif(nr * nsim, min = rep(object[, 1L], nsim), max = rep(object[, 2L], nsim)),
        nrow = nr, ncol = nsim, dimnames = list(rownames(object), paste("r", 1L:nsim, sep = "_"))
      )
    } else {
      ## FIXME: (Z) probably needs to be done on the transformed rather than the uniform scale...
      ## TODO: (ML) Otherwise akward features for heavily skewed distributions: observational vs. probability scale
      nam <- rownames(object)
      
      ## FIXME: (ML) Alternative computation, which is test below.
      object2 <- sapply(prob, function(x) qunif(x, min = object[, 1L], max = object[, 2L]))

      object <- object[, 1L]  %*% t(1 - prob) + object[, 2L] %*% t(prob)

      stopifnot(all.equal(object, object2))  # FIXME: (ML) Test alternative computation

      dimnames(object) <- list(nam, paste("q", prob, sep = "_"))
    }
    nc <- NCOL(object)
  }
  if(!is.null(dim(object)) & nc == 1L) object <- drop(as.matrix(object))  
  # FIXME: (ML) object can be a data.frame, so make sure drop works by converting to matrix

  ## compute quantile residuals  
  if(!is.null(trafo)) object <- trafo(object)  # TODO: (ML) Why on the normal scale? Common behaviour compared to traditional diagnostics.
  return(object)
}
