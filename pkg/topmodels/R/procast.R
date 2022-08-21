#' Procast: Probabilistic Forecasting
#' 
#' Generic function and methods for computing various kinds of probabilistic
#' forecasts from (regression) models.
#' 
#' The function \code{procast} provides a unified framework for probabilistic
#' forcasting (or procasting, for short) based on probabilistic (regression)
#' models, also known as distributional regression approaches. Typical types
#' of predictions include quantiles, probabilities, (conditional) expectations,
#' variances, and (log-)densities. Internally, \code{procast} methods typically
#' compute the predicted parameters for each observation and then compute the
#' desired outcome for the distributions with the respective parameters.
#'
#' Some quantities, e.g., the moments of the distribution (like mean or variance),
#' can be computed directly from the predicted parameters of the
#' distribution while others require an additional argument \code{at} which the
#' distribution is evaluated (e.g., the probability of a quantile or an
#' observation of the response.
#' 
#' The default \code{procast} method leverages the S3 classes and methods for
#' probability distributions from the \pkg{distributions3} package. It proceeds
#' in two steps: First, \code{\link[distributions3]{prodist}} is used to obtain the
#' predicted probability distribution object. Second, the extractor methods such
#' as \code{quantile}, \code{\link[distributions3]{cdf}}, etc. are used to
#' compute quantiles, probabilities, etc. from the distribution objects.
#'
#' Therefore, to enable \code{procast} for a certain type of model object, the
#' recommended approach is to implement a \code{prodist} method which can then
#' be leveraged. However, if the \pkg{distributions3} package does not support
#' the necessary probability distribution, then it may also be necessary to
#' implement a new distribution objects, see \code{\link[distributions3]{apply_dpqr}}.
#'
#' Before adopting the \pkg{distributions3} framework as the recommended workflow
#' for procasting, the package had taken a different approach which is described
#' in the following. Note, however, that this will be discontinued when we have
#' converted all procasting methods to the new workflow.
#'
#' Old workflow: The function \code{procast_setup} is a convenience wrapper that makes setting
#' up \code{procast} methods easier for package developers. It takes a data
#' frame of predicted parameters \code{pars} and a function \code{FUN} which is
#' to be evaluated at the parameters. This can either have the interface
#' \code{FUN(pars, \dots)} when the desired quantity can be predicted
#' directly from the predicted parameters -- or the interface \code{FUN(at,
#' pars, \dots)} if an additional argument \code{at} is needed.
#' \code{procast_setup} takes care of suitable expanding \code{at} to the
#' dimensions of \code{pars}.
#' 
#' @aliases procast procast_setup
#' @param object a fitted model object. For the \code{default} method this
#' needs to have a \code{\link[distributions3]{prodist}} method.
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param na.action function determining what should be done with missing
#' values in \code{newdata}.  The default is to employ \code{NA}.
#' @param type character specifying the type of probabilistic forecast to
#' compute.  In \code{procast_setup} the \code{type} is only used for nice
#' labels of the returned data frame. 
#' @param at specification of values at which the forecasts should be
#' evaluated, typically a numeric vector but possibly also a matrix or data
#' frame.  Additionally, \code{at} can be the character string
#' \code{"function"} or \code{"list"}, see details below.
#' @param drop logical. Should forecasts be returned in a data frame (default)
#' or (if possible) dropped to a vector, see details below.
#' @param FUN function to be used for forecasts. Either of type \code{FUN(pars,
#' \dots)} or \code{FUN(at, pars, \dots)}, see details below.
#' @param pars a data frame of predicted distribution parameters.
#' @param \dots further parameters passed to methods.
#' @param elementwise logical. Should each element of distribution only be evaluated at the
#' corresponding element of \code{at} (\code{elementwise = TRUE}) or at all elements
#' in \code{at} (\code{elementwise = FALSE}). Elementwise evaluation is only possible
#' if the number of observations is length of \code{at} are the same and in that case a vector of
#' the same length is returned. Otherwise a matrix is returned. The default is to use
#' \code{elementwise = TRUE} if possible, and otherwise \code{elementwise = FALSE}.
#' @param drop logical. Should the result be simplified to a vector if possible (by
#' dropping the dimension attribute)? If \code{FALSE} a matrix is always returned.
#' @return Either a \code{data.frame} of predictions (in case of multivariate
#' forecasts, or if \code{drop = FALSE}, default) or a vector (in case of a
#' univariate forecast and additionally \code{drop = TRUE}). Unless \code{at}
#' is the character string \code{"function"} or \code{"list"} in which case a
#' (list of) function(s) is returned.
#' @keywords regression
#' @examples
#' ## linear regression models (homoscedastic Gaussian response)
#' m <- lm(dist ~ speed, data = cars)
#' 
#' ## medians on observed data
#' procast(m)
#' procast(m, drop = TRUE)
#' 
#' ## probability integral transform (PIT) on observed data
#' procast(m, type = "probability", at = cars$dist)
#' 
#' ## log-likelihood contributions
#' procast(m, type = "density", at = cars$dist, log = TRUE)
#' 
#' ## log-likelihood sum
#' sum(procast(m, type = "density", at = cars$dist, log = TRUE))
#' logLik(m)
#' 
#' 
#' ## medians on new data
#' nd <- data.frame(speed = c(10, 15, 20))
#' procast(m, newdata = nd)
#' 
#' ## different quantile for each observation
#' procast(m, newdata = nd, at = c(0.25, 0.5, 0.75), elementwise = TRUE)
#' 
#' ## all combinations of quantiles and observations
#' procast(m, newdata = nd, at = c(0.25, 0.5, 0.75), elementwise = FALSE)
#' 
#' @export 
procast <- function(object, newdata = NULL, na.action = na.pass, type = "quantile", at = 0.5, drop = FALSE, ...) {
  UseMethod("procast")
}

#' @rdname procast
#' @importFrom distributions3 prodist cdf pdf log_pdf variance skewness kurtosis
#' @export 
procast.default <- function(object, newdata = NULL, na.action = na.pass,
  type = c("quantile", "mean", "variance", "probability", "density", "loglikelihood", "distribution", "parameters", "kurtosis", "skewness"),
  at = 0.5, drop = FALSE, ...)
{
  ## match type
  type <- match.arg(type[1L], c(
    "quantile", "mean", "variance",
    "probability", "cdf",
    "density", "pdf", "pmf",
    "loglikelihood", "log_pdf",
    "distribution", "parameters",
    "kurtosis", "skewness"))
  if(type == "cdf") type <- "probability"
  if(type %in% c("pdf", "pmf")) type <- "density"
  if(type == "log_pdf") type <- "loglikelihood"

  ## FIXME: how to handle 'size' in binomial family?
  ## extract probability distribution object
  pd <- if(is.null(newdata)) {
    distributions3::prodist(object)
  } else {
    distributions3::prodist(object, newdata = newdata, na.action = na.action)
  }
  
  ## evaluate type of procast
  pc <- switch(type,
    "quantile"      = quantile(pd, at, ...),
    "mean"          = mean(pd),
    "variance"      = distributions3::variance(pd),
    "probability"   = distributions3::cdf(pd, at, ...),
    "density"       = distributions3::pdf(pd, at, ...),
    "loglikelihood" = distributions3::log_pdf(pd, at, ...),
    "distribution"  = pd,
    "parameters"    = as.matrix(pd),
    "skewness"      = distributions3::skewness(pd, at, ...),
    "kurtosis"      = distributions3::kurtosis(pd, at, ...)
  )
  
  ## convert to data frame if drop = FALSE
  if(drop) {
    if(!is.null(dim(pc)) && NCOL(pc) == 1L) pc <- drop(pc)
  } else {
    if(inherits(pc, "distribution")) {
      pc <- as.data.frame(pc)
      colnames(pc) <- type
    }
    if(is.null(dim(pc))) {
      pc <- as.matrix(pc)
      if(ncol(pc) == 1L) colnames(pc) <- type
    }
    if(!inherits(pc, "data.frame")) pc <- as.data.frame(pc)
  }
  
  return(pc)
}

#' @rdname procast
#' @export 
procast_setup <- function(pars,
                          FUN,
                          at = NULL,
                          drop = FALSE,
                          type = "procast",
                          elementwise = NULL,
                          ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## * Is 'at' some kind of 'data'
  ## * or the special type 'list' or 'function'
  ## * or missing altogether ('none')
  attype <- if (is.null(at) || names(formals(FUN))[1L] == "pars") {
    "none"
  } else if (is.character(at)) {
    match.arg(at, c("function", "list"))
  } else {
    "data"
  }

  # -------------------------------------------------------------------
  # PREPARE OUTPUT CONDITIONAL ON `attype`
  # -------------------------------------------------------------------
  if (attype == "none") {

    ## If 'at' is missing: prediction is just a transformation of the parameters
    rval <- FUN(pars, ...)
    if (is.null(dim(rval))) names(rval) <- rownames(pars)
  } else {

    ## If 'at' is 'list':
    ## prediction is a list of functions with predicted parameters as default
    if (attype == "list") {
      FUN2 <- function(at, pars = pars, ...) NULL
      body(FUN2) <- call("FUN",
        at = as.name("at"),
        pars = as.call(structure(sapply(c("data.frame", names(pars)), as.name), .Names = c("", names(pars)))),
        as.name("...")
      )
      ff <- formals(FUN2)

      rval <- lapply(1L:nrow(pars), function(i) {
        FUN2i <- FUN2
        formals(FUN2i) <- as.pairlist(c(ff[1L], as.pairlist(as.list(pars[i, ])), ff[3L]))
        FUN2i
      })
      return(rval)

      ## Otherwise ('at' is 'data' or 'function'):
      ## set up a function that suitably expands 'at' (if necessary)
      ## and then either evaluates it at the predicted parameters ('data')
      ## or sets up a user interface with predicted parameters as default ('function')
      ## FIXME: (ML) added drop = FALSE below, to work for glm w/ family = poisson. check!
    } else {
      FUN2a <- function(at, pars, ...) {
        n <- NROW(pars)
        if (is.null(dim(at))) {
          if (length(at) == 1L) at <- rep.int(as.vector(at), n)
          if (length(at) == n) {
            if(is.null(elementwise)) elementwise <- TRUE
          } else {
            if(isTRUE(elementwise)) stop(
              sprintf("length of 'at' differs from number of observations: %s != %s", length(at), n)
            )
            elementwise <- FALSE
          }
        } else {
          if(NROW(at) != n) stop(
            sprintf("number of rows of 'at' differs from number of observations: %s != %s", NROW(at), n)
          )
          at <- as.matrix(at)
        }
        if (elementwise && is.null(dim(at))) {
          rv <- FUN(at, pars = pars, ...)
          names(rv) <- rownames(pars)
        } else {
          if(is.null(dim(at))) {
            at <- matrix(rep(at, each = n), nrow = n)
            cnam <- paste(substr(type, 1L, 1L),
              round(at[1L, ], digits = pmax(3L, getOption("digits") - 3L)),
              sep = "_")
          } else {
            cnam <- paste(substr(type, 1L, 1L), 1L:NCOL(at), sep = "_")
          }
          rv <- FUN(as.vector(at), pars = pars[rep(1L:n, ncol(at)), , drop = FALSE], ...)
          rv <- matrix(rv, nrow = n)
          rownames(rv) <- rownames(pars)
          colnames(rv) <- cnam
        }
        return(rv)
      }

      if (attype == "function") {
        rval <- function(at, pars = pars, ...) NULL
        body(rval) <- call("FUN2a",
          at = as.name("at"),
          pars = as.call(structure(sapply(c("data.frame", names(pars)), as.name), .Names = c("", names(pars)))),
          as.name("...")
        )
        ff <- formals(FUN2a)
        formals(rval) <- as.pairlist(c(ff[1L], as.pairlist(structure(
          lapply(names(pars), function(n) call("$", as.name("pars"), n)),
          .Names = names(pars)
        )), ff[3L]))
        return(rval)
      } else {
        rval <- FUN2a(at, pars = pars, ...)
      }
    }
  }

  # -------------------------------------------------------------------
  # RETURN
  # -------------------------------------------------------------------
  ## Default: return a data.frame (drop=FALSE) but optionally this can
  ## be dropped to a vector if possible (drop=TRUE)
  if (drop) {
    if (!is.null(dim(rval)) && NCOL(rval) == 1L) {
      # TODO: (ML) How can this condition ever be fulfilled? Compare code coverage.
      rval <- drop(rval)
    }
  } else {
    if (is.null(dim(rval))) {
      rval <- as.matrix(rval)
      if (ncol(rval) == 1L) colnames(rval) <- type
    }
    if (!inherits(rval, "data.frame")) rval <- as.data.frame(rval)
  }

  return(rval)
}


#' @rdname procast
#' @method procast disttree
#' @param use_distfamily For intern use only, will not be supported in the future.
#' @export
procast.disttree <- function(object,
                             newdata = NULL,
                             na.action = na.pass, # FIXME: (ML) Currently not supported
                             type = c("quantile", "location", "scale", "parameter", "density", "probability"),
                             at = 0.5,
                             drop = FALSE,
                             use_distfamily = TRUE,
                             ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES AND NECESSARY ARGUMENTS
  # -------------------------------------------------------------------
  ## First if working for `NO()`
  ## FIXME: (ML) Extend for all families
  #if (object$info$family$family.name != "Normal Distribution") {
  #  warning("So far not tested for specified family")
  #}

  ## Get family
  family <- object$info$family

  ## Predicted means
  ## FIXME: (ML) Fix `na.action`, not supported by predict.distforest()
  # pars <- predict(object, newdata = newdata, type = "parameter", na.action = na.action)
  pars <- predict(object, newdata = newdata, type = "parameter")

  ## Types of predictions
  type <- match.arg(type)

  ## Set up function that computes prediction from model parameters
  if (!use_distfamily) {
    warning("For `use_distfamily = FALSE`, `qnorm()`, `dnorm()`, and `pnorm()` are used")
    FUN <- switch(type,
      "quantile" = function(at, pars, ...) qnorm(at, mean = pars$mu, sd = pars$sigma, ...),
      "location" = function(pars) pars$mu,
      "scale" = function(pars) pars$sigma,
      "parameter" = function(pars) pars,
      "density" = function(at, pars, ...) dnorm(at, mean = pars$mu, sd = pars$sigma, ...),
      "probability" = function(at, pars, ...) pnorm(at, mean = pars$mu, sd = pars$sigma, ...)
    )
  } else {
    ## FIXME: (ML) `disttree` unfortunately does not support a vector of parameters!
    ## Here ugly workaround, which must be improved (probabily straight in `disttree`).
    FUN <- switch(type,
      "quantile" = function(at, pars, ...) {
        object <- cbind(at = at, pars)
        rval <- sapply(
          1:NROW(object),
          function(idx) {
            family$qdist(
              p = as.numeric(object[idx, grepl("at", names(object))]),
              eta = as.numeric(family$linkfun(object[idx, !grepl("at", names(object))])),
              ...
            )
          }
        ) 
        rval <- if (NCOL(at) > 1) t(rval) else rval
        return(rval)
      },
      "quantile" = function(at, pars) {
        sapply(
          1:NROW(pars),
          function(i) {
            family$qdist(
              p = at[i], eta = as.numeric(family$linkfun(pars[i, ])),
              lower.tail = TRUE, log.p = FALSE
            )
          }
        )
      },
      "location" = function(pars) pars$mu,
      "scale" = function(pars) pars$sigma,
      "parameter" = function(pars) pars,
      "density" = function(at, pars, ...) {
        object <- cbind(at = at, pars)
        rval <- sapply(
          1:NROW(object),
          function(idx) {
            family$ddist(
              y = as.numeric(object[idx, grepl("at", names(object))]),
              eta = as.numeric(family$linkfun(object[idx, !grepl("at", names(object))])),
              ...
            )
          }
        ) 
        rval <- if (NCOL(at) > 1) t(rval) else rval
        return(rval)
      },
      "probability" = function(at, pars, ...) {
        object <- cbind(at = at, pars)
        rval <- sapply(
          1:NROW(object),
          function(idx) {
            family$pdist(
              q = as.numeric(object[idx, grepl("at", names(object))]),
              eta = as.numeric(family$linkfun(object[idx, !grepl("at", names(object))])),
              ...
            )
          }
        ) 
        rval <- if (NCOL(at) > 1) t(rval) else rval
        return(rval)
      }
    )
  }

  # -------------------------------------------------------------------
  # CALL WOKRHORSE AND RETURN
  # -------------------------------------------------------------------
  procast_setup(pars, FUN = FUN, at = at, drop = drop, type = type, ...)
}

#' @rdname procast
#' @method procast distforest
#' @export
procast.distforest <- procast.disttree
