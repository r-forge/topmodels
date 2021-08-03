# -------------------------------------------------------------------
# Programming outline: `procast()`
# -------------------------------------------------------------------
# - For scores we need at = data.frame(y = response, x = model.matrix)
# -------------------------------------------------------------------


#' Procast: Probabilistic Forecasting
#' 
#' Generic function and methods for computing various kinds of forecasts from
#' probabilistic (regression) models (probabilistic 4-casting).
#' 
#' The function \code{procast} provides a unified framework for probabilistic
#' 4-casting based on probabilistic (regression) models, also known as
#' distributional regression approaches. Typical types of predictions include
#' quantiles, probabilities, (conditional) expectations, variances,
#' (log-)densities, or scores. Internally, \code{procast} methods typically
#' compute the predicted parameters for each observation and then transform
#' these to the desired outcome. Some quantities (e.g., expectations or
#' variances) can be computed directly from the predicted parameters of the
#' distribution while others require an additional argument \code{at} which the
#' distribution is evaluated (e.g., the probability of a quantile or an
#' observation of the response. The argument \code{at} can also be the
#' character \code{"function"} or \code{"list"} so that a single function or
#' list of functions is set up that can be evaluated \code{at} different values
#' later on.
#' 
#' The function \code{procast_setup} is a convenience wrapper that make setting
#' up \code{procast} methods easier for package developers. It takes a data
#' frame of predicted parameters \code{pars} and a function \code{FUN} which is
#' to be evaluated at the parameters. This can either have the interface
#' \code{FUN(pars, \dots)} when the desired quantity can be predicted
#' directly from the predicted parameters -- or the interface \code{FUN(at,
#' pars, \dots)} if an additional argument \code{at} is needed.
#' \code{procast_setup} takes care of suitable expanding \code{at} to the
#' dimensions of \code{pars} and optionally setting up a (list of) function(s)
#' to be returned.
#' 
#' @aliases procast procast.lm procast_setup
#' @param object a fitted model object. For the \code{default} method this
#' needs to be \code{formula}-based so that \code{\link[stats]{model.frame}}
#' can be used to extract the response from the original data the model was
#' fitted to or \code{\link[stats]{terms}} can be used to set up the response
#' on \code{newdata}.
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
#' @return Either a \code{data.frame} of predictions (in case of multivariate
#' forecasts, or if \code{drop = FALSE}, default) or a vector (in case of a
#' univariate forecast and additionally \code{drop = TRUE}). Unless \code{at}
#' is the character string \code{"function"} or \code{"list"} in which case a
#' (list of) function(s) is returned.
#' @keywords regression
#' @examples
#' 
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
#' procast(m, newdata = nd, at = c(0.25, 0.5, 0.75))
#' 
#' ## all combinations of quantiles and observations
#' procast(m, newdata = nd, at = rbind(c(0.25, 0.5, 0.75)))
#' 
#' ## function for computing quantiles (vectorized)
#' qnt1 <- procast(m, newdata = nd, at = "function")
#' ## as before
#' qnt1(0.5)
#' qnt1(c(0.25, 0.5, 0.75))
#' qnt1(rbind(c(0.25, 0.5, 0.75)))
#' 
#' ## list of functions
#' qnt2 <- procast(m, newdata = nd, at = "list")
#' qnt2[[1]]
#' qnt2[[1]](0.5)
#' qnt2[[1]](c(0.25, 0.5, 0.75))
#' 
#' 
#' @export 
procast <- function(object, newdata = NULL, na.action = na.pass, type = "quantile", at = 0.5, ...) {
  UseMethod("procast")
}

#' @rdname procast
#' @export 
procast_setup <- function(pars,
                          FUN,
                          at = NULL,
                          drop = FALSE,
                          type = "procast",
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
        if (!is.data.frame(at)) {
          if (length(at) == 1L) at <- rep.int(as.vector(at), n)
          if (length(at) != n) at <- rbind(at)
        }
        if (is.matrix(at) && NROW(at) == 1L) {
          at <- matrix(rep(at, each = n), nrow = n)
          rv <- FUN(as.vector(at), pars = pars[rep(1L:n, ncol(at)), , drop = FALSE], ...)
          rv <- matrix(rv, nrow = n)
          rownames(rv) <- rownames(pars)
          colnames(rv) <- paste(substr(type, 1L, 1L),
            round(at[1L, ], digits = pmax(3L, getOption("digits") - 3L)),
            sep = "_"
          )
        } else {
          rv <- FUN(at, pars = pars, ...)
          names(rv) <- rownames(pars)
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
#' @method procast lm
#' @export
procast.lm <- function(object,
                       newdata = NULL,
                       na.action = na.pass,
                       type = c(
                         "quantile", "location", "scale", "parameter",
                         "density", "probability", "score"
                       ),
                       at = 0.5,
                       drop = FALSE,
                       ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES AND NECESSARY ARGUMENTS
  # -------------------------------------------------------------------
  ## Predicted means
  pars <- if (missing(newdata) || is.null(newdata)) {
    object$fitted.values
  } else {
    predict(object, newdata = newdata, na.action = na.action)
  }
  pars <- data.frame(mu = pars)

  ## Add maximum likelihood estimator of constant varians
  pars$sigma <- summary(object)$sigma * sqrt(df.residual(object) / nobs(object))

  ## Types of predictions
  type <- match.arg(type)

  ## Set up function that computes prediction from model parameters
  FUN <- switch(type,
    "quantile" = function(at, pars, ...) qnorm(at, mean = pars$mu, sd = pars$sigma, ...),
    "location" = function(pars) pars$mu,
    "scale" = function(pars) pars$sigma,
    "parameter" = function(pars) pars,
    "density" = function(at, pars, ...) dnorm(at, mean = pars$mu, sd = pars$sigma, ...),
    "probability" = function(at, pars, ...) pnorm(at, mean = pars$mu, sd = pars$sigma, ...),
    "score" = function(at, pars, ...) (at$y - pars$mu)^2 * at$x
  )

  # -------------------------------------------------------------------
  # CALL WOKRHORSE AND RETURN
  # -------------------------------------------------------------------
  procast_setup(pars, FUN = FUN, at = at, drop = drop, type = type, ...)
}


#' @rdname procast
#' @method procast crch
#' @export
procast.crch <- function(object,
                         newdata = NULL,
                         na.action = na.pass,
                         type = c(
                           "quantile", "location", "scale", "parameter",
                           "density", "probability", "score"
                         ),
                         at = 0.5,
                         drop = FALSE,
                         ...) { # FIXME: (ML) Additional parameters currently not used
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES AND NECESSARY ARGUMENTS
  # -------------------------------------------------------------------
  ## Call
  cl <- match.call() # FIXME: (ML) Change to normal function call, do not use `eval(cl, parent.frame())`
  cl[[1]] <- quote(predict)
  cl$drop <- NULL

  ## Check possible types of predictions
  type <- match.arg(type)

  ## FIXME: (ML) Implement score (estfun)
  if (type == "score") {
    stop("`type=score` not supported yet.")
  }

  ## * Is 'at' some kind of 'data'
  ## * or the special type 'list' or 'function'
  ## * or missing altogether ('none')
  attype <- if (is.null(at) || type %in% c("location", "scale", "parameter")) {
    "none"
  } else if (is.character(at)) {
    match.arg(at, c("function", "list"))
  } else {
    "data"
  }

  ## Set up internal function that computes prediction
  cl$type <- switch(type,
    "quantile" = "quantile",
    "location" = "location",
    "scale" = "scale",
    "parameter" = "parameter",
    "density" = "density",
    "probability" = "probability"
  )

  # -------------------------------------------------------------------
  # PREPARE OUTPUT CONDITIONAL ON `attype`
  # -------------------------------------------------------------------
  if (attype == "function") {
    rval <- eval(cl, parent.frame())
  } else if (attype == "list") {
    ## FIXME: (ML) Argument `type' == `list' is not yet fully supported in predict.crch()
    stop("Argument `type' == `list' is not yet supported for crch model classes.")
  } else {
    rval <- eval(cl, parent.frame())

    ## Default: return a data.frame (drop=FALSE) but optionally this can
    ## be dropped to a vector if possible (drop=TRUE)
    if (drop) {
      if (!is.null(dim(rval)) && NCOL(rval) == 1L) {
        # NOTE: (ML) How can condition be fulfilled? Compare code coverage.
        rval <- drop(rval)
      }
    } else {
      if (is.null(dim(rval))) {
        rval <- as.matrix(rval)
        if (ncol(rval) == 1L) colnames(rval) <- type
      }
      if (!inherits(rval, "data.frame")) rval <- as.data.frame(rval)
    }
  }

  # -------------------------------------------------------------------
  # RETURN
  # -------------------------------------------------------------------
  return(rval)
}

#' @rdname procast
#' @method procast disttree
#' @param use_distfamily For intern use only, will not be supported in the future.
#' @export
procast.disttree <- function(object,
                             newdata = NULL,
                             na.action = na.pass, # FIXME: (ML) Currently not supported
                             type = c(
                               "quantile", "location", "scale", "parameter",
                               "density", "probability", "score"
                             ),
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

  ## FIXME: (ML) Implement score (estfun)
  if (type == "score") {
    stop("`type=score` not supported yet.")
  }

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
    ## Here ugly workaround, which must be improved (probabily straight in `disstree`).
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
#' @method procast glm
#' @export
procast.glm <- function(object,
                        newdata = NULL,
                        na.action = na.pass,
                        type = c(
                          "quantile", "location", "scale", "parameter",
                          "density", "probability", "score"
                        ),
                        at = 0.5,
                        drop = FALSE,
                        ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES AND NECESSARY ARGUMENTS
  # -------------------------------------------------------------------
  ## Get family
  family <- substr(family(object)$family, 1L, 17L)

  ## Get weight, nobs, n if available
  weights <- if (is.null(attr(at, "weights"))) weights(object) else attr(at, "weights")
  nobs <- if (is.null(attr(at, "nobs"))) nobs(object) else attr(at, "nobs")
  n <- if (is.null(attr(at, "n"))) NULL else attr(at, "n")

  if (is.null(n)) {
    ## FIXME: (ML) This is probably not correct, as we need initializing for `at`
    ## NOTE: (ML) weights and nobs property of the (new) response?
    ## NOTE: (ML) n property of the model?
    y <- newresponse(object, newdata = newdata)

    weights <- attr(y, "weights")
    nobs <- attr(y, "nobs")
    n <- attr(y, "n")
  }

  ## Predicted means
  ## FIXME: (ML) Is it really always only one predicted parameter
  pars <- if (missing(newdata) || is.null(newdata)) {
    fitted(object)
  } else {
    predict(object, newdata = newdata, na.action = na.action)
  }
  pars <- data.frame(mu = pars)

  ## Types of predictions
  type <- match.arg(type)

  ## FIXME: (ML) Implement score (estfun)
  if (family == "gaussian") {

    ## Add maximum likelihood estimator of constant varians
    ## FIXME: (ML) do we need second part (compare countreg)
    pars$sigma <- sqrt(summary(object)$dispersion) * sqrt(df.residual(object) / nobs(object))

    FUN <- switch(type,
      "quantile" = function(at, pars, ...) qnorm(at, mean = pars$mu, sd = pars$sigma, ...),
      "location" = function(pars) pars$mu,
      "scale" = function(pars) pars$sigma,
      "parameter" = function(pars) pars,
      "density" = function(at, pars, ...) dnorm(at, mean = pars$mu, sd = pars$sigma, ...),
      "probability" = function(at, pars, ...) pnorm(at, mean = pars$mu, sd = pars$sigma, ...),
      "score" = function(at, pars, ...) (at$y - pars$mu)^2 * at$x
    )
  } else if (family == "poisson") {
    FUN <- switch(type,
      "quantile" = function(at, pars, ...) qpois(at, lambda = pars$mu, ...),
      "location" = stop("not yet implemented"),
      "scale" = stop("not yet implemented"),
      "parameter" = stop("not yet implemented"),
      "density" = function(at, pars, ...) dpois(at, lambda = pars$mu, ...),
      "probability" = function(at, pars, ...) ppois(at, lambda = pars$mu, ...),
      "score" =  stop("not yet implemented")
    )
  } else if (family == "Negative Binomial") {
    ## FIXME: (ML) check implementation and compare w/ countreg -> different number of bins
    pars$theta <- object$theta
    if(is.null(pars$theta)) pars$theta <- get(".Theta", environment(family(object)$variance))

    FUN <- switch(type,
      "quantile" = function(at, pars, ...) qnbinom(at, mu = pars$mu, size = pars$theta, ...),
      "location" = stop("not yet implemented"),
      "scale" = stop("not yet implemented"),
      "parameter" = stop("not yet implemented"),
      "density" = function(at, pars, ...) dnbinom(at, mu = pars$mu,, size = pars$theta, ...),
      "probability" = function(at, pars, ...) pnbinom(at, mu = pars$mu,, size = pars$theta, ...),
      "score" =  stop("not yet implemented")
    )
  } else if (family == "binomial") {
    FUN <- switch(type,
      "quantile" = function(at, pars, ...) qbinom(at, size = n, prob = pars$mu, ...),
      "location" = stop("not yet implemented"),
      "scale" = stop("not yet implemented"),
      "parameter" = stop("not yet implemented"),
      "density" = function(at, pars, ...) {
        ## FIXME: (ML) copied from countreg, is this needed (1st part why, 2nd part handled in `dbinom`)
        # at <- y * weights / n  # TODO: (ML) Here `at <- at * weights / n` ?
        # if (!isTRUE(all.equal(as.numeric(at), as.numeric(round(at))))) {
        #  stop("binomial quantile residuals require integer response")
        # }
        # at <- round(at)
        dbinom(at, size = n, prob = pars$mu, ...)
      },
      "probability" = function(at, pars, ...) {
        ## FIXME: (ML) copied from countreg, is this needed (1st part why, 2nd part no error in `pbinom`)
        # at <- y * weights / n  # TODO: (ML) Here `at <- at * weights / n` ?
        # if (!isTRUE(all.equal(as.numeric(at), as.numeric(round(at))))) {
        #  stop("binomial quantile residuals require integer response")
        # }
        # at <- round(at)
        pbinom(at, size = n, prob = pars$mu, ...)
      },
      "score" = stop("not yet implemented"),
    )
  } else {
    stop(sprintf("family %s not implemented yet", family))
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
