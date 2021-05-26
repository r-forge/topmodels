# --------------------------------------------------------------------
# PROBABILISTIC FORECASTING METHOD: `procast()`
# --------------------------------------------------------------------
# TODO: (ML) Do we need/want support weights -> `newresponse()`?
# TODO: (ML) Which `type`s should be supported?
# FIXME: (ML) Support types mean/expectation/response and variance/dispersion instead of
#   location and scale.
# NOTE: (Z) For scores we need at = data.frame(y = response, x = model.matrix).
# TODO: (ML) In case `at = function`, it always returns a named vector. Do we want that?

# --------------------------------------------------------------------
# Set up generic function
# --------------------------------------------------------------------
procast <- function(object, newdata = NULL, na.action = na.pass, type = "quantile", at = 0.5, ...) {
  UseMethod("procast")
}


# --------------------------------------------------------------------
# Helper function for nice formatting
# --------------------------------------------------------------------
procast_setup <- function(pars, FUN, at = NULL, drop = FALSE, type = "procast", ...) {

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
            round(at[1L, ], digits = pmax(3L, getOption("digits") - 3L)), sep = "_")
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
          .Names = names(pars))), ff[3L]))
        return(rval)
      } else {
        rval <- FUN2a(at, pars = pars, ...)
      }
    }
  }

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


# --------------------------------------------------------------------
# `lm` class
# --------------------------------------------------------------------
procast.lm <- function(object, 
                       newdata = NULL, 
                       na.action = na.pass,
                       type = c("quantile", "location", "scale", "parameter", 
                         "density", "probability", "score"),
                       at = 0.5,
                       drop = FALSE, 
                       ...) {
  
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

  procast_setup(pars, FUN = FUN, at = at, drop = drop, type = type, ...)
}


# --------------------------------------------------------------------
# `crch` class
# --------------------------------------------------------------------
procast.crch <- function(object,
                         newdata = NULL,
                         na.action = na.pass,
                         type = c("quantile", "location", "scale", "parameter",
                           "density", "probability", "score"),
                         at = 0.5,
                         drop = FALSE,
                         ...) { # FIXME: (ML) Additional parameters currently not used

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

  return(rval)
}


# --------------------------------------------------------------------
# `disttree` and `distforest` class
# --------------------------------------------------------------------
procast.disttree <- procast.distforest <- function(object,
                                                   newdata = NULL,
                                                   na.action = na.pass, # FIXME: (ML) Currently not supported
                                                   type = c("quantile", "location", "scale", "parameter",
                                                     "density", "probability", "score"),
                                                   at = 0.5,
                                                   drop = FALSE,
                                                   use_internals = TRUE, # FIXME: (ML) Just for development
                                                   ...) { # FIXME: (ML) Additional parameters currently not always supported

  ## First if working for `NO()`
  ## FIXME: (ML) Extend for all families
  if (object$info$family$family.name != "Normal Distribution") {
    warning("So far not tested for specified family")
  }

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
  if (use_internals == FALSE) {
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
      "density" = function(at, pars) {
        sapply(
          1:NROW(pars),
          function(i) {
            family$ddist(
              y = at[i], eta = as.numeric(family$linkfun(pars[i, ])),
              log = FALSE, weights = NULL, sum = FALSE
            )
          }
        )
      },
      "probability" = function(at, pars) {
        sapply(
          1:NROW(pars),
          function(i) {
            family$pdist(
              q = at[i], eta = as.numeric(family$linkfun(pars[i, ])),
              lower.tail = TRUE, log.p = FALSE
            )
          }
        )
      }
    )
  }

  procast_setup(pars, FUN = FUN, at = at, drop = drop, type = type, ...)
}


# --------------------------------------------------------------------
# `glm` class
# --------------------------------------------------------------------
procast.glm <- function(object,
                        newdata = NULL,
                        na.action = na.pass,
                        type = c("quantile", "location", "scale", "parameter",
                          "density", "probability", "score"),
                        at = 0.5,
                        drop = FALSE,
                        ...) {

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
  if (object$family$family == "gaussian") {

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
  } else if (object$family$family == "poisson") {
    FUN <- switch(type,
      "quantile" = function(at, pars, ...) qpois(at, lambda = pars$mu, ...),
      "location" = stop("not yet implemented"),
      "scale" = stop("not yet implemented"),
      "parameter" = stop("not yet implemented"),
      "density" = function(at, pars, ...) dpois(at, lambda = pars$mu, ...),
      "probability" = function(at, pars, ...) ppois(at, lambda = pars$mu, ...),
      "score" =  stop("not yet implemented")
    )
  } else if (object$family$family == "binomial") {

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
    stop(sprintf("family %s not implemented yet", object$family$family))
  }

  ## set up function that computes prediction from model parameters
  procast_setup(pars, FUN = FUN, at = at, drop = drop, type = type, ...)
}
