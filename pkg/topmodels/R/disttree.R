## old procast_disttree/distforest
## still relying on procast_setup (included at the end)
## ignored for now
##
## TODO: write prodist.disttree instead
## - could either use the pdist, ddist, etc. from within the disttree
## - or resort to GAMLSS() family instead, provided this is used 

procast_disttree <- function(object,
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

procast_distforest <- procast_disttree

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

