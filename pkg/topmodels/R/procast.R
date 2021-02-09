## probabilistic forecasting method: procast


procast <- function(object, newdata = NULL, na.action = na.pass, type = "quantile", at = 0.5, ...)
{
  UseMethod("procast")
}

procast.lm <- function(object, newdata = NULL, na.action = na.pass,
                  type = c("quantile", "mean", "variance", "parameter", "density", "probability", "score"),
                  at = 0.5, drop = FALSE, ...)
{
  ## predicted means
  pars <- if(missing(newdata) || is.null(newdata)) {
    object$fitted.values
  } else {
    predict(object, newdata = newdata, na.action = na.action)
  }
  pars <- data.frame(mu = pars)

  ## add maximum likelihood estimator of constant varians
  pars$sigma <- summary(object)$sigma * sqrt(df.residual(object)/nobs(object))
  ## NOTE: (ML) check if constant sigma in case of `na.action' is what we want

  ## types of predictions
  type <- match.arg(type, c("quantile", "mean", "variance", "parameter", "density", "probability", "score"))

  ## set up function that computes prediction from model parameters
  FUN <- switch(type,
    "quantile" = function(at, pars, ...) qnorm(at, mean = pars$mu, sd = pars$sigma, ...),
    "mean" = function(pars) pars$mu,
    "variance" = function(pars) pars$sigma^2,
    "parameter" = function(pars) pars,
    "density" = function(at, pars, ...) dnorm(at, mean = pars$mu, sd = pars$sigma, ...),
    "probability" = function(at, pars, ...) pnorm(at, mean = pars$mu, sd = pars$sigma, ...),
    "score" = function(at, pars, ...) (at$y - pars$mu)^2 * at$x
  )

  ## NOTE: for scores we need at = data.frame(y = response, x = model.matrix)

  procast_setup(pars, FUN = FUN, at = at, drop = drop, type = type, ...)
}

procast_setup <- function(pars, FUN, at = NULL, drop = FALSE, type = "procast", ...)
{
  ## - is 'at' some kind of 'data'
  ## - or the special type 'list' or 'function'
  ## - or missing altogether ('none')
  attype <- if(is.null(at) || names(formals(FUN))[1L] == "pars") {
    "none"
  } else if(is.character(at)) {
    match.arg(at, c("function", "list"))
  } else {
    "data"
  }

  if(attype == "none") {

    ## if 'at' is missing: prediction is just a transformation of the parameters
    rval <- FUN(pars, ...)
    if(is.null(dim(rval))) names(rval) <- rownames(pars)

  } else {

    ## if 'at' is 'list':
    ## prediction is a list of functions with predicted parameters as default
    if(attype == "list") {

      FUN2 <- function(at, pars = pars, ...) NULL
      body(FUN2) <- call("FUN",
        at = as.name("at"),
        pars = as.call(structure(sapply(c("data.frame", names(pars)), as.name), .Names = c("", names(pars)))),
        as.name("..."))
      ff <- formals(FUN2)

      rval <- lapply(1L:nrow(pars), function(i) {
        FUN2i <- FUN2
        formals(FUN2i) <- as.pairlist(c(ff[1L], as.pairlist(as.list(pars[i, ])), ff[3L]))
        FUN2i
      })
      return(rval)
    
    ## otherwise ('at' is 'data' or 'function'):
    ## set up a function that suitably expands 'at' (if necessary)
    ## and then either evaluates it at the predicted parameters ('data')
    ## or sets up a user interface with predicted parameters as default ('function')
    } else {
      FUN2a <- function(at, pars, ...) {
    	  n <- NROW(pars)
	      if(!is.data.frame(at)) {
          if(length(at) == 1L) at <- rep.int(as.vector(at), n)
          if(length(at) != n) at <- rbind(at)
	      }
    	  if(is.matrix(at) && NROW(at) == 1L) {
    	    at <- matrix(rep(at, each = n), nrow = n)
    	    rv <- FUN(as.vector(at), pars = pars[rep(1L:n, ncol(at)), ], ...)
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
      
      if(attype == "function") {
        rval <- function(at, pars = pars, ...) NULL
        body(rval) <- call("FUN2a",
          at = as.name("at"),
          pars = as.call(structure(sapply(c("data.frame", names(pars)), as.name), .Names = c("", names(pars)))),
          as.name("..."))
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

  ## default: return a data.frame (drop=FALSE) but optionally this can
  ## be dropped to a vector if possible (drop=TRUE)
  if(drop) {
    if (!is.null(dim(rval)) && NCOL(rval) == 1L) {
      # NOTE: (ML) How can condition be fulfilled? Compare code coverage.
      rval <- drop(rval)
    }
  } else {
    if(is.null(dim(rval))) {
      rval <- as.matrix(rval)
      if(ncol(rval) == 1L) colnames(rval) <- type
    }
    if(!inherits(rval, "data.frame")) rval <- as.data.frame(rval)
  }

  return(rval)
}
