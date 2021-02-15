qresiduals <- function(object, ...) {
  UseMethod("qresiduals")
}

qresiduals.default <- function(object, newdata = NULL, trafo = qnorm, type = c("random", "quantile"), nsim = 1L, prob = 0.5, ...)
{
  ## type of residual for discrete distribution (if any)
  type <- match.arg(type)

  ## if 'object' is not a vector/matrix, apply procast(..., type = "probability") method
  if(is.object(object) | !is.numeric(object)) {
    y <- newresponse(object, newdata = newdata)
    object <- procast(object, newdata = newdata, at = cbind(y - .Machine$double.eps^0.8, y), type = "probability")
    if(inherits(object, "try-error")) stop("could not obtain probability integral transform from 'object'")
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
      ## FIXME: probably needs to be done on the transformed rather than the uniform scale...
      nam <- rownames(object)
      object <- object[, 1L]  %*% t(1 - prob) + object[, 2L] %*% t(prob)
      dimnames(object) <- list(nam, paste("q", prob, sep = "_"))
    }
    nc <- NCOL(object)
  }
  if(!is.null(dim(object)) & nc == 1L) object <- drop(object)

  ## compute quantile residuals  
  if(!is.null(trafo)) object <- trafo(object)
  return(object)
}


qresiduals.crch <- function(object, 
                            newdata = NULL, 
                            trafo = qnorm, 
                            type = c("random", "quantile"), 
                            nsim = 1L, 
                            prob = 0.5, 
                            mass_redist = c("quantile", "random", "none"),
                            ...) {
# TODO: (ML) 
# * Better argument name for `mass_redist' and for its types. 
# * Is really an own S3 method necessary or would be an extension of the default method 
#   be more appropriate? 
# * Or even adapt in procast(), but here problem with missing y. 
# * Maybe best to save truncation (e.g., crch, gamlss) as attribute in return value of `procast()' 
#   and check in `qresiduals.default()' for truncation?!

  ## type of residual for discrete distribution (if any)
  type <- match.arg(type)
  mass_redist <- match.arg(mass_redist)

  ## if 'object' is not a vector/matrix, apply procast(..., type = "probability") method
  if(is.object(object) | !is.numeric(object)) {
    y <- newresponse(object, newdata = newdata)

    object <- procast(object, newdata = newdata, 
      at = cbind(y - .Machine$double.eps^0.8, y), type = "probability")

    if (inherits(object, "try-error")) {
      stop("could not obtain probability integral transform from 'object'")
    }

    # If `object' is censored, point mass must be randomly or evenly redistributed
    if ((mass_redist != "none") & any(sapply(attr(object, "cens"), is.finite))) {
      left <- attr(object, "cens")$left 
      right <- attr(object, "cens")$right
      
      idx <- which(y <= left | y >= right)
      if (mass_redist == "quantile") {
        # TODO: (ML) Check which approach is correct

        # Get quantiles equally distributed between 0 and 1
        #p <- seq(length.out = length(idx)) / (length(idx) + 1)  

        # Get quantiles equally distributed within potential bins
        p <- seq(0, 1, length.out = length(idx) + 1L)
        p <- filter(p, rep(1 / 2, 2), sides = 2)[-length(p)]
        object[idx, ] <- sample(qunif(p))
      } else {
        object[idx, ] <- runif(length(idx))
      }
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
      ## FIXME: probably needs to be done on the transformed rather than the uniform scale...
      nam <- rownames(object)
      object <- object[, 1L]  %*% t(1 - prob) + object[, 2L] %*% t(prob)
      dimnames(object) <- list(nam, paste("q", prob, sep = "_"))
    }
    nc <- NCOL(object)
  }
  if(!is.null(dim(object)) & nc == 1L) object <- drop(object)

  ## compute quantile residuals  
  if(!is.null(trafo)) object <- trafo(object)
  return(object)
}
