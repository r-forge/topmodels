qresiduals <- function(object, ...) {
  UseMethod("qresiduals")
}

qresiduals.default <- function(object, 
                               newdata = NULL, 
                               trafo = qnorm, 
                               type = c("random", "quantile"), 
                               nsim = 1L, 
                               prob = 0.5, 
                               ...) {
# TODO: (ML) 
# * Better argument name for `mass_redist' and for its types. 
# * Or is an own S3 method necessary be more appropriate (previous version)?
# * Or even adapt in procast(), but here problem with missing y. 
# * Now save truncation (e.g., crch, gamlss) as attribute in return value of `procast()' 
#   and check in `qresiduals.default()' for censor point?!

  ## type of residual for discrete distribution (if any)
  type <- match.arg(type)

  ## if 'object' is not a vector/matrix, apply procast(..., type = "probability") method
  if(is.object(object) | !is.numeric(object)) {
    y <- newresponse(object, newdata = newdata)

    object <- procast(object, newdata = newdata, 
      at = cbind(y - .Machine$double.eps^0.8, y), type = "probability")

    # TODO: (ML) There is no `try()` environment, which errors can be caught
    if (inherits(object, "try-error")) {
      stop("could not obtain probability integral transform from 'object'")
    }

    # TODO: (ML) This is actually not needed, but must we make sure that not values
    # below/above censoring point are within y?!
    ## If `object' is censored, point mass must be randomly or evenly redistributed
    #if (mass_redist & any(sapply(attr(object, "cens"), is.finite))) {
    #  left <- attr(object, "cens")$left 
    #  right <- attr(object, "cens")$right
    # 
    #  idx <- which(y <= left | y >= right)
    #  object[idx, ] <- object[idx, ] * runif(length(idx))
    #} 
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
