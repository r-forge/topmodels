qresiduals <- function(object, ...) {
  UseMethod("qresiduals")
}

qresiduals.default <- function(object, newdata = NULL, trafo = qnorm, type = c("random", "quantile"), nsim = 1L, prob = 0.5, ...)
{
  ## type of residual for discrete distribution (if any)
  type <- match.arg(type)

  ## if 'object' is not a vector/matrix, apply p4(..., type = "probability") method
  if(is.object(object) | !is.numeric(object)) {
    y <- newresponse(object, newdata = newdata)
    object <- p4(object, newdata = newdata, at = cbind(y - .Machine$double.eps^0.8, y), type = "probability")
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
