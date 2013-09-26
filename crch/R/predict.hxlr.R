predict.hxlr <-
function(object, newdata = object$data, type = c("class", "probability", "cumprob", "location", "scale"), thresholds = object$thresholds, ...) {
  type <- match.arg(type)
  
  ## get coefficients
  gamma <- object$coefficients2$location
  delta <- object$coefficients2$scale

  ## matrices for location and scale model from newdata
  x <- model.matrix(delete.response(terms(object$formula)), newdata)
  z <- model.matrix(object$scale, newdata)
  
  ## location and scale 
  mu <- drop(x %*% gamma)
  sigma <- exp(drop(z %*% delta))

  ## cumulative probabilities P(y<threshold[i]|x)
  cumprob <- NULL
  for(i in 1:length(thresholds)) { 
    cumprob <- cbind(cumprob, plogis((thresholds[i] - mu)/sigma))
  }
  
  ## category probabilities P(threshold[i-1]<=y<threshold[i]|x)
  prob <- t(apply(cbind(0, cumprob, 1), 1, diff))

  rval <- switch(type,
      "location" = mu,
      "scale" = sigma,
      "cumprob" = cumprob,
      "probability" = prob,
      "class" = apply(prob, 1, which.max)
    )

  return(rval)
}
