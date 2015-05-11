## density
dtlogis <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  x <- data.frame(x = x, mean, sd)$x
  denom <- plogis((right - mean)/sd) - plogis((left-mean)/sd) 
  ifelse(x < left, dlogis(-Inf, log = log), 
  ifelse(x > right, dlogis(Inf, log = log), 
  dlogis((x - mean)/sd, log = log)/(sd*denom)^(1 - log) - log(sd*denom)*log))
}

## distribution function
ptlogis <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  q <- data.frame(q = q, mean, sd)$q
  denom <- plogis((right - mean)/sd) - plogis((left - mean)/sd)
  lower <- if(lower.tail) left else right
  qtmp <- (plogis((lower - mean)/sd) - plogis((q - mean)/sd)) *
    (-1)^lower.tail
  ifelse(q < left, plogis(-Inf, lower.tail = lower.tail, log.p = log.p),
  ifelse(q > right, plogis(Inf, lower.tail = lower.tail, log.p = log.p),
  if(log.p) log(qtmp) - log(denom) else qtmp/denom ))
}

## quantiles
qtlogis <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  upper <- if(lower.tail) right else left
  p <- plogis((lower-mean)/sd, lower.tail = lower.tail) * (1 - p) + 
    p*plogis((upper - mean)/sd, lower.tail = lower.tail)
  qlogis(p, lower.tail = lower.tail)*sd + mean
}

## random numbers
rtlogis <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  qtlogis(runif(n), mean, sd, left = left, right = right)
}

## scores
stlogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  x <- data.frame(x = x, mean, sd)$x
  dxm <- x - mean
  dlm <- left - mean
  drm <- right - mean
  dlm2 <- if(is.finite(left)) dlm else 0
  drm2 <- if(is.finite(right)) drm else 0
  sd2 <- sd^2
  
  sc <- sclogis(x, mean, sd, left = -Inf, right = Inf)

  denom <-  (plogis(drm/sd) - plogis(dlm/sd))
  enum1 <- (dlogis(drm/sd) - dlogis(dlm/sd))/sd
  enum2 <- (drm2*dlogis(drm/sd) - dlm2*dlogis(dlm/sd))/sd2
  
  for(w in which) {
    if(w == "mu")
      score2 <- ifelse(x < left, 0,
        ifelse(x > right, 0,
        sc[,"dmu"] + enum1/denom))
    if(w == "sigma")
      score2 <- ifelse(x < left, 0,
        ifelse(x > right, 0,
        sc[,"dsigma"] + enum2/denom))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
htlogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  x <- data.frame(x = x, mean, sd)$x
  sd2 <- sd^2
  dlm <- left - mean
  drm <- right - mean
  dxm <- x - mean
  
  scorel <- if(is.finite(left)) sclogis(left, mean, sd, left = -Inf, right = Inf) else 
    data.frame(dmu = 0, dsigma = 0)
  scorer <- if(is.finite(right)) sclogis(right, mean, sd, left = -Inf, right = Inf) else 
    data.frame(dmu = 0, dsigma = 0)
  
  hc <- hclogis(x, mean, sd, which = c("mu", "sigma", "mu.sigma", "sigma.mu"), 
    left = -Inf, right = Inf)

  denom <-  (plogis(drm/sd) - plogis(dlm/sd))
  enum1 <- (dlogis(drm/sd) - dlogis(dlm/sd))/sd
  drm2 <- if(is.finite(right)) drm else 0
  dlm2 <- if(is.finite(left)) dlm else 0
  enum2 <- (drm2*dlogis(drm/sd) - dlm2*dlogis(dlm/sd))/sd2
  enum3 <- (scorer[,"dmu"]*dlogis(drm/sd)/sd - scorel[,"dmu"]*dlogis(dlm/sd)/sd)
  enum4 <- drm2/sd2*dlogis(drm/sd)*(scorer[,"dsigma"] - 1/sd) - 
             dlm2/sd2*dlogis(dlm/sd)*(scorel[,"dsigma"] - 1/sd)
  enum5 <- (scorer[,"dsigma"]*dlogis(drm/sd)/sd - scorel[,"dsigma"]*dlogis(dlm/sd)/sd)
  
  

  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- 
        ifelse(  x < left,  0,
          ifelse(x > right, 0,
                            hc[,"d2mu"] + enum1^2/denom^2 + enum3/denom))

    if(w == "sigma")
      hess[[w]] <-           
        ifelse(x < left,    0, 
          ifelse(x > right, 0,
                            hc[,"d2sigma"] + enum2^2/denom^2 + enum4/denom)) 

    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- 
        ifelse(x < left,    0,
          ifelse(x > right, 0,
                            hc[,"dmu.dsigma"] +  enum5/denom + enum1*enum2/denom^2))
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}
