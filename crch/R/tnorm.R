## density
dtnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  x <- data.frame(x = x, mean, sd)$x
  denom <- pnorm((right - mean)/sd) - pnorm((left - mean)/sd) 
  ifelse(x < left, dnorm(-Inf, log = log), 
  ifelse(x > right, dnorm(Inf, log = log), 
  dnorm((x - mean)/sd, log = log)/(sd*denom)^(1 - log) - log(sd*denom)*log))
}

## distribution function
ptnorm <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  q <- data.frame(q = q, mean, sd)$q
  denom <- pnorm((right - mean)/sd) - pnorm((left - mean)/sd)
  lower <- if(lower.tail) left else right
  qtmp <- (pnorm((lower - mean)/sd) - pnorm((q - mean)/sd)) *
    (-1)^lower.tail
  ifelse(q < left, pnorm(-Inf, lower.tail = lower.tail, log.p = log.p),
  ifelse(q > right, pnorm(Inf, lower.tail = lower.tail, log.p = log.p),
  if(log.p) log(qtmp) - log(denom) else qtmp/denom ))
}

## quantiles
qtnorm <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  upper <- if(lower.tail) right else left
  p <- pnorm((lower-mean)/sd, lower.tail = lower.tail) * (1 - p) + 
    p*pnorm((upper - mean)/sd, lower.tail = lower.tail)
  qnorm(p, lower.tail = lower.tail)*sd + mean
}

## random numbers
rtnorm <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  qtnorm(runif(n), mean, sd, left = left, right = right)
}

## scores
stnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
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
  
  sc <- scnorm(x, mean, sd, left = -Inf, right = Inf)

  denom <-  (pnorm(drm/sd) - pnorm(dlm/sd))
  enum1 <- (dnorm(drm/sd) - dnorm(dlm/sd))/sd
  enum2 <- (drm2*dnorm(drm/sd) - dlm2*dnorm(dlm/sd))/sd2
  
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
htnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"),
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
  
  scorel <- if(is.finite(left)) scnorm(left, mean, sd, left = -Inf, right = Inf) else 
    data.frame(dmu = 0, dsigma = 0)
  scorer <- if(is.finite(right)) scnorm(right, mean, sd, left = -Inf, right = Inf) else 
    data.frame(dmu = 0, dsigma = 0)
  
  hc <- hcnorm(x, mean, sd, which = c("mu", "sigma", "mu.sigma", "sigma.mu"), 
    left = -Inf, right = Inf)

  denom <-  (pnorm(drm/sd) - pnorm(dlm/sd))
  enum1 <- (dnorm(drm/sd) - dnorm(dlm/sd))/sd
  drm2 <- if(is.finite(right)) drm else 0
  dlm2 <- if(is.finite(left)) dlm else 0
  enum2 <- (drm2*dnorm(drm/sd) - dlm2*dnorm(dlm/sd))/sd2
  enum3 <- (scorer[,"dmu"]*dnorm(drm/sd)/sd - scorel[,"dmu"]*dnorm(dlm/sd)/sd)
  enum4 <- drm2/sd2*dnorm(drm/sd)*(scorer[,"dsigma"] - 1/sd) - 
             dlm2/sd2*dnorm(dlm/sd)*(scorel[,"dsigma"] - 1/sd)
  enum5 <- (scorer[,"dsigma"]*dnorm(drm/sd)/sd - scorel[,"dsigma"]*dnorm(dlm/sd)/sd)
  
  

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
