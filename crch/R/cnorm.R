## density
dcnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  x <- data.frame(x = x, mean, sd)$x
  ifelse(x <= left, pnorm((left-mean)/sd, log.p = log), 
  ifelse(x >= right, pnorm((right-mean)/sd, log.p = log, lower.tail = FALSE), 
  dnorm((x-mean)/sd, log = log)/sd^(1 - log) - log(sd) * log))
}

## distribution function
pcnorm <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  q <- data.frame(q = q, mean, sd)$q
  if(lower.tail){ 
    ifelse(q < left, 0, 
    ifelse(q >= right, 1, 
    pnorm((q-mean)/sd, log.p = log.p)))
  } else {
    ifelse(q <= left, 1, 
    ifelse(q > right, 0, 
    pnorm((q-mean)/sd, lower.tail = FALSE, log.p = log.p)))
  }
}

## random numbers
rcnorm <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  rval <- rnorm(n) * sd + mean
  pmax(pmin(rval, right), left)
}

## quantiles
qcnorm <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  rval <- qnorm(p, lower.tail = lower.tail, log.p = log.p) * sd + mean
  pmax(pmin(rval, right), left)
}

## scores
scnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
  left = -Inf, right = Inf) {
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  x <- data.frame(x = x, mean, sd)$x
  dxm <- x - mean
  dlm <- left - mean
  drm <- right - mean
  sd2 <- sd^2
  millsl <- dnorm(left, mean, sd)/pnorm(left, mean, sd)
  millsr <- dnorm(right, mean, sd)/pnorm(right, mean, sd, lower.tail = FALSE)
  for(w in which) {
    if(w == "mu")
      score2 <- ifelse(x <= left, - millsl,
        ifelse(x >= right, millsr,
        dxm / sd2))
    if(w == "sigma")
      score2 <- ifelse(x <= left, - millsl * dlm/sd,
        ifelse(x >= right, millsr * drm/sd,
        (dxm^2 - sd2) / sd^3))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
hcnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), 
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
  millsl <- dnorm(left, mean, sd)/pnorm(left, mean, sd)
  millsr <- dnorm(right, mean, sd)/pnorm(right, mean, sd, lower.tail = FALSE)
  scorel <- scnorm(left, mean, sd, which = "mu", left = -Inf, right = Inf)
  scorer <- scnorm(right, mean, sd, which = "mu", left = -Inf, right = Inf)

  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- 
          ifelse(  x <= left, - scorel * millsl - millsl^2,
            ifelse(x >= right,  scorer * millsr - millsr^2,
                                -1 / sd2))
    
    if(w == "sigma")
      hess[[w]] <-           
        ifelse(x <= left,    ( 2 * dlm/sd2 - dlm^2/sd2*scorel)*millsl-
          millsl^2*dlm^2/sd2, 
          ifelse(x >= right, (- 2 * drm/sd2 + drm^2/sd2*scorer)*millsr-
            millsr^2*drm^2/sd2, (sd2 - 3 * (x - mean)^2) / sd2^2)) 

    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- 
        ifelse(x <= left,    (  1/sd - dlm/sd*scorel) * millsl - 
          dlm/sd * millsl^2,
          ifelse(x >= right, (- 1/sd + drm/sd*scorer) * millsr - 
            drm/sd * millsr^2, -2 * dxm / sd^3))
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}
