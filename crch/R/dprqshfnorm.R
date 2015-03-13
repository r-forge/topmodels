dcnorm <- function(x, ..., left = -Inf, right = Inf, log = FALSE) {
  ifelse(x <= left, pnorm(left, ..., log.p = log), 
  ifelse(x >= right, pnorm(right,  ..., log.p = log, lower.tail = FALSE), 
  dnorm(x, ..., log = log)))
}

pcnorm <- function(q, ..., left = -Inf, right = Inf) {
  ifelse(q < left, 0, 
  ifelse(q >= right, 1, 
  pnorm(q, ..., log = log)))
}

rcnorm <- function(n, ..., left = -Inf, right = Inf) {
  rval <- rnorm(n, ...)
  pmax(pmin(rval, right), left)
}

qcnorm <- function(p, ..., left = -Inf, right = Inf) {
  rval <- qnorm(p, ...)
  pmax(pmin(rval, right), left)
}


scnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = -Inf, right = Inf) {
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


hcnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = -Inf, right = Inf)
{
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
        ifelse(x <= left,    ( 2 * dlm/sd2 - dlm^2/sd2*scorel)*millsl-millsl^2*dlm^2/sd2, 
          ifelse(x >= right, (- 2 * drm/sd2 + drm^2/sd2*scorer)*millsr-millsr^2*drm^2/sd2,
                             (sd2 - 3 * (x - mean)^2) / sd2^2)) 

    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- 
        ifelse(x <= left,    (  1/sd - dlm/sd*scorel) * millsl - dlm/sd * millsl^2,
          ifelse(x >= right, (- 1/sd + drm/sd*scorer) * millsr - drm/sd * millsr^2,
                             -2 * dxm / sd^3))
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}


## TODO: check if fcens correct!

fcnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = -Inf, right = Inf)
{
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  x <- data.frame(x = x, mean, sd)$x
  fish <- list()
  sd2 <- sd^2
  dlm <- left - mean
  drm <- right - mean
  millsl <- dnorm(left, mean, sd)/pnorm(left, mean, sd)
  millsr <- dnorm(right, mean, sd)/pnorm(right, mean, sd, lower.tail = FALSE)
  for(w in which) {
    if(w == "mu")
      fish[[w]] <- 
        ifelse(  x <= left,    dlm/sd2 * millsl + millsl^2,
          ifelse(x >= right, - drm/sd2 * millsr + millsr^2,
                               1 / sd2))
    if(w == "sigma")
      fish[[w]] <- 
        ifelse(x <= left,    (  2 * dlm/sd2 - dlm^3/sd^3) * millsl - dlm^2/sd2 * millsl^2,
          ifelse(x >= right, (- 2 * drm/sd2 + drm^3/sd^3) * millsr - drm^2/sd2 * millsr^2,
                              2 / sd2))
    if(w %in% c("mu.sigma", "sigma.mu"))
      fish[[w]] <- 
        ifelse(x <= left,    (  1/sd - dlm^2/sd^3) * millsl - dlm/sd * millsl^2,
          ifelse(x >= right, (- 1/sd + drm^2/sd^3) * millsr - drm/sd * millsr^2,
                             0))
  }
  fish <- do.call("cbind", fish)
  colnames(fish) <- gsub("mu", "dmu", colnames(fish))
  colnames(fish) <- gsub("sigma", "dsigma", colnames(fish))
  colnames(fish)[colnames(fish) == "dmu"] <- "d2mu"
  colnames(fish)[colnames(fish) == "dsigma"] <- "d2sigma"
  fish
}


