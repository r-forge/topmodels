dclogis <- function(x, mean, sd, left = -Inf, right = Inf, log = FALSE) {
  ifelse(x <= left, plogis((left - mean)/sd, log.p = log), 
  ifelse(x >= right, plogis((right - mean)/sd, log.p = log, lower.tail = FALSE), 
  dlogis((x - mean)/sd, log = log)/sd^(1 - log) - log(sd) * log))
}

pclogis <- function(q, mean, sd, left = -Inf, right = Inf) {
  ifelse(q < left, 0, 
  ifelse(q >= right, 1, 
  plogis((q - mean)/sd, log = log)))
}

rclogis <- function(n, mean, sd, left = -Inf, right = Inf) {
  rval <- rlogis(n) * sd + mean
  pmax(pmin(rval, right), left)
}

qclogis <- function(p, mean, sd, left = -Inf, right = Inf) {
  rval <- qlogis(p) * sd + mean
  pmax(pmin(rval, right), left)
}


sclogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = -Inf, right = Inf) {
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  x <- data.frame(x = x, mean, sd)$x
  dxm <- x - mean
  dlm <- left - mean
  drm <- right - mean 
  sd2 <- sd^2
  millsl <- dlogis(dlm/sd)/sd/plogis(dlm/sd)
  millsr <- dlogis(drm/sd)/sd/plogis(drm/sd, lower.tail = FALSE)
  for(w in which) {
    if(w == "mu")
      score2 <- ifelse(x <= left, - millsl,
        ifelse(x >= right, millsr,
        (1 - 2 * plogis(-dxm/sd))/sd))
    if(w == "sigma")
      score2 <- ifelse(x <= left, - millsl * dlm/sd,
        ifelse(x >= right, millsr * drm/sd,
        (1 - 2 * plogis(-dxm/sd))*dxm/sd2 - 1/sd))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}


hclogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = -Inf, right = Inf)
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
  millsl <- dlogis(dlm/sd)/sd/plogis(dlm/sd)
  millsr <- dlogis(drm/sd)/sd/plogis(drm/sd, lower.tail = FALSE)
  scorel <- sclogis(left, mean, sd, which = "mu", left = -Inf, right = Inf)
  scorer <- sclogis(right, mean, sd, which = "mu", left = -Inf, right = Inf)
  score  <- sclogis(x, mean, sd, left = -Inf, right = Inf)
  for(w in which) {
    if(w == "mu")
      hess[[w]] <- 
          ifelse(  x <= left, - scorel * millsl - millsl^2,
            ifelse(x >= right,  scorer * millsr - millsr^2,
                                - 2/sd2 * dlogis(dxm/sd)))
    if(w == "sigma")
      hess[[w]] <- 
        ifelse(x <= left,    (  2 * dlm/sd2 - dlm^2/sd2*scorel)*millsl-millsl^2*dlm^2/sd2,
          ifelse(x >= right, (- 2 * drm/sd2 + drm^2/sd2*scorer)*millsr-millsr^2*drm^2/sd2,
                             - score[,"dmu"]*dxm/sd2 - 2 * dxm^2/sd2^2 * dlogis(dxm/sd) - 
                               score[, "dsigma"]/sd))

    if(w %in% c("mu.sigma", "sigma.mu"))## not correct
      hess[[w]] <- 
        ifelse(x <= left,    (  1/sd - dlm/sd*scorel) * millsl - dlm/sd * millsl^2,
          ifelse(x >= right, (- 1/sd + drm/sd*scorer) * millsr - drm/sd * millsr^2,
                             -score[, "dmu"]/sd - 2*dxm/sd^3*dlogis(dxm/sd)))
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}


## TODO: check if fcens correct!

fclogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = -Inf, right = Inf)
{
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  x <- data.frame(x = x, mean, sd)$x
  fish <- list()
  sd2 <- sd^2
  dlm <- left - mean
  drm <- right - mean
  dxm <- x - mean
  millsl <- d(left, mean, sd)/p(left, mean, sd)
  millsr <- d(right, mean, sd)/p(right, mean, sd, lower.tail = FALSE)
  for(w in which) {
    if(w == "mu")
      fish[[w]] <- 
        ifelse(  x <= left,    dlm/sd2 * millsl + millsl^2,
          ifelse(x >= right, - drm/sd2 * millsr + millsr^2,
                               - 2/sd/6 ))
    if(w == "sigma")
      fish[[w]] <- 
        ifelse(x <= left,    (  2 * dlm/sd2 - dlm^3/sd^3) * millsl - dlm^2/sd2 * millsl^2,
          ifelse(x >= right, (- 2 * drm/sd2 + drm^3/sd^3) * millsr - drm^2/sd2 * millsr^2,
                               - 2/sd/6 - 1/sd2))
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


