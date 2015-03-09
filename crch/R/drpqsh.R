dcnorm <- function(x, mean = 0, sd = 1, left = 0, right = Inf, log = FALSE) {
  ifelse(x <= left, pnorm(left, mean, sd, log.p = log), 
    ifelse(x >= right, pnorm(right, mean, sd, log.p = log, lower.tail = FALSE), 
    dnorm(x, mean , sd, log = log)))
}

pcnorm <- function(q, mean = 0, sd = 1, left = 0, right = Inf, lower.tail = TRUE, 
  log.p = FALSE) {
  ifelse(x <= left, pnorm(left, mean, sd, log.p = FALSE, lower.tail = lower.tail),
    ifelse(x >= right, pnorm(right, mean, sd, log.p = FALSE, lower.tail = lower.tail),
    pnorm(x, mean, sd, log.p = FALSE, lower.tail = lower.tail)))
}

rcnorm <- function(n, mean = 0, sd = 1, left = 0, right = Inf) {
  rval <- rnorm(n, mean, sd)
  pmax(pmin(rval, right), left)
}

qcnorm <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, left = 0, right = Inf) {
 rval <- qnorm(p, mean, sd, lower.tail = lower.tail, log.p = log.p)
 pmax(pmin(rval, right), left)
}


scnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = 0, right = Inf)
{
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
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

hcnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = 0, right = Inf)
{
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  sd2 <- sd^2
  dlm <- left - mean
  drm <- right - mean
  millsl <- dnorm(left, mean, sd)/pnorm(left, mean, sd)
  millsr <- dnorm(right, mean, sd)/pnorm(right, mean, sd, lower.tail = FALSE)
  for(w in which) {
    if(w == "mu")
      hess[[w]] <- 
          ifelse(  x <= left, - dlm/sd2 * millsl - millsl^2,
            ifelse(x >= right,  drm/sd2 * millsr - millsr^2,
                                -1 / sd2))
    if(w == "sigma")
      hess[[w]] <- 
        ifelse(x <= left,    (  2 * dlm/sd2 - dlm^3/sd^3) * millsl - dlm^2/sd2 * millsl^2,
          ifelse(x >= right, (- 2 * drm/sd2 + drm^3/sd^3) * millsr - drm^2/sd2 * millsr^2,
                             (sd2 - 3 * (x - mean)^2) / sd2^2))
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- 
        ifelse(x <= left,    (  1/sd - dlm^2/sd^3) * millsl - dlm/sd * millsl^2,
          ifelse(x >= right, (- 1/sd + drm^2/sd^3) * millsr - drm/sd * millsr^2,
                             -2 * (x - mean) / sd^3))
  }
  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}

## TODO: check if fcnorm correct!
fcnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = 0, right = Inf)
{
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  n <- length(x)
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
                             2 / sd2^2))
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





###############################################################################
### Censored logistic distribution ############################################
###############################################################################


dclogis <- function(x, mean = 0, sd = 1, left = 0, right = Inf, log = FALSE) {
  ifelse(x <= left, plogis(left, mean, sd, log.p = log), 
    ifelse(x >= right, plogis(right, mean, sd, log.p = log, lower.tail = FALSE), 
    dlogis(x, mean , sd, log = log)))
}

pclogis <- function(q, mean = 0, sd = 1, left = 0, right = Inf, lower.tail = TRUE, 
  log.p = FALSE) {
  ifelse(x <= left, plogis(left, mean, sd, log.p = FALSE, lower.tail = lower.tail),
    ifelse(x >= right, plogis(right, mean, sd, log.p = FALSE, lower.tail = lower.tail),
    plogis(x, mean, sd, log.p = FALSE, lower.tail = lower.tail)))
}

rclogis <- function(n, mean = 0, sd = 1, left = 0, right = Inf) {
  rval <- rlogis(n, mean, sd)
  plogis(pmin(rval, right), left)
}

qclogis <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, left = 0, right = Inf) {
 rval <- qlogis(p, mean, sd, lower.tail = lower.tail, log.p = log.p)
 pmax(pmin(rval, right), left)
}


sclogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = 0, right = Inf)
{
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  data <- data.frame(x, mean, sd)
  x <- data$x
  mean <- data$mean
  sd <- data$sd


  score <- NULL
  dxm <- x - mean
  dlm <- left - mean
  drm <- right - mean
  sd2 <- sd^2
  millsl <- dlogis(left, mean, sd)/plogis(left, mean, sd)
  millsr <- dlogis(right, mean, sd)/plogis(right, mean, sd, lower.tail = FALSE)
  for(w in which) {
    if(w == "mu")
      score2 <- ifelse(x <= left, - millsl,
        ifelse(x >= right, millsr,
        (1 - 2 * plogis(-x, - mean, sd))/sd))
    if(w == "sigma")
      score2 <- ifelse(x <= left, - millsl * dlm/sd,
        ifelse(x >= right, millsr * drm/sd,
        (1 - 2 * plogis(-x, - mean, sd))*dxm/sd2 - 1/sd))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

hclogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = 0, right = Inf)
{
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  sd2 <- sd^2
  dlm <- left - mean
  drm <- right - mean
  dxm <- x - mean
  millsl <- dlogis(left, mean, sd)/plogis(left, mean, sd)
  millsr <- dlogis(right, mean, sd)/plogis(right, mean, sd, lower.tail = FALSE)
  scorel <- sclogis(left, mean, sd, which = "mu", left = -Inf, right = Inf)
  scorer <- sclogis(right, mean, sd, which = "mu", left = -Inf, right = Inf)
  score <- sclogis(x, mean, sd, left = left, right = right)
  for(w in which) {
    if(w == "mu")
      hess[[w]] <- 
          ifelse(  x <= left, - scorel * millsl - millsl^2,
            ifelse(x >= right,  scorer * millsr - millsr^2,
                                - 2/sd * dlogis(x, mean, sd)))
    if(w == "sigma")
      hess[[w]] <- 
        ifelse(x <= left,    (  2 * dlm/sd2 - dlm^2/sd*scorel)*millsl-millsl^2*dlm^2/sd2,
          ifelse(x >= right, (- 2 * drm/sd2 + drm^2/sd*scorer)*millsr-millsr^2*drm^2/sd2,
                             - score[,"dmu"]*dxm/sd2 - 2 * dxm^2/sd^3 * dlogis(x,mean,sd) - 
                               score[, "dsigma"]/sd)) 

    if(w %in% c("mu.sigma", "sigma.mu"))
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






g <- function(par) {
  fit <- fitfun(par)
  grad <- with(fit, sclogis(y, mu, sigma, left = left, right = right))
  grad <- cbind(grad[,1]*x, grad[,2] * fit$sigma * z)
  return(-colSums(weights * grad))
}

h <- function(par) {
  fit <- fitfun(par)
  hess <- with(fit, hclogis(y, mu, sigma, left = left, right = right,
    which = c("mu", "sigma", "mu.sigma", "sigma.mu")))
  grad <- with(fit, sclogis(y, mu, sigma, left = left, right = right, which = "sigma"))

  hess[, "d2sigma"] <- (hess[, "d2sigma"] + grad/fit$sigma)*fit$sigma^2
  hess[, "dmu.dsigma"] <- hess[, "dsigma.dmu"] <- hess[, "dmu.dsigma"]*fit$sigma
  hessmu <- crossprod(hess[,"d2mu"]*x, x)
  hessmusigma <- crossprod(hess[,"dmu.dsigma"]*x, z)
  hesssigmamu <- crossprod(hess[,"dsigma.dmu"]*z, x)
  hesssigma <- crossprod(hess[,"d2sigma"]*z, z)
  -cbind(rbind(hessmu, hesssigmamu), rbind(hessmusigma, hesssigma))
}



###############################################################################
### Censored student-t distribution ############################################
###############################################################################


dct <- function(x, df, ncp, left = 0, right = Inf, log = FALSE) {
  ifelse(x <= left, pt(left, df, ncp, log.p = log), 
    ifelse(x >= right, pt(left, df, ncp, log.p = log, lower.tail = FALSE), 
    dt(x, df, ncp, log = log)))
}

pct <- function(q, df, ncp, left = 0, right = Inf, lower.tail = TRUE, 
  log.p = FALSE) {
  ifelse(x <= left, pt(left, df, ncp, log.p = FALSE, lower.tail = lower.tail),
    ifelse(x >= right, pt(right, df, ncp, log.p = FALSE, lower.tail = lower.tail),
    pt(x, df, ncp, log.p = FALSE, lower.tail = lower.tail)))
}

rct <- function(n, df, ncp, left = 0, right = Inf) {
  rval <- rt(n, df, ncp)
  pt(pmin(rval, right), left)
}

qct <- function(p, df, ncp, lower.tail = TRUE, log.p = FALSE, left = 0, right = Inf) {
 rval <- qt(p, df, ncp, lower.tail = lower.tail, log.p = log.p)
 pmax(pmin(rval, right), left)
}


sct <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = 0, right = Inf)
{
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  data <- data.frame(x, mean, sd)
  x <- data$x
  mean <- data$mean
  sd <- data$sd


  score <- NULL
  dxm <- x - mean
  dlm <- left - mean
  drm <- right - mean
  sd2 <- sd^2
  millsl <- dlogis(left, mean, sd)/plogis(left, mean, sd)
  millsr <- dlogis(right, mean, sd)/plogis(right, mean, sd, lower.tail = FALSE)
  for(w in which) {
    if(w == "mu")
      score2 <- ifelse(x <= left, - millsl,
        ifelse(x >= right, millsr,
        dxm/sd2 * (df + 1) / (df + dxm^2/sd2)))
    if(w == "sigma")
      score2 <- ifelse(x <= left, - millsl * dlm/sd,
        ifelse(x >= right, millsr * drm/sd,
        (1 - 2 * plogis(-x, - mean, sd))*dxm/sd2 - 1/sd))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

hclogis <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"), left = 0, right = Inf)
{
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  sd2 <- sd^2
  dlm <- left - mean
  drm <- right - mean
  dxm <- x - mean
  millsl <- dlogis(left, mean, sd)/plogis(left, mean, sd)
  millsr <- dlogis(right, mean, sd)/plogis(right, mean, sd, lower.tail = FALSE)
  scorel <- sclogis(left, mean, sd, which = "mu", left = -Inf, right = Inf)
  scorer <- sclogis(right, mean, sd, which = "mu", left = -Inf, right = Inf)
  score <- sclogis(x, mean, sd, left = left, right = right)
  for(w in which) {
    if(w == "mu")
      hess[[w]] <- 
          ifelse(  x <= left, - scorel * millsl - millsl^2,
            ifelse(x >= right,  scorer * millsr - millsr^2,
                                - 2/sd * dlogis(x, mean, sd)))
    if(w == "sigma")
      hess[[w]] <- 
        ifelse(x <= left,    (  2 * dlm/sd2 - dlm^2/sd*scorel)*millsl-millsl^2*dlm^2/sd2,
          ifelse(x >= right, (- 2 * drm/sd2 + drm^2/sd*scorer)*millsr-millsr^2*drm^2/sd2,
                             - score[,"dmu"]*dxm/sd2 - 2 * dxm^2/sd^3 * dlogis(x,mean,sd) - 
                               score[, "dsigma"]/sd)) 

    if(w %in% c("mu.sigma", "sigma.mu"))
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






g <- function(par) {
  fit <- fitfun(par)
  grad <- with(fit, sclogis(y, mu, sigma, left = left, right = right))
  grad <- cbind(grad[,1]*x, grad[,2] * fit$sigma * z)
  return(-colSums(weights * grad))
}

h <- function(par) {
  fit <- fitfun(par)
  hess <- with(fit, hclogis(y, mu, sigma, left = left, right = right,
    which = c("mu", "sigma", "mu.sigma", "sigma.mu")))
  grad <- with(fit, sclogis(y, mu, sigma, left = left, right = right, which = "sigma"))

  hess[, "d2sigma"] <- (hess[, "d2sigma"] + grad/fit$sigma)*fit$sigma^2
  hess[, "dmu.dsigma"] <- hess[, "dsigma.dmu"] <- hess[, "dmu.dsigma"]*fit$sigma
  hessmu <- crossprod(hess[,"d2mu"]*x, x)
  hessmusigma <- crossprod(hess[,"dmu.dsigma"]*x, z)
  hesssigmamu <- crossprod(hess[,"dsigma.dmu"]*z, x)
  hesssigma <- crossprod(hess[,"d2sigma"]*z, z)
  -cbind(rbind(hessmu, hesssigmamu), rbind(hessmusigma, hesssigma))
}



