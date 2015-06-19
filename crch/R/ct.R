## density
dct <- function(x, mean = 0, sd = 1, df, left = -Inf, right = Inf, log = FALSE) {
  ifelse(x <= left, pt((left - mean)/sd, df = df, log.p = log), 
  ifelse(x >= right, pt((right - mean)/sd, df = df, log.p = log, lower.tail = FALSE), 
  dt((x - mean)/sd, df = df, log = log)/sd^(1 - log) - log(sd) * log))
}

## distribution function
pct <- function(q, mean = 0, sd = 1,  df, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  q <- data.frame(q = q, mean, sd)$q
  if(lower.tail){
    ifelse(q < left, 0, 
    ifelse(q >= right, 1, 
    pt((q - mean)/sd, df = df, log.p = log.p)))
  } else {
    ifelse(q <= left, 1, 
    ifelse(q > right, 0, 
    pt((q - mean)/sd, df = df, lower.tail = FALSE, log.p = log.p)))
  }
}

## random numbers
rct <- function(n, mean = 0, sd = 1, df, left = -Inf, right = Inf) {
  rval <- rt(n, df = df) * sd + mean
  pmax(pmin(rval, right), left)
}

## quantiles
qct <- function(p, mean = 0, sd = 1, df, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  rval <- qt(p, df = df, lower.tail = lower.tail, log.p = log.p) * sd + mean
  pmax(pmin(rval, right), left)
}



## scores
sct <- function(x, mean = 0, sd = 1, df, which = c("mu", "sigma"), 
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
  millsl <- dt(dlm/sd, df = df)/sd/pt(dlm/sd, df = df)
  millsr <- dt(drm/sd, df = df)/sd/pt(drm/sd, df = df, lower.tail = FALSE)
  for(w in which) {
    if(w == "mu")
      score2 <- ifelse(x <= left, - millsl,
        ifelse(x >= right, millsr,
        dxm/sd2 * (df + 1) / (df + dxm^2/sd2)))
    if(w == "sigma")
      score2 <- ifelse(x <= left, - millsl * dlm/sd,
        ifelse(x >= right, millsr * drm/sd,
        dxm^2/sd^3 * (df + 1) / (df + dxm^2/sd2) - 1/sd))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
hct <- function(x, mean = 0, sd = 1, df, which = c("mu", "sigma"), 
  left = -Inf, right = Inf)
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
  millsl <- dt(dlm/sd, df = df)/sd/pt(dlm/sd, df = df)
  millsr <- dt(drm/sd, df = df)/sd/pt(drm/sd, df = df, lower.tail = FALSE)
  scorel <- sct(left, mean, sd, df, which = "mu", left = -Inf, right = Inf)
  scorer <- sct(right, mean, sd, df, which = "mu", left = -Inf, right = Inf)
  score  <- sct(x, mean, sd, df, left = -Inf, right = Inf)
  for(w in which) {
    if(w == "mu")
      hess[[w]] <- 
          ifelse(  x <= left, - scorel * millsl - millsl^2,
            ifelse(x >= right,  scorer * millsr - millsr^2,
                                (df + 1)*(dxm^2 - df*sd2) / (df*sd2 + dxm^2)^2))
    if(w == "sigma")
      hess[[w]] <- 
        ifelse(x <= left,    (  2 * dlm/sd2 - dlm^2/sd2*scorel)*millsl -
          millsl^2*dlm^2/sd2,
          ifelse(x >= right, (- 2 * drm/sd2 + drm^2/sd2*scorer)*millsr -
            millsr^2*drm^2/sd2,
                             dxm^2 * (df + 1) * (-3*sd2*df -dxm^2) / 
                               (sd2 * (df*sd2 + dxm^2)^2) + 1/sd2))

    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- 
        ifelse(x <= left,    (  1/sd - dlm/sd*scorel) * millsl - 
                              dlm/sd * millsl^2,
          ifelse(x >= right, (- 1/sd + drm/sd*scorer) * millsr - 
                              drm/sd * millsr^2,
                             - 2* dxm * (df + 1) *sd * df / (df*sd2 + dxm^2)^2))
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}
