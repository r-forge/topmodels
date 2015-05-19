## density
dtt <- function(x, mean = 0, sd = 1, df, left = -Inf, right = Inf, log = FALSE) {
  x <- data.frame(x = x, mean, sd)$x
  denom <- pt((right - mean)/sd, df = df) - pt((left-mean)/sd, df = df) 
  ifelse(x < left, dt(-Inf, df = df, log = log), 
  ifelse(x > right, dt(Inf, df = df, log = log), 
  dt((x - mean)/sd, df = df, log = log)/(sd*denom)^(1 - log) - log(sd*denom)*log))
}

## distribution function
ptt <- function(q, mean = 0, sd = 1, df, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  q <- data.frame(q = q, mean, sd)$q
  denom <- pt((right - mean)/sd, df = df) - pt((left - mean)/sd, df = df)
  lower <- if(lower.tail) left else right
  qtmp <- (pt((lower - mean)/sd, df = df) - pt((q - mean)/sd, df = df)) *
    (-1)^lower.tail
  ifelse(q < left, pt(-Inf, df = df, lower.tail = lower.tail, log.p = log.p),
  ifelse(q > right, pt(Inf, df = df, lower.tail = lower.tail, log.p = log.p),
  if(log.p) log(qtmp) - log(denom) else qtmp/denom ))
}

## quantiles
qtt <- function(p, mean = 0, sd = 1, df, lower.tail = TRUE, log.p = FALSE, 
  left = -Inf, right = Inf) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  p <- pt((lower-mean)/sd, df = df) * (1 - p) + p*pt((right - mean)/sd, df = df)
  qt((p - mean)/sd, df = df, lower.tail = lower.tail)
}

## random numbers
rtt <- function(n, mean = 0, sd = 1, df, left = -Inf, right = Inf) {
  qtt((runif(n) - mean)/sd, df = df, left = left, right = right)
}

## scores
stt <- function(x, mean = 0, sd = 1, df, which = c("mu", "sigma"), 
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
  
  sc <- sct(x, mean, sd, df = df, left = -Inf, right = Inf)

  denom <-  (pt(drm/sd, df = df) - pt(dlm/sd, df = df))
  enum1 <- (dt(drm/sd, df = df) - dt(dlm/sd, df = df))/sd
  enum2 <- (drm2*dt(drm/sd, df = df) - dlm2*dt(dlm/sd, df = df))/sd2
  
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
htt <- function(x, mean = 0, sd = 1, df, which = c("mu", "sigma"), 
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
  
  scorel <- if(is.finite(left)) sct(left, mean, sd, df, left = -Inf, right = Inf) else 
    data.frame(dmu = 0, dsigma = 0)
  scorer <- if(is.finite(right)) sct(right, mean, sd, df, left = -Inf, right = Inf) else 
    data.frame(dmu = 0, dsigma = 0)
  
  hc <- hct(x, mean, sd, df, which = c("mu", "sigma", "mu.sigma", "sigma.mu"), 
    left = -Inf, right = Inf)

  denom <-  (pt(drm/sd, df = df) - pt(dlm/sd, df = df))
  enum1 <- (dt(drm/sd, df = df) - dt(dlm/sd, df = df))/sd
  drm2 <- if(is.finite(right)) drm else 0
  dlm2 <- if(is.finite(left)) dlm else 0
  enum2 <- (drm2*dt(drm/sd, df = df) - dlm2*dt(dlm/sd, df = df))/sd2
  enum3 <- (scorer[,"dmu"]*dt(drm/sd, df = df)/sd - scorel[,"dmu"]*dt(dlm/sd, df = df)/sd)
  enum4 <- drm2/sd2*dt(drm/sd, df = df)*(scorer[,"dsigma"] - 1/sd) - 
             dlm2/sd2*dt(dlm/sd, df = df)*(scorel[,"dsigma"] - 1/sd)
  enum5 <- (scorer[,"dsigma"]*dt(drm/sd, df = df)/sd - scorel[,"dsigma"]*dt(dlm/sd, df = df)/sd)
  
  

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

