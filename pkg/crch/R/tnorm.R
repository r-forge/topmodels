## density
dtnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cdtnorm", x, mean, sd, left, right, log))
  if(is.matrix(x)) {
    rval <- matrix(rval, ncol = ncol(x), nrow = nrow(x))
    colnames(rval) <- colnames(x)
    rownames(rval) <- rownames(x)
  }
  return(rval)
}


## distribution function
ptnorm <- function(q, mean = 0, sd = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  input <- data.frame(q = as.numeric(q), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("cptnorm", q, mean, sd, left, right, lower.tail, log.p))
  if(is.matrix(q)) {
    rval <- matrix(rval, ncol = ncol(q), nrow = nrow(q))
    colnames(rval) <- colnames(q)
    rownames(rval) <- rownames(q)
  }
  return(rval)
}

## quantiles
qtnorm <- function(p, mean = 0, sd = 1, left = -Inf, right = Inf,
  lower.tail = TRUE, log.p = FALSE) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  upper <- if(lower.tail) right else left
  p <- pnorm((lower-mean)/sd, lower.tail = lower.tail) * (1 - p) + 
    p*pnorm((upper - mean)/sd, lower.tail = lower.tail)
  rval <- qnorm(p, lower.tail = lower.tail)*sd + mean
  if(is.matrix(p)) {
    rval <- matrix(rval, ncol = ncol(p), nrow = nrow(p))
    colnames(rval) <- colnames(p)
    rownames(rval) <- rownames(p)
  }
  return(rval)
}

## random numbers
rtnorm <- function(n, mean = 0, sd = 1, left = -Inf, right = Inf) {
  qtnorm(runif(n), mean, sd, left = left, right = right)
}



## scores
stnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf,
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  
  for(w in which) {
    if(w == "mu")
      score2 <- with(input, .Call("stnorm_mu", x, mean, sd, left, right))
    if(w == "sigma")
      score2 <- with(input, .Call("stnorm_sigma", x, mean, sd, left, right))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
htnorm <- function(x, mean = 0, sd = 1, left = -Inf, right = Inf,
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- with(input, .Call("htnorm_mu", x, mean, sd, left, right))  
    if(w == "sigma")
      hess[[w]] <- with(input, .Call("htnorm_sigma", x, mean, sd, left, right))  
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- with(input, .Call("htnorm_musigma", x, mean, sd, left, right))  
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}

.erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
.erfc <- function(x) 2 * pnorm(x * sqrt(2), lower.tail = FALSE)
.erfcx <- function(x) 2 * pnorm(x * sqrt(2), lower.tail = FALSE) * exp(x^2)
.F1 <- function(x, y) {
    delta <- exp(x^2 - y^2)
    fx <- is.finite(x)
    fy <- is.finite(y)
    sx <- sign(x)
    sy <- sign(y)
    ifelse(fx & !fy, sy / .erfcx(sy * x),
    ifelse(!fx & fy, sx / .erfcx(sx * y),
    ifelse(abs(x) > y & y >= 0, (exp(-y^2) - exp(-x^2)) / (.erf(x) - .erf(y)),
    ifelse(x < 0 & y < 0, (1 - delta) / (delta * .erfcx(-y) - .erfcx(-x)),
    ifelse(x > 0 & y > 0, (1 - delta) / (.erfcx(x) - delta * .erfcx(y)),
    (1 - delta) * exp(-x^2) / (.erf(y) - .erf(x)))))))
}
.F2 <- function(x, y) {
    delta <- exp(x^2 - y^2)
    fx <- is.finite(x)
    fy <- is.finite(y)
    sx <- sign(x)
    sy <- sign(y)
    ifelse(fx & !fy, sy * x / .erfcx(sy * x),
    ifelse(!fx & fy, sx * y / .erfcx(sx * y),
    ifelse(abs(x) > y & y >= 0, (y * exp(-y^2) - x * exp(-x^2)) / (.erf(x) - .erf(y)),
    ifelse(x < 0 & y < 0, (x - y * delta) / (delta * .erfcx(-y) - .erfcx(-x)),
    ifelse(x > 0 & y > 0, (x - y * delta) / (.erfcx(x) - delta * .erfcx(y)),
    (x - y * delta) * exp(-x^2) / (.erf(y) - .erf(x)))))))
}


## Using the expressions in
## https://github.com/cossio/TruncatedNormal.jl/blob/fc904152f2da11a257e3ccdd3e49ef118b81d437/notes/normal.pdf
## to avoid catastrophic cancellation

etnorm <- function (mean = 0, sd = 1, left = -Inf, right = Inf) {
    rmm <- (right - mean) / sd / sqrt(2)
    lmm <- (left - mean) / sd / sqrt(2)
    ifelse(rmm == Inf & lmm == -Inf, mean,
           mean + sqrt(2 / pi) * .F1(lmm, rmm) * sd)
}

sdtnorm <- function (mean = 0, sd = 1, left = -Inf, right = Inf) {
    rmm <- (right - mean) / sd / sqrt(2)
    lmm <- (left - mean) / sd / sqrt(2)
    ifelse(rmm == Inf & lmm == -Inf, sd,
           sd * sqrt(1 + 2 / sqrt(pi) * .F2(lmm, rmm) - 2 / pi * (.F1(lmm, rmm))^2))
}

