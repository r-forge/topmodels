## density
dtt <- function(x, mean = 0, sd = 1, df, left = -Inf, right = Inf, log = FALSE) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    df = as.numeric(df), left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("dtt", x, mean, sd, df, left, right, log))
  if(is.matrix(x)) {
    rval <- matrix(rval, ncol = ncol(x), nrow = nrow(x))
    colnames(rval) <- colnames(x)
    rownames(rval) <- rownames(x)
  }
  return(rval)
}


## distribution function
ptt <- function(q, mean = 0, sd = 1, df, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  input <- data.frame(q = as.numeric(q), mean = as.numeric(mean), sd = as.numeric(sd), 
    df = as.numeric(df), left = as.numeric(left), right = as.numeric(right))
  rval <- with(input, .Call("ptt", q, mean, sd, df, left, right, lower.tail, log.p))
  if(is.matrix(q)) {
    rval <- matrix(rval, ncol = ncol(q), nrow = nrow(q))
    colnames(rval) <- colnames(q)
    rownames(rval) <- rownames(q)
  }
  return(rval)
}

## quantiles
qtt <- function(p, mean = 0, sd = 1, df,  left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE) {
  if(log.p) p <- exp(p) 
  lower <- if(lower.tail) left else right
  upper <- if(lower.tail) right else left
  p <- pt((lower-mean)/sd, df = df) * (1 - p) + p*pt((upper - mean)/sd, df = df)
  rval <- qt(p, df = df, lower.tail = lower.tail)*sd + mean
  if(is.matrix(p)) {
    rval <- matrix(rval, ncol = ncol(p), nrow = nrow(p))
    colnames(rval) <- colnames(p)
    rownames(rval) <- rownames(p)
  }
  return(rval)
}

## random numbers
rtt <- function(n, mean = 0, sd = 1, df, left = -Inf, right = Inf) {
  qtt(runif(n), mean = mean, sd = sd, df = df, left = left, right = right)
}

## scores
stt <- function(x, mean = 0, sd = 1, df, left = -Inf, right = Inf, 
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    df = as.numeric(df), left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  
  for(w in which) {
    if(w == "mu")
      score2 <- with(input, .Call("stt_mu", x, mean, sd, df, left, right))
    if(w == "sigma")
      score2 <- with(input, .Call("stt_sigma", x, mean, sd, df, left, right))
    score <- cbind(score, score2)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

## Hessian
htt <- function(x, mean = 0, sd = 1, df, left = -Inf, right = Inf, 
  which = c("mu", "sigma")) {
  input <- data.frame(x = as.numeric(x), mean = as.numeric(mean), sd = as.numeric(sd), 
    df = as.numeric(df), left = as.numeric(left), right = as.numeric(right))
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  hess <- list()
  for(w in which) {       
    if(w == "mu")         
      hess[[w]] <- with(input, .Call("htt_mu", x, mean, sd, df, left, right))  
    if(w == "sigma")
      hess[[w]] <- with(input, .Call("htt_sigma", x, mean, sd, df, left, right))  
    if(w %in% c("mu.sigma", "sigma.mu"))
      hess[[w]] <- with(input, .Call("htt_musigma", x, mean, sd, df, left, right))  
  }

  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}

## Expectation
ett <- function(mean = 0, sd = 1, df, left = -Inf, right = Inf) {
  rmm <- (right-mean)/sd
  lmm <- (left-mean)/sd
  pncens <- pt(rmm, df = df)-pt(lmm, df = df)
  rval <- mean + sd*df*((dt(rmm, df = df) * (1+rmm^2/df)^(is.finite(rmm))-
    dt(lmm, df = df)*(1+lmm^2/df)^(is.finite(lmm)))/(1-df)) /pncens
  rval
}

