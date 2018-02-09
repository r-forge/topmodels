library("crch")
library("zoo")


load("RainIbk2.rda")
data <- na.omit(RainIbk)
for(i in 1:ncol(data)) data[,i] <- (data[,i]- mean(data[,i]))/sd(data[,i])

source("../../R/crch.glmnet.R")
source("../../R/cnorm.R")


cvind <- sort(sample(1:10, replace = TRUE, nrow(data)))
pred <- data.frame(rain = data$rain, NR = 0, lasso = 0, boost = 0)

for(i in 1:10) {
  train <- data[cvind != i,]
  test <- data[cvind == i,]
  NR <- crch(rain~tp_ensmean|tp_enssd, train, left = 0)
  lasso <- crch(rain~.|., train, control = crch.glmnet(reltol = 1E-3, maxit = 100, lambda.min.ratio = 0.01, mstop = "bic"), left = 0)
  boost <- crch(rain ~ . | ., train, method = "boosting", maxit = 500, mstop = "bic", left = 0)

  pred[cvind == i,]$NR <- predict(NR, newdata = test, type = "density", at = pred$rain[cvind == i], log = TRUE)
  pred[cvind == i,]$boost <- predict(boost, newdata = test, type = "density", at = pred$rain[cvind == i], log = TRUE)
  pred[cvind == i,]$lasso <- predict(lasso, newdata = test, type = "density", at = pred$rain[cvind == i], log = TRUE)
}

boot <- NULL
for(i in 1:100) {
  bootind <- sample(nrow(data), replace = TRUE)
  boot <- rbind(boot, colSums(pred[bootind,-1]))
}

save(boot, file = "bootresults.rda")


t <- NULL
nvar <- seq(5, ncol(data), 10)
for(i in nvar) {
  lassotime <- system.time(crch(rain~.|., data[,1:i], control = crch.glmnet(reltol = 1E-2, maxit = 100, lambda.min.ratio = 0.01, mstop = "cv"), left = 0))
  boosttime <- system.time(crch(rain ~ . | ., data[,1:i], method = "boosting", maxit = 500, mstop = "cv", left = 0))
  t <- rbind(t, c(boosttime[3], lassotime[3]))
}
save(t, nvar, file = "timings.rda")
  


boost <- list()
pred <- NULL
files <- list.files("data/", pattern = "rainLondon*")
for(leadtime in seq(30, 174, 24)) {
  file <- paste0("rainLondon_", leadtime, ".rda")
  load(paste0("data/", file))
  data <- rainLondon
  #for(i in 2:ncol(data)) data[,i] <- (data[,i]- mean(data[,i], na.rm = TRUE))/sd(data[,i], na.rm = TRUE)
  data2 <- data 
  data <- na.omit(rainLondon)


  boost[[as.character(leadtime)]] <- crch(rain ~ . | ., data, method = "boosting", maxit = 500, mstop = "bic", left = 0)

  prob <- predict(boost[[as.character(leadtime)]], newdata = tail(data2, 1), type = "probability", at = 0)
  expe <- predict(boost[[as.character(leadtime)]], newdata = tail(data2, 1), type = "response")
  pred <- rbind(pred, data.frame(prob, expe))
}
pred <- zoo(pred, order.by = as.POSIXct(rownames(pred)))
save(pred, file = "londonpred.rda")


  
