library("zoo")
source("~/Documents/meteoR/crch/R/crch.glmnet.R")
source("~/Documents/meteoR/crch/R/cnorm.R")

source("~/Documents/meteoR/varselect/results.R")
dyn.load("~/Documents/meteoR/crch/src/crch.glmnetold.so")



#load("~/Documents/meteoR/varselect/data/Tmax11036_66.rda")
#data <- data[, - grep("sd_", names(data))]
#anom <- anomalies(data)

load("~/Documents/meteoR/crch/inst/EGU2017/anom.rda")


boost <- crch(obs ~ . | ., anom$data, method = "boosting", maxit = 100, mstop = "cv")

lasso <- crch(obs~.|., anom$data, control = crch.glmnet(reltol = 1E-3, maxit = 100, lambda.min.ratio = 0.01))

plot(lasso)


ngr <- crch(obs~t2m_dmin_mean|t2m_dmin_sd, train)


predboost <- predict(boost, newdata = test, type = "response")
predlasso <- predict(lasso, newdata = test, type = "response")
predngr <- predict(ngr, newdata = test)

err <- as.data.frame(cbind((predngr - test$obs)^2, (predboost - test$obs)^2, (predlasso - test$obs)^2))

library("scoringRules")
err <- as.data.frame(cbind(
NGR = crps(y=as.numeric(test$obs), mean=predngr, sd=predict(ngr, newdata = test, type = "scale"), family = "normal"),
boost = crps(y=as.numeric(test$obs), mean=predboost, sd=predict(boost, newdata = test, type = "scale"), family = "normal"),
LASSO = crps(y=as.numeric(test$obs), mean=predlasso, sd=predict(lasso, newdata = test, type = "scale"), family = "normal")))






booterr <- NULL
for(b in 1:250) {
  bootind <- sample(nrow(err), replace = TRUE)
  booterr <- rbind(booterr, sqrt(colMeans(err[bootind,])))
}
