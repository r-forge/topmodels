library("crch")
library("STAGE")
library("zoo")
source("~/Dropbox/meteoR/crch/R/crch.glmnet.R")
source("~/Dropbox/meteoR/crch/R/cnorm.R")

source("~/Dropbox/meteoR/varselect/results.R")
dyn.load("~/Dropbox/meteoR/crch/src/crch.glmnet.so")


## create RRR10444_*.rda files
if(FALSE) {
station <- 10444
obs <- sybase("synop", station, parameter = "RRR", begin = as.POSIXct("2010-01-01"))
obs2 <- obs

for(lt in seq(30, 138,12)) {
  nwp2 <- NULL
  for(lt2 in seq(-9,0,3)) {
    nwp2 <- rbind(nwp2, dmo("ecmwf", station, step = lt + lt2, 
      deaccumulate = 3, begin = as.POSIXct("2011-01-01"), end = as.POSIXct("2015-12-31")))
  } 
  nwp2$tp[nwp2$tp<0] <- 0
#  nwp2$tp <- sqrt(nwp2$tp)
  for(i in 2:ncol(nwp2)) {
    nwp2 <- cbind(nwp2, nwp2$tp * nwp2[,i])
    names(nwp2)[ncol(nwp2)] <- paste0("tp", names(nwp2)[i])
  }
  max <- !lt %in% seq(30, 126, 24)
  nwp2 <- aggregate(nwp2, as.Date(index(nwp2)+!max*5*60*60), sum)
#  nwp2$date <- as.zoo(as.POSIXct(index(nwp2))+6*60*60, order.by = index(nwp2))
  index(nwp2) <- as.POSIXct(index(nwp2)) + (max*18 + (!max)*6)*60*60

  data <- merge(obs=obs, nwp2, all = FALSE)
  names(data)[1] <- "obs"
  data <- data[,colSums(is.na(data)) < 0.15*nrow(data)]
  data <- na.omit(data)
  save(data, file = paste0("RRR", station, "_", lt, ".rda"))
}
}

## boosting model to predict Goettingen precipitation
boost <- anom2 <- list(NULL)
for(lt in seq(30, 138,12)) {
  load(paste0("RRR10444_", lt, ".rda"))
  anom <- anomalies(data, simple = TRUE)
  left <- -anom$mall[,1]/anom$sall[,1]
  boost[[(lt-30)/12+1]] <- crch(obs~.|., anom$data, control = crch.boost(maxit = 400, mstop = "bic"), left = left)
  anom2[[(lt-30)/12+1]] <- anom
}

## get latest prediction
nwp <- dmo.latest(station = 10444,model = "ecmwf", deaccumulate = 3, data = TRUE, runhour = 0)
save(nwp, file = "latest.rda")

## fit Innsbruck boosting and LASSO
load("RRR11120_30.rda")
data$obs[data$obs<0] <- 0

## subset needed to avoid LASSO paths to diverge
data <- subset(data, select = c('obs','tp_mean','r700_mean','sund_mean','w700_mean','r1000_mean',
'str_mean','ttr_mean','tpz850_mean','tpz1000_mean','tpu700_mean','tpv700_mean','persistence',
'tpcape_mean','tpt2m_mean','tpu850_mean','tpv850_mean','tpw850_mean','tp_sd','d2m_sd','tpcape_sd',
'tpu700_sd','tpv700_sd','t500_sd','q1000_sd','v850_sd','fg10m_mean','fg10m_sd','tpvo700_mean', 'tpq700_mean','tptp_mean','t1000_sd','bld_sd','t700_mean','cape_mean','z500_mean','q700_mean'))
anom <- anomalies(data, simple = TRUE)


data2 <- as.data.frame(anom$data)
left <- -anom$mall[,1]/anom$sall[,1]

formula <- "obs ~tp_mean + r700_mean + sund_mean + w700_mean + r1000_mean + str_mean + ttr_mean +tpz850_mean + tpz1000_mean + tpu700_mean + tpv700_mean + persistence + tpcape_mean + tpt2m_mean + tpu850_mean + tpv850_mean + tpw850_mean+tp_sd + d2m_sd + tpcape_sd + tpu700_sd + tpv700_sd + t500_sd + q1000_sd + v850_sd + fg10m_mean + fg10m_sd|tp_sd + t700_mean + cape_mean + z500_mean + q700_mean + tpu700_sd + tpv700_sd + tpvo700_mean + tpq700_mean + tpu850_mean + tpv850_mean + tptp_mean + tp_mean + t1000_sd + bld_sd"

## fit models
lasso <- crch(obs~.|., data2, control = crch.glmnet(reltol = 1E-3, maxit = 10000, mstop = "cv"), left = left)
boost <- crch(obs~.|., data2, control = crch.boost(maxit = 3000, nu = 0.1, mstop = "cv"), left = left)
save(lasso, boost, file = "fitted.rda")

## bootstrap out of bag predictive performance
bootfun <- function(i) {
  train <- sample(nrow(data), replace = TRUE)
  test <- which(!c(1:nrow(data)) %in% unique(train))
  lasso <- crch(formula, data2[train,], control = crch.glmnet(reltol = 1E-3, maxit = 10000, mstop = "cv"), left = left[train])
  boost <- crch(formula, data2[train,], control = crch.boost(maxit = 1500, nu = 0.15, mstop = "cv"), left = left[train])
  print(boost$mstopopt)
  mod <- crch(obs~tp_mean|tp_sd, data2[train,], left = left[train])
  fc_mu <- predict(lasso, newdata = anom$data[test,])
  fc_sigma <- predict(lasso, newdata = anom$data[test,], type = "scale")
  ll_lasso <- sum(dcnorm(anom$data$obs[test], fc_mu, fc_sigma, left = left[test], log = TRUE))
  fc_mu <- predict(boost, newdata = anom$data[test,])
  fc_sigma <- predict(boost, newdata = anom$data[test,], type = "scale")
  ll_boost <- sum(dcnorm(anom$data$obs[test], fc_mu, fc_sigma, left = left[test], log = TRUE))
  fc_mu <- predict(mod, newdata = anom$data[test,])
  fc_sigma <- predict(mod, newdata = anom$data[test,], type = "scale")
  ll_mod <- sum(dcnorm(anom$data$obs[test], fc_mu, fc_sigma, left = left[test], log = TRUE))
  c(ll_mod, ll_boost, ll_lasso)
}

boot <- NULL
for(i in 1:100) boot <- rbind(boot, bootfun(i))
save(boot, file = "bootresults.rda")

