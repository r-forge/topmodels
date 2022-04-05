# ------------------------------------------------------------
# Creates the graphs needed for the abstract.
#
# The data set as well as the PDFs will be checked in to 
# the repository for reproduceability.
#
# Author: Reto, 2022-04-05
# ------------------------------------------------------------


rm(list = objects())
rdsfile <- "../data_IWSM_day1to3.rda"
stopifnot(file.exists(rdsfile))
load(rdsfile, verbose = TRUE)


# ------------------------------------------------------------
# Prepare datasets by appending the ensemble mean and ensemble
# standard deviation. The transformer is used to transform
# the data (tp_[0-9]) before calculating the ensemble mean
# and standard deviation. We are using
# - train_id:  identity
# - train_pow: power transformed (using 1/1.5; better than sqrt for this data set)
#
# - ensmean: mean(<transformed tp_[0-9]>)
# - enssd:   sd(<transformed tp_[0-9]>)
# - enswet:  Number of wet members in ensemble (a number 0 to 5); wet is identified
#            if the forecast is larger than transformer(<wet_threshold>)
# ------------------------------------------------------------
prepare_training_data <- function(x, transformer = identity, wet_treshold = 1) {
    stopifnot(is.data.frame(x))
    stopifnot(is.function(transformer))
    idx <- grep("^tp_[0-9]$", names(x))
    for (i in idx) x[, i] <- transformer(x[, i])
    x <- transform(x, observation = transformer(observation))
    # Appending ensemble mean and standard deviation
    x <- transform(x, ensmean = apply(x[, idx], 1, mean),
                      enssd   = apply(x[, idx], 1, sd),
                      enswet  = apply(x[, idx], 1, function(x) sum(x > transformer(wet_treshold))))
    # Appending additional information
    x <- transform(x, year = as.integer(format(x$datetime, "%Y")),
                      yday = as.integer(format(x$datetime, "%j")))
    # Remoe
    return(x)
}

train_id  <- prepare_training_data(data, identity, 1)
train_pow <- prepare_training_data(data, function(x) x^(1/1.5), 1)


# ------------------------------------------------------------
# Estimating the models we would like to show
# ------------------------------------------------------------

library("crch")

# Homoscedastic non-censored Gaussian model (identity)
m_lm_id <- lm(observation ~ ensmean, data = train_id)

# Heteroscedastic censored Gaussian model (identity)
m_chgauss_id <- crch(observation ~ ensmean | log(enssd), data = train_id, left = 0, dist = "gaussian")

# Heteroscedastic censored logistic response model (power transformed data)
m_chlogis_pow <- crch(observation ~ ensmean | log(enssd), data = train_pow, left = 0, dist = "logistic")

models <- list(m_lm_id = m_lm_id, m_chgauss_id = m_chgauss_id, m_chlogis_pow =  m_chlogis_pow)
sapply(models, logLik)
sapply(models, BIC)


# ------------------------------------------------------------
# Model analysis
# ------------------------------------------------------------
library("topmodels")
library("ggplot2")

###### This works somewhat better but is less well calibrated
#####f <- formula(observation ~ ensmean*enswet + sin(yday/365*2*pi) + cos(yday/365*2*pi) | log(ensmean)*enswet + sin(yday/365*2*pi) + cos(yday/365*2*pi))
#####xx <- crch(f, data = train_pow, left = 0, dist = "logistic")
#####summary(xx)
#####logLik(xx)
#####AIC(xx, m_chlogis_pow)
#####BIC(xx, m_chlogis_pow)
#####autoplot(c(pithist(xx, plot = F), pithist(m_chlogis_pow, plot = F)), single_graph = TRUE, col = 2:3)
#####autoplot(c(rootogram(xx, plot = F), rootogram(m_chlogis_pow, plot = F)), single_graph = TRUE, col = 2:3)
#####autoplot(c(wormplot(xx, plot = F), wormplot(m_chlogis_pow, plot = F)), single_graph = TRUE, col = 2:3)
#####logLik(m_chlogis_pow)


###############
# Rootogram (hanging)
library("gridExtra")
bk <- 20
g1 <- autoplot(rootogram(m_lm_id, breaks = bk, plot = FALSE)) +
        ggtitle("Heteroscedastic Gaussian")
g2 <- autoplot(rootogram(m_chgauss_id, breaks = bk, plot = FALSE)) +
        ggtitle("Homoscedastic Censored Logistic")
ggsave(file = "Stauffer-rootograms.pdf", grid.arrange(g1, g2, ncol = 2),
       width = 7, height = 2.5)



###############
# QQ-R-Plot plus wormplot
models <- list(m_lm_id = m_lm_id, m_chgauss_id = m_chgauss_id, m_chlogis_pow = m_chlogis_pow)
g1 <- autoplot(do.call(c, lapply(models, function(m) qqrplot(m, plot = FALSE))),
               single_graph = TRUE, col = seq_along(models) + 1) +
        ggtitle("Q-Q Residuals")

g2 <- autoplot(do.call(c, lapply(models, function(m) wormplot(m, plot = FALSE))),
               single_graph = TRUE, col = seq_along(models) + 1) +
        ggtitle("Wormplot") + ylim(-1, 1)

ggsave(file = "Stauffer-qqresiduals.pdf", grid.arrange(g1, g2, ncol = 2, widths = c(1, 1.5)),
       width = 7, height = 2.5)




















