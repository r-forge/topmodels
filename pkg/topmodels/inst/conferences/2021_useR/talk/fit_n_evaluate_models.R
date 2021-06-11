# -------------------------------------------------------------------
# - NAME:   fit_n_evaluate_models.R
# - AUTHOR: Moritz N. Lang
# - DATE:   2021-06-10
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2021-06-11 on thinkmoritz
# -------------------------------------------------------------------

## clean up
rm(list = objects())

# -------------------------------------------------------------------
# PRELIMINARIES
# -------------------------------------------------------------------
library("topmodels")
library("crch")
library("disttree")
library("scoringRules")
library("ggplot2")

# -------------------------------------------------------------------
# LOAD DATA
# -------------------------------------------------------------------
data("RainAxams")

# learning data: 24 years (1985 - 2008, both inlcuded)
# testing data: 4 successive years (2009, 2010, 2011, 2012)
train <- RainAxams[RainAxams$year < 2009,]
test <- RainAxams[RainAxams$year %in% c(2009, 2010, 2011, 2012),]

# -------------------------------------------------------------------
# FIT MODELS
# -------------------------------------------------------------------
# tree and forest formula
dt_formula <- df_formula <- formula(
  robs ~ tppow_mean + tppow_sprd + tppow_min + tppow_max +
  tppow_mean0612 + tppow_mean1218 + tppow_mean1824 + tppow_mean2430 +
  tppow_sprd0612 + tppow_sprd1218 + tppow_sprd1824 + tppow_sprd2430 +
  capepow_mean + capepow_sprd + capepow_min + capepow_max +
  capepow_mean0612 + capepow_mean1218 + capepow_mean1224 + capepow_mean1230 +
  capepow_sprd0612 + capepow_sprd1218 + capepow_sprd1224 + capepow_sprd1230 +
  dswrf_mean_mean + dswrf_mean_max +
  dswrf_sprd_mean + dswrf_sprd_max +
  msl_mean_mean + msl_mean_min + msl_mean_max +
  msl_sprd_mean + msl_sprd_min + msl_sprd_max +
  pwat_mean_mean + pwat_mean_min + pwat_mean_max +
  pwat_sprd_mean + pwat_sprd_min + pwat_sprd_max +
  tmax_mean_mean + tmax_mean_min + tmax_mean_max +
  tmax_sprd_mean + tmax_sprd_min + tmax_sprd_max +
  tcolc_mean_mean + tcolc_mean_min + tcolc_mean_max +
  tcolc_sprd_mean + tcolc_sprd_min + tcolc_sprd_max +
  t500_mean_mean + t500_mean_min + t500_mean_max +
  t700_mean_mean + t700_mean_min + t700_mean_max +
  t850_mean_mean + t850_mean_min + t850_mean_max +
  t500_sprd_mean + t500_sprd_min + t500_sprd_max +
  t700_sprd_mean + t700_sprd_min + t700_sprd_max +
  t850_sprd_mean + t850_sprd_min + t850_sprd_max +
  tdiff500850_mean + tdiff500850_min + tdiff500850_max +
  tdiff700850_mean + tdiff700850_min + tdiff700850_max +
  tdiff500700_mean + tdiff500700_min + tdiff500700_max +
  msl_diff
)

# baseline crch formula
formula_crch1 <- formula(robs ~ tppow_mean | log(tppow_sprd + 0.001))

# heuristic crch formula
formula_crch2 <- formula(
  robs ~ tppow_mean + tppow_mean1218 * capepow_mean1218 + tppow_max +
  dswrf_mean_mean + tcolc_mean_mean + msl_diff + pwat_mean_mean + tdiff500850_mean |
  tppow_sprd + tppow_sprd1218 * capepow_mean1218 + dswrf_sprd_mean + tcolc_sprd_mean +
  tdiff500850_mean
)

## fit all models
models <- list()

m_dt1 <- disttree(dt_formula, data = train, family = dist_list_cens_normal,
                 censtype = "left", censpoint = 0,
                 control = disttree_control(teststat = "quad", testtype = "Univ",
                                            type.tree = "ctree",
                                            intersplit = TRUE,
                                            mincriterion = 0, minsplit = 50,
                                            minbucket = 20))

m_df1 <- distforest(df_formula, data = train, family = dist_list_cens_normal,
                 ntree = 100, censtype = "left", censpoint = 0,
                 control = disttree_control(teststat = "quad", testtype = "Univ",
                                            type.tree = "ctree",
                                            intersplit = TRUE,
                                            mincriterion = 0, minsplit = 50,
                                            minbucket = 20), mtry = 27)
m_crch1 <- crch(formula_crch1, censored = TRUE, left = 0, data = train)
m_crch2 <- crch(formula_crch2, left = 0, control = crch.control(method = "boosting"), data = train)

m_lm1 <- lm(robs ~ tppow_mean, data = train)

## predict all models
pred <- lapply(list("tree" = m_dt1, "forest" = m_df1, "tobit1" = m_crch1, "tobit2" = m_crch2, 
  "lm" = m_lm1), 
  function(x) procast(x, type = "parameter", newdata = test))

# -------------------------------------------------------------------
# SCORING (OUT-OF-SAMPLE)
# -------------------------------------------------------------------
scoring_crps <- list()
scoring_crps <- lapply(pred, function(x) {
  crps_cnorm(
    test$robs, 
    location = x[, 1], 
    scale = x[, 2],
    lower = 0, 
    upper = Inf
    )
  }
)

scores <- list("crps" = data.frame(do.call("cbind", scoring_crps)))

scores_boot <- lapply(scores, function(x) {
  kboot <- 500
  x_boot <- data.frame(matrix(NA, ncol = ncol(x), nrow = kboot, dimnames = list(NULL, names(x))))
  for (i in 1:kboot) {
     s <- sample(1:nrow(x), nrow(x), replace = TRUE)
     x_boot[i,] <- apply(x[s, ], 2, mean)
  }
  x_boot
})

par(mfrow = c(1, 2))
boxplot(scores[["crps"]], main = "Raw CRPS")
boxplot(scores_boot[["crps"]], main = "Bootstrapped CRPS")

# -------------------------------------------------------------------
# GRAPHICAL EVALUATION (IN-SAMPLE)
# -------------------------------------------------------------------
## set theme for all plots
theme_set(theme_minimal())

## rootogram
root1_lm1 <- rootogram(m_lm1, plot = FALSE, breaks = -1:14 + 1e-12)
root1_df1 <- rootogram(m_df1, plot = FALSE, breaks = -1:14 + 1e-12)

autoplot(c("lm" = root1_lm1, "random_forest" = root1_df1), col = 1:2)

## reliagram
rel1_lm1 <- reliagram(m_lm1, thresholds = 0.1, plot = FALSE)
rel1_df1 <- reliagram(m_df1, thresholds = 0.1, plot = FALSE)

autoplot(c("lm" = rel1_lm1, "random_forest" = rel1_df1), 
  col = 1:2, fill = 1:2, single_graph = TRUE, legend = TRUE)

## pithist
pit1_lm1 <- pithist(m_lm1, plot = FALSE)
pit1_df1 <- pithist(m_df1, plot = FALSE)

autoplot(c("lm" = pit1_lm1, "random_forest" = pit1_df1))

## qqrplot
qq1_lm1 <- qqrplot(m_lm1, plot = FALSE)
qq1_df1 <- qqrplot(m_df1, plot = FALSE)

autoplot(c("lm" = qq1_lm1, "random_forest" = qq1_df1),
  col = 1:2, fill = 1:2, single_graph = TRUE, legend = TRUE)

## qqrplot
worm1_lm1 <- wormplot(m_lm1, plot = FALSE)
worm1_df1 <- wormplot(m_df1, plot = FALSE)

autoplot(c("lm" = worm1_lm1, "random_forest" = worm1_df1),
  col = 1:2, fill = 1:2, single_graph = TRUE, legend = TRUE)

# -------------------------------------------------------------------
# GRAPHICAL EVALUATION (OUT-OF-SAMPLE)
# -------------------------------------------------------------------
## set theme for all plots
theme_set(theme_minimal())

## rootogram
root2_lm1 <- rootogram(m_lm1, newdata = test, plot = FALSE, breaks = -1:14 + 1e-12)
root2_df1 <- rootogram(m_df1, newdata = test, plot = FALSE, breaks = -1:14 + 1e-12)

autoplot(c("lm" = root2_lm1, "random_forest" = root2_df1), col = 1:2)

## reliagram
rel2_lm1 <- reliagram(m_lm1, newdata = test, thresholds = 0.1, plot = FALSE)
rel2_df1 <- reliagram(m_df1, newdata = test, thresholds = 0.1, plot = FALSE)

autoplot(c("lm" = rel2_lm1, "random_forest" = rel2_df1), 
  col = 1:2, fill = 1:2, single_graph = TRUE, legend = TRUE)

## pithist
pit2_lm1 <- pithist(m_lm1, newdata = test, plot = FALSE)
pit2_df1 <- pithist(m_df1, newdata = test, plot = FALSE)

autoplot(c("lm" = pit2_lm1, "random_forest" = pit2_df1))

## qqrplot
qq2_lm1 <- qqrplot(m_lm1, newdata = test, plot = FALSE)
qq2_df1 <- qqrplot(m_df1, newdata = test, plot = FALSE)

autoplot(c("lm" = qq2_lm1, "random_forest" = qq2_df1),
  col = 1:2, fill = 1:2, single_graph = TRUE, legend = TRUE)

## qqrplot
worm2_lm1 <- wormplot(m_lm1, newdata = test, plot = FALSE)
worm2_df1 <- wormplot(m_df1, newdata = test, plot = FALSE)

autoplot(c("lm" = worm2_lm1, "random_forest" = worm2_df1),
  col = 1:2, fill = 1:2, single_graph = TRUE, legend = TRUE)

save(root1_lm1, root1_df1, rel1_lm1, rel1_df1, 
  pit1_lm1, pit1_df1, qq1_lm1, qq1_df1, worm1_lm1, worm1_df1, 
  root2_lm1, root2_df1, rel2_lm1, rel2_df1, 
  pit2_lm1, pit2_df1, qq2_lm1, qq2_df1, worm2_lm1, worm2_df1, 
  file = "Data/topmodels_axams.rda")
