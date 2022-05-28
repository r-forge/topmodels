# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `topmodels.R`
# --------------------------------------------------------------------
## FIXME: (ML) Implement more systematic testing

# --------------------------------------------------------------------
# Run some tests 
# --------------------------------------------------------------------
suppressPackageStartupMessages(require("crch"))
suppressPackageStartupMessages(require("ggplot2"))

data("CrabSatellites", package = "countreg")
CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]

m1 <- lm(dist ~ speed, data = cars)
m2 <- crch(dist ~ speed | speed, left = 3, data = cars)
m3 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
m4 <- disttree(dist ~ speed, data = cars)
m5 <- distforest(dist ~ speed, data = cars)

# FIXME: (ML) These all fail as for `do.call()` w/i `topmodels()` model object is not found
expect_silent( 
  rootogram(m1, plot = "base")
)

expect_silent( 
  topmodels(m1, plot = "base", col = "black", confint_level = 0.5, single_page = TRUE)
)

## FIXME: (ML) First line includes `Inf`, so there is a warning in `ggplot2`: Why is `Inf` there? Better to remove NAs with `na.rm = TRUE` in ggplot2?
expect_warning( 
  topmodels(m2, plot = "ggplot2")
)

expect_silent( 
  topmodels(m3, plot = "ggplot2", colour = "black", nsim = 10L)
)

expect_silent( 
  topmodels(m4, plot = "base")
)

expect_silent( 
  topmodels(m5, plot = "ggplot2", colour = "black", nsim = 10L)
)

