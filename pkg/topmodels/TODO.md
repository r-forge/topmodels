# TODOs: Status quo and wishlist for `topmodels`

## 1. General
### 1.1 Description file
* Fix naming and description of package [**done**]
* Check authors' contacts [**done**]

### 1.2 Outlook
* Provide `geom_()` for `ggplot2` instead/additionally to `autoplot()`
* Provide infrastucture for `scoringRules`: in-sample vs. out-of-sample  and aggregated vs. observation-wise contributions
* Implement discretized log-score

## 2. Functions

### 2.1 Functions for forecasting: Summary

Function name | S3 classes supported | S3 classes planned | TODOs
--- | --- | --- | ---
`procast()` | `lm`, `crch`, (`disttree`) | `glm`, `countreg`, `betareg`, | yes
 | | (`glm`)  | `gam`, `gamlss`, `bamlss` | 
`procast_setup()` | none | none | none
`newresponse()` | `default` | no | yes
`qresiduals()` | `default` | no | yes 

### 2.2 Functions for forecasting: TODOs
* `procast()` 
    * Improve S3 method for `crch`: Use functions instead of `eval()` and `disttree`
    * Improve S3 method for `disttree`: Check compatability w/ forests and other families, and make usage of (not yet implemented) vectorized family functions
    * Extend S3 classes (maybe w/ more flexible default method)
    * Implement `at = data.frame(y = response, x = model.matrix)` for scores
    * Implement scores.
* `newresponse()`
    * Implement (correctly) weights
    * Use `expand.model.frame()` 
    * How to handle `glm` objects
* `qresiduals()`
    * Fix `type = "quantile"`: Probably needs to be done on the transformed scale
    * Check if correctly working for discrete variables (no?!)

### 2.3 Functions for graphical model assessment: Summary

Class | S3 classes | `c()` | `plot()` | `lines()`/`points()` | `autoplot()` | TODOs
--- | --- | --- | --- | --- | --- | ---
`reliagram` | `default` | yes | yes | yes | yes | few
`pithist` | `default` | yes | yes | yes | yes | few
`qqrplot` | `default` | yes | yes | yes | yes | few 
`rootogram` | `default`| yes | yes | - | yes | few
`wormplot` | yes | yes | yes | yes | yes | few

### 2.4 Functions for graphical model assessment: TODOs
* `pithist()`
    * Streamline code and optional plotting arguments for the user [DONE]
    * Implement `type = "proportional"` 
    * Check again if we need two different CI computations
* `qqrplot()`
    * Setup as a generic function and include S3 methods [DONE]
* `reliagram()`
    * Improve default method, get rid of `verification` package? [DONE]
    * Implement generic functions [DONE]
    * Get newest fancy version of Reto
* `rootogram()`
    * Move away from `countreg` default: Make argument `object` mandatory [DONE]
    * Merge `rootogram_procast()` and `rootogram_glm` [DONE]
* `wormplot()`
    * Start from scratch: Needs to be setup [DONE]

# Update 2021-04-16
* How to do perform initilization correctly for `glm` objects?
* Should we increase the difference of the evaluation `at` in `qqresiduals()` to work for dicscrete distribution (binom, pois etc.) which leads to a wrong results for qnorm (compare tests). [SOLVED]
* How should we handle values outside (truncation)/censoring points. These lead to skewed PIT histograms.
* Is the graphical representation of the range in Q-Q plot correct?
* What types should `procast()` support: location vs. mean/expectation/response, scale vs. variance/dispersion?
* Why is the range not working in `qqrplot`.

## Errors in countreg
* Wrongly set argument `range = TRUE` leads to an error when randomized quantile residuals make no sense, which seems not be caught.
```
m_gauss <- glm(dist ~ speed, data = cars, family = gaussian)
qqrplot(m_gauss, range = TRUE)
```

# Update 2021-05-18

## Demonstration
```
devtools::load_all()
data("CrabSatellites", package = "countreg")
CrabSatellites2 <- CrabSatellites[CrabSatellites$satellites <= 1, ]
m1 <- glm(satellites ~ width + color, data = CrabSatellites, family = poisson)
m2 <- glm(satellites ~ width + color, data = CrabSatellites2, family = binomial)
m3 <- lm(dist ~ speed, data = cars)

rootogram(m1)
rootogram(m2)
rootogram(m3)

r1 <- reliagram(m1)
r1b <- reliagram(m1, minimum = 30)
r2 <- reliagram(m2)
r3 <- reliagram(m3)
plot(c(r1, r2))
plot(c(r1, r2), col = c(1, 2))
plot(c(r1, r2), confint = c(FALSE, TRUE))
plot(c(r1, r2), confint = c(1, 2), single_graph = TRUE)
lines(r3, col = 2)

autoplot(r1)
reliagram(m1, flavor = "tidyverse")

topmodels(m1, pages = 1)
topmodels(m1, pages = 1, fill = 2)
topmodels(m1, pages = 1, nsim = 10)
topmodels(m1, pages = 1, nsim = 10, flavor = "t")

```

## General Points
* Concept for plotting functions: `ref = TRUE`, `ref = c(col, col)`, `...`
* topmodels()
* flavor "base" vs. "tidyverse"
* autoplot: 
    * same arguments than base
    * same concept with groups as for base plotting functions 
    * `col = colour`, `pch = symbol`, `lty = linetype`, `cex = lwd = size`, ...
    * not all finished implemented, especially gropus can be difficult
 
## Current TODOs:
* Several bugfixes, especially in `autoplot()` with `single_graph = TRUE`
* `pithist()` proportional methodos
* legend()

## Outlook (until conference)
* Further streamlining incl. FIXMEs, clean up code and comments
* Check naming of in/out arguments
* Implement proper testing
* Adapt manuals (also for plotting funs)
* Maybe simple homepage
* `procast()`
* `newresponse()`
* `qresiduals()`

## Questions for Achim
* `na.action`: see `newresponse()`
* `response_type` + `newresponse()` in general
* rootogram: breaks no special options, `+.rootogram()`
* procast: weights, named data.frames as rval, scores (estfun), which models supported
* reliagram: extend to left/right border, minimum, new confint representation, `reliagram.crch()`
* topmodels: arguments compare `bamlss()`
* wormplot: `ciplot()`
