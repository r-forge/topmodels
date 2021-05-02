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
`procast()` | `lm`, `crch`, (`disttree`) | `glm`, `countreg`, `betareg`, | many
 | | (`glm`)  | `gam`, `gamlss`, `bamlss` | 
`procast_setup()` | none | none | none
`newresponse()` | `default` | no | few
`qresiduals()` | `default` | no | many 

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
`reliagram` | `default` | yes | yes | yes | no | few
`pithist` | `default` | yes | yes | yes | yes | few
`qqrplot` | `default` | yes | yes | yes | no | many 
`rootogram` | `default`| yes | yes | no | yes | few
`wormplot` | no | no | no | no | no | setup

### 2.4 Functions for graphical model assessment: TODOs
* `pithist()`
    * Streamline code and optional plotting arguments for the user
    * Implement `type = "proportional"`
    * Check again if we need two different CI computations
* `qqrplot()`
    * Setup as a generic function and include S3 methods. 
* `reliagram()`
    * Improve default method, get rid of `verification` package?
    * Implement generic functions
    * Get newest fancy version of Reto
* `rootogram()`
    * Move away from `countreg` default: Make argument `object` mandatory
    * Merge `rootogram_procast()` and `rootogram_glm`
* `wormplot()`
    * Start from scratch: Needs to be setup

## Current question/problems
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

