# TODOs: Status quo and wishlist for `topmodels`
## Updated 2022-04-16 by Moritz N. Lang

## 1. General
### 1.1 Description file
* Fix description of package.

### 1.2 `pkgdown` homepage
* Update and extend vignettes.
* Write working papers.

### 1.3 Outlook
* Provide infrastucture for `scoringRules`: in-sample vs. out-of-sample  and aggregated vs. observation-wise contributions?
* Implement discretized log-score?
* S3 methods for topmodels: `plot()`, `autoplot()`?
* Employ `distributions3` in `procast()`.
* Support more model classes.
* Allow varying colors, etc. for additional graphical components in the `autoplot()`.

## 2. Package

### 2.1 Functions for forecasting: Summary

Function name | S3 classes supported | S3 classes planned | TODOs
--- | --- | --- | ---
`procast()` | `lm`, `crch`, (`disttree`) | `glm`, `countreg`, `betareg`, | yes
 | | (`glm`)  | `gam`, `gamlss`, `bamlss` | 
`procast_setup()` | none | none | none
`newresponse()` | `default` | no | yes
`pitresiduals()` | `default` | no | yes 

### 2.2 Functions for forecasting: TODOs
* In general
    * Correctly support `na.action`.
* `procast()`
    * Implement `countreg` and maybe `pscl` soon (update: tried, not that easy).
    * Which `type`s should be supported.
    * Improve S3 method for `crch`: Use functions instead of `eval()`.
    * Improve S3 method for `disttree`: Check compatability w/ forests and other families, and make usage of (not yet implemented) vectorized family functions.
    * In case `at = function`, it always returns a named vector. Do we want that?
    * Extend S3 classes (maybe w/ more flexible default method).
    * Implement `at = data.frame(y = response, x = model.matrix)` for scores.
    * Implement scores.
    * Employ `distributions3`!
* `newresponse()`
    * Check implementation of weights.
    * Rethink how and where to define `response_type` for `rootogram()`.
* `qresiduals()`
    * Fix `type = "quantile"`: Compare order statistics!!

### 2.3 Functions for graphical model assessment: Summary

Class | S3 classes | `c()` | `plot()` | `lines()`/`points()` | `autoplot()` | TODOs
--- | --- | --- | --- | --- | --- | ---
`rootogram` | `default`| yes | yes | - | yes | few
`reliagram` | `default` | yes | yes | yes | yes | few
`pithist` | `default` | yes | yes | yes | yes | few
`qqrplot` | `default` | yes | yes | yes | yes | few 
`wormplot` | `default` | yes | yes | yes | yes | few

### 2.4 Functions for graphical model assessment: TODOs
* In general
    * Check order of different layers.
    * Check names and consistency of all (plotting) arguments.
    * Rewrite examples w/i manuals.
    * Extend tests!
* In detail
    * Check again initilization for `glm` objects?
    * How should we handle values outside (truncation)/censoring points - these can lead to skewed PIT histograms.
    * Should a `+.object` S3-method be implemented?
    * Support more `trafo`s.
    * Change argument `trafo` to `scale`?
    * Allow varying graphical parameters in `ggplot2` facets. E.g., varying colors for CI intervals.
* `rootogram()`
    * Improve `breaks` and `width`
        * Get rid of argument `response_type`.
        * Extend breaks outside the range of observed frequencies in case expected probability mass exists.
        * Improve workaround in case of negative expecteded frequencies.
        * Improve workaround in case os point masses.
    * Implement empirical CIs based on random drawns from distribution?
    * Allow different `na.action`s?
    * Plot solely an expected line w/o points for many bars.
* `pithist()`
    * Implement proportional method (`type = "proportional"`).
    * Check if we really need two different CI computations.
    * Allow CI and ref to vary for non-equidstant breaks [done] -> improve for autoplot().
    * Implement CI and ref for trafo.
    * Fix ordering of graphics for `autoplot`.
* `qqrplot()`
    * Implement additional CI intervals (same for wormplot)? Compare R package qqplotr and Aldor-Noiman et al. (2013) [done].
* `wormplot()`
    * ciplot (Burren and Frederiks) for other `trafo`: Compare order statistics!
    * Fix CI polygon: CI polygon seems sometimes to be mirrored?! [done]
* `reliagram()`
    * Implement additionally to absolute histogram a frquency histogram (for base and ggplot2).
    * Write `geom`s and adapt `autoplot` accordingly.

## 3. Random notes
* Errors in countreg
    * Wrongly set argument `range = TRUE` leads to an error when randomized quantile residuals make no sense, which seems not be caught.
        * `m_gauss <- glm(dist ~ speed, data = cars, family = gaussian)`
        *  `qqrplot(m_gauss, range = TRUE)`
* Old notes by Achim/Jakob
    * Z: create R-Forge repos "topmodels" [done]
    * Z: svndump "crch" -> import to "topmodels" on R-Forge [done]
    * Z: package skeleton "topmodels" on R-Forge [done]
    * Z: procast() generic plus flexible procast_setup() [done]
    * JM: procast.crch() with new procast_setup() infrastructure.
    * JM: reliagram() prototype
    * JM: port pithist() from "countreg" + add type = "proportional" + Manu style
    * Z: port qqrplot() and rootogram() from "countreg" plus adaptations
    * JM/Z/.../CK/IK/.../SW/MS/GS/???: procast methods for betareg, countreg, glm, gamlss, mgcv::gam, ...
    * scoring rules for fitted model objects:
    * in-sample vs. out-of-sample / aggregated vs. observation-wise contributions
    * out-of-sample logLik()/logs(), crps(), ..., discretized log-score
* Other stuff
    * `distributions3`: Finish manual
    * Master student: Discretized log-score? Responde Alex.

## 4. Current notes (to be sorted):


