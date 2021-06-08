# TODOs: Status quo and wishlist for `topmodels`

## 1. General
### 1.1 Description file
* Fix naming and description of package
* Check authors' contacts [**done**]

### 1.2 `pkgdown` homepage
* Set up homepage
* One vignette 'get started'
* One vignette 'technical specifications/framework'
* One vignette per each graphical evaluation method

### 1.3 Long term outlook
* Provide `geom_()` for `ggplot2` instead/additionally to `autoplot()`
* Provide infrastucture for `scoringRules`: in-sample vs. out-of-sample  and aggregated vs. observation-wise contributions
* Implement discretized log-score
* S3 methods for topmodels: plot, autoplot?

## 2. Package

### 2.1 Functions for forecasting: Summary

Function name | S3 classes supported | S3 classes planned | TODOs
--- | --- | --- | ---
`procast()` | `lm`, `crch`, (`disttree`) | `glm`, `countreg`, `betareg`, | yes
 | | (`glm`)  | `gam`, `gamlss`, `bamlss` | 
`procast_setup()` | none | none | none
`newresponse()` | `default` | no | yes
`qresiduals()` | `default` | no | yes 

### 2.2 Functions for forecasting: TODOs
* In general
    * Correctly support `na.action`
* `procast()` 
    * Which `type`s should be supported
    * Improve S3 method for `crch`: Use functions instead of `eval()`
    * Improve S3 method for `disttree`: Check compatability w/ forests and other families, and make usage of (not yet implemented) vectorized family functions
    * In case `at = function`, it always returns a named vector. Do we want that?
    * Extend S3 classes (maybe w/ more flexible default method)
    * Implement `at = data.frame(y = response, x = model.matrix)` for scores
    * Implement scores
* `newresponse()`
    * Check implementation of weights
    * Rethink how and where to define `response_type` for `rootogram()`
* `qresiduals()`
    * Fix `type = "quantile"`: Compare order statistics!!

### 2.3 Functions for graphical model assessment: Summary

Class | S3 classes | `c()` | `plot()` | `lines()`/`points()` | `autoplot()` | TODOs
--- | --- | --- | --- | --- | --- | ---
`reliagram` | `default` | yes | yes | yes | yes | few
`pithist` | `default` | yes | yes | yes | yes | few
`qqrplot` | `default` | yes | yes | yes | yes | few 
`rootogram` | `default`| yes | yes | - | yes | few
`wormplot` | `default` | yes | yes | yes | yes | few

### 2.4 Functions for graphical model assessment: TODOs
* In general
    * Check order of different layers
    * Check consistency of `xlim` and `ylim`
    * Check names and consistency of all (plotting) arguments
* In detail
    * Check again initilization for `glm` objects?
    * How should we handle values outside (truncation)/censoring points - these can lead to skewed PIT histograms.
* `reliagram()`
    * Implement additionally to absolute histogram a frquency histogram (for base and ggplot2).
* `pithist()`
    * Implement proportional method (`type = "proportional"`)
    * Check again if we need two different CI computations
* `wormplot()`
    * ciplot (Burren and Frederiks) for other `trafo`: Compare order statistics!! 

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

## 4. Current Questions
* Package
    * rootogram: breaks no special options, `+.rootogram()`
    * procast: named data.frames as rval, scores (estfun), support which types
* Check layers in plots
* Check naming of arguments
* Organization of manuals
* Help to setup pkgdown 
* Content of presentation use-R
    * Kurz distributional modeling motivieren
    * Idee mit Baeumen und Waeldern anreissen, verfuegbar in disttree
    * Beispiel auf Datensatz zeigen
    * Frage: Wie weiss man, ob das Modell ok funktioniert oder nicht?
    * Scoring rules aber zugehoerige Grafiken
    * An Beispiel herzeigen
    * Ausblick: Allgemeine Toolbox fuer generelle Modelle mit Base-Grafiken und ggplot2.

## 5. Next steps
* Streamline `procast()`
* Check naming of in/out arguments
* Adapt manuals (also for plotting funs)
* Implement simple homepage
* Some testing by working group (Reto)

