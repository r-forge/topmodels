# crch 1.2-0

* Improved implementations of expectation and variance of censored and truncated
  normal distributions (contributed by Ioannis Kosmidis).

* Turned `vignette("crch", package = "crch")`  from Sweave into Quarto
  vignettes. Some improvements/updates in the text.

* New package web page (via `altdoc`/`quarto`) at
  <https://topmodels.R-Forge.R-project.org/crch/>

* Achim Zeileis takes over maintenance from Jakob W. Messner. Added Ioannis Kosmidis
  and Georg J. Mayr as contributors.


# crch 1.1-2

* Fixed name of first argument in `crps()` method to be consistent with the generic.


# crch 1.1-1

* Added methods for `is_discrete()` and `is_continuous()` for the new distributions3
  objects.

* Replaced deprecated C function `finite()` with `isfinite()` (again).


# crch 1.1-0

* Added support for [distributions3](https://CRAN.R-project.org/package=distributions3)
  workflows for censored and truncated normal, logistic, and Student's t
  distributions: `CensoredNormal()`, `TruncatedNormal()`, `CensoredLogistic()`,
  `TruncatedLogistic()`, `CensoredStudentsT()`, `TruncatedStudentsT()`.
  See the corresponding manual pages for examples illustrating how to
  work with the distributions in practice, for computing moments, probabilities,
  densities, simulating random values, etc.
  
* Added `prodist()` method for extracting the `distributions3` objects above
  from fitted `crch` objects, either in-sample or out-of-sample.
  
* Bug fix in the computation of the mean of censored or truncated logistic
  distributions with large (or infinite) censoring/truncation points.


# crch 1.0-4

* Replaced deprecated C function `finite()` with `isfinite()`.


# crch 1.0-3

* Updated contact information.


# crch 1.0-1

* Added argument `type` to `crch()` which can be set to `"crps"` for parameter
  estimation with minimum CRPS instead of maximum likelihood. 

* Added S3 method for `crps()` from `scoringRules` for `crch` objects.

* Improvements for the `predict()` method:

  - New types `"parameter"`, `"density"`, `"probability"`, and `"crps"`.
  - With `type = "response"` now the expected value and not the location 
    parameter is returned (not equal for censored and truncated 
    distributions). For better backward compatibility, the default type is
    set to `"location"`.

* Added `pit()`, `rootogram()`, and `simulate()` methods for `crch` objects.

* Changed argument names `mean` and `sd` to `location` and `scale` in logistic and
  Student's t distribution functions 

* Added new function `crch.stabsel()` for stability selection based on    
  `crch.boost.fit()`. Some S3 methods for the returned class `stabsel.crch`
  are also provided.


# crch 1.0-0

* New release accompanying the R Journal paper: [Heteroscedastic Censored 
  and Truncated Regression with crch](https://doi.org/10.32614/RJ-2016-012)
  by Messner, Mayr, and Zeileis. See also `citation("crch")`. 
  
* Added `estfun()` method for `crch` objects


# crch 0.9-2

* The `crch()` function now supports coefficient optimization by
  boosting to automatically select the most relevant input variables 
  in high-dimensional data settings. Extractor and plotting functions 
  for corresponding `crch.boost` objects are also available.

* Transferred functions to estimate density, distribution, score, and 
  Hessian matrices to C code to accelerate coefficient optimization.

* Added option to `crch()` to avoid computation of covariance matrix.

* Added `left` and `right` arguments to `predict.crch()` and 
  `predict.crch.boost()` to allow quantile predictions for non-constant 
  censoring or truncation points.


# crch 0.9-1

* Added `model.matrix()` and `model.frame()` methods for `crch` objects

* Bug fix in `predict.crch()`: In previous versions predictions for 
  models with other link functions than the log gave wrong results


# crch 0.9-0

* Added vignette to introduce the `crch()` function with some
  theoretical background and an illustrating example:
  `vignette("crch", package = "crch")`

* The `crch()` function now also supports truncated responses.
  Furthermore added a wrapper function `trch()` to fit truncated
  regression models.

* `crch()`: Analytical gradients and Hessian matrices are provided for most 
  models to speed up maximum likelihood optimization
  (not available for Student's t distribution with degrees of freedom
  estimation).  

* `crch()`: For the scale model a link function can now be specified
  (`"log"`, `"identity"`, or `"quadratic"`). In previous version only the log
  was supported.

* Added functions for probability density, cumulative distribution, 
  random numbers, and quantiles for censored and truncated normal, 
  logistic, and Student's t distributions.

* The `residuals()` method for `crch` objects now also provides quantile 
  residuals (Dunn and Smyth 1996).

* Added `update()` method for `crch` objects.


# crch 0.1-0

* First official release of the package on CRAN. See `citation("crch")`
  for the accompanying manuscripts. Note that the interface of both
  `crch()` and `hxlr()` is still under development and might change in
  future versions of the package.
