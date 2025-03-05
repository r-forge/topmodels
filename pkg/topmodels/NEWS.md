# topmodels 0.3-0

* Entire package now leverages `distributions3` for object-oriented computations
  on distributions fitted/predicted by various kinds of models.
  
* In particular, `procast()` first obtains `prodist()` (probability distribution)
  and then applies the standard methods for computing densities, probabilities,
  quantiles, and moments.
  
* Similarly, `proresiduals()` obtains the predicted distributions and compares
  with the `newresponse()` to obtain (randomized) quantile residuals by default.
  Alternatively, PIT residuals or Pearson residuals as well as raw response
  residuals are available. The function `proresiduals()` also replaces both
  `pitresiduals()` and `qresiduals()` which were provided by earler versions of
  `topmodels`.
  
* Via the same approach `proscore()` implements various kinds of scoring rules,
  in particular log-score (negative log-likelihood), (continuous) ranked
  probability score (CRPS), mean absolute error (MAE), mean squared error (MSE),
  and the Dawid-Sebastiani-Score (DSS). The standard log-likelihood (without
  sign change) is also available.
  
* For the CRPS one can either leverage the functions from the `scoringRules`
  package (if available) or the new `crps.distribution()` method for numeric
  approximation/numeric integration to calculate the CRPS for univariate
  distributions. This is also used when no analytic solution is available in the
  `scoringRules` package.

* The graphical functions `rootogram()`, `pithist()`, `qqrplot()`, `wormplot()`,
  and `reliagram()` are also switched to the new infrastructure based on
  `distributions3`, notably via `procast()` and `proresiduals()`.

* The pointwise and simultaneous confidence intervals in `rootogram()` now rely
  on the exact `PoissonBinomial()` distribution (now available in `distributions3`)
  rather than its binomial approximation.
  
* In addition to pointwise and simultaneous confidence intervals for `rootogram()`,
  `"tukey"` confidence intervals are now available (and the default) which simply correspond to
  limits of -1 and 1 for hanging or suspended rootograms. For other flavors of
  rootograms these limits are transformed correspondingly.

* For Q-Q residual plots in `qqrplot()` there are two methods for pointwise
  confidence intervals (`"beta"` and `"normal"`) and two for simultaneous intervals
  (`"ks"` and `"ell"`, the latter requiring package `qqconf`). Additionally,
  `"pointwise"` and `"simultaneous"` are simply aliases for the preferred
  corresponding methods `"beta"` and `"ell"`, respectively.

* New distribution/model interfaces were added first in `topmodels` but some
  subsequently moved to other packages: `GAMLSS()` is now in `gamlss.dist`,
  `BAMLSS()` is now in `bamlss`, and `Empirical()` is still in `topmodels` for
  now.

* New wrapper function `promodel()` that adds the class `"promodel"` (for
  probabilistic model) to an existing model object so that `predict()`
  dispatches to `procast()` and `residuals()` dispatches to `proresiduals()`.
  This facilitates using model functionality based on the standard `predict()`
  and `residuals()` methods like `marginaleffects`.


# topmodels 0.2-0

* New version, presented at DAGStat 2022 and at useR! 2022 (together with
  `distributions3`).

* Some conceptual changes in the generation of graphical evaluation tools for
  both base R and `ggplot2` style graphics.

* `autoplot()` builds now on newly written `geom_*()` and `stat_*()` functions.


# topmodels 0.1-0

* First version, presented at useR! 2021.

* Diagnostic graphics for Q-Q plots of randomized residuals, PIT (probability
  integral transform) histograms, reliability diagrams, wormplots, and
  rootograms. All graphical evaluations can be rendered both in base R graphics
  and `ggplot2`.

* Basic probabilistic forecasting infrastructure for `lm`, `crch`, `disttree`
  and `glm` model classes. Not all families and forecasting types are fully
  supported yet.
