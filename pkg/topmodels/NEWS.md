# topmodels 0.3-0

* Added `crps.distribution` method for numeric approximation/numeric integration to calculate
  the CRPS for univariate distributions where no analytic solution is available in the
  `scoringRules` package.


# topmodels 0.2-0

* Version, presented at DAGStat 2022.

* Some conceptual changes in the generation of graphical evaluation tools for both base *R* and
  `ggplot2` style graphics.

* `autoplot()` builds now on newly written `geom_*()` and `stat_*()` functions. 


# topmodels 0.1-0

* First version, presented at useR! 2021.

* Diagnostic graphics for Q-Q plots of randomized residuals, PIT (probability integral transform)
  histograms, reliability diagrams, wormplots, and rootograms. All graphical evaluations can be
  rendered both in base *R* graphics and `ggplot2`.

* Basic probabilistic forecasting infrastructure for `lm`, `crch`, `disttree` and `glm` model
  classes. Not all families and forecasting types are fully supported yet.
