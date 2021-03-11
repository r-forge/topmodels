---
title: |
  Probability Distribution Forecasts: Learning with Random Forests and Graphical Assessment 
author:
- name: Moritz N. Lang
  affiliation: Department of Statistics, Faculty of Economics and Statistics, Universität Innsbruck
  email: Moritz.Lang@uibk.ac.at
  number: 1
- name: Reto Stauffer
  affiliation: Digital Science Center, Universität Innsbruck
  number: 1,2
- name: Achim Zeileis
  number: 1 
output: pdf_document
bibliography: ref.bib
---

**Keywords**: Distributional regression, probabilistic forecasts, regression trees, random forests, graphical model assessment.

**Webpages**: R-Forge repositories for the *R* packages [`disttree`](https://R-Forge.R-project.org/projects/partykit/pkg/disttree/) and [`topmodels`](https://R-Forge.R-project.org/projects/topmodels/pkg/topmodels/).

Forecasts in terms of entire probability distributions (often called "probabilistic
forecasts" for short) - as opposed to predictions of only the mean of these
distributions - are of prime importance in many different disciplines
from natural sciences to social sciences and beyond. Hence, distributional
regression models have been receiving increasing interest over the last decade.
Here, we make contributions to two common challenges in distributional
regression modeling:

1. Obtaining sufficiently flexible regression models that can capture complex
patterns in a data-driven way.
2. Assessing the goodness-of-fit of distributional models both in-sample and
out-of-sample using visualizations that bring out potential deficits of these
models

Distributional trees and forests [@Schlosser+Hothorn+Stauffer+Zeileis:2019]
blend the recursive partitioning strategy of decision trees with the concept of
distributional (regression) modeling. The resulting tree-based models can capture non-linear 
effects and interactions, and automatically select the relevant
covariates that are associated with changes in the respective
distribution parameters. To evaluate the goodness-of-fit of the distributional 
random forest models used in this study, we further present a toolbox that
provides an unifying infrastructure to obtain
predictions of probabilities, densities, scores, and Hessian values for
probabilistic models. The unifying prediction method provides users with many
different graphical evaluation tools, such as reliability diagrams, PIT
histograms, rootograms [@Kleiber+Zeileis:2016], and randomized Q-Q plots
[@Dunn+Smyth:1996].

The distributional random forests are provided in the *R*
package **disttree**, available on R-Forge [@Schlosser+Lang+Zeileis:2019]. The
unifying toolbox for inference and forecasting of probabilistic
(distributional) models is provided in the *R* package **topmodels**, also
available on R-Forge [@Zeileis+Kleiber+Kosmidis:2018]. The package includes the
generic function `procast` to compute various types of predictions (e.g.,
density functions, probabilities, quantiles), which so far supports the model
classes `lm`, `crch`, and `disttree`. For these, **topmodels** provides
routines to easily graphically assess and compare different probabilistic
models and model types using `ggplot2` [@Wickham:2016] and base *R* graphics. 


## References
