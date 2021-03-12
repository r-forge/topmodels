---
title: |
  Probability Distribution Forecasts: Learning with Random Forests and Graphical Assessment 
author:
- name: Moritz N. Lang
  affiliation: Faculty of Economics and Statistics, Universität Innsbruck
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

**Webpages**: R-Forge repositories for the *R* packages [**disttree**](https://R-Forge.R-project.org/projects/partykit/pkg/disttree/) and [**topmodels**](https://R-Forge.R-project.org/projects/topmodels/pkg/topmodels/).

Forecasts in terms of entire probability distributions (often called
"probabilistic forecasts" for short) -- as opposed to predictions of only the
mean of these distributions -- are of prime importance in many different
disciplines from natural sciences to social sciences and beyond. Hence,
distributional regression models have been receiving increasing interest over
the last decade. Here, we make contributions to two common challenges in
distributional regression modeling:

1. Obtaining sufficiently flexible regression models that can capture complex
   patterns in a data-driven way.
2. Assessing the goodness-of-fit of distributional models both in-sample and
   out-of-sample using visualizations that bring out potential deficits of these
   models.

Regarding challenge 1, we present the *R* package **disttree**
[@disttree], that implements distributional trees and forests
[@Schlosser+Hothorn+Stauffer+Zeileis:2019]. These blend the recursive
partitioning strategy of classical regression trees and random forests with
distributional modeling. The resulting tree-based models can capture nonlinear
effects and interactions and automatically select the relevant covariates that
determine differences in the underlying distributional parameters.

For graphically evaluating the goodness-of-fit of the resulting probabilistic
forecasts (challenge 2), the *R* package **topmodels**
[@topmodels] is introduced, providing extensible
probabilistic forecasting infrastructure and corresponding diagnostic graphics
such as Q-Q plots of randomized residuals, PIT (probability integral transform)
histograms, reliability diagrams, and rootograms. In addition to distributional
trees and forests other models can be plugged into these displays, which can be
rendered both in base *R* graphics and **ggplot2** [@Wickham:2016].

## References
