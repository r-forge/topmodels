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

---

As a consequence of the growing importance of probabilistic predictions in
various application fields due to a necessary functional risk management and
strategy, there is an increasing demand for appropriate probabilistic model
evaluation. Besides proper probabilistic scores [@Gneiting+Raftery:2007], which
evaluate not only the expectation but the entire predictive distribution,
graphical assessment methods are particularly advantageous to diagnose possible
model misspecification problems.  

Probabilistic predictions are often based on distributional regression models,
whereby the computation of predictive distributions, probabilities, and
quantiles is generally dependent on the *R* package being used.  Therefore,
routines to graphically evaluate probabilistic models are not always available
and if so then only for specific types of models and distributions provided by
the corresponding package. An easy to use unified infrastructure to graphical
assess and compare different probabilistic model types does not yet exist.

Using distribution based random forest models
[@Schlosser+Hothorn+Stauffer+Zeileis:2019] as an example, we present a toolbox
providing an unifiying infrastructure to gain predictions of probabilities,
densities, scores, and Hessian values for probabilistic models. The unifying
prediction method provides users with many different graphical evaluation
tools, such as reliability diagrams, PIT histograms, rootograms
[@Kleiber+Zeileis:2016], and randomized Q-Q plots [@Dunn+Smyth:1996].

The distributional random forests used for illustration are provided in the *R*
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
