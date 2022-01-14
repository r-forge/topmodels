---
title: |
  Graphical Model Assessment of Probabilistic Forecasts
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

**Keywords**: Distributional regression, ensemble forecasts, graphical model assessment.


As a consequence of the growing importance of probabilistic predictions in
various application fields due to a necessary functional risk management and
strategy, there is an increasing demand for appropriate probabilistic model
evaluation. Besides proper scoring rules, which can evaluate not only the
expectation but the entire predictive distribution, graphical assessment
methods are particularly advantageous to diagnose possible model
misspecifications.

Probabilistic forecasts are often based on distributional regression models,
whereby the computation of predictive distributions, probabilities, and
quantiles is generally dependent on the software (package) being used.
Therefore, routines to graphically evaluate probabilistic models are not always
available and if so then only for specific types of models and distributions
provided by the corresponding package. An easy to use unified infrastructure to
graphical assess and compare different probabilistic model types does not yet
exist. Trying to fill that gap, we present a common conceptual framework
accompanied by a flexible and object-oriented software implementation in the
*R* package **topmodels** (https://topmodels.R-Forge.R-project.org/). 

The package includes visualizations for PIT (probability integral transform)
histograms, Q-Q (quantile-quantile) plots of (randomized) quantile residuals,
rootograms, reliability diagrams, and worm plots. All displays can be rendered
in base *R* as well as in **ggplot2** and provide different options for, e.g.,
computing confidence intervals, scaling or setting graphical parameters.
Using examples of post-processing precipitation ensemble forecasts, we further
discuss how all theses types of graphics can be compared to each other and
which types of displays are particularly useful for bringing out which types of
model deficiencies.


