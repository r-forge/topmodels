---
title: |
  Visualizing Goodness-of-Fit of Probabilistic Regression Models 
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

**Keywords**: Graphical model assessment, distributional regression, probabilistic forecasting.

Forecasts in terms of entire probability distributions (often called
"probabilistic forecasts" for short) -- as opposed to predictions of only the
mean of these distributions -- are of prime importance in many different
disciplines from natural sciences to social sciences and beyond. Hence,
distributional regression models have been receiving increasing interest over
the last decade.  In order to evaluate the goodness of fit of these models and
to identify possible model deficits, besides proper scoring rules, graphical
methods are particularly beneficial.

Probabilistic predictions are often based on distributional regression models,
whereby the computation of predictive distributions, probabilities, and
quantiles is generally dependent on the *R* package being used. Therefore,
routines to graphically evaluate probabilistic models are not always available
and if so then only for specific types of models and distributions provided by
the corresponding package. An easy to use unified infrastructure to graphical
assess and compare different probabilistic model types does not yet exist.

The *R* package **topmodels** (https://topmodels.r-forge.r-project.org) aims to
fill that gap, providing extensible probabilistic forecasting infrastructure
and corresponding diagnostic graphics to evaluate both the marginal and
conditional calibration of such probabilistic models. The graphical displays
include Q-Q plots of randomized residuals, worm plots, PIT (probability
integral transform) histograms, reliability diagrams, and rootograms; while the
modular object-oriented implementation supports various model objects such as
`lm`, `glm`, `crch`, `disttree`, and more to come.  All graphical displays can
be rendered both in base *R* graphics and **ggplot2**.
