---
title: |
  Visualizing Goodness of Fit of Probabilistic Regression Models 
author:
- name: Moritz N. Lang
  affiliation: Department of Statistics, Universität Innsbruck
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

Modeling and forecasting in terms of entire probability distributions -- as
opposed to capturing only the mean of a certain variable -- is of prime
importance in many different disciplines from natural sciences to social
sciences and beyond. Hence, likelihood-based models and distributional
regression have been receiving increasing interest over the last decades. In
order to assess the goodness of fit along with potential deficits of such models,
graphical evaluations in terms of visualizations are an important complement to
numerical evaluations (e.g., based on proper scoring rules).

Various types of diagnostic graphics have been suggested in the literature for this
purpose, e.g., PIT (probability integral transform) histograms, Q-Q
(quantile-quantile) plots of (randomized) quantile residuals, rootograms,
reliability diagrams, and worm plots, among others. However, discussion in the
literature and usage in practice is somewhat scattered across different
scientific communities, different types of variables (continuous vs. discrete
vs. categorical), and different software packages. To overcome these problems,
we present the visualizations above in a common conceptual framework accompanied
by a flexible and object-oriented software implementation in the R package
_topmodels_ (https://topmodels.R-Forge.R-project.org/). It is discussed how
all these types of graphics can be understood as different strategies for visualizing
either the _marginal_ fit (on the original scale of the observations) or the
_conditional_ fit (on a probability or quantile scale). Using several empirical
examples, it is highlighted which types of displays are particularly useful
for bringing out which kind of model deficits.
