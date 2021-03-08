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

**Keywords**: Distributional regression, probabilistic forecasts, random forest, graphical model assessment.

**Webpages**: *R* repositories [disttree](https://R-Forge.R-project.org/projects/partykit/pkg/disttree/) and [topmodels](https://R-Forge.R-project.org/projects/topmodels/pkg/topmodels/) on R-Forge.

As a consequence of the growing importance of probabilistic predictions in
various application areas due to a necessary functional risk management and
strategy, there is an increased demand for appropriate probabilistic model
evaluation methods. Besides proper probabilistic scores
[@Gneiting+Raftery:2007], which evaluate not only the expectation but the
entire predictive distribution, simple graphical assessment methods are
particularly advantageous to diagnose possible model misspecification problems.  

Probabilistic predictions are mostly based on distributional (regression)
models, whereby the computation of predictive distributions, probabilities and
quantiles is generally *R* package dependent. Therefore, unified graphical
evaluation routines for probabilistic models are limited and graphical model
assessment are at best provided by the respective package or have to be written
by the user himself.

Using distribution based random forest models
[@Schlosser+Hothorn+Stauffer+Zeileis:2019], we present a planned unifying
infrastructure that provides predictions of probabilities, densities, scores,
and Hessian values for various probabilistic models and distributional
regression models. By means of the unifying prediction S3 method, many different
graphical assessment tools are available to the users out of the box, such as
reliability diagrams, PIT histograms, rootograms, and randomized QQ plots.

The random forest models used for illustration are provided in the *R* package
**disttree** which is available on R-Forge and includes routines for estimating
distributional trees and forests [@Schlosser+Lang+Zeileis:2019]. The proposed
unifying infrastructure for inference and forecasting in probabilistic models
is provided in the *R* package **topmodels**, also available on R-Forge
[@Zeileis+Kleiber+Kosmidis:2018]. The actual workhorse of the package is the S3
function `procast`, a generic function for computing various types of
predictions (e.g., density functions, probabilities, quantiles), which is
currently implemented for `lm`, `crch`, and `disttree` model classes. For
these, **topmodels** provides various graphical assessment routines, e.g.,
`reliagram`, `pithist` and `rootogram` for probabilistic (distributional)
forecasts.

## References
