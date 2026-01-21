# Getting started

[![R-CMD-check](https://github.com/Deltares-research/weathergenr/actions/workflows/R-CMD-check.yaml/badge.svg?branch=master)](https://github.com/Deltares-research/weathergenr/actions/workflows/R-CMD-check.yaml?query=branch%3Amaster)
[![CRAN
status](https://www.r-pkg.org/badges/version/weathergenr.png)](https://CRAN.R-project.org/package=weathergenr)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://deltares-research.github.io/weathergenr/LICENSE)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)

# Overview

**weathergenr** is an R package that implements a semiparametric,
multivariate, multisite stochastic weather generator designed for
climate risk and stress-testing applications. The approach is
conceptually based on the framework of Steinschneider & Brown (2013) and
is adapted here for gridded datasets and modern NetCDF-based workflows.

The package is intended for workflows such as:

- climate risk and stress testing
- hydrological and water-resources modelling
- scenario analysis for climate adaptation studies

## Methodological framework

The generator represents climate variability across multiple time scales
by coupling low-frequency climate dynamics with realistic daily weather
sequences. It consists of three components:

The generator combines low-frequency climate variability with realistic
daily weather sequences using a three-stage approach:

**1. Low-frequency climate variability (WARM)** Interannual to decadal
variability is modeled with wavelet autoregressive methods applied to
annual climate aggregates. This preserves persistence and spectral
structure and defines annual climate states that condition daily
weather.

**2. Daily weather generation (Markov chain + KNN)** Wet-dry persistence
is simulated with a multi-state Markov chain, while daily precipitation
and temperature values are generated via K-nearest-neighbour resampling.
This maintains seasonality, cross-variable dependence, and spatial
coherence.

**3. Climate perturbation and stress testing** Quantile-based
perturbations impose controlled changes in means, variability, and
extremes, enabling systematic climate stress testing while preserving
internal consistency.

## Intended use

The resulting framework is intended for *bottom-up climate vulnerability
assessments* to explore system response under a wide range of plausible
climate conditions.

## Key features

- netcdf i/o for convenient coupling with other models
- efficient generation, filtering, and evaluation of large stochastic
  ensembles
- integrated diagnostics and visualization for validation and analysis

# Installation

Install the latest version from GitHub:

``` r

# install.packages("devtools")
# devtools::install_github("Deltares-research/weathergenr")
```

# Getting started

A quick tutorial is available here:
<https://deltares-research.github.io/weathergenr/articles/getting_started.html>

# References

> Steinschneider, S., & Brown, C. (2013). *A semiparametric
> multivariate, multisite weather generator with low-frequency
> variability for use in climate risk assessments.* **Water Resources
> Research**, 49(11), 7205-7220. (<https://doi.org/10.1002/wrcr.20528>)

# License

Licensed under the MIT License. See `LICENSE` for details.
