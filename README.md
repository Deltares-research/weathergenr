# Getting started


[![R-CMD-check](https://github.com/Deltares-research/weathergenr/actions/workflows/R-CMD-check.yaml/badge.svg?branch=master)](https://github.com/Deltares-research/weathergenr/actions/workflows/R-CMD-check.yaml?query=branch%3Amaster)
[![CRAN
status](https://www.r-pkg.org/badges/version/weathergenr.png)](https://CRAN.R-project.org/package=weathergenr)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
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

1.  *Modeling Low-frequency variability*: Interannual to decadal
    variability is simulated using wavelet autoregressive modeling
    (WARM) applied to annual climate aggregates (typically
    precipitation). Wavelet decomposition isolates low-frequency
    components, which are stochastically simulated and recombined to
    produce synthetic annual sequences that preserve observed
    persistence and spectral structure. These annual states condition
    the daily generator.

2.  *Daily weather generation (Markov chain + KNN)*: Daily precipitation
    occurrence is modeled using a multi-state Markov chain to reproduce
    wet and dry spell persistence, while precipitation and temperature
    amounts are generated via k-nearest-neighbour (KNN) resampling. This
    preserves seasonality, cross-variable dependence, and spatial
    correlations, with historical analogues selected conditional on the
    simulated annual climate state.

3.  *Climate perturbations and quantile mapping*: Quantile-based
    perturbation methods (e.g. quantile mapping) are used to impose
    controlled changes in temperature and precipitation distributions,
    enabling systematic stress-testing of shifts in means, variability,
    and extremes while retaining internal temporal and spatial
    consistency.

## Intended use

The resulting framework is intended for *bottom-up climate vulnerability
assessments* to explore system response under a wide range of plausible
climate conditions.

## Key features

- support for **gridded NetCDF inputs** and outputs
- multivariate and multisite weather generation
- explicit representation of low-frequency climate variability
- efficient filtering and selection of stochastic realizations
- designed for reproducible, scenario-based workflows

# Installation

Install the latest version from GitHub:

``` r
# install.packages("devtools")
# devtools::install_github("Deltares-research/weathergenr")
```

# Getting started

A worked end-to-end tutorial is available here:
https://deltares-research.github.io/weathergenr/articles/getting_started.html

# References

> Steinschneider, S., & Brown, C. (2013). *A semiparametric
> multivariate, multisite weather generator with low-frequency
> variability for use in climate risk assessments.* **Water Resources
> Research**, 49(11), 7205-7220. (https://doi.org/10.1002/wrcr.20528)

# License

Licensed under the MIT License. See `LICENSE` for details.
