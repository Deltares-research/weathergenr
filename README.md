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
is adapted here for gridded datasets and netcdf-based workflows.

The package is intended for workflows such as:

- climate risk and stress testing
- hydrological and water-resources modelling
- scenario analysis for climate adaptation studies

## Methodological framework

The generator represents climate variability across multiple time scales
by coupling low-frequency climate dynamics with realistic daily weather
sequences. It consists of three components. The first two are coupled
inside `generate_weather()` — the annual states from (1) condition the
daily resampling in (2). The third is a separate stage you apply to the
generated series, not a step inside the generator; see *Composing the
stages* below.

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
internal consistency. Applied by `apply_climate_perturbations()` to an
already-generated series.

## Composing the stages

`run_weather_generator()` runs stages 1–2 and evaluates the result
against the observed record. Perturbation is deliberately outside that
pipeline: evaluation asks whether the synthetic series reproduces the
observed climate, whereas a perturbed series is *meant* to depart from
it. So stress scenarios are built by generating first, then perturbing:

``` r
gen <- generate_weather(
  obs_data = ncdata$data, obs_grid = ncdata$grid, obs_dates = ncdata$date,
  vars = c("precip", "temp", "temp_min", "temp_max"), ...
)

# Turn resampled dates into daily values
eval_data <- prepare_evaluation_data(
  gen_output = gen, obs_data = ncdata$data, obs_dates = ncdata$date,
  grid_ids = seq_along(ncdata$data),
  variables = c("precip", "temp", "temp_min", "temp_max")
)
sim <- eval_data$sim_data[[1]]

future <- apply_climate_perturbations(
  data = sim, grid = ncdata$grid,
  date = sim[[1]]$date,
  precip_mean_factor = rep(0.7, 12),
  temp_delta = rep(2, 12)
)
```

Note `vars` must include `temp_min` and `temp_max`, which
`apply_climate_perturbations()` requires, and that the date vector comes
from `sim`, not from `gen$dates`. The
[Climate Perturbations](https://deltares-research.github.io/weathergenr/articles/climate_perturbation.html)
article covers both points.

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
https://deltares-research.github.io/weathergenr/articles/getting_started.html

## Styling and exporting diagnostics

Use `config$return_plots = TRUE` with `config$save_plots = FALSE` to collect
all diagnostics without writing figures. `run_weather_generator()` returns
four generation figures in `result$gen_output$plots` and the evaluation plots
in `result$evaluation`. Direct `generate_weather()` calls accept the same
`return_plots`, `save_plots` and `plot_style` arguments.

```r
config$save_plots <- FALSE
config$return_plots <- TRUE
config$plot_style <- list(
  theme = ggplot2::theme(text = ggplot2::element_text(size = 10, family = "sans")),
  colors = c(observed = "#0072B2", simulated = "#666666", reference = "#D55E00"),
  ensemble_palette = c("#B3CDE3", "#005B96"),
  power_palette = c("#FFFFCC", "#253494")
)
result <- run_weather_generator(
  obs_data = ncdata$data, obs_grid = ncdata$grid, obs_dates = ncdata$date,
  config = config, out_dir = "output"
)
export_weather_plot(
  result$gen_output$plots$obs_power_spectra, "output/spectra.png",
  width = 180, height = 90, dpi = 300, bg = "white"
)
```

`plot_style` is an ordinary R list; callers own YAML parsing or token translation.
Omitted settings retain existing defaults. Partial themes preserve local layout
unless explicitly overridden, and explicit font sizes survive automatic export
scaling. See `?export_weather_plot` for color-role and palette semantics.


For the optional compact diagnostic profile, restyle returned objects explicitly:

```r
p <- style_weather_plot(result$evaluation$daily_mean, "daily_mean",
  plot_style = config$plot_style,
  panel_layout = list(panel_width_mm = 70, panel_height_mm = 70,
    preferred_columns = c(3, 2), max_width_mm = 180))
export_weather_plot(p, "output/daily_mean.png", dpi = 600)
```

This profile uses gray stochastic traces, black central summaries, blue observations
and red significance guides. Comparison circles are ensemble medians; boxplot
diamonds are means. Wavelet composites and the narrow WARM summary retain their
local layouts. `arrange_weather_panels()` also arranges other single-panel ggplots
with fixed physical plotting areas, measured label space and row-major order.
It never shrinks text or marks to fit. See `?style_weather_plot` and
`?arrange_weather_panels` for profile semantics and layout settings.

# References

> Steinschneider, S., & Brown, C. (2013). *A semiparametric
> multivariate, multisite weather generator with low-frequency
> variability for use in climate risk assessments.* **Water Resources
> Research**, 49(11), 7205-7220. (https://doi.org/10.1002/wrcr.20528)

# License

Licensed under the MIT License. See `LICENSE` for details.
