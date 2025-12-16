# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

weathergenr is an R package for generating synthetic gridded daily weather series using a stochastic weather generator. It implements the methodology from Steinschneider et al. (2013), combining wavelet-based annual simulation with Markov chain daily disaggregation.

## Build and Development Commands

```r
# Document package (regenerate NAMESPACE and .Rd files)
devtools::document()

# Run all tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-generateWeatherSeries.R")

# Run package checks (full R CMD check)
devtools::check()

# Build and install locally
devtools::build()
devtools::install()
```

## Architecture

### Core Simulation Pipeline

The main entry point is `generateWeatherSeries()` which orchestrates:

1. **Preprocessing**: Removes leap days, computes water years, aggregates grid cells to domain averages
2. **Annual simulation (WARM)**: `wavelet_arima()` generates synthetic annual precipitation using wavelet decomposition + ARIMA modeling. Traces are filtered by `filter_warm_simulations()` based on statistical criteria.
3. **Daily disaggregation**: `resample_weather_dates()` uses:
   - Annual KNN (`knn_sample()`) to select historical years matching simulated annual totals
   - Three-state Markov chain (`markov_next_state()`, `monthly_markov_probs()`) for wet/dry/extreme spell persistence
   - Daily KNN resampling for precipitation and temperature anomalies

### Climate Change Application

`imposeClimateChanges()` applies climate perturbations to generated weather:
- Precipitation: quantile mapping with mean/variance change factors
- Temperature: additive delta factors (supports transient/step changes)
- PET: recalculated via `pet_hargreaves()`

### Key Data Structures

- **weather.data**: Named list of data frames, one per grid cell, with columns for each climate variable
- **weather.grid**: Data frame with `id`, `xind`, `yind`, `x`, `y` for spatial metadata
- **weather.date**: Date vector aligned with weather.data rows

### NetCDF I/O

- `read_netcdf()`: Reads gridded NetCDF into the weather.data/grid/date structure
- `write_netcdf()`: Writes simulation results back to NetCDF

## Key Parameters

- `month.start`: 1 for calendar year, other values (e.g., 10) for water year simulation
- `warm.signif.level`: Significance level for wavelet component detection (default 0.90)
- `mc.wet.quantile` / `mc.extreme.quantile`: Define three precipitation states for Markov chain
- `dry.spell.change` / `wet.spell.change`: Monthly vectors for adjusting spell persistence

## Testing

Tests use testthat (edition 3). Manual/integration tests are in `tests/tests_manual/` and are not run automatically.
