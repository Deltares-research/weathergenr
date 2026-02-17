# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**weathergenr** is an R package implementing a gridded, semiparametric, multivariate, multisite stochastic weather generator based on the Steinschneider & Brown (2013) framework. It is used for climate risk assessment, hydrological modeling, and bottom-up climate vulnerability assessments.

## Build, Test, and Lint Commands

### Testing
```r
devtools::test()              # Run all tests
testthat::test_file("tests/testthat/test-calendar.R")  # Run single test file
```

### Linting
```bash
Rscript tools/lint.R              # Lint whole package
Rscript tools/lint.R --changed    # Lint only changed files
Rscript tools/lint.R R/file.R     # Lint specific file
```

### Development Workflow (from tools/dev_workflow.R)
```r
source("tools/dev_workflow.R")
publish_docs()           # Complete 9-step publish workflow
quick_site(open=TRUE)    # Fast site preview (no checks)
build_vignettes_only()   # Build Quarto vignettes only
check_only()             # R CMD check only
deep_clean()             # Remove all build artifacts
```

### Building and Installing
```r
devtools::document()     # Generate/update documentation
devtools::build()        # Build package
devtools::install()      # Install locally
devtools::check()        # Full R CMD check
```

## Architecture

The package couples three statistical approaches:

1. **WARM (Wavelet Autoregressive Modeling)**: Low-frequency climate variability simulation applied to annual aggregates, preserving interannual to decadal patterns
2. **Daily Weather Generation**: Multi-state Markov chain (dry/wet/extreme) with K-nearest neighbor resampling for daily values
3. **Climate Perturbation**: Quantile-based modifications for controlled climate stress testing

### Core Module Organization (R/)

| Module | Purpose |
|--------|---------|
| `generator.R` | Main workflow combining WARM, Markov chains, and daily resampling |
| `warm_filtering.R` | Filters WARM realizations using distributional/tail/spectral criteria |
| `resample.R` | Daily weather disaggregation and KNN-based resampling |
| `wavelet_warm.R` | Wavelet analysis for WARM annual series simulation |
| `wavelet_cwt.R` | Continuous wavelet transform (Morlet) implementations |
| `wavelet_modwt.R` | Maximum overlap discrete wavelet transform (wrapper around waveslim) |
| `quantile_mapping.R` | Gamma quantile mapping for precipitation with tail adjustments |
| `climate_perturbations.R` | Apply controlled climate perturbations to generated weather |
| `evaluate_generator.R` | Metrics for evaluating synthetic weather quality |
| `io_netcdf.R` | Read/write NetCDF gridded data with CF conventions |
| `io_zarr.R` | Read/write Zarr format (cloud-optimized arrays) |
| `calendar.R` | Calendar utilities (water years, leap day handling, 365-day normalization) |
| `pet.R` | Potential evapotranspiration (Hargreaves method) |

### Key Design Patterns

- **365-day calendar**: The package normalizes all dates to a 365-day calendar by removing leap days
- **Gridded workflow**: Operates on spatial grids; functions support parallel processing over grid cells
- **Parallel processing**: Uses PSOCK clusters with `doParallel` for grid cell iteration
- **Seed control**: RNG seeds passed through workflow for reproducibility

### Main Entry Points

- `run_weather_generator()`: End-to-end workflow combining all components
- `generate_weather()`: Core generation function
- `simulate_warm()`: WARM simulation for annual traces
- `filter_warm_pool()`: Adaptive filtering of WARM realizations
- `apply_climate_perturbations()`: Apply scenario changes to generated weather
- `evaluate_weather_generator()`: Comprehensive quality metrics

## Testing

13 test files in `tests/testthat/` covering:
- Calendar conversions and water year handling
- WARM simulation and filtering
- Wavelet analysis (CWT, MODWT)
- Daily resampling and KNN logic
- Quantile mapping and climate perturbations
- NetCDF I/O
- Integration tests for the full generator

Tests use testthat v3.0+ with `Edition: 3` configuration.

## CI/CD

GitHub Actions workflow (`.github/workflows/R-CMD-check.yaml`) runs on:
- macOS, Windows, Ubuntu
- R versions: release, devel, oldrel-1
- Includes Quarto for vignette rendering

## Documentation

- Vignettes are Quarto format (`.qmd`) in `vignettes/`
- pkgdown site configuration in `_pkgdown.yml`
- Documentation generated via roxygen2 (`devtools::document()`)
