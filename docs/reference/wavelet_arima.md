# Wavelet Autoregressive Modeling (WARM)

Simulates synthetic time series by modeling each wavelet component
(signal or noise) with an ARIMA model, and then summing the simulated
components. Optionally enforces variance matching to preserve
statistical properties of the original components.

## Usage

``` r
wavelet_arima(
  wavelet.components = NULL,
  sim.year.num = NULL,
  sim.num = 1000,
  seed = NULL,
  match.variance = TRUE,
  variance.tolerance = 0.1,
  check.diagnostics = FALSE,
  verbose = TRUE
)
```

## Arguments

- wavelet.components:

  A list or matrix where each column (or list element) is a numeric
  vector corresponding to a wavelet component (low-frequency signal or
  noise).

- sim.year.num:

  Integer. Desired length (number of years or timesteps) of each
  simulated series.

- sim.num:

  Integer. Number of synthetic series to produce. Default: 1000.

- seed:

  Optional. Integer random seed for reproducibility.

- match.variance:

  Logical. If TRUE, rescale simulated components to match the variance
  of the original components (default: TRUE). Recommended for preserving
  statistical properties.

- variance.tolerance:

  Numeric. Relative tolerance for variance matching (default: 0.1 =
  10%). Only applies if match.variance = TRUE.

- check.diagnostics:

  Logical. If TRUE, perform basic ARIMA model diagnostics and issue
  warnings if models appear inadequate (default: FALSE).

- verbose:

  Logical. If TRUE, emit informative logger::log_info() messages
  (default: TRUE). If FALSE, suppresses all logger::log_info() output
  from this function.

## Value

A matrix of dimension `sim.year.num` x `sim.num`, where each column is a
synthetic time series realization.
