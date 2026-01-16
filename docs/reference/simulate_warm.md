# Wavelet Autoregressive Modeling (WARM)

Simulates synthetic series by modeling each wavelet component with an
ARIMA model, simulating each component forward, and summing across
components. Optionally rescales simulated components to match the
observed component variance.

## Usage

``` r
simulate_warm(
  components = NULL,
  n = NULL,
  n_sim = 1000,
  seed = NULL,
  match_variance = TRUE,
  var_tol = 0.1,
  check_diagnostics = FALSE,
  verbose = TRUE
)
```

## Arguments

- components:

  Matrix, data.frame, or list of numeric vectors. Wavelet components,
  typically produced by
  [`extract_wavelet_components`](https://deltares-research.github.io/weathergenr/reference/extract_wavelet_components.md).

- n:

  Integer. Length of each simulated series.

- n_sim:

  Integer. Number of realizations to generate.

- seed:

  Optional integer. Base RNG seed for reproducibility.

- match_variance:

  Logical. If `TRUE`, rescales each simulated component to match the
  observed component standard deviation.

- var_tol:

  Numeric in \[0, 1\]. Relative tolerance used to trigger variance
  rescaling.

- check_diagnostics:

  Logical. If `TRUE`, runs simple ARIMA diagnostics and warns on issues.

- verbose:

  Logical. If `TRUE`, emits informational logs (via
  [`logger::log_info`](https://daroczig.github.io/logger/reference/log_level.html)
  when available).

## Value

Numeric matrix of dimension `n x n_sim`. Each column is a simulated
realization.
