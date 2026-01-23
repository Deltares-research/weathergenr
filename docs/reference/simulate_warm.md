# Simulate synthetic series with Wavelet Autoregressive Modeling (WARM)

Generates n_sim synthetic realizations of length n.

Two operating modes are supported.

Component mode (n \>= bypass_n) - Each wavelet component is centered. -
A stationary non seasonal ARMA(p,q) model is fit once per component
using a grid search over p in 0..max_p and q in 0..max_q with AIC
selection. - Simulations are generated from the selected ARMA model
using a fixed simulator (no S3 simulate dispatch). - Component means are
re added and components are summed. - If fitting fails for a component,
a block bootstrap fallback is used for that component only.

Bypass mode (n \< bypass_n) - Component simulation is skipped. - A
stationary non seasonal ARMA(p,q) model is fit to the original series
and simulated directly. - If fitting fails, a block bootstrap fallback
is used.

## Usage

``` r
simulate_warm(
  components = NULL,
  n = NULL,
  n_sim = 1000,
  seed = NULL,
  series_obs = NULL,
  bypass_n = 15L,
  match_variance = TRUE,
  var_tol = 0.1,
  check_diagnostics = FALSE,
  verbose = TRUE
)
```

## Arguments

- components:

  Matrix, data.frame, or list of numeric vectors. Wavelet components. In
  component mode, each component must have length n. In bypass mode,
  components are only used to reconstruct the original series if
  series_obs is not provided.

- n:

  Integer. Length of each simulated series.

- n_sim:

  Integer. Number of realizations to generate. Default is 1000.

- seed:

  Optional integer. Base RNG seed for reproducibility.

- series_obs:

  Optional numeric vector. Observed original series. Provide this in
  bypass mode to avoid reconstructing the series from components.

- bypass_n:

  Integer. Threshold for bypass mode. Default is 25.

- match_variance:

  Logical. If TRUE, rescale simulated outputs to match observed standard
  deviation when relative mismatch exceeds var_tol.

- var_tol:

  Numeric in \[0, 1\]. Relative standard deviation tolerance.

- check_diagnostics:

  Logical. If TRUE, runs a Ljung Box test on ARMA residuals and warns
  when residual autocorrelation remains.

- verbose:

  Logical. If TRUE, prints informational messages.

## Value

Numeric matrix with dimension n by n_sim. Each column is a simulated
realization.

## Details

ARMA configuration - Non seasonal: seasonal terms are not considered. -
Stationary: fits enforce stationarity constraints via parameter
transforms. - No drift and no mean: include.mean is FALSE. The code
centers series before fitting and re adds the observed mean after
simulation.

Fallback logic - Fallback is applied per component in component mode.
Components that fit successfully use ARMA simulation; failed components
use block bootstrap. - In bypass mode, fallback is applied to the whole
series.
