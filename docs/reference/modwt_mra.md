# MODWT multiresolution analysis (MRA)

Extracts additive time-domain components from a MODWT decomposition
using multiresolution analysis. Each component captures variance within
a dyadic scale band. The components sum (approximately, under periodic
boundary handling) to the original series.

## Usage

``` r
modwt_mra(
  x,
  filter = c("la8", "haar", "d4", "d6", "d8", "la16"),
  n_levels = NULL,
  boundary = "periodic",
  include_smooth = TRUE,
  max_period_frac = 1/3
)
```

## Arguments

- x:

  Numeric vector. Input time series. Must contain no missing values.

- filter:

  Character. Wavelet filter name. One of `"la8"` (default), `"haar"`,
  `"d4"`, `"d6"`, `"d8"`, `"la16"`.

- n_levels:

  Integer scalar or NULL. Number of decomposition levels (J). If NULL, a
  conservative default is selected subject to stability and max-period
  caps.

- boundary:

  Character. Boundary handling method. Only `"periodic"` is supported.

- include_smooth:

  Logical. If TRUE, includes the smooth (approximation) component at the
  coarsest level (SJ).

- max_period_frac:

  Numeric scalar in (0, 1\]. Maximum represented period as a fraction of
  record length. Default is `1/3`, meaning no structure beyond
  approximately n/3.

## Value

A list with additive components and summary diagnostics. The list
contains: - `components`: numeric matrix (n x n_components) with columns
D1..DJ and optionally SJ - `periods`: numeric vector of representative
periods per component - `variance`: component variances -
`variance_fraction`: component variance divided by total variance of x -
`total_variance`: variance of x - `sum_component_variance`: sum of
component variances - `cross_covariance_sum`: sum of cross-covariances
(pairwise, without doubles) - `variance_identity_error`: total_var -
(sum_var + 2 \* cross_cov_sum) - `reconstruction_error`: max absolute
error between x and reconstructed series - `n_levels`, `filter`,
`filter_length`, `boundary`, `include_smooth`, `max_period_frac`, and
`is_additive`
