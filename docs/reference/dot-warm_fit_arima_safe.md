# Fit ARIMA for WARM with safe fallback

Fit ARIMA for WARM with safe fallback

## Usage

``` r
.warm_fit_arima_safe(
  x,
  max_p,
  max_q,
  stationary = TRUE,
  include_mean = FALSE,
  allow_drift = FALSE
)
```

## Arguments

- x:

  Numeric vector (centered).

- max_p:

  Integer. Maximum AR order.

- max_q:

  Integer. Maximum MA order.

- stationary:

  Logical. Whether to enforce stationarity.

- include_mean:

  Logical. Whether to include mean.

- allow_drift:

  Logical. Whether to allow drift.

## Value

List with model or parameter entries, or NULL.
