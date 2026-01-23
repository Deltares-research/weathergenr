# Fit ARIMA for bypass mode

Fit ARIMA for bypass mode

## Usage

``` r
.fit_warm_arima_forecast(
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
