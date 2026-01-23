# Simulate from ARIMA fit for WARM

Simulate from ARIMA fit for WARM

## Usage

``` r
.warm_simulate_from_fit(fit, n, n_sim)
```

## Arguments

- fit:

  List. Output from .warm_fit_arima_safe or .fit_warm_arima_forecast.

- n:

  Integer. Simulation length.

- n_sim:

  Integer. Number of realizations.

## Value

Numeric matrix n by n_sim.
