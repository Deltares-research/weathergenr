# Fit stationary non seasonal ARMA model by AIC grid search

Fits ARMA(p,q) models using stats::arima with include.mean = FALSE and
selects the model with minimum AIC.

## Usage

``` r
.fit_arma_aic(x, max_p = 2L, max_q = 2L)
```

## Arguments

- x:

  Numeric vector. Should be centered.

- max_p:

  Integer. Maximum AR order.

- max_q:

  Integer. Maximum MA order.

## Value

List with ar, ma, sigma2, residuals, and order. Returns NULL if no
candidate fit succeeds.
