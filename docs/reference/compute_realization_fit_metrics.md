# Compute fit metrics for each realization

Computes Mean Absolute Error (MAE) for all key metrics. MAE measures the
average magnitude of errors between simulated and observed values. Lower
MAE indicates better fit.

## Usage

``` r
compute_realization_fit_metrics(obs.results, sim.results, variables)
```

## Arguments

- obs.results:

  Observed data results from process_observed_data

- sim.results:

  Simulated data results from process_simulated_data

- variables:

  Character vector of variable names
