# Compute fit metrics for each realization

Computes Mean Absolute Error (MAE) for all key metrics. MAE measures the
average magnitude of errors between simulated and observed values. Lower
MAE indicates better fit.

## Usage

``` r
.summarize_realization_fit(obs_results, sim_results, variables)
```

## Arguments

- obs_results:

  Observed data results from .summarize_observed_data

- sim_results:

  Simulated data results from .summarize_simulated_data

- variables:

  Character vector of variable names
