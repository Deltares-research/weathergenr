# Compute correlation matrix on monthly anomalies

Compute correlation matrix on monthly anomalies

## Usage

``` r
.compute_anomaly_correlations(data, variables)
```

## Arguments

- data:

  Data frame with \`date\`, \`id\`, and variables.

- variables:

  Character vector of variable names to correlate.

## Value

Long-format data frame of upper-triangle pairwise correlations.
