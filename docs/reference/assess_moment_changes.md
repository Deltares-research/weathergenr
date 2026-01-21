# Assign qualitative assessments for moment preservation

Classifies moment changes based on absolute percent change thresholds
that differ by metric.

## Usage

``` r
assess_moment_changes(moments_df)
```

## Arguments

- moments_df:

  Data.frame. Output from
  [`compute_moment_diagnostics`](https://deltares-research.github.io/weathergenr/reference/compute_moment_diagnostics.md)
  containing at least `metric` and `pct_change`.

## Value

Character vector of length `nrow(moments_df)` with assessment labels.

## Details

The function applies metric-specific thresholds:

- mean/variance: excellent \< 5, good \< 15, else poor

- cv: excellent \< 10, good \< 20, else poor

- others (sd, skewness, kurtosis): good \< 15, acceptable \< 30, else
  poor
