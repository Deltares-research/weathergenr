# Safely Compute Correlation with Pairwise Completeness Check

Computes a correlation coefficient between two numeric vectors using
only finite paired observations. If the number of valid pairs is below a
required minimum, the correlation is not computed.

## Usage

``` r
.safe_cor(x, y, min_pairs = 3L, method = "pearson")
```

## Arguments

- x, y:

  Numeric vectors of equal length. Only pairs where both values are
  finite are used.

- min_pairs:

  Integer \>= 1. Minimum number of paired observations required to
  compute the correlation.

- method:

  Character. Correlation method passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). One of
  `"pearson"`, `"spearman"`, or `"kendall"`.

## Value

Named numeric vector with elements:

- `value`: correlation coefficient, or `NA_real_`,

- `n`: number of finite paired observations.
