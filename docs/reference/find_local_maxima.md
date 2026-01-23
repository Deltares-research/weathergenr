# Find local maxima indices in a numeric vector

Scans a numeric vector and returns indices of local maxima. A local
maximum is defined by comparing each element to its immediate neighbors.

Uses vectorized operations for O(n) complexity instead of loop-based
O(n^2).

## Usage

``` r
find_local_maxima(x, strict = TRUE)
```

## Arguments

- x:

  Numeric vector. Values to scan for local maxima.

- strict:

  Logical scalar. If TRUE, requires strict inequality on both sides. If
  FALSE, allows ties but still requires at least one strict inequality
  so flat plateaus do not produce multiple peaks.

## Value

Integer vector of peak indices. Returns an empty integer vector if fewer
than three values are available.
