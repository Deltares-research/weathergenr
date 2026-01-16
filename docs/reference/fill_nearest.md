# Fill NA Values with the Nearest Non-NA Neighbor

Fills missing values in a numeric vector by propagating the nearest
available non-missing value. This is a deterministic forward-fill
followed by a backward-fill. It is mainly used to stabilize curves after
interpolation (e.g., regridded spectra).

## Usage

``` r
fill_nearest(x)
```

## Arguments

- x:

  Numeric vector possibly containing `NA`.

## Value

Numeric vector with missing values filled. If `x` is all `NA`, it is
returned unchanged.
