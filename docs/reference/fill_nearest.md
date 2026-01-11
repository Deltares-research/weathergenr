# Fill NA values with nearest non-NA value

Forward and backward fills NA values in a vector using the nearest
non-NA value. Useful for handling edge effects in interpolated data.

## Usage

``` r
fill_nearest(v)
```

## Arguments

- v:

  Numeric vector possibly containing NAs

## Value

Numeric vector with NAs filled
