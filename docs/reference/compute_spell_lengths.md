# Compute run-lengths (spell lengths) for dry or wet states

Computes consecutive run lengths of dry or wet states relative to a
threshold.

## Usage

``` r
compute_spell_lengths(x, threshold, below = TRUE)
```

## Arguments

- x:

  Numeric vector. Precipitation values.

- threshold:

  Numeric scalar. Wet-day threshold.

- below:

  Logical. If `TRUE`, defines spells where `x < threshold` (dry spells).
  If `FALSE`, defines spells where `x >= threshold` (wet spells).

## Value

Numeric vector of positive integers giving spell lengths. May be length
zero.

## Details

Uses [`rle()`](https://rdrr.io/r/base/rle.html) on the logical state
series to identify consecutive runs. If no spells of the requested type
occur, the function returns `numeric(0)`.
