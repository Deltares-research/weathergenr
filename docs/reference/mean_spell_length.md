# Mean Spell Length Above or Below a Threshold

Computes the mean length of contiguous runs ("spells") in a numeric time
series after threshold-based classification.

The input vector is first converted to a binary occurrence series using
the supplied threshold. Consecutive days in the same state are grouped
into spells, and the mean spell length is calculated for either the
below-threshold or above-threshold state.

This function is typically used to diagnose wet or dry spell persistence
in daily precipitation or other hydro-meteorological time series and is
suitable for validating stochastic weather generators and Markov-chain
occurrence models.

## Usage

``` r
mean_spell_length(x, threshold = 0, below = TRUE)
```

## Arguments

- x:

  Numeric vector. Daily values of a weather variable (for example
  precipitation or temperature). Missing values are not explicitly
  handled and should be removed beforehand if present.

- threshold:

  Numeric scalar. Threshold used to classify days into two states.
  Values less than or equal to the threshold are considered
  "below-threshold"; values above the threshold are considered
  "above-threshold".

- below:

  Logical. If TRUE (default), the mean length of below-threshold spells
  (for example dry spells) is returned. If FALSE, the mean length of
  above-threshold spells (for example wet spells) is returned.

## Value

A single numeric value giving the mean spell length (in days) for the
selected state. Returns:

- `NA_real_` if `x` has zero length,

- `0` if no spells of the selected type are present.

## Details

A spell is defined as a maximal sequence of consecutive days belonging
to the same threshold-defined state. Spell lengths are computed from
transitions in the binary occurrence series derived from `x`.

The function treats values exactly equal to the threshold as
below-threshold. This convention should be kept consistent with other
occurrence or Markov-state definitions used in the analysis.

## See also

[`markov_next_state`](https://deltares-research.github.io/weathergenr/reference/markov_next_state.md)
