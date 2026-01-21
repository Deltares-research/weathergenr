# Get Positions of a Specific Occurrence State Transition Within Candidate Indices

Returns the positions within \`day0_idx\` where a specified state
transition occurs, based on wet/dry/extreme day thresholds.

## Usage

``` r
match_transition_positions(
  state_from,
  state_to,
  precip_vec,
  day0_idx,
  wet_threshold,
  extreme_threshold
)
```

## Arguments

- state_from:

  Integer. The current occurrence state (0 = dry, 1 = wet, 2 = extreme).

- state_to:

  Integer. The next occurrence state (0 = dry, 1 = wet, 2 = extreme).

- precip_vec:

  Numeric vector of precipitation values.

- day0_idx:

  Integer vector of indices representing candidate "day 0" positions in
  the time series.

- wet_threshold:

  Numeric. Threshold separating dry and wet days.

- extreme_threshold:

  Numeric. Threshold above which days are considered extreme.

## Value

Integer vector of positions within \`day0_idx\` where the transition
from \`state_from\` to \`state_to\` occurs.

## Details

Occurrence states are defined as: - 0: Dry day (precipitation \<=
wet_threshold) - 1: Wet day (precipitation \> wet_threshold and \<=
extreme_threshold) - 2: Extreme wet day (precipitation \>
extreme_threshold)

## Examples

``` r
precip_vec <- c(0, 0, 5, 15, 30, 0, 2, 25, 40, 0)
day0_idx <- 1:(length(precip_vec) - 1)
wet_threshold <- 1
extreme_threshold <- 20

match_transition_positions(0, 1, precip_vec, day0_idx, wet_threshold, extreme_threshold)
#> [1] 2 6
match_transition_positions(1, 2, precip_vec, day0_idx, wet_threshold, extreme_threshold)
#> [1] 4 7
match_transition_positions(2, 0, precip_vec, day0_idx, wet_threshold, extreme_threshold)
#> [1] 5 9
```
