# Estimate Monthly Occurrence-State Markov Transition Probabilities

Estimates month-specific transition probabilities for a three-state
precipitation occurrence Markov chain (dry, wet, very wet) from observed
daily precipitation data, and maps these probabilities onto a simulated
time axis.

Transition probabilities are inferred from observed lag-1 to lag-0
precipitation pairs using month-dependent wet and extreme thresholds. A
Dirichlet-type smoothing is applied to avoid zero-probability
transitions, with the effective smoothing strength decreasing
automatically with the number of observed transitions per month.

The function supports both calendar-year and water-year simulations. The
regime is inferred from month.list: if month.list\[1\] == 1,
calendar-year logic is used; otherwise, water-year logic is assumed.

Optional spell-persistence modifiers are applied to dry and wet states
to influence transition persistence. After these adjustments, transition
probabilities are explicitly normalized to ensure stochastic validity
(non-negativity and rows summing to one).

The estimated monthly transition probabilities are assigned to all
simulated days that fall within the target simulation year and
corresponding calendar month.

## Usage

``` r
monthly_markov_probs(
  precip.lag0,
  precip.lag1,
  month.lag0,
  month.lag1,
  year.lag0 = NULL,
  year.lag1 = NULL,
  wet.threshold,
  extreme.threshold,
  month.list,
  sim.months,
  sim.water.years,
  year.idx,
  sim.start.year,
  dry.spell.change,
  wet.spell.change,
  sim.length,
  alpha = 1
)
```

## Arguments

- precip.lag0:

  Numeric vector. Observed daily precipitation at lag 0 (current day)
  used to estimate state transitions.

- precip.lag1:

  Numeric vector. Observed daily precipitation at lag 1 (previous day)
  used to estimate state transitions.

- month.lag0:

  Integer vector (1-12). Calendar month corresponding to precip.lag0.

- month.lag1:

  Integer vector (1-12). Calendar month corresponding to precip.lag1.

- year.lag0:

  Optional integer vector. Calendar year or water year corresponding to
  precip.lag0. If provided together with year.lag1, transitions are
  restricted to same-year pairs.

- year.lag1:

  Optional integer vector. Calendar year or water year corresponding to
  precip.lag1.

- wet.threshold:

  Numeric vector of length 12. Monthly precipitation thresholds
  separating dry and wet states, aligned to month.list.

- extreme.threshold:

  Numeric vector of length 12. Monthly precipitation thresholds
  separating wet and very wet states, aligned to month.list.

- month.list:

  Integer vector of length 12 defining the simulated ordering of months
  (for example 1:12 for calendar years or c(10:12, 1:9) for
  October-start water years).

- sim.months:

  Integer vector. Calendar month for each simulated day.

- sim.water.years:

  Integer vector. Simulation year (calendar or water year, depending on
  regime) for each simulated day.

- year.idx:

  Integer. Index of the current simulated year (1-based).

- sim.start.year:

  Integer. First year of the simulation, interpreted as a calendar year
  or water year depending on the simulation regime.

- dry.spell.change:

  Numeric vector of length 12. Monthly multiplicative factors
  controlling dry-state persistence.

- wet.spell.change:

  Numeric vector of length 12. Monthly multiplicative factors
  controlling wet-state persistence.

- sim.length:

  Integer. Total number of simulated days. Determines the length of the
  returned probability vectors.

- alpha:

  Numeric. Base smoothing strength for transition-probability
  estimation. The effective smoothing applied in each month is alpha /
  sqrt(N_m), where N_m is the number of observed transitions in that
  month.

## Value

A named list of nine numeric vectors, each of length sim.length,
representing time-varying transition probabilities:

- p00_final: P(dry to dry)

- p01_final: P(dry to wet)

- p02_final: P(dry to very wet)

- p10_final: P(wet to dry)

- p11_final: P(wet to wet)

- p12_final: P(wet to very wet)

- p20_final: P(very wet to dry)

- p21_final: P(very wet to wet)

- p22_final: P(very wet to very wet)

## Details

Precipitation occurrence states are defined as:

- Dry:

  Precipitation less than or equal to the monthly wet threshold

- Wet:

  Precipitation greater than the wet threshold and less than or equal to
  the extreme threshold

- Very wet:

  Precipitation greater than the extreme threshold

All transition-probability rows are explicitly normalized after
smoothing and spell adjustments to ensure valid Markov-chain behavior.
