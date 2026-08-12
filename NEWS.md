# weathergenr (development version)

## Performance

* Calendar fields are extracted with a single `as.POSIXlt()` conversion instead
  of `as.integer(format(date, "%Y"))`, which routed every value through string
  formatting and an integer re-parse. On a 27,375-day series
  `compute_water_year()` drops from 536 ms to 7.5 ms. Applied across
  `compute_water_year()`, `find_leap_day_indices()`, `build_historical_dates()`,
  `read_netcdf()`, the evaluation summaries, and the perturbation time indices.
* `resample_weather_dates()` accumulates resampled dates in a plain double
  vector rather than a `Date`-classed one. Assigning into a classed vector
  dispatches `[<-.Date` and copies the whole vector on every assignment, making
  the daily loop quadratic in `n_years`; the class is restored once on exit.
  The saving grows with simulation length, from about 8% at 20 years to about
  a third at 80 years (14.7 s to 9.6 s).
* `read_netcdf()` parses the NetCDF time axis once instead of once per variable
  when dropping Feb 29.
* `resample_weather_dates()` no longer keeps a per-year cache of the drawn
  observed subset. Its key was built from `annual_knn_n` order-dependent draws
  with replacement, so it never hit, while retaining every subset for the
  lifetime of the call in each parallel worker.
* `.align_obs_sim_periods()` derives the year-filter row index once per side
  rather than once per grid per realization.
* The three-state precipitation classification uses a sum of two comparisons
  instead of nested `ifelse()`, which allocated a logical mask and both branch
  vectors on every call. The encoding is monotone in precipitation, so the
  result is identical; the operation itself is about 4.5x faster.
* `estimate_monthly_markov_probs()` no longer allocates nine vectors of length
  `365 * n_years` per simulated year to fill 365 rows of each -- quadratic in
  `n_years` in aggregate. The computation now returns one row per month
  internally, and the daily resampling loop indexes it by month. The exported
  function's return value is unchanged: it still broadcasts those rows across
  the simulated time axis, now with one vectorised assignment instead of twelve
  full-length `which()` scans.
* End to end on the bundled ntoum fixture (30 years, 3 realizations,
  `save_plots = FALSE`): generation 9.3 s to 5.7 s, evaluation 3.6 s to 1.1 s,
  together 12.9 s to 6.8 s (-47%). Timings are the minimum of four runs in one
  session, since machine load shifts absolute figures between sessions.
  **Outputs are unchanged**: for a fixed `seed` these changes are bit-identical,
  verified by the end-to-end baseline gate over both water-year and
  calendar-year scenarios.

## Breaking changes

* `estimate_monthly_markov_probs()` and `match_transition_positions()` now error
  when `wet_threshold` exceeds `extreme_threshold`. The three-state encoding has
  always assumed the thresholds are ordered; an inverted pair previously passed
  silently and produced a "wet" state no observation could occupy.
  `generate_weather()` already validated `extreme_q > wet_q`, so this only
  affects direct calls to the two exported functions.

## Bug fixes

* `evaluate_weather_generator()` no longer emits tidyselect deprecation warnings
  ("Use of `.data` in tidyselect expressions was deprecated in tidyselect
  1.2.0"). Thirty `.data$` references inside `dplyr::select()`,
  `dplyr::rename()`, `tidyr::pivot_longer()` and `tidyr::pivot_wider()` now pass
  column names as strings. Results are unchanged; the `.data` pronoun is still
  used in data-masking verbs, where it remains correct.

## Internal

* `DESCRIPTION` declares authorship with `Authors@R` instead of the deprecated
  `Author`/`Maintainer` pair.
* The abandoned Zarr I/O prototypes moved from `inst/experiments/` to
  `dev/scripts/zarr-prototype/`, so no part of the development tree sits under
  `inst/` any more.
* Added `testthat` coverage for `R/qm_diagnostics.R` and the wavelet plotting
  functions, including regression tests for the spell-length and dry-day
  diagnostics' `NA` handling.
* Added pkgdown and lint GitHub Actions workflows, plus a `.lintr` baseline
  that scopes linting to semantic findings.

# weathergenr 1.2.0

## Bug fixes

* Fixed warnings raised by the evaluation and diagnostic plotting functions.
* Fixed netCDF template propagation when writing gridded output.

## Internal

* `inst/experiments/` and `inst/extcode/` are no longer installed with the
  package. Both held development scratch material — including abandoned Zarr
  I/O prototypes — that was previously shipping into user libraries.
* Consolidated the two `utils::globalVariables()` declarations into
  `R/globals.R` and removed `R/zzz_globals.R`. Of the 42 names the removed
  file declared, 28 were unused. `colorRampPalette` was among them: declaring
  a function there also silences R CMD check's "no visible global function
  definition", which would have masked the missing `grDevices` import had a
  call been reintroduced.
* Dropped the `LazyData` field. The package ships no `data/` directory, so it
  produced an R CMD check note without effect.
* `AGENTS.md` is now the canonical agent-instruction file, with `CLAUDE.md`
  importing it.

Versions before 1.2.0 predate this changelog; see the commit history.
