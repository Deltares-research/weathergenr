# weathergenr (development version)

## New features

* `theme_weathergenr()` is a new exported theme, applied to every figure the
  package draws. Before it, a single run wrote PNGs in three different themes:
  the evaluation diagnostics used `theme_bw(base_size = 12)`, the wavelet and
  WARM panels used `theme_light()` — except one WARM panel on `theme_bw()` — and
  `annual_precip.png` applied no theme at all, so it rendered in ggplot2's
  default grey. It is exported because `evaluate_weather_generator()` returns its
  diagnostics as ggplot objects: without it, a caller re-rendering one of those
  plots cannot reproduce the package's styling.

  It wraps `ggplot2::theme_light()` at `base_size = 12` — the size the evaluation
  pipeline already produced — with title and subtitle sizes given as `rel()` so
  they track `base_size` rather than being stranded at the fixed 14 and 10 points
  they used to be.

  When the package writes a figure it scales that base size with the figure's
  own dimensions, clamped to 9–16 pt. A fixed point size makes text physically
  identical everywhere, which is not what a reader sees: these figures are
  viewed fit-to-window or embedded at a common width, so each is scaled by its
  own size first. At a fixed 12 pt the apparent size varied about twofold across
  a run — large on the 8 × 4.5 in WARM figures, small on the 12 × 13.1 in
  conditional-correlation grid. Scaling brings that spread under 1.2. Calling
  `theme_weathergenr()` directly is unaffected and still gives 12 pt.

* `generate_weather()` gains `plot_dpi` (default 300) and `plot_device`
  (default `NULL`), forwarded from `run_weather_generator()`'s `config`. The four
  figures it writes — `obs_power_spectra.png`, `warm_annual_precip.png`,
  `warm_annual_stats.png` and `warm_annual_wavelet.png` — previously ignored both
  settings and used `ggsave()`'s own defaults, so `plot_dpi` governed only 12 of
  the 17 PNGs a run produces. Both arguments are additive and default to the
  previous behaviour.

## Bug fixes

* `annual_precip.png` rendered in ggplot2's default grey theme and at a hardcoded
  300 dpi, silently ignoring `plot_dpi` and `plot_device`. It bypassed the shared
  export path entirely; it no longer does.

* `precip_conditional_correlations.png` was written at a fixed 8 × 8.5 in
  regardless of how many panels it actually drew. Its size was read from fields
  that only `facet_wrap()` populates, and that figure uses `facet_grid()`, so its
  three regime rows were being squeezed into a two-row canvas. Figure size is now
  derived from the panel grid a plot actually renders.

* The significance threshold in `warm_annual_wavelet.png` is now drawn across the
  full period axis, matching `obs_power_spectra.png`. The panel was being given
  the cone-of-influence-masked threshold, which is `NA` at every scale the record
  is too short to test — on a 20-year annual series that is everything beyond
  roughly 3.5 years, so the red dashed line appeared as a stub over the shortest
  periods and then vanished. It now receives the unmasked curve, as the observed
  power-spectrum plot already did.

  Filtering is unchanged: the masked curve remains the one
  `identify_significant_peaks()` and the testable-scale count read, because there
  an `NA` means "not testable" rather than "not significant". Note that the
  plotted line is therefore the permissive threshold, drawn over scales the cone
  of influence cannot actually test — the same convention `obs_power_spectra.png`
  uses.

## Other changes

* Figure dimensions are derived from the panel grid for every output rather than
  set per call site. Sizes remain internal — there are no width or height
  arguments. Most geometry is unchanged: the family used by the faceted
  observed-vs-simulated diagnostics reproduces the previous arithmetic exactly.
  What does change: figures carrying a title are 0.6 in taller, since the old
  allowance was sized for neither the title nor the two-line subtitles three of
  them use; the three WARM figures and `warm_annual_stats.png` are 8 × 4.5 in
  (from 8 × 5); `obs_power_spectra.png` is 9 × 4.5 in (from 8 × 5); and
  `annual_precip.png` is 8 × 4.5 in (from 8 × 4).

  A faceted figure is also now sized for the panels it actually draws rather than
  the grid that was requested. With the usual four variables these agree; with
  fewer, a plot that previously got a full-size canvas holding one panel and
  three empty slots now gets a canvas that fits it.

* Every write now sets `units`, `dpi` and background colour explicitly instead of
  leaving them to `ggsave()`'s defaults.

* Line weights, point sizes and alpha values are shared constants grouped by geom
  family, each role carrying its full specification. The observed series is now
  drawn at one weight in every figure — it was previously 0.9 in one module, 1.25
  in another, and the device default in a third — and is deliberately the
  heaviest line in any plot that draws one, so it reads clearly against the
  ensemble behind it. In `annual_precip.png` and `warm_annual_precip.png` it
  previously competed with the simulated traces rather than standing out from
  them.

  Ensemble traces are split into two roles by how they are coloured, because the
  alpha that reads correctly depends on it: monochrome grey traces are drawn
  lighter (0.30) than hue-scaled per-realization traces (0.70). A single value
  cannot serve both — grey at the hue scale's alpha reads as a dark mass, and
  hue-scaled lines at grey's alpha are close to invisible.

* The monthly-pattern legend used `legend.position = c(1, 0)`, superseded in
  ggplot2 3.5.0 by `legend.position = "inside"`.

# weathergenr 1.4.0

## New features

* `apply_climate_perturbations()` gains `temp_range_factor`, a monthly
  multiplier on the diurnal temperature range applied about its own midpoint, so
  `temp_min` and `temp_max` move apart or together without shifting the daily
  mean. It accepts a length-12 vector or an `n_years x 12` matrix like the other
  factors, and ramps with `temp_transient`. `NULL` (the default) reproduces the
  previous behaviour exactly.

  The gap it fills: `temp_delta` is added to `temp`, `temp_min` and `temp_max`
  alike, so warming on its own left the diurnal range untouched — and the
  Hargreaves PET this package computes goes as the square root of that range.
  Warming alone therefore moves PET by only about 2.3% per degree, and the range
  decides how much more. At +4 degC, PET rises 3.7% with the range a tenth
  smaller, 9.3% with it fixed, and 14.7% with it a tenth larger. Holding it
  fixed was a choice for the middle of a fourfold span, not a neutral default,
  and it matters wherever PET drives the water balance.

  There is no defensible non-null default, because observed and projected
  diurnal-range trends vary by region and season; supply values for the region
  being studied, or vary the factor as a stress-test dimension in its own right.

## Other changes

* Documented why the WARM pool's peak-matching criterion often constrains
  nothing, following the significance fix above. It is length-adaptive by
  design: the cone of influence leaves only the shorter periods testable, so on
  a 20-year annual series no observed peak qualifies and every candidate passes,
  while a longer record activates it. `peak_match_frac_min = 1.0` therefore
  reads as maximally strict but is only as strict as the record allows — check
  `diagnostics$spectral_diag$n_sig_peaks_found` before concluding it did
  anything.

  Matching the strongest peaks regardless of significance was considered and
  rejected. It would not be redundant — measured against the spectral
  correlation it discriminates on something largely independent (r = 0.11) — but
  on the packaged 20-year record the two strongest peaks sit at 5.8 and 16.5
  years, and requiring simulated traces to reproduce a 16.5-year cycle seen in
  20 years of data selects for a sampling artefact.

## Bug fixes

* `dry_spell_factor` and `wet_spell_factor` now mean what their names say: each
  is a multiplier on the **mean spell length**, so 2 doubles the mean run of dry
  (or wet) days. Previously they divided individual off-diagonal transition
  probabilities and renormalised, which moved spell length in the right
  direction by an unpredictable amount — on the packaged record a
  `dry_spell_factor` of 3 lengthened mean dry spells by 1.46x, and 1.5 by 1.11x.

  `wet_spell_factor` also now slows the return to dry from **both** non-dry
  states. Previously it adjusted only the wet row, so an extreme day ended a wet
  spell at the historical rate while a merely wet one did not, and the factor
  did progressively less the wetter the month.

  Both default to 1, which is a no-op, so runs that do not set them are
  unaffected. Runs that do set them will change, and a given factor now produces
  a larger effect than before. Note the factors are not independent of wet-day
  frequency: lengthening dry spells lowers the wet-day fraction, which is a
  property of the Markov chain rather than an artefact — adjust both together to
  hold it roughly fixed.

* `apply_climate_perturbations(seed = )` no longer leaves the caller's random
  number stream reseeded. It was the only seeded entry point in the package that
  did not restore `.Random.seed` on exit, so calling it silently changed every
  random draw the caller made afterwards. If the caller had never drawn a random
  number, no `.Random.seed` is left behind either.

* `apply_climate_perturbations()`'s validation errors named arguments the
  function does not have — `'climate.data'`, `'sim.dates'`,
  `'change.factor.precip.mean'` and others, left over from before a rename. They
  now name the actual arguments (`data`, `date`, `precip_mean_factor`, …). The
  latitude error also names both accepted column names, `lat` and `y`, rather
  than only `y`.

* `generate_weather()` could fail with `set.seed(NA)` partway through a run.
  Derived seeds were built by integer addition onto a base drawn from
  `sample.int(.Machine$integer.max, 1L)`, which overflows to `NA` when the base
  falls within the offset of the maximum. The arithmetic now wraps instead. Every
  seed that did not overflow is unchanged, so simulated output is identical.
  `simulate_warm()` had the same pattern and is fixed with it.

* `generate_weather(save_plots = FALSE)` still built the three WARM filtering
  diagnostic panels over every selected realization, skipping only the write.
  `save_plots` is now passed through, so the plots are not built either.

* One `.log()` call in `generate_weather()` ignored `verbose`, so a failure to
  write the WARM filtering plots was reported even in a quiet run.

* `simulate_warm()` rejected list-form `components` whenever `n` differed from
  the observed series length — that is, whenever the simulation length differed
  from the record, which is the normal case. The list branch validated component
  lengths against `n` (the simulated length) where the matrix and data.frame
  branches correctly used the observed length. Matrix input returning a result
  where the identical list input errored is now fixed, and the error message
  names the length it actually requires.

* `simulate_warm()`'s `bypass_n` documentation said the default was 25; it is
  15.

* `estimate_monthly_markov_probs()` gains `wet_q` and `extreme_q`. Its fallback
  for a month whose supplied wet or extreme threshold is not finite derived a
  pooled threshold at fixed quantiles of 0.2 and 0.8, regardless of the wet and
  extreme definitions the rest of the run was using, so a degenerate month could
  silently switch definition. The new arguments default to the previous fixed
  values, and `generate_weather()` now threads its own `wet_q` / `extreme_q`
  through. Behaviour is unchanged for the package's own pipeline, which fills
  such thresholds with the configured quantiles before this point — the fallback
  matters to direct callers of the exported function.

* `filter_warm_bounds_defaults()`'s `tail_tol_log` default changes from
  `log(1.03)` to `log(1.25)`. The 3% figure matched the `mean` and `sd` bounds,
  but a mean and a standard deviation are stable statistics whereas tail mass —
  computed over the handful of points beyond a quantile of a ~20-year record —
  is not, and its natural spread across candidates drawn from the fitted WARM
  model is tens of percent.

  At `log(1.03)` the tail criterion admitted 0.4% of candidates. It was
  therefore always the lowest pass rate, always the criterion the relaxation
  loop selected, and because that loop stops as soon as the pool reaches
  `n_select`, the loop rather than the default was setting the operative bound.
  At `generate_weather()`'s own defaults that left a pool of 24 out of 5,000 for
  20 realizations, so the documented filter-then-rank design was in practice
  relax-until-just-enough and the ranking had almost nothing to choose between.

  `log(1.25)` still rejects roughly 97% of candidates, so it remains a real
  constraint, and it needs no relaxation on the packaged example — pool 157 of
  5,000. Realization selection changes, so simulated output changes with it.
  On the packaged fixture 26 of 40 evaluation error metrics improved and 14
  worsened, median −7.9%, consistent with the ranking having a genuine pool to
  choose from; with two realizations per scenario that is directional evidence
  rather than a significant result.

  Pass `filter_bounds = list(tail_tol_log = log(1.03))` to `filter_warm_pool()`,
  or the same entry via `warm_filter_bounds` to `generate_weather()`, to restore
  the previous value.

* WARM's variance matching pinned the mean of every realization it corrected.
  Rescaling a simulated trace onto the observed standard deviation also
  re-centred it on the observed *mean*, replacing the trace's own mean rather
  than preserving it. On the packaged fixture this affected 60% of a 3,000-trace
  pool, putting a point mass exactly on the observed mean: those realizations
  carried no spread in their own 20-year mean, so `filter_warm_pool()`'s mean
  criterion was structurally unable to reject them. The practical effect on
  filtering was small — that criterion admitted 99.5% of candidates before the
  fix and 98.4% after — but the ensemble was misrepresenting its own
  variability. Traces are now rescaled about their own mean, restoring roughly
  58% more spread in the simulated multi-decadal mean, which is the
  low-frequency variability WARM exists to reproduce rather than noise to be
  removed.

* WARM compared a *population* standard deviation against `stats::sd()` targets.
  The internal column-sd helper omitted the Bessel correction, so every
  corrected trace landed at `target_sd * sqrt(n / (n - 1))` rather than on the
  target — a one-sided +1.3% at `n = 40`, against a `filter_warm_pool()` sd
  tolerance of 3%. Both sides now use the sample standard deviation.

  Both fixes change simulated output: realization selection shifts, so the daily
  analogue dates and every evaluation statistic downstream shift with it. A run
  reproduced from a given seed will not match one from an earlier version. Pool
  sizes and relaxation behaviour are unaffected (the sd criterion moves 72.4% to
  71.2% on the packaged fixture). Evaluation metrics move in both directions;
  `mae_mean_precip` in particular worsens, because the previous ensemble matched
  the observed mean artificially well.

* `read_netcdf()` ignored the CF `calendar` attribute and always decoded the
  time axis as a proleptic Gregorian calendar. Files on a `noleap` / `365_day`
  calendar — including every file `write_netcdf()` produces with its default
  `calendar = "noleap"` — therefore read back one day short per Gregorian leap
  year crossed: a run written as 2020-01-01 to 2021-12-31 returned
  2020-01-01 to 2021-12-30, drifting 25 days over a century. `read_netcdf()`
  now decodes `noleap` / `365_day` axes on a 365-day calendar, so a
  `write_netcdf()` / `read_netcdf()` round trip is exact.

  CF-aware readers such as `xarray` and the Wflow toolchain already honoured
  the attribute, so files written by earlier versions were correct on disk and
  need no migration — only `read_netcdf()` was misreading them. Any workflow
  that compensated for the old off-by-one when reading `noleap` files must drop
  that correction. Files on `standard`, `gregorian`, `proleptic_gregorian`, or
  `julian` calendars, and files with no `calendar` attribute, decode exactly as
  before.

* `read_netcdf()` now raises an informative error for `360_day`, `all_leap`, and
  `366_day` axes, whose dates (30 February; 29 February in a common year) have
  no representation in R's `Date` class and were previously decoded to silently
  wrong dates. An unrecognised calendar name warns and falls back to the
  standard decode rather than failing the read.

## Other changes

* `write_netcdf()` gains an optional `dates` argument. When supplied it is
  checked against the number of time steps, `origin_date`, and `calendar`, and
  the time coordinate is computed from it rather than assumed contiguous. This
  catches the case where `origin_date` is given as a conventional epoch such as
  `"1970-01-01"` instead of the first date of the series — the time coordinate
  is written as `0:(nt - 1)`, so that silently relocated the whole series.
  The argument is optional and the default behaviour is unchanged.

* `write_netcdf()` now rejects a `calendar` value it cannot write correctly,
  instead of accepting arbitrary text. The `calendar` attribute is also no
  longer written inside `try()`, so failing to record it is an error rather
  than a silently mislabelled file.

* `write_netcdf()`'s `origin_date` documentation now states that it must be the
  date of the first time step. The previous example suggested `"1970-01-01"`,
  which is wrong for any series not starting on that date.

* `generate_weather()`'s `warm_signif` was documented as the "wavelet
  significance level used to retain low-frequency components in WARM". It does
  not retain anything: every MODWT MRA component is modelled, significant or
  not, and the `cwt_to_modwt_map` a significance-based selection would use is
  computed for diagnostics and never consumed. Anyone who tuned `warm_signif`
  expecting it to control which components are simulated was misled. Its two
  real effects are now documented — the peak-significance threshold in
  `filter_warm_pool()`, and a conditional MODWT level bump that rarely fires on
  20-40 year records. Behaviour is unchanged.

  `generate_weather()` also gains a Details section setting out how the annual
  generator relates to the wavelet autoregressive model of Kwon et al. (2007)
  and Steinschneider & Brown (2013). The original retains globally-significant
  CWT scales as signals and models the remainder as a residual term; this
  package uses an exactly-additive MODWT multiresolution analysis, which leaves
  no residual to relegate a non-significant band to, and on records of this
  length the red-noise significance test detects nothing at any conventional
  level. Both departures are deliberate and are now stated rather than implied.

# weathergenr 1.3.1

## Other changes

* Console log lines now read `HH:MM:SS - tag - message` (e.g.
  `14:51:42 - resample - Processing realization: 1/2`) instead of
  `[YYYY-MM-DD HH:MM:SS] [TAG] message`. The date is dropped and the phase tag
  is lower-cased, so `weathergenr` output interleaves cleanly with the
  `blueearth_cst` toolkit's logging. Message text is unchanged; anything parsing
  the old prefix needs updating.

# weathergenr 1.3.0

## Bug fixes

* `generate_weather()`'s `relax_priority` argument had no effect. It was
  accepted, validated and documented as being forwarded to
  [filter_warm_pool()] as `relax_order`, but the call passed a hard-coded
  ordering and discarded the argument. It is now forwarded as documented.
  `relax_priority` breaks ties between WARM filters that share the lowest pass
  rate, so a run where one filter is strictly most restrictive is unaffected.
* `run_weather_generator()` now forwards `config$relax_priority`, which was
  missing from the arguments it passes on to `generate_weather()`.
* `run_weather_generator()` no longer fails when `config` omits
  `warm_filter_bounds` or `relax_priority`. An absent entry arrives as `NULL`
  rather than as a missing argument, so the formal defaults never applied and
  validation rejected the `NULL`; both now fall back to their documented
  defaults.
* `generate_weather()` validates `relax_priority` against its documented
  contract (each of `mean`, `sd`, `tail_low`, `tail_high`, `wavelet` exactly
  once) on entry. `filter_warm_pool()` already enforced this, but only after the
  wavelet analysis and pool simulation had run.

* `evaluate_weather_generator()`'s fit metrics were computed against the wrong
  observed statistics. `.summarize_realization_fit()` filtered the *simulated*
  side to one statistic but joined the *observed* side unfiltered, so `stat` and
  `type` were missing from the join keys. Observed `stats.season` holds a `mean`,
  an `sd` and a `skewness` row per grid-month-variable, and observed `wetdry`
  holds a `days` and a `spells` row per grid-month-stat -- so every simulated
  mean was differenced against all three observed statistics, and every wet-day
  count against spell lengths too, with the results averaged together.

  Four of the six metrics in `fit_summary` were affected, and with them
  `overall_score` and `rank`. The errors were large and inflated the apparent
  error: on the bundled fixture `mae_mean_precip` falls from 2.38 to 0.25,
  `mae_days_Wet` from 9.32 to 0.52, and `mae_mean_temp` from 17.11 to 0.04
  (temperature means had been compared against temperature skewness). Realization
  *ranking* also changes in some cases. `mae_cor_crossgrid` and
  `mae_cor_intervariable` were already correct and are unchanged.

  **Any previously recorded fit metrics or realization rankings should be
  regenerated.** Simulated output itself is unaffected: `generate_weather()` is
  byte-identical, as the end-to-end baseline confirms -- only the evaluation
  summary changes.

* `evaluate_weather_generator()` no longer emits tidyselect deprecation warnings
  from `tidyr::separate()`, `tidyr::unite()` and `dplyr::mutate(.after =)`. Six
  remaining `.data$` references in tidyselect arguments now pass column names as
  strings, completing the change begun for `select()`/`rename()`/`pivot_*()`.

## New features

* `evaluate_weather_generator()` gains `plot_dpi` (default `300`) and
  `plot_device` (default `NULL`), forwarded to [ggplot2::ggsave()]. Writing the
  diagnostic PNGs dominates an evaluation run that saves plots — 4.7 s of 5.7 s
  on the bundled fixture — and that cost scales with resolution: `plot_dpi = 150`
  takes 5.0 s, `plot_dpi = 96` takes 4.3 s, against 1.1 s for
  `save_plots = FALSE`. `plot_device` accepts a faster device such as
  `ragg::agg_png`; `ragg` is not a dependency, so pass the function only if it is
  installed. Both are reachable through `run_weather_generator()` as
  `config$plot_dpi` and `config$plot_device`. Defaults reproduce the previous
  output exactly.

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
* `generate_weather(parallel = TRUE)` no longer starts more workers than there
  are realizations. Realizations are the only unit of parallelism in the daily
  resampling — it carries state across days and years, so it cannot be split
  within a realization — and every extra worker idled while still paying PSOCK
  startup and memory. With `n_cores = 10` and `n_realizations = 2` a run drops
  from 21.9 s to 12.0 s, output unchanged. `filter_warm_pool()` still uses the
  full `n_cores`, since its work does divide further. When the cap leaves a
  single worker the daily loop runs sequentially, which produces the same
  output and skips cluster setup entirely.
* The daily resampling loop no longer rebuilds a `"month.day"` string key per
  simulated day, nor coerces the month to character to read monthly means. Both
  are now integer lookups.
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
  `save_plots = FALSE`): generation 8.7 s to 4.9 s, evaluation 3.7 s to 1.2 s,
  together 12.4 s to 6.1 s (-51%). Timings are the minimum of four runs in one
  session, since machine load shifts absolute figures between sessions.
  **Outputs are unchanged**: for a fixed `seed` these changes are bit-identical,
  verified by the end-to-end baseline gate over both water-year and
  calendar-year scenarios.

## Bug fixes (parallel execution)

* `generate_weather(parallel = TRUE)` no longer holds an idle PSOCK cluster
  while `filter_warm_pool()` runs. `filter_warm_pool()` spawns its own cluster,
  once per relaxation iteration, so the two overlapped: a run with
  `n_cores = 3` peaked at 8 R processes where 4 were expected. The cluster is
  now created immediately before the daily disaggregation that uses it. Output
  is unchanged, bit for bit. On a machine with cores to spare this does not
  change runtime; it matters at the default
  `n_cores = parallel::detectCores() - 1`, where the overlap oversubscribed
  every core.

## Breaking changes

* `generate_weather(parallel = TRUE)` now reproduces `parallel = FALSE` for the
  same `seed`. It previously did not. `parallel::clusterSetRNGStream()` switches
  workers to the L'Ecuyer-CMRG generator, and `resample_weather_dates()` called
  `set.seed()` without naming a generator, so the seed was applied to whichever
  one happened to be active: the same integer produced one stream in the master
  and a different one on a worker. Seeding now pins Mersenne-Twister, so a seed
  denotes one stream everywhere.

  **Sequential output is unchanged** -- the master already used
  Mersenne-Twister, so the pin is a no-op there, as the end-to-end baseline
  confirms. **Output from previous parallel runs is not reproducible under this
  version**; re-run rather than mixing old and new parallel ensembles. Parallel
  results were already stable across worker counts and remain so.

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
