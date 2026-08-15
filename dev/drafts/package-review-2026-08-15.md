# weathergenr 1.3.1 — architecture, consistency, scientific coherence review

Backing document for the `review-2026-08-15` board origin. Every finding below
is tracked; this draft is shared by all of them, so it is promoted on closure
rather than deleted with any one item.

| Finding | Item | Area | Queue |
|---|---|---|---|
| C5 NetCDF noleap round trip | `t2608151253` | io | 1 |
| C9 + C10 `.variance_match_matrix()` | `t2608151253a` | warm | 2 |
| C1 `warm_signif` component selection | `t2608151253b` | warm | 3 |
| A3 equivalence tests for the two forks | `t2608151253c` | tests | 4 |
| A1 perturbation not wired into the pipeline | `t2608151254` | api | 5 |
| C3 random observed benchmark window | `t2608151254a` | warm | 6 |
| C4 PET diurnal-range response | `t2608151254b` | perturbations | 7 |
| C2 Cholesky ordering dependence | `t2608151254c` | warm | 8 |
| C6 Markov spell-length factors | `t2608151254d` | resample | 9 |
| C7 hardcoded thresholds, `dirichlet_alpha` | `t2608151254e` | resample | 10 |
| C8 + B9 `simulate_warm()` argument defects | `t2608151254f` | warm | 11 |
| B3 + B4 + B5 `generate_weather()` defects | `t2608151254g` | api | 12 |
| B1 + B7 `apply_climate_perturbations()` | `t2608151254h` | perturbations | 13 |
| B11 `.Rbuildignore` for `vignettes/.quarto` | `t2608151254i` | ci | 14 |
| B8 duplicated `test-wavelet.R` | `t2608151254j` | tests | 15 |
| A4 export surface | `t2608151254k` | api | 16 |
| A2 god modules | `t2608151255` | architecture | *watch* |
| B2 `set.seed()` discipline | `t2608151255a` | warm | *watch* |
| B6 argument-name drift | `t2608151255b` | api | *watch* |
| A5 internal naming conventions | `t2608151255c` | architecture | *watch* |
| B10 lint / glue opacity | recorded on `t2608061641` | ci | — |

B10 has no item of its own. The pre-existing lint item covers clearing the
findings; the observation that `.log()`'s glue interpolation makes lintr
structurally unable to see into `log_filter_iteration()` and
`log_final_summary()` is recorded as a correction on that note, because it
changes what the note's checklist should do rather than adding work.

Evidence base: full read of `generator.R`, `io_netcdf.R`, `wavelet_warm.R`,
`climate_perturbations.R`, `pet.R`; targeted reads of `resample.R`,
`warm_filtering.R`, `wavelet_modwt.R`, `quantile_mapping.R`. Plus
`Rscript tools/lint.R` (20 warnings), `check_only()` (0 errors, 2 warnings,
0 notes, 3m34s), and four empirical probes run against the working tree.

---

## Verdict

The numerical core is careful work. The Cholesky whitening algebra in
`simulate_warm()` is correct, the FAO-56 Hargreaves implementation is correct
(including the 0.408 MJ→mm conversion), the Markov lag guards are right, the
performance comments explain real decisions, and `R CMD check` is clean of
errors and notes.

The problems are not in the mathematics. They are in three places: the
pipeline stops one stage short of the package's stated purpose; several
documented knobs do something different from what they say; and two numerical
routines have been hand-forked for speed with nothing testing that the forks
still agree with their originals.

Five findings are worth acting on: C5, C9/C10, C1, A3 and A1 below.

---

## A. Architecture

### A1. The third component is not wired into the pipeline

`README.md` describes three coupled components. Components 1 and 2 (WARM,
Markov+KNN) are coupled inside `generate_weather()`. Component 3 — climate
perturbation, the reason for a stress-testing package — is not reachable from
either entry point. `run_weather_generator()` chains generate → prepare-eval →
evaluate. `apply_climate_perturbations()`, `adjust_precipitation_qm()`,
`calculate_monthly_pet()` and `write_netcdf()` appear nowhere in it.

This is more than a missing call. The interfaces do not meet:

- `generate_weather()` returns resampled *dates* and a date axis. No values.
- `apply_climate_perturbations()` requires per-cell data frames with `precip`,
  `temp`, `temp_min`, `temp_max` (`climate_perturbations.R:150`).
- `generate_weather()` validates and simulates only `precip` and `temp`
  (`generator.R:208`).

So `temp_min`/`temp_max` — required by stage 3 and by Hargreaves PET — are
never produced by stage 2. Every user must hand-assemble the join, as
`vignettes/climate_perturbation.qmd` does. That is the largest structural gap
in the package.

### A2. Two god modules

`warm_filtering.R` — 1903 lines, 14 functions. `evaluate_generator.R` — 1836
lines, 17 functions, with `evaluate_weather_generator()` alone spanning
lines 102–360. Together 28% of the package. Both are past the point where the
`*_plots.R` separation convention (which is otherwise held consistently) is
enough on its own.

### A3. Two hand-forked numeric routines, no equivalence test

- `.compute_gws_batch()` (`warm_filtering.R:1012`) reimplements the CWT from
  `analyze_wavelet_spectrum()`. Its comments assert the equivalence —
  *"CWT grid (identical to analyze_wavelet_spectrum internals)"*,
  *"mirrors analyze_wavelet_spectrum exactly"* — but nothing enforces it.
- `.knn_draw_one_rank()` (`resample.R:~355`) reimplements `knn_sample()`'s
  distance computation, partial-sort neighbour selection and rank
  probabilities.

Both forks exist for speed, and both are reasonable optimisations. But each
fork *is* the scientific criterion in its path: the spectral acceptance test
for WARM traces, and daily analogue selection. Neither has a test comparing it
against the canonical routine — I grepped `tests/testthat/` for both. The
failure mode is silent: someone tunes `dj`, `s0`, `k0`, the padding, or the
detrend in `analyze_wavelet_spectrum()`, the fork stays behind, the numbers
move, and the full suite stays green.

This is cheap to close: one `test_that()` per fork asserting
`expect_equal(fork(x), canonical(x))` on a small fixture.

### A4. 48 exports for 2 entry points

`NAMESPACE` exports helpers with no plausible external audience:
`criteria_string_compact`, `generate_symmetric_dummy_points`,
`get_result_index`, `match_transition_positions`. Several carry
`@keywords internal` and are exported anyway. (The wavelet utilities —
`fill_nearest`, `gws_regrid`, `extract_signif_curve` — are general-purpose
enough to defend; the four above are not.)

This matters more here than in most packages: per `AGENTS.md`, `blueearth_cst`
pins weathergenr by git tag, so every export is a contract you cannot quietly
withdraw.

### A5. Three competing internal-naming conventions

`.`-prefixed (`.markov_month_probs`, `.date_parts`), bare-unexported
(`fit_monthly_distributions`, `compute_target_parameters`,
`apply_quantile_mapping`, `find_local_maxima`, `log_filter_iteration`), and
`@keywords internal` + exported (`criteria_string_compact`,
`get_result_index`). All three appear within single files.

---

## B. Consistency

### B1. Four RNG save/restore idioms, one of them absent

| Site | Behaviour |
|---|---|
| `quantile_mapping.R:260` | saves, restores, and `rm()`s if `.Random.seed` did not exist — the correct one |
| `generator.R:248`, `wavelet_warm.R:143`, `warm_filtering.R:195`, `evaluate_generator.R:156`, `wavelet_cwt.R:576` | save/restore, but leak a `.Random.seed` into a session that had none |
| `climate_perturbations.R:193` | `set.seed()` with **no restore** — permanently reseeds the caller's session |

`apply_climate_perturbations()` is a user-facing entry point. Calling it with a
seed silently hijacks the caller's RNG stream for everything that follows.

### B2. `.set_seed_fixed_kind()` is used in one file

`AGENTS.md` marks it IMPORTANT: never a bare `set.seed()` in code that can run
on a worker. It is called only from `resample.R` (2 sites). Six other
`set.seed()` calls use the bare form.

None of the six currently runs on a worker — I checked `.compute_gws_batch()`'s
`parLapplyLB` body, which is fully deterministic — so this is not a live bug.
But `simulate_warm()` is one `parallel = TRUE` refactor away from being a
worker call, and the rule is currently held by convention rather than by
construction.

### B3. Seed arithmetic can overflow to `NA`

`warm_seed`/`daily_seed` are drawn from `sample.int(.Machine$integer.max, 1L)`
(`generator.R:261-262`), then incremented: `warm_seed + 1L` (`:434`),
`daily_seed + n` (`:508`, `:536`), `seed + k * 1000L`
(`wavelet_warm.R:400`). Integer overflow yields `NA`; `set.seed(NA)` errors.
Probability is ~1e-5 per run, but it fails hard, deep into a long run, and only
under a specific seed — the worst kind to debug. `seed %% 1e6L + k` fixes it.

### B4. `save_plots = FALSE` still builds the plots

`generator.R:437` passes `make_plots = TRUE` to `filter_warm_pool()`
unconditionally. Only the `ggsave()` is skipped. The three diagnostic panels
over the selected realizations are constructed regardless.

### B5. `.log()` at `generator.R:457` omits `verbose = verbose`

Its sibling at `:401` includes it. The WARM-plot failure message therefore
ignores the verbosity setting.

### B6. Argument names drift across the internal boundary

`out_dir` (generate) vs `output_dir` (evaluate); `warm_filter_bounds` →
`filter_bounds`; `relax_priority` → `relax_order`. `run_weather_generator()`
absorbs the renames, so anyone calling the pieces directly hits them.

### B7. Error messages name arguments that do not exist

`apply_climate_perturbations()` reports `'climate.data'`, `'sim.dates'`,
`'change.factor.precip.mean'`, `'change.factor.precip.variance'`,
`'change.factor.temp.mean'` (`:123-130`) for arguments actually named `data`,
`date`, `precip_mean_factor`, `precip_var_factor`, `temp_delta`. Line 120
marks this deliberate ("legacy messages preserved"), but a user who passes
`data = NULL` is told to fix `'climate.data'`, which is not in the signature.

Related: the latitude check (`:145-148`) accepts `lat` or `y`, but the message
names only `'y'`.

### B8. `test-wavelet.R` is a strict subset of `test-warm.R`

All 8 `test_that()` titles in `test-wavelet.R` are duplicated verbatim among
`test-warm.R`'s 24. Its header comment cites `R/wavelet.R`, a file that no
longer exists — leftover from the `wavelet_cwt.R`/`wavelet_warm.R` split.

### B9. `simulate_warm()` doc says `bypass_n` default is 25; the formal is `15L`

`wavelet_warm.R:55` vs `:103`.

### B10. Static analysis is degraded by `.log()`'s glue interpolation

Lint is 20 warnings, nothing severe. But the composition is telling: 12 are
`.Random.seed` binding noise (one `utils::globalVariables()` entry in
`globals.R` clears them all), and 6 are *false* "unused local" reports, all in
`log_filter_iteration()` and `log_final_summary()` — the variables are used, but
only inside `.log("{n_pass}")` glue strings, which lintr cannot see through.
Narrow in scope (two functions), but it means genuine dead locals in those two
would not be flagged either.

### B11. `check_only()` currently exits non-zero on master

0 errors, 0 notes, 2 warnings — both because `vignettes/.quarto/_freeze/` (18
cached artifacts) ships into the built package and `inst/doc` does not exist.
Per `AGENTS.md` the release gate is `check_only(build_vignettes = TRUE)`, and
the wrapper errors on WARNINGs. Adding `^vignettes/\.quarto$` to
`.Rbuildignore` clears both.

---

## C. Scientific coherence

### C1. `warm_signif` does not select components — the doc says it does

`generator.R:70` documents it as *"Wavelet significance level used to retain
low-frequency components in WARM."* It does not retain anything.

Tracing it: `analyze_wavelet_additive()` runs the CWT purely to choose
`n_levels` (`wavelet_modwt.R:527-546`). It then computes `cwt_to_modwt_map` —
which MODWT levels the significant CWT periods imply (`:557-577`) — returns it,
and **nothing consumes it**. `simulate_warm()` receives
`warm_analysis$components`, the full MODWT MRA basis, and fits an ARMA to every
component, significant or not.

This departs from Steinschneider & Brown (2013) and Kwon et al. (2007), where
AR models are fit to the *significant* wavelet bands and the remainder is
treated as noise. Here every band is modelled the same way.

I have not measured what that costs in practice, and I would not assume it is
large — an ARMA fit to a weak, non-significant MRA band reproduces that band's
own small variance, so it may contribute little rather than something spurious.
What is verified is narrower and enough on its own: the documented behaviour
does not exist, the computed mapping is dead code, and the method differs from
the paper the README cites in a way no documentation discloses.

Compounding it, `warm_signif` is reused for a third job: `filter_warm_pool`'s
`wavelet_args$signif_level` (`generator.R:439`). One knob, three meanings.

Fix either way, but pick one: gate component selection on `cwt_to_modwt_map`,
or correct the documentation and say plainly that the significance test is
diagnostic. Given the README cites S&B 2013 as the method, this deserves an
explicit statement.

### C2. Cholesky whitening is order-dependent and destroys scale separation

The algebra is right — I checked it. `Z = X R⁻¹` with `t(R) %*% R = Σ` gives
`cov(Z) = I`; resumming via `rowSums(R)` correctly reconstitutes
`rowSums(Z_sim %*% R)` (`wavelet_warm.R:496-500`).

The concern is methodological. `chol()` is triangular, so whitening is
*sequential*: Z₁ ∝ D₁, Z₂ mixes D₁ and D₂, …, Z_J mixes everything. The fitted
ARMAs are therefore no longer scale-specific — which is the entire premise of
WARM — and the result depends on component ordering (D1…DJ, S), a convention
rather than a physical fact. The code half-acknowledges this at `:407-412`
("no longer smooth regardless of its name").

A symmetric whitening (`Σ^{-1/2}` via eigendecomposition) removes the ordering
dependence at the same cost. At minimum, the header comment at `:17-20` — which
still says independent modelling *does not* preserve cross-component covariance
"unless you implement coupling externally" — now contradicts the Cholesky path
and should be updated.

### C3. The WARM filter benchmarks against a single random window of the observations

`filter_warm_pool()` (`warm_filtering.R:209-214`): when the observed record is
longer than the simulated length, one random contiguous window is drawn, and
*all* acceptance statistics — mean, sd, tail mass, global wavelet spectrum —
are computed against that window rather than the full record.

Symmetrically (`:222-228`), when simulations are longer, each candidate trace is
scored on its own random window — but `selected` returns the **full** trace
(`:561`).

Two consequences:

1. The acceptance criterion moves with the seed. Re-running with a different
   seed changes which traces are admissible, not just which are drawn.
2. A trace admitted on a random 40-year slice is delivered as a 100-year trace
   whose full-length statistics were never tested.

The intent — compare like lengths — is sound. The single-draw implementation is
not. A median over several windows would fix (1); (2) needs at least a
documented statement in `@details`, since nothing in "filtering thresholds"
prepares a user for it.

### C4. Diurnal range is invariant under warming, so PET barely responds

`apply_climate_perturbations()` adds the same monthly `temp_delta` to `temp`,
`temp_min` and `temp_max` (`:465-467`). So `temp_max - temp_min` is unchanged by
construction, and Hargreaves —

    PET = 0.0023 · Ra · (T + 17.8) · √(ΔT)

— responds to warming only through `(T + 17.8)`. At T = 15 °C, +4 °C raises PET
by ~12%.

For a package whose stated purpose is bottom-up water-resources vulnerability
assessment, PET drives the drought limb of the response surface, and this
systematically damps it. Holding DTR fixed is a defensible choice (DTR change
*is* genuinely uncertain), but it is undocumented and there is no knob to
perturb it. A `temp_range_factor` argument, or an explicit `@details`
statement, would close the gap.

The Hargreaves implementation itself is correct FAO-56 — I verified the
`0.408` MJ m⁻² d⁻¹ → mm d⁻¹ conversion and the Ra formulation.

### C5. weathergenr cannot read back the NetCDF it writes

`write_netcdf()` writes `vals = 0:(nt-1)` against `units = "days since
<origin_date>"` and stamps `calendar = "noleap"` (`io_netcdf.R:546-553, 650`).
`read_netcdf()` ignores the calendar attribute entirely and decodes with
proleptic-Gregorian arithmetic (`.parse_time_to_date`, `:140-149`).

Verified round trip — 730 noleap days from a 2020-01-01 origin:

```
intended:              2020-01-01 .. 2021-12-31
read_netcdf() returns: 2020-01-01 .. 2021-12-30
```

One day lost per leap year — roughly 24 days over a century-long stress-test
run. CF-compliant consumers (xarray, Wflow) that *do* honour `calendar` get the
intended dates, so weathergenr's own reader is the outlier; but either way the
two halves of the package's I/O disagree, and a user who writes then reads
their own output gets silently shifted dates.

Two further gaps in the same function:

- `write_netcdf()` takes no date vector, so it *cannot* check that
  `origin_date` is the series start. The roxygen example suggests
  `"1970-01-01"` (`:391`), which against `vals = 0:(nt-1)` would silently place
  a 2020-start run in 1970.
- `calendar` is unvalidated free text. Passing `"standard"` mislabels noleap
  data with no complaint.

### C6. Markov spell factors are not spell-length multipliers, and are asymmetric

`dry_spell_factor[m]` divides p01 and p02 then renormalises
(`resample.R:1073-1078`); `wet_spell_factor[m]` divides p10 only (`:1080-1085`).

- They change the unconditional wet-day fraction as well as spell persistence —
  the two knobs are not orthogonal, and `dry_spell_factor = 1.5` does not mean
  50% longer dry spells.
- The extreme-state row (p20/p21/p22) is never touched, so transitions *out of*
  the extreme state keep historical persistence regardless of either factor.

The docs ("monthly dry-spell adjustment factors applied in the Markov-chain
persistence logic") are vague enough not to be wrong, but give the user no way
to calibrate. State the transformation explicitly, or convert the arguments to
target spell lengths and solve for the multiplier.

### C7. Two smaller method issues in `.markov_month_probs()`

- **Fallback ignores the configured thresholds** (`resample.R:1005-1009`). When
  a monthly threshold is non-finite, the fallback hardcodes quantiles 0.2 and
  0.8 rather than the user's `wet_q`/`extreme_q`. Rarely reached, but it makes
  a degenerate month silently disagree with the run's own definition of
  wet/extreme.
- **Dirichlet smoothing decays faster than standard** (`:1036`):
  `dirichlet_alpha / sqrt(n_transitions_m)` makes the prior *weight* shrink
  with sample size, so influence decays as ~n^-1.5 instead of the usual ~n^-1.
  With ~30 years of daily data (n ≈ 900 per month) the effective alpha is
  ~0.03 — essentially no smoothing exactly where you would want it, on rare
  extreme→extreme transitions. Undocumented, and not exposed through
  `generate_weather()`.

### C8. `simulate_warm()` rejects list-form components whenever `n != n_fit`

`wavelet_warm.R:297` tests `any(lens != n)` — the *simulated* length. It should
test `n_fit`, the observed length, as the matrix and data.frame branches
correctly do (`:275`, `:284`).

Verified against the working tree, 40 observed years, `n = 25`:

```
matrix input: 25 x 3
list input:   ERROR: All component vectors must have length 'n'.
```

Not reached from `generate_weather()` (it passes a matrix), and the existing
list-input test uses `n == n_fit`, which is why it has gone unnoticed.

### C9. `.variance_match_matrix()` silently pins the mean, making the `mean` filter inert for half the pool

The function is documented as variance matching. It also replaces the column
*mean* with `target_mean`, unconditionally, for every column it corrects
(`wavelet_warm.R:938-941`).

In `simulate_warm()`'s Tier-1 call, `target_mean = mean(obs_total)` where
`obs_total` is the anomaly reconstruction. So every corrected trace has its mean
forced to one value; `generate_weather.R:428` then adds
`warm_series_obs_mean` back, yielding traces whose mean is exactly the observed
mean. Verified — 40 years, 2000 traces:

```
traces with mean pinned:              1048 / 2000  (52%)
sd of trace means, pinned vs not:     NA (zero spread)  vs  1.7615
```

`filter_warm_pool`'s `mean` criterion (3% relative tolerance) therefore cannot
reject any of that 52%. The filter stage is provably inert for half the pool.

One scope caveat: `filter_warm_pool` scores `obs_use`/`sim_use`, which may be
random windows (C3). A window of a mean-pinned trace has a window mean that is
*not* exactly `obs_mean`, so the inertness holds cleanly only when
`n_obs0 == n_sim0` and no windowing occurs. That is the default path —
`n_years` defaults to the observed year count (`generator.R:330`) — so it is the
common case, not a corner case.

### C10. `.variance_match_matrix()` mixes population and sample sd, biasing the sd filter

`.fast_col_sd()` (`wavelet_warm.R:904`) computes a **population** sd —
`sqrt(m2 - m1²)`, no n/(n-1) correction. `.variance_match_matrix()` compares it
against, and rescales it to, targets computed with `stats::sd()`, which is the
**sample** sd. So every corrected trace lands at a sample sd of
`target × sqrt(n/(n-1))`.

Verified — 40 years, 2000 traces, `match_variance = TRUE`:

```
target (stats::sd of obs):        11.561910
sample sd of every clamped trace: 11.709206   (= target × sqrt(40/39))
traces clamped:                   1048 / 2000  (52%)
systematic bias:                  +1.274%
```

Every clamped trace lands on that value *exactly* — the correction creates a
point mass, not a distribution. Two consequences:

1. `filter_warm_pool`'s default `sd` tolerance is **3%**
   (`filter_warm_bounds_defaults()`). A systematic +1.27% bias consumes 42% of
   that budget, one-sidedly, and biases which traces survive.
2. Half the ensemble carries an annual standard deviation pinned to a single
   value. Trace-to-trace spread in low-frequency variance is part of the
   uncertainty a stress test is meant to sample; clamping removes it for the
   affected subset.

The author clearly knows the distinction — `.compute_gws_batch()` applies
`* (n1/(n1-1))` explicitly and comments it "sample variance"
(`warm_filtering.R:1060`). `.fast_col_sd()` just omits it.

Note this is a numeric-output change: it needs the baseline gate run before and
after, per `AGENTS.md`. Same for C9 — they are one function and should be
fixed together.

---

## Suggested order

1. **C5** — NetCDF calendar round trip. Silent, cumulative, crosses the
   `blueearth_cst` boundary.
2. **C9 + C10** — `.variance_match_matrix()`. One function, two defects: an
   inert `mean` filter for half the pool, and a systematic +1.27% sd bias
   against a 3% threshold. Baseline gate before and after.
3. **C1** — `warm_signif`. Either wire up `cwt_to_modwt_map` or correct the
   docs; a method claim tied to a cited paper should not be ambiguous.
4. **A3** — equivalence tests for the two forks. Cheap, and it protects
   everything above from silent regression.
5. **A1** — decide whether stage 3 belongs in `run_weather_generator()` or the
   README should stop calling the components coupled.

B11 (`.Rbuildignore`) is a one-line fix that turns the release gate green.
