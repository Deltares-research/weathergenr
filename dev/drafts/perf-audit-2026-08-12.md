---
title: Performance audit — profile of the generate/evaluate pipeline
type: draft
area: performance
status: closed
created: 2026-08-12
updated: 2026-08-12
---

# weathergenr — performance audit

> [!note] Outcome — all ten findings closed, 2026-08-12
> Implemented in `63e7adc`, `d868dc4`, `d99e2f7`, `205d0f5`, `7ac48cb`,
> `5205d9e`, `400c9ca`, `d2ade8a`, guarded by the baseline gate added in
> `694ddb7`. Every change verified bit-identical for a fixed seed except where a
> deliberate behaviour change is recorded in `NEWS.md`.
>
> **Cumulative, bundled fixture, 30 years, 3 realizations, `save_plots = FALSE`,
> minimum of four runs in one session:** generation 8.67 → 4.92 s (−43%),
> evaluation 3.73 → 1.19 s (−68%), together 12.40 → 6.11 s (**−51%**).
>
> Where the audit was wrong, and worth remembering:
> - **#3 was a near-null result.** Ranked third on projected runtime; delivered
>   ~0.02 s, because #1 had already made `compute_water_year()` ~70× faster.
>   Fixing the hottest thing first re-orders everything under it.
> - **#8 was underrated.** Filed as minor polish; delivered −20% on generation.
> - **#10's premise was wrong.** The finding assumed idle cores needed a second
>   axis of parallelism. There isn't one — the daily loop carries state across
>   days and years. The real defect was the converse: more workers were started
>   than could ever receive work.
> - **#5 mattered less for speed than for what it uncovered** — investigating it
>   surfaced that `parallel = TRUE` never reproduced `parallel = FALSE`
>   (`clusterSetRNGStream` switches the generator; `set.seed()` does not pin it).
> - **The gate found a correctness bug the audit never looked for**: fit metrics
>   were joined without `stat`/`type`, so every simulated mean was compared
>   against observed mean, sd *and* skewness. `mae_mean_temp` was overstated 430×.
>   Fixed in `34e252f`.

Ten ranked findings from a profiled end-to-end run, each labelled with its
evidence basis and whether it
changes RNG behaviour.

## How this was measured

- **Profile**: `Rprof(interval = 0.01, line.profiling = TRUE)` over
  `run_weather_generator()` on the packaged fixture
  `inst/extdata/ntoum_era5_data.nc`, with
  `vars = c("precip","temp")`, `year_start_month = 10`, `n_years = NULL`
  (all complete years), `warm_pool_size = 2000`, `n_realizations = 3`,
  `annual_knn_n = 100`, `parallel = FALSE`, `save_plots = TRUE`, `seed = 100`.
  Total wall clock **24.15 s**. R 4.6.0, Windows.
- **Microbenchmarks**: candidate replacements timed in isolation on a
  27,375-day (75-year, no-leap) date vector, with equality of results asserted
  via `identical()`.
- The profiling and benchmark scripts were one-offs in session scratch and are
  not committed. The parameters above are sufficient to reproduce the run; if
  the numbers need re-checking after a change, re-create the harness under
  `dev/scripts/`.

Every finding is labelled **RNG-neutral** (does not change the number or order of
random draws, so a fixed `seed` reproduces bit-identical output) or
**RNG-affecting**. RNG-affecting items are ranked lower regardless of speedup,
because a downstream consumer pins this package by Git tag.

## Profile summary (share of the 24.15 s run)

| Component | total % |
| --- | --- |
| `resample_weather_dates` (daily KNN + Markov loop) | 41.9 |
| `evaluate_weather_generator` | 41.9 |
| ├ `create_all_diagnostic_plots` (ggplot render + `ggsave`) | 24.7 |
| └ `.align_obs_sim_periods` | 11.3 |
| `compute_water_year` | 12.7 |
| `format.POSIXlt` (self time, all callers) | 14.0 |
| `match_transition_positions` | 13.6 |
| `filter_warm_pool` / WARM wavelet pool | < 1 |

The WARM pool spectral code is **not** a bottleneck: at `warm_pool_size = 2000`
it barely registers (0.14 s). It is already vectorised (`mvfft` over chunks,
Morlet daughters precomputed once, no per-column recomputation of scales). Do
not spend effort there; the big number in `warm_pool_size` is misleading.

---

## Findings, ranked

### 1. Date fields are extracted by string formatting — ~70–90× slower than needed

**Basis: measured. RNG-neutral. Highest value/effort ratio in the package.**

`as.integer(format(date, "%Y"))` converts `Date` → `POSIXlt` → formatted string →
integer re-parse. There are ~25 such sites across `R/`. `format.POSIXlt` alone is
**14.0 % of total self time** in the profile.

On 27,375 dates:

| | current | `as.POSIXlt` once | speedup |
| --- | --- | --- | --- |
| `format(d,"%Y")` + `format(d,"%m")` | 540 ms | 6.0 ms | 90× |
| `compute_water_year(d, 10)` | 536 ms | 7.5 ms | 71× |
| `which(format(d,"%m-%d")=="02-29")` | 184 ms | 8.0 ms | 23× |

Results verified `identical()`.

`calendar.R:102-103`:

```r
lt        <- as.POSIXlt(date)
cal_year  <- lt$year + 1900L
cal_month <- lt$mon + 1L
```

Same treatment for `find_leap_day_indices()` (`calendar.R:138` →
`which(lt$mon == 1L & lt$mday == 29L)`), `build_historical_dates()`
(`calendar.R:223-226`, one `as.POSIXlt` feeding year/month/day),
`generator.R:323,327-330`, `evaluate_generator.R:533,587,770,778-779,830,838-839,1725`,
`climate_perturbations.R:200-202`, `io_netcdf.R:237,313`.

Cheapest form: one internal helper, e.g. `.date_parts(date)` returning
`list(year, month, day)` from a single `as.POSIXlt()`, and route every site
through it.

### 2. `sim_obs_date` is a `Date`-classed accumulator — O(n²) in the daily loop

**Basis: measured. RNG-neutral. One-line fix.**

`resample.R:319` initialises `sim_obs_date <- as.Date(rep(NA, n_sim_day))`, then
assigns into it element-by-element (`resample.R:555, 601, 611, 648, 669, 712`).
Every `x[i] <- v` on a classed vector dispatches `[<-.Date` → `NextMethod()`,
which defeats R's in-place-modification optimisation and **copies the whole
vector on each assignment**.

30,000 sequential assignments into a length-30,000 vector:

| | time |
| --- | --- |
| `x <- as.Date(rep(NA, n))` (current) | 3,587 ms |
| `x <- numeric(n)`, re-class at end | below timer resolution |

In the profile this is 3.5 % self time at `resample.R:712` — but that was
`n_years ≈ 30`. It grows quadratically: at `n_years = 75` it is ~6× worse per
realization, at 120 ~16×.

The function **already re-classes on exit** (`resample.R:717`,
`class(sim_obs_date) <- "Date"`), so the fix is just the initialiser:

```r
sim_obs_date <- rep(NA_real_, n_sim_day)   # was: as.Date(rep(NA, n_sim_day))
```

Verified: assigning `Date` values into a plain double vector and re-classing at
the end yields a result `identical()` to the current code, NAs included.
Note `NA_real_`, not `NA` — plain `NA` gives a logical vector and coerces on
first assignment.

### 3. `.align_obs_sim_periods` recomputes the same water-year vector per grid, per realization

**Basis: measured (11.3 % of the run). RNG-neutral. Compounds with #1.**

`.filter_grid_years()` (`evaluate_generator.R:582-592`) calls
`compute_water_year(df$date, ...)` on every grid data frame. It is invoked
`n_grids` times for obs (`:645`) and `n_grids × n_realizations` times for sim
(`:649-653`) — with `eval_max_grids = 25` and 3 realizations that is **100 calls**
on the *same* date vector. The function's own comment at `:612` already states
the assumption that all grids share a date vector.

Fix: compute the year vector once per side, derive `keep_idx <- yrs %in% years_keep`
once, and apply `df[keep_idx, , drop = FALSE]` across grids. Same pattern in
`.summarize_observed_data` (`:770-779`) and `.summarize_simulated_data`
(`:830-839`), which rebuild month/day/year per grid.

`.pick_contiguous_year_window()` draws from `sample.int`, but hoisting the year
vector does not move that draw, so the change is RNG-neutral.

### 4. `year_subset_cache` never hits — pure overhead plus a memory leak

**Basis: measured (0 % hit rate). RNG-neutral.**

`resample.R:335,400-401` keys the cache on
`paste(obs_year_draw, collapse = ",")`, where `obs_year_draw` is
`annual_knn_n` (100–120) **order-dependent draws with replacement**. Two
simulated years colliding is astronomically unlikely. Simulated:

| n_years | annual_knn_n | unique keys | hits |
| --- | --- | --- | --- |
| 30 | 100 | 30 / 30 | 0 |
| 30 | 120 | 30 / 30 | 0 |
| 75 | 100 | 75 / 75 | 0 |
| 75 | 120 | 75 / 75 | 0 |

So it never saves work, and it is worse than neutral: each entry retains ~7
vectors of length `annual_knn_n × 365 ≈ 36,500–43,800` plus a `split()` lookup
list, and entries are **never evicted** — the environment grows to `n_years`
entries and is held for the whole call, replicated in every PSOCK worker. At
`n_years = 75`, `annual_knn_n = 120` that is on the order of 100+ MB of live
memory per worker, for zero benefit.

Options: (a) delete the cache and build the subset inline each year — simplest,
strictly faster, big memory win; or (b) if caching is wanted, key on
`sort(unique(obs_year_draw))` and cap the environment size. (b) only helps if
duplicate year *sets* actually occur, which the numbers above suggest they do not.

### 5. Nested clusters: `filter_warm_pool` oversubscribes cores

**Basis: inspection (not exercised — the pipeline runs `parallel = FALSE`).
Correctness-adjacent, not just speed.**

`generate_weather()` creates and registers a cluster of `n_cores` workers
(`generator.R:255-256`) and keeps it alive via `on.exit` (`:259`) for the whole
call. It then calls `filter_warm_pool(parallel = parallel, n_cores = n_cores)`
(`:411-426`), which creates **its own** `parallel::makeCluster(use_cores)` at
`warm_filtering.R:1098` and again at `:1390`.

Result: `2 × n_cores` R processes while the first cluster sits idle, each holding
a copy of `sim_padded`. On a machine where `n_cores = detectCores() - 1` this
oversubscribes every physical core. Worth resolving either way (defer the
`generate_weather` cluster until after WARM filtering, or have `filter_warm_pool`
reuse the registered backend).

Note `evaluate_weather_generator`'s cluster (`evaluate_generator.R:865`) is fine —
`generate_weather` has already returned and stopped its cluster by then.

### 6. Nested `ifelse()` for three-state classification — 4.5× slower than arithmetic

**Basis: measured. RNG-neutral.**

`match_transition_positions` is **13.6 % of the run**; `ifelse` is 12.1 % of self
time. `ifelse()` allocates a logical mask, two branch vectors, and does the
subsetting in R. The state encoding is monotone in precipitation, so it is a sum
of two comparisons.

On 1e5 values:

| | time |
| --- | --- |
| `ifelse(p <= wt, 0L, ifelse(p <= et, 1L, 2L))` | 6.8 ms |
| `(p > wt) + (p > et)` | 1.5 ms |

Results verified `identical()` after `as.integer()`. Note `(p > wt) + (p > et)`
returns integer already when both operands are logical.

Sites: `resample.R:1204-1208` (`match_transition_positions`) and
`resample.R:919-923` (`estimate_monthly_markov_probs`, over the full observed
lag vectors once per simulated year).

**Guard required.** The two forms agree only when
`wet_threshold < extreme_threshold`. If the thresholds are inverted, `ifelse`
yields 0/1 while the sum yields 0/2. `generate_weather` enforces
`extreme_q > wet_q` (`generator.R:203`), but `match_transition_positions` is
**exported** and takes thresholds directly, so a caller can violate it. Add the
check to the exported function if this substitution is applied there.
(`NA` inputs behave identically under both forms — no issue.)

### 7. `estimate_monthly_markov_probs` allocates the full simulation length per year

**Basis: asymptotic — the allocation is quadratic by construction. Measured
timing does *not* yet show it dominating at `n_years <= 120`. RNG-neutral.**

Called once per simulated year from the `y` loop (`resample.R:511`), it allocates
**nine** vectors of length `n_days_sim = 365 × n_years` (`:925-927`) and fills
only the 365 rows belonging to the current year. Total allocation is
`9 × 365 × n_years²` doubles. It also runs 12 full-length
`which(sim_month == mm & sim_wyear == ...)` scans per call (`:932`) —
also quadratic in aggregate.

Measured, holding the observed record fixed and varying `n_years` (single
unreplicated `system.time` per row — treat as indicative only):

| n_years | cost of the full year loop | per call |
| --- | --- | --- |
| 30 | 0.74 s | 24.7 ms |
| 60 | 1.20 s | 20.0 ms |
| 120 | 3.30 s | 27.5 ms |

Per-call cost is essentially **flat** over this range, so the quadratic term is
not yet dominating at `n_years <= 120` — the per-call cost is still carried by
the O(n_obs) state classification of the observed lag vectors. The finding rests
on the allocation arithmetic instead, which needs no measurement: at
`n_years = 120` the nine vectors are ~394k doubles per call × 120 calls
≈ 380 MB of allocation churn, plus 1,440 full-length `which()` scans, all to
fill 365 rows per call.

Fix: return the nine 12-element monthly probability vectors (or a 3×3×12 array)
and index them by month at use time in `markov_next_state`, instead of
broadcasting onto the full simulated time axis. That also removes the 12
full-length `which()` scans. No draws occur in this function, so it is RNG-neutral.

### 8. Per-day allocations inside the `j in 2:365` loop

**Basis: measured (small individually). RNG-neutral. Do these alongside #2/#6.**

The inner loop runs `365 × n_years` times per realization:

- `key <- paste(cur_month, cur_day, sep = ".")` (`resample.R:594`) — `paste` is
  1.6 % self time. The 365-day calendar is fixed, so the 365 keys (or better,
  integer `month * 100L + day` indices into a preallocated lookup) can be built
  once outside the year loop. Same at `:536`.
- `obs_month_mean_precip[as.character(cur_month)]` (`:679-680, 682-683`) — four
  `as.character()` calls per day to index a `setNames(..., 1:12)` vector that
  indexes correctly by integer. `as.character` is 0.75 % self time.
- `as.matrix(cbind(precip_day0_anom, temp_day0_anom))` inside
  `.knn_draw_one_rank` (`:338, 689`) — `p == 2` is known at this call site;
  the `as.matrix` and the `ncol` dispatch can be skipped by passing the two
  vectors directly.
- `expand_indices` (`:752-756`) is 3.0 % total; it allocates a
  `length(base_idx) × length(offset_vec)` vector twice per day (once for the
  ±3 window, again for the ±30 fallback). Minor, but the `rep()` index vectors
  are constant and could be hoisted.

### 9. Plot rendering is ~25 % of the run and is a knob, not a bug

**Basis: measured. RNG-neutral.**

`create_all_diagnostic_plots` is 24.7 % of total; `ggsave` alone is 22.7 %,
with `.export_multipanel_plot` (`evaluate_generator_plots.R:156-182`) the single
largest self-time line in the profile (20.4 %). `dpi = 300` is hard-coded at
`:177`, and `generator.R:381,431-436` write four more PNGs.

Not a defect — but worth knowing that `save_plots = FALSE` removes roughly a
quarter of the wall clock, and that passing `device = ragg::agg_png` (or
lowering `dpi` for iterative runs) typically halves the remaining render cost.
Exposing `dpi` / `device` through `plot_config` would make that a user choice.
`ragg` would be a new dependency, so weigh it against the deliberately small
`Imports` set.

### 10. Parallel topology leaves cores idle

**Basis: inspection. RNG-affecting if the unit of parallelism changes — rank last.**

`generate_weather` parallelises over realizations only
(`foreach(n = seq_len(n_realizations))`, `generator.R:452`), and
`.summarize_simulated_data` caps at `min(n_cores, n_realizations)`
(`evaluate_generator.R:861`). Typical `n_realizations` is 3–20; on a 16-core
machine with `n_realizations = 3`, 13 cores idle.

The current path is reproducibility-safe: `resample_weather_dates` receives an
explicit `daily_seed + n` per realization (`generator.R:468`), so per-realization
output is **invariant to worker count** — `clusterSetRNGStream` only seeds the
workers, it does not determine the result. That is the property to preserve.

So the cost of #10 is not that the existing path is fragile; it is that any *new*
unit of parallelism (chunking years, parallelising grids in evaluation) needs its
own explicit, index-derived seed rather than inheriting a worker stream —
otherwise results start depending on how work was distributed. That is a design
decision, not a mechanical change. Findings 1–8 are the cheaper path to the same
wall clock.

---

## Documentation drift found along the way

`AGENTS.md` states: *"Grid-cell iteration runs in a PSOCK cluster
(`parallel::makeCluster` + `doParallel::registerDoParallel`, `R/generator.R`)."*
The cluster in `R/generator.R` iterates over **realizations**, not grid cells
(`generator.R:452`). Grid cells are iterated sequentially inside
`compute_area_averages` and the evaluation summaries. Worth correcting so future
work does not optimise against the wrong model.

## Suggested order of work

1. #1 (date parts) and #2 (Date accumulator) — largest measured wins, both
   mechanical, both bit-identical for a fixed seed. Do these first.
2. #3 and #4 — remove redundant recomputation and dead cache; #4 is also the
   main memory win.
3. #6 and #7 — inner-loop arithmetic and the Markov allocation shape.
4. #5 — resolve the nested-cluster oversubscription (correctness-adjacent).
5. #8, #9 — polish and knobs.
6. #10 — only with an explicit reproducibility decision.

Each of #1–#8 should be verifiable by asserting that
`generate_weather(..., seed = 100)` produces output `identical()` to the
pre-change result on the packaged fixture. That is worth writing as a regression
test before touching anything.

## Refs

- `R/resample.R` — findings 2, 4, 6, 7, 8.
- `R/calendar.R`, `R/evaluate_generator.R` — findings 1, 3.
- `R/warm_filtering.R:1098,1390` — finding 5.
- `AGENTS.md` "Conventions" — the stale grid-cell parallelism claim.
