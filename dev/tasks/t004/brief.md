# Brief: journal figures for the perturbed climate

Task `t004` (formerly `t2611161814`). Written 2026-11-16 for an implementing agent with no prior context.

## Outcome

Running one dev script from the repository root produces six publication-quality figures that show what
`apply_climate_perturbations()` does to temperature and precipitation, drawn from a whole ensemble of
realizations rather than one cell of one realization. A short findings file reports how much interannual
precipitation variability the perturbation removes.

The figures are about the *imposed climate change*. Figures of stochastic variability alone already exist
(`R/evaluate_generator_plots.R`) and are out of scope.

## Context

- The existing `plot_climate_perturbations()` (`R/climate_perturbations_plots.R`, uncommitted owner work)
  plots target against outcome per month and year. Its output is in `output/perturbation-preview/`. It
  verifies the implementation; it does not describe the perturbed climate. Leave it unchanged.
- Perturbation mechanics to read first: `R/climate_perturbations.R` and `apply_quantile_mapping()` in
  `R/quantile_mapping.R`. Temperature gets an additive monthly delta. Wet-day precipitation is mapped from
  a fitted monthly Gamma to a target Gamma whose mean and variance are scaled. A transient factor `f` ramps
  linearly from 1 to `2f - 1` over the years (temperature from 0 to `2 * delta`), so it averages to the
  stated value.
- `dev/scripts/manual-weather-run.R` shows the full generate → prepare → perturb call sequence on the
  bundled Ntoum fixture (six cells). Reuse that sequence.

### Suspected defect to quantify, not fix

With `enforce_target_mean = TRUE` (the default), `apply_quantile_mapping()` rescales every month-year group
so its wet-day mean equals *monthly climatology × factor* (`R/quantile_mapping.R`, the `enforce_target_mean`
block). Each year's wet-day mean intensity is therefore identical, and interannual variability of monthly
and annual totals survives only through wet-day counts. This is inferred from the code and from
`output/perturbation-preview/precipitation_targets-transient.png`; it has not been measured.

## Scope and constraints

Explicit restrictions:

- No package behaviour change. Do not edit `R/`, `tests/`, `man/`, `NAMESPACE`, `DESCRIPTION`, `NEWS.md`,
  or `.github/`. If a figure seems to need a package change, stop and report it.
- Do not fix the `enforce_target_mean` behaviour. Changing it moves numeric output and is the owner's
  decision.
- The primary checkout holds uncommitted owner work. Do not commit, stash, reset, switch branches, or push.
  Touch only the files this task owns.
- Do not start other agents. Never `source("tools/dev_workflow.R")`.
- Add no dependency to `DESCRIPTION`. The script may use a layout package (`patchwork` or `cowplot`) only if
  it is already installed; otherwise compose with facets.

Task-owned surfaces (starting points, not an allowlist):

- `dev/scripts/perturbation-figures.R` — new; may be split into a second helper file beside it.
- `dev/tasks/t004/findings.md` — new; variance diagnostic and anything the owner must decide.
- `dev/README.md` — add a short usage paragraph next to the manual-run section. The file already has
  uncommitted owner edits; add your hunk without rewriting theirs.
- `output/perturbation-figures/` — generated results (git-ignored).

Unknowns to resolve by inspection:

- Runtime of a 20-realization, 30-year run on this machine. Follow the global rule for long processes
  (visible Herdr `watch` tab when `HC_TERMINAL_ID` is set).
- The wet-day fraction of the fixture at thresholds 0 and 1 mm/day. The preview showed about 100% wet days
  at threshold 0. Report both in the findings.
- Whether the helper statistics in `tests/testthat/test-workflow-e2e.R` (spell lengths, anomaly
  correlations) can be copied into the script. Copy definitions; do not source test files.

## Experiment design

Script settings live in one editable `config` list, with a reduced "smoke" configuration selectable by a
command-line argument so the figures can be rebuilt quickly.

| Setting | Full | Smoke |
|---|---|---|
| Realizations | 20 | 4 |
| Years | 30 calendar years from 2030 (`year_start_month = 1`) | 16 |
| Variables | `precip`, `temp`, `temp_min`, `temp_max` | same |
| Seeds | fixed generator and perturbation seeds | same |

Calendar years keep the simulation-year index equal to the calendar year. Keep `n_years >= 16` (wavelet
requirement). Cache the generated and prepared ensemble as RDS in the output directory and reuse it when
present, so figure work does not rerun the generator.

Scenarios, all with the same factor in every month:

| Label | Temperature | Wet-day mean | Other | Used in |
|---|---|---|---|---|
| S0 baseline | none | none | unperturbed ensemble | all |
| S1 warm-dry | +2 °C | × 0.85 | defaults | 1, 2, 4, 6 marker |
| S2 hot-wet | +4 °C | × 1.15 | defaults | 1, 2, 4, 5, 6 marker |
| M1 mean | 0 | × 1.15 | variance scales with mean | 1, 3 |
| M2 variance | 0 | × 1.00 | variance × 1.5, `scale_var_with_mean = FALSE` | 1, 3 |
| M3 tail | 0 | × 1.15 | `exaggerate_extremes = TRUE`, threshold 0.95, `extreme_k = 1.2` | 1, 3 |

S1 and S2 run in step mode (`*_transient = FALSE`) for figures 2, 5 and 6, and in transient mode for
figures 1 and 4. M1–M3 are step runs. Leave `precip_intensity_threshold` at its default.

Conventions for every figure:

- Wet day for *statistics*: `precip >= 1` mm/day.
- Area mean: unweighted mean over the six cells. Used for figures 2, 4 and 6.
- Representative cell: cell 1, named in the caption. Used for figures 1 and 3, because quantile mapping
  acts per cell and an area mean smooths the tail.
- Ensemble summaries: median line with 5–95% ribbon across realizations (figure 4 adds 25–75%).
- Compare perturbed with baseline *within the same realization* wherever a ratio or difference is shown.
  Pairing removes sampling noise; unpaired pooled curves were the weakness of the preview.

## Work

### Phase 0 — variance diagnostic (do first)

For S2 in step mode, run the perturbation twice: `enforce_target_mean = TRUE` and `FALSE`. For each
realization and cell compute, for baseline and both perturbed series:

- interannual SD and CV of annual precipitation totals;
- SD across years of each month's wet-day mean (threshold as used by the mapping), averaged over months;
- lag-1 autocorrelation of annual totals;
- correlation between baseline and perturbed annual-total anomalies.

Report the perturbed/baseline CV ratio as median and range over realization × cell, for each setting, in
`findings.md` with a three-line conclusion. A pure rescaling gives a CV ratio of 1.

Falsifier: if the CV ratio with `enforce_target_mean = TRUE` is 0.9 or higher, the suspected defect is not
material; say so plainly. Either way continue with the figures using package defaults.

### Phase 1 — figures

| # | File stem | Content |
|---|---|---|
| 1 | `fig1-mechanics` | How the perturbation transforms values; four panels |
| 2 | `fig2-seasonal-cycle` | Monthly climatology shift against ensemble spread |
| 3 | `fig3-daily-distribution` | Wet-day exceedance curves and paired quantile ratios |
| 4 | `fig4-transient` | Annual trajectories under the transient ramp |
| 5 | `fig5-preserved-structure` | Statistics the perturbation should leave unchanged |
| 6 | `fig6-response-surfaces` | Climate indices across the stress-test grid |

**Figure 1 — mechanics** (double column, 2 × 2).

- (a) Cell 1, the month with the largest baseline wet-day mean. Baseline empirical CDF as points, fitted
  baseline Gamma and M1 target Gamma as lines (parameters from the perturbation diagnostic's `base_gamma`
  and `target_gamma`). Trace one value, the baseline 90th percentile: value → probability → mapped value,
  with arrows.
- (b) Transfer function, mapped against original wet-day precipitation, for M1, M2 and M3 as analytic
  curves (`qgamma(pgamma(x))` with the tail transform for M3), plus a 1:1 line. Underlay the actual paired
  daily values for M1 as faint points so the spread added by the per-year mean rescale is visible.
- (c) Realized annual-mean warming (perturbed minus baseline) against year for S1 and S2, step and
  transient.
- (d) Specified precipitation-factor path against year for S1 and S2, step and transient, as lines, with
  the realized annual area-mean precipitation ratio as points.

**Figure 2 — seasonal cycle** (double column, three panels: precipitation mm/month, mean temperature °C,
PET mm/month). Month on x. One climatology per realization; ribbon and median across realizations for S0,
S1, S2 in step mode. Add the observed climatology as black points.

**Figure 3 — daily distribution** (single or 1.5 column, two stacked panels sharing x). Cell 1, wet days,
pooled over years within each realization; quantiles on a fixed probability grid from 0.5 to 0.999 so
bands are defined.

- Top: precipitation against wet-day exceedance probability (log scale, rare events to the right), for S0
  and M1–M3.
- Bottom: paired quantile ratio, perturbed/baseline, on the same x. Reference lines at ratio 1, at the
  mean factor 1.15, and a vertical line at exceedance 0.05 where the tail transform starts.

**Figure 4 — transient trajectories** (double column; rows: annual area-mean precipitation mm/yr, annual
mean temperature °C; columns: S1, S2). Fan charts for baseline (grey) and perturbed (colour). Overlay the
imposed path as a dashed line: baseline ensemble mean × factor path, and baseline ensemble mean + delta
path. If phase 0 confirms the variance loss, also write `fig4b-transient-enforce-off` with a third fan for
`enforce_target_mean = FALSE` in the precipitation row.

**Figure 5 — preserved structure** (double column, small multiples of perturbed against baseline with a
1:1 line; S2 step against S0; equal x and y range within each panel). One point per cell × realization
unless stated:

- wet-day fraction, per month;
- mean wet-spell length; mean dry-spell length;
- lag-1 autocorrelation of daily precipitation;
- cross-cell correlation of daily precipitation, per cell pair;
- correlation of daily precipitation and temperature anomalies about calendar-month means;
- interannual CV of annual precipitation.

**Figure 6 — response surfaces** (double column). Step runs on the grid ΔT ∈ {0, 1, 2, 3, 4, 5} °C ×
wet-day mean factor ∈ {0.7, 0.8, …, 1.3}. Ensemble-mean, area-mean percentage change from S0, as tiles
with labelled contours and a diverging scale centred on zero:

- (a) annual precipitation; (b) Rx1day, the mean annual maximum daily precipitation, computed per cell
  then averaged; (c) annual PET; (d) aridity index, precipitation over PET.
- (e) PET change against ΔT as lines for `temp_range_factor` 0.9, 1.0 and 1.1.

Mark S1 and S2 on panels (a)–(d). The realization count for this figure may be reduced if runtime
requires; state the count used in the caption.

### Presentation

- `theme_weathergenr()` as the base; Okabe–Ito colours, one fixed colour per scenario across all figures.
- Final-size design: 90 mm single column, 180 mm double column, text 7–9 pt at that size.
- No title inside the image. Panel tags `(a)`, `(b)`, …. Axis labels carry units.
- Write each figure as PNG at 300 dpi and as vector PDF, plus its plotted data as CSV.
- Write `captions.md` with a publication-style caption per figure: what is shown, scenario settings,
  ensemble size, cell or area mean, wet-day threshold, and band definition.

### Phase 2 — documentation

Add the usage paragraph to `dev/README.md`. Update the Progress checklist in
`dev/tasks/t004/t004-perturbed-climate-journal-figures.md`, ticking steps with what each produced.

## Acceptance and verification

- `Rscript dev/scripts/perturbation-figures.R` (full) and its smoke variant exit 0 from the repository
  root and write six PNG, six PDF, six CSV, and `captions.md` under `output/perturbation-figures/`.
- `findings.md` exists with the phase 0 table, the wet-day fractions, and a stated conclusion.
- Each PNG has been opened and inspected at final size: no overlapping or clipped text, legends readable,
  colours consistent between figures. Fix and re-render before reporting.
- Numeric sanity checks, reported with their values:
  - figure 2: perturbed minus baseline monthly temperature equals the scenario delta;
  - figure 3: the M1 paired ratio is flat near 1.15 below the tail, and M3 departs from it above
    exceedance 0.05;
  - figure 6: the (ΔT = 0, factor = 1.0) node is 0% in every panel, and annual precipitation change tracks
    the factor.
  A failed check is a finding to report, not something to tune away.
- `Rscript tools/lint.R dev/scripts/perturbation-figures.R` reports no lints, or the remaining ones are
  listed with a reason.
- `git status --short` shows changes only in the task-owned surfaces.

The package test suite and the numeric baseline gate are not required: no package code changes.

## Return

Report to the coordinator: file paths, the phase 0 conclusion in two sentences, check results with values,
any deviation from this brief with its reason, and open decisions for the owner. Put detail in
`findings.md`, not in the message.
