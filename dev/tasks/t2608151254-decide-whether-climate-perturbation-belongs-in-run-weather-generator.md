---
title: Decide whether climate perturbation belongs in run_weather_generator()
type: todo-item
status: backlog
effort: 2
area: api
origin: review-2026-08-15
queue: 1
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — apply_climate_perturbations() is unreachable from either entry point, and the interfaces do not meet: it needs temp_min/temp_max per cell, while generate_weather() validates only precip/temp and returns dates rather than values.
> **Why** — The README calls the three components coupled; component 3 is hand-assembled by every user. Either wire it up or stop claiming the coupling.
> **Effort** — Large if wired up, since `generate_weather()` would have to return values rather than dates and carry temp_min/temp_max; trivial if the answer is to correct the README instead. Decide that first.

*Deferred 2026-08-15 at the owner's request — moved to the back of the queue,
not blocked. Nothing external is waiting on it; the scope question below is a
design decision the owner wants to take separately.*

## Progress

- [x] **The note's premise was wrong on two counts, and both made the job look bigger than it is.** `generate_weather()` does not validate "only precip/temp" — `R/generator.R:289` requires `vars` to *include* them, so `vars = c("precip","temp","temp_min","temp_max")` is accepted and every column rides along regardless, since the generator returns resampled dates rather than values. And the date-to-values join is not work to be done: `prepare_evaluation_data()` (exported, `R/evaluate_generator.R:1626`) already performs it, because `run_weather_generator()` needs it to feed the evaluator. Its `sim_data` is a list of per-cell data.frames with `date` first — exactly `apply_climate_perturbations()`'s `data` argument. The extra `date` column is tolerated: the required-columns check is a `setdiff()`.
- [x] Ran the chain end to end before deciding, rather than reasoning from shapes. 20 years, 6 cells, 2 realizations, no new code: precip mean ratio 0.700 against a 0.70 target, temp delta 2.000 against 2.0. Scripts under `.tmp/scratchpad/2026-08-17_1430/`.
- [x] **Scope chosen: documentation, no new code.** Not the middle option — a thin `apply_perturbations_to_realization()` would hide three lines behind a permanent export contract, which is the wrong trade for a package whose export surface was being trimmed in the same release.
- [x] Full integration rejected on a principled ground, not a cost one. `run_weather_generator()` evaluates its output *against the observed record*; a perturbed series is meant to depart from observations, so folding perturbation in would break what the evaluation means. `blueearth_cst` already composes the two calls itself — the only real consumer had voted.
- [x] `vignettes/climate_perturbation.qmd` perturbs `ncdata$data`, the **observed** record, not realizations. That was the actual docs defect. It keeps the observed-data demo as an introduction and gains a section doing it on a realization.
- [x] `README.md` corrected: components 1–2 are coupled inside `generate_weather()`, component 3 is a stage applied afterwards. `NEWS.md` entry written under *Other changes*.
- [x] Three traps found by running it, now documented and pinned by tests in `test-climate_perturbations.R`: `vars` must include `temp_min`/`temp_max`; the date vector must come from the returned frames, not `gen_output$dates`, which can be longer because incomplete years are dropped; and `n_years x 12` matrices are sized by **calendar** years spanned, so a 20-water-year series needs 21 rows.
- [x] That last trap was also a live documentation bug: `apply_climate_perturbations()`'s "Year indexing convention (critical)" section called `year_idx` a simulation-year index running `1..n_years`, but `R/climate_perturbations.R:247` derives it as `calendar_year - min(calendar_year) + 1`. The two agree only for calendar-year runs. Corrected; no behavior changed.

## Resolution

Closed as **documentation**, landed in `b73938c`. The README's "three coupled
components" claim is now true as written, and the coupling it describes needed
no code because it was already there.

Verified: `test-climate_perturbations.R` 81 pass; `devtools::test()` 1275 pass,
0 fail; baseline gate 8 pass, 0 fail (numeric output unmoved);
`check_only(build_vignettes = TRUE)` 0/0/0, which rebuilt the vignette.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § A1.
- Interface gap: `R/climate_perturbations.R:150` (required columns) versus
  `R/generator.R:208` (validated vars) and `:576` (return shape).
