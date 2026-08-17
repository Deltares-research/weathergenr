---
title: Trim the export surface
type: todo-item
status: backlog
effort: 2
area: api
origin: review-2026-08-15
queue: 1
branch: refactor/trim-exports
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — NAMESPACE exports criteria_string_compact, generate_symmetric_dummy_points, get_result_index and match_transition_positions, which have no plausible external audience; several carry @keywords internal and are exported anyway.
> **Why** — blueearth_cst pins this package by Git tag, so every export is a contract that cannot be quietly withdrawn. Deciding now is cheaper than after another consumer appears.
> **Effort** — Large not in code but in coordination: unexporting is breaking, the downstream consumer pins by Git tag, and the release has to be sequenced so nothing is pulled from under it.

## Progress

- [x] Checked the consumer against every export rather than the shortlist. `blueearth_cst` references **four**: `run_weather_generator()`, `read_netcdf()`, `write_netcdf()`, `apply_climate_perturbations()`. The other 42 are unreferenced there — but that is not the removal criterion, since the repo is public and `generate_weather()` is obviously legitimate API despite being unused downstream.
- [x] Found the objective signal the note was reaching for: **14** functions are exported *and* declare `@keywords internal`. Milder than it sounds, though — the keyword only hides the `.Rd` from the reference index, it does not affect export, so these are usable-but-unadvertised, which is a recognised idiom rather than a bug.
- [x] Shortlist confirmed and extended by one: `criteria_string_compact`, `generate_symmetric_dummy_points`, `get_result_index`, `match_transition_positions`, `relax_bounds_one_filter`. All five are internal-use only, referenced in tests but in no vignette, and none is in `_pkgdown.yml`.
- [x] Owner chose the narrow scope over unexporting all 14 or curating to ~15. `filter_warm_bounds_defaults()` in particular stays exported despite its internal keyword: users need it to inspect the bounds they are overriding, so it is a case for dropping the keyword rather than the export.
- [x] Implemented on `refactor/trim-exports`: `@export` removed, `@noRd` added so no `.Rd` is generated — matching the fourteen internal helpers that already work that way. That also resolved a failure the first attempt hit: three of the five had examples calling the function unqualified, which `R CMD check` runs where only exports are visible. Dropping the `.Rd` removes the problem rather than papering over it with `:::`.
- [x] Exports 46 → 41. `check_only(build_vignettes = TRUE)` 0/0/0; `WEATHERGENR_BASELINE=1 devtools::test()` 1044 pass, 0 fail, 0 skip.
- [x] `NEWS.md` entry written under a new development heading, naming the five and telling anyone affected to use `weathergenr:::` or open an issue.
- [x] ~~**Held, not merged.**~~ Bundled with `t2608151254` so consumers absorb one break rather than two. **The hold dissolved when that item resolved to documentation only** — with no second break to bundle with, there was nothing left to wait for.
- [x] **Re-applied on `master` rather than rebased**, in `e47d2d1`. The rebase conflicted on the whole of `evaluate_generator_plots.R`, `resample.R` and `warm_filtering.R` — every line, markers at line 1. Cause: the original `f278dba` also flipped those three files from LF to CRLF, so with `core.autocrlf=false` and no `.gitattributes`, git saw a total rewrite. The semantic change is five roxygen edits (`@export` → `@keywords internal` + `@noRd`); those were reproduced directly, leaving line endings untouched. Resolving a whole-file line-ending conflict by hand would have risked silently reverting master's plot-export work in the same files.
- [x] Exports 47 → 42, not the 46 → 41 recorded above: `master` gained `theme_weathergenr()` after the branch was cut.
- [x] Bumped to **2.0.0** in `efc071d`. The tag is deliberately not created — `auto_push: false`, and `blueearth_cst` pins by Git tag, so tagging and pushing is an owner-run step.

## Closing state

Landed on `master`. `refactor/trim-exports` is now redundant and can be deleted;
it is kept for the moment only as the record of the original attempt.

Verified after re-application: `devtools::test()` 1275 pass, 0 fail;
`check_only(build_vignettes = TRUE, document = FALSE)` 0 errors, 0 warnings,
0 notes. Neither `_pkgdown.yml` nor any vignette referenced the five.

The other twelve `@keywords internal`-but-exported functions were left alone
and are a separate decision if it is ever wanted: `compute_spectral_metrics`,
`compute_tailmass_metrics`, `plot_filter_diagnostics`, `extract_signif_curve`,
`extract_wavelet_components`, `fill_nearest`, `gws_regrid`,
`morlet_parameters`, `morlet_wavelet`, `modwt_decompose`, `modwt_reconstruct`,
and `filter_warm_bounds_defaults` (which should keep its export and lose the
keyword).

## Refs

- `dev/drafts/package-review-2026-08-15.md` § A4.
- Related watch-item `t2608151255c` (three competing internal-naming
  conventions) — settling that would make the export list self-evident.
