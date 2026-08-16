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
- [ ] **Held, not merged.** Owner chose to bundle this with the deferred `t2608151254` so consumers absorb one break rather than two.

## Held for 2.0.0

The work is complete and verified on `refactor/trim-exports`; the branch is
deliberately unmerged and `master` is untouched. Nothing here is urgent — the
exports do no harm while they exist, and the cost of the break is paid once,
when it lands.

To finish: decide `t2608151254` (whether climate perturbation belongs in
`run_weather_generator()`), land both, then bump to 2.0.0 and tag. Rebase this
branch onto `master` first — `master` has moved since it was cut, and the
branch is local and unpushed, so a rebase is the right resync.

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
