---
title: LOG
type: log
---

# LOG — closed work

Terminal ledger for items closed off the board — one line per closed item,
newest first. A browsable index *over* Git history, not a copy of it: full task
bodies stay recoverable with `git log --diff-filter=D -- dev/tasks/`.

Rows are appended by `todoboard done`, never hand-edited.

| Closed     | ID           | Item                                     | Area         |
| ---------- | ------------ | ---------------------------------------- | ------------ |
| 2026-08-15 | t2608151254c | Resolve the order-dependent Cholesky whitening in simulate_warm() | warm         |
| 2026-08-15 | t2608151451  | Decide whether the WARM peak criterion should be a similarity filter | warm         |
| 2026-08-15 | t2608151254d | Document or reshape the Markov spell-length factors | resample     |
| 2026-08-15 | t2608151254h | Fix the RNG restore and legacy error messages in apply_climate_perturbations() | perturbations |
| 2026-08-15 | t2608151254g | Fix the three consistency defects in generate_weather() | api          |
| 2026-08-15 | t2608151254f | Fix the list-components length check and bypass_n doc in simulate_warm() | warm         |
| 2026-08-15 | t2608151254e | Make .markov_month_probs() respect the configured thresholds | resample     |
| 2026-08-15 | t2608151522  | Tail criterion rejects the whole pool at its default bound | warm         |
| 2026-08-15 | t2608151507  | Score the non-spectral WARM criteria on full series, not a window [dropped] | warm         |
| 2026-08-15 | t2608151254a | Make the WARM filter's observed benchmark deterministic | warm         |
| 2026-08-15 | t2608151353  | Stop the WARM pool filter deciding from plotting-only spectra | warm         |
| 2026-08-15 | t2608151254j | Remove the duplicated test-wavelet.R     | tests        |
| 2026-08-15 | t2608151254i | Unignore vignettes/.quarto so check_only() passes | ci           |
| 2026-08-15 | t2608151253c | Add equivalence tests for the two forked numeric routines | tests        |
| 2026-08-15 | t2608151253b | Reconcile warm_signif with the documented WARM component selection | warm         |
| 2026-08-15 | t2608151253a | Fix the mean pinning and sample-sd bias in .variance_match_matrix() | warm         |
| 2026-08-15 | t2608151253  | Fix the NetCDF noleap calendar round trip | io           |
