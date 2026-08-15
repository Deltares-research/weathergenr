---
title: Clear lintr findings and make the lint job blocking
type: todo-item
status: backlog
effort: 2
area: ci
queue:
created: 2026-08-06
updated: 2026-08-06
---

> [!note] Overview
> **What** — Fix the 28 findings Rscript tools/lint.R reports under the .lintr baseline, then remove continue-on-error from .github/workflows/lint.yaml so the job gates.
> **Why** — The lint job is advisory, so a red result is invisible and the debt grows unchecked.
> **Effort** — large

## Progress

- [ ] Decide the `.Random.seed` save/restore fix (10 findings across `evaluate_generator.R`,
      `generator.R`, `resample.R`, `warm_filtering.R`): either `# nolint` per site, or
      switch to `get()`/`assign(envir = globalenv())`, which also removes the `<<-`.
      The second touches RNG-critical code — re-run the full suite if taken.
- [ ] Remove 16 unused locals (`grid_count_original`, `metrics_list`, `out`, `n_grids`,
      `n_vars`, `start_time`, `progress_str`, `nt`, `pool_pct`, `n_pass`, `rate`, `crit`,
      `status_icon`, `status_txt`, `pct`) — check each is genuinely dead and not a
      dropped verbose-message input before deleting.
- [ ] Replace 2 `seq_len(length(x))` with `seq_along(x)` in `warm_filtering_plots.R`.
- [ ] Confirm `Rscript tools/lint.R` exits 0, then drop `continue-on-error` from
      `.github/workflows/lint.yaml`.

## Refs

- `.lintr` — the baseline config; its header explains which style linters are frozen
  and how to re-enable one at a time.
- Findings were 1,844 before the baseline and 28 after.

### Correction — 2026-08-15

Re-measured during the `review-2026-08-15` sweep: `Rscript tools/lint.R` now
reports **20**, not 28. Two counts in the checklist above are wrong, and one of
them would cause damage:

- **12** are `.Random.seed` bindings, not 10 — `wavelet_cwt.R:578`/`:584` and
  `wavelet_warm.R:150`/`:156` are also in the set.
- **6**, not 16, are unused locals, and all six are **false positives**:
  `pool_pct`, `n_pass`, `rate`, `crit`, `status_icon`, `status_txt` and `pct`
  in `log_filter_iteration()`/`log_final_summary()` are used inside `.log()`
  glue strings (`"{n_pass}"`), which lintr cannot see through. Deleting them
  breaks the filter logging. Suppress or restructure; do not remove.

The `seq_len(length(x))` pair in `warm_filtering_plots.R` still stands.
See `dev/drafts/package-review-2026-08-15.md` § B10.
