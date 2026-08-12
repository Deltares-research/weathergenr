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
