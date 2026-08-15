---
title: Make the WARM filter's observed benchmark deterministic
type: todo-item
status: backlog
effort: 1
area: warm
origin: review-2026-08-15
queue: 5
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — filter_warm_pool() scores candidates against a single random contiguous window of the observations, and returns full-length traces that were scored on a window.
> **Why** — Acceptance thresholds move with the seed rather than only the draw, and a trace admitted on a 40-year slice ships as a 100-year trace whose full-length statistics were never tested.
> **Effort** — Small code change, but it moves numeric output and the right statistic (median over windows? full record?) is a judgment call, not a lookup.

## Progress

- [ ] Baseline gate before starting — this moves numeric output.
- [ ] Replace the single observed window (`R/warm_filtering.R:209-214`) with something stable: the full record, or a median of the statistic across several windows. The full record is simplest but reintroduces the length mismatch the windowing was there to avoid.
- [ ] Decide what to do about the per-trace simulation windows (`:222-228`) when traces are longer than the observed record — scoring on a window and shipping the full trace (`:561`) is the part that most needs either a fix or an explicit statement.
- [ ] Document the chosen behaviour in `filter_warm_pool()`'s `@details`; nothing in "filtering thresholds" currently prepares a user for windowing at all.
- [ ] Baseline gate after; re-record once the delta is reviewed.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C3.
- Note this interacts with `t2608151253a`: the mean-filter inertness measured
  there holds cleanly only on the unwindowed path, which is the default since
  `n_years` defaults to the observed year count (`R/generator.R:330`).
