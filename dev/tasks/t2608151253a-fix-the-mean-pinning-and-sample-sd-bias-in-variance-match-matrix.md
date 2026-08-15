---
title: Fix the mean pinning and sample-sd bias in .variance_match_matrix()
type: todo-item
status: backlog
effort: 1
area: warm
origin: review-2026-08-15
queue: 1
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — The helper replaces the column mean with target_mean for every corrected column, not just the variance, and compares a population sd from .fast_col_sd() against stats::sd sample targets.
> **Why** — 52 percent of traces get their mean pinned to one value, making filter_warm_pool's 3 percent mean criterion inert for half the pool; the sd path adds a systematic +1.27 percent against a 3 percent threshold.
> **Effort** — Two small edits to one function, but both move numeric output, so the cost is in the baseline gate and in reviewing the delta rather than in the code.

## Progress

- [ ] Run the end-to-end baseline gate **before** touching anything and keep the output — this is a numeric-output change (`AGENTS.md` Workflow).
- [ ] Decide whether `.variance_match_matrix()` should touch the mean at all. Rescaling about the column's own mean preserves it; the current `+ target_mean` replaces it (`R/wavelet_warm.R:938-941`).
- [ ] Fix `.fast_col_sd()` (`:904`) to return a sample sd — `* n/(n-1)` on the variance — or switch the targets to a population sd. One or the other, not a mix.
- [ ] Re-check `filter_warm_pool`'s `mean` and `sd` pass rates afterwards; both criteria were being fed biased inputs, so the default 3 percent tolerances may now need revisiting.
- [ ] Run the baseline gate again, review the named keys that moved, and re-record only once the delta is understood.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C9 and § C10 — the measured
  numbers (1048/2000 traces mean-pinned; every clamped trace at exactly
  `target * sqrt(40/39)`, +1.274 percent).
- `.compute_gws_batch()` (`R/warm_filtering.R:1060`) already applies
  `* (n1/(n1-1))` and calls it "sample variance" — the convention to match.
