---
title: Tail criterion rejects the whole pool at its default bound
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
> **What** — compute_tailmass_metrics' logdiff against the default tail_tol_log = log(1.03) admits under 1 percent of candidates at every trace length tested, so the pool only survives because the relaxation loop then loosens the bound.
> **Why** — A default that rejects 99 percent of candidates is not a filter, it is a trigger for relaxation; the effective bound is whatever the loop happens to reach rather than the documented one.
> **Effort** — small

## Progress

- [ ] Reproduce independently of the measurement that surfaced this. Traces from `simulate_warm(match_variance = TRUE)` against the packaged driver gave a tail pass rate of 0.4-0.8% at every length from 21 to 120 years, against `tail_tol_log = log(1.03)`.
- [ ] Work out whether the bound or the metric is wrong. `log(1.03)` is about 0.03, while observed `logdiff_low` values sit around 2.6 — two orders of magnitude apart, which suggests the tolerance was chosen for a differently-scaled quantity rather than merely set too tight.
- [ ] Check the lower and upper tails separately. In the measurement `logdiff_low` was large and stable while `logdiff_high` shrank with trace length, so they may not want a shared bound at all.
- [ ] Trace what this does to the relaxation loop. If the tail filter always has the lowest pass rate it is always the one `which.min(rates)` picks, so it may be absorbing relaxation budget that other criteria need — check `relax_log` on a real run.
- [ ] Whatever changes, baseline gate before and after.

## Refs

- Surfaced while measuring `t2608151507`, which was dropped; the tail finding
  is independent of that item's question and of trace length.
- `compute_tailmass_metrics()` normalises each side by its own length as of
  `t2608151254a`, so this is not a length-mismatch artifact.
- `filter_warm_bounds_defaults()` sets `tail_tol_log = log(1.03)` and
  `relax_tail_tol_log_max = log(2.0)`.
