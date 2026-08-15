---
title: Score the non-spectral WARM criteria on full series, not a window
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
> **What** — Only the spectral criterion needs matched lengths. Mean, sd and tail mass are length-robust, so compute them on the full observed record and the full simulated trace instead of on the harmonised window both currently share.
> **Why** — A 60-year trace is judged on 20 of its years for three of the four criteria, discarding two thirds of the evidence about the thing being filtered.
> **Effort** — Small once the groundwork from `t2608151254a` is in: the tail metric is already normalised per side, so the change is which series each criterion is handed. It moves numeric output, so the cost is the baseline gate and the delta review.

## Progress

- [ ] Split the inputs in `filter_warm_pool()`: pass `obs_series` / `sim_series` to the mean, sd and tail-mass criteria, and keep the harmonised `obs_use` / `sim_use` for `compute_spectral_metrics()` only.
- [ ] Confirm the premise holds for tail mass in practice, not just in principle. It is now normalised per side (`compute_tailmass_metrics()`), but `thr_low`/`thr_high` are still observed quantiles at 0.20/0.80, and a quantile from a 21-point sample carries real variance — check that comparing a 60-point simulated mass against a 21-point observed threshold behaves sanely before relying on it.
- [ ] Decide whether the observed side should use the full record even when it is *longer* than the simulation. Using all 40 observed years to judge a 20-year trace is more data but changes what the criterion means.
- [ ] Baseline gate before and after; re-record only once the delta is reviewed.
- [ ] Document the final behaviour in `filter_warm_pool()`'s `@details` — deferred from `t2608151254a` because what to write depends on this outcome. Nothing there currently prepares a user for harmonisation at all.

## Refs

- Closed sibling `t2608151254a` made the harmonisation deterministic and
  measured why it matters: with 60-year traces scored on 20-year windows, no
  trace passed in more than 52% of its 40 possible windows.
- `compute_tailmass_metrics()` used to divide both sides by the observed
  length; that was fixed there, which is what makes this item cheap.
