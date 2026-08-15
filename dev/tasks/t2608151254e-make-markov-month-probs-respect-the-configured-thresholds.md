---
title: Make .markov_month_probs() respect the configured thresholds
type: todo-item
status: backlog
effort: 1
area: resample
origin: review-2026-08-15
queue: 6
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — The non-finite-threshold fallback hardcodes quantiles 0.2 and 0.8 instead of the run's wet_q/extreme_q, and dirichlet_alpha divided by sqrt(n) is neither documented nor reachable from generate_weather().
> **Why** — A degenerate month silently disagrees with the run's own definition of wet and extreme, and the smoothing decays as n^-1.5 so it vanishes exactly where rare extreme-to-extreme transitions would want it.
> **Effort** — Small. The fallback fix is mechanical; the `dirichlet_alpha` half is a short judgment about whether the 1/sqrt(n) scaling was intended.

## Progress

- [ ] Thread `wet_q`/`extreme_q` into the non-finite-threshold fallback (`R/resample.R:1005-1009`) in place of the hardcoded 0.2 and 0.8.
- [ ] Add a test that reaches the fallback — a month with a degenerate threshold — since nothing currently exercises it.
- [ ] Decide whether `dirichlet_alpha / sqrt(n_transitions_m)` (`:1036`) is intended. Standard Dirichlet smoothing holds alpha fixed and lets its influence decay as ~1/n; dividing by sqrt(n) makes it decay as ~n^-1.5, so with ~900 transitions per month the effective alpha is ~0.03.
- [ ] If intended, document the reasoning at the call site. If not, drop the sqrt and pick a fixed alpha.
- [ ] Decide whether `dirichlet_alpha` should be reachable from `generate_weather()` at all — today it is hardcoded at its 1.0 default.
- [ ] Baseline gate if the smoothing changes.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C7.
