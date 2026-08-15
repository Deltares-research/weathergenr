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

- [x] Reproduced: 0.4% tail pass rate against `log(1.03)` on 2,000 candidates from the fitted model.
- [x] **The metric is fine; the bound is wrong.** My opening claim that `logdiff` sat around 2.6 — two orders off the tolerance — was an artifact of the unrepresentative probe that surfaced this item. Measured properly, `M_obs_low` is 0.083 against a simulated median of 0.081, and `M_obs_high` 0.039 against 0.051: the masses agree well. `tail_eps = 1e-5` is four orders below them, so it distorts nothing. The real medians are `logdiff_low` 0.238 and `logdiff_high` 0.416, against a bound of 0.0296.
- [x] Diagnosis: 3% is the right shape of tolerance for a mean or a standard deviation, which are stable statistics. Tail mass over the four-or-so points beyond a quantile of a 21-year record is not, and its natural spread across candidates from the fitted model is tens of percent. The bound asks for a precision the statistic cannot deliver.
- [x] Traced the relaxation loop, which is where this actually bites. The tail filter always has the lowest pass rate, so `which.min(rates)` always picks it, and the loop stops the moment the pool reaches `n_select`. **The loop, not the default, was setting the operative bound.** At `generate_weather()`'s own defaults (5,000 candidates, 20 realizations) that left a pool of **24**, after 4 relaxation iterations — so the documented filter-then-rank design was in practice relax-until-just-enough, with the ranking choosing 20 from 24.
- [x] Measured the alternatives at production scale, pool size and relaxation iterations: `log(1.03)` → 24 after 4 iterations; `log(1.25)` → 157, none; `log(1.50)` → 297, none; `log(2.00)` → 495, none.
- [x] Changed the default to `log(1.25)`. It still rejects about 97% of candidates so it stays a real constraint, needs no relaxation, and sits well inside `relax_tail_tol_log_max = log(2.0)` — a ceiling the author already treated as acceptable.
- [x] Checked lower and upper tails separately, as the step asked: `logdiff_low` is the larger and more stable of the two, `logdiff_high` shrinks with trace length. A shared bound is defensible at `log(1.25)`, where each admits a comparable share; splitting them would be a separate change and is not needed to fix this.
- [ ] **Owner decision pending** — baseline delta reviewed but *not* re-recorded. See below.

## Baseline delta — awaiting acceptance

**60 of 106 keys differ.** Realization selection changes, so the analogue dates
and every downstream statistic change with it. All `config.*` and structural
keys hold.

Unusually for this sweep, the delta has a direction. Across the 40 evaluation
error metrics: **26 improved, 14 worsened, median relative change −7.9%**, and
every metric family improved in at least half its keys (most in 3 of 4). That is
consistent with the ranking finally having a genuine pool to select from rather
than taking whatever survived. With two realizations per scenario it is
directional evidence, not a significant result — a 26/14 split is p ≈ 0.08 —
so it should be read as "no sign of degradation, some sign of improvement".

Re-record with `Rscript tools/record_baseline.R --force` once accepted.

## Refs

- Surfaced while measuring `t2608151507`, which was dropped; the tail finding
  is independent of that item's question and of trace length.
- `compute_tailmass_metrics()` normalises each side by its own length as of
  `t2608151254a`, so this is not a length-mismatch artifact.
- `filter_warm_bounds_defaults()` sets `tail_tol_log = log(1.03)` and
  `relax_tail_tol_log_max = log(2.0)`.
