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

- [x] Implemented the split, then **reverted it**. The premise — that mean, sd and tail mass are length-robust so unwindowing is a free improvement — is false, and the measurement says so clearly.
- [x] Measured with realistic candidates: the actual driver, and traces from `simulate_warm(match_variance = TRUE)`, which is what the pool is really made of. Pass rates against the default bounds, same fitted model and same observed record in every row, varying only trace length:

  | trace length | mean pass | sd pass | tail pass |
  |---|---|---|---|
  | 21 | 94.9% | 70.0% | 0.4% |
  | 40 | 99.1% | 64.4% | 0.4% |
  | 60 | 100.0% | 59.1% | 0.5% |
  | 120 | 100.0% | 50.2% | 0.8% |

- [x] The criteria move in **opposite directions**, so unwindowing is not a neutral use of more evidence. The mean criterion loosens to vacuity (100% by 60 years) because `simulate_warm()` centres traces on the observed mean and a longer trace converges to it. The sd criterion tightens sharply (70% → 50%) because a longer trace converges to the *model's* sd while the benchmark stays a noisy 21-point sample estimate. The tolerances were implicitly calibrated against matched lengths; changing the lengths without recalibrating them changes what the filter selects for, in a direction nobody chose.
- [x] Caught one of my own measurements being unrepresentative before drawing from it. The first pass used raw `arima.sim()` traces, which have no variance matching and gave a 1.2% sd pass rate — an artifact of the fixture, not the package. Redone with real pool traces.
- [x] Kept the groundwork from `t2608151254a`: `compute_tailmass_metrics()` normalising each side by its own length is correct regardless, and stays.
- [x] Documented harmonisation in `filter_warm_pool()`'s `@details` — the step deferred from `t2608151254a`. Doing it here rather than opening another item, since the behaviour to document is the existing one and it is now settled.

## Outcome — dropped

Not the right change. Landing it would have loosened one criterion to vacuity
and tightened another, purely as a function of how long the simulation happens
to be, under tolerances tuned for matched lengths.

The underlying tension is real and unresolved: harmonisation means a 60-year
trace is judged on 20 of its years, which does discard evidence. But the fix is
not "use the full series"; it would be recalibrating the bounds as a function of
length, or expressing the criteria as standardised distances rather than raw
relative differences. That is a larger piece of work and needs its own item if
it is ever wanted.

## Separate finding — see `t2608151522`

The tail criterion passes **under 1% of candidates at every length tested**
(0.4% to 0.8%) against its default `tail_tol_log = log(1.03)`. It is rejecting
essentially the whole pool, and the pool only survives because the relaxation
loop then loosens the bound. That is independent of this item and of trace
length.

## Refs

- Closed sibling `t2608151254a` made the harmonisation deterministic and
  measured why it matters: with 60-year traces scored on 20-year windows, no
  trace passed in more than 52% of its 40 possible windows.
- `compute_tailmass_metrics()` used to divide both sides by the observed
  length; that was fixed there, which is what makes this item cheap.
