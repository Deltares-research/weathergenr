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

- [x] Baseline captured before any edit — gate green at `5e76d51`, CSV copied aside for delta comparison.
- [x] Decided: `.variance_match_matrix()` must **not** touch the mean. It now rescales each column about its own mean, and the `target_mean` parameter is gone (internal, 3 call sites, all updated). The pinning was not load-bearing — it was the inconsistency: corrected columns were pinned to exactly the target while uncorrected ones kept their own mean, so the two paths disagreed about how the mean is carried. With preservation they agree.
- [x] `.fast_col_sd()` now returns the sample sd (`* n/(n-1)`, `NA` for `n < 2`), matching the `stats::sd()` targets it is compared against. Kept the `m2 - m1^2` form: cancellation costs ~2 of ~16 significant digits on annual totals, nowhere near the tolerances in play.
- [x] Pass rates re-checked on the packaged fixture, 3,000 candidates: mean 99.5% → 98.4%, sd 72.4% → 71.2%, both 72.4% → 70.0%. The pool does not starve and the default 3% bounds need no revisiting.
- [x] Verified the per-component branch separately. Both the packaged fixture and my probe take the Tier 2 (Cholesky) path, which **skips** per-component matching — so neither exercised `:474-484`. Forced it with a single variable component: the total mean reconstitutes to within 0.004%, no double-counting, no NaN.
- [ ] **Owner decision pending** — baseline delta reviewed but *not* re-recorded. See below.

## Baseline delta — awaiting acceptance

`Rscript tools/record_baseline.R --dry-run`: **60 of 106 keys differ.**

What did *not* move is the reassuring part: every `config.*` key, `gen.n_days`,
`gen.n_realizations`, `gen.hash_sim_dates`, `gen.rlzN.n_na`, `eval.fit_columns`
and `eval.fit_n_rows`. Structure and shape are intact; only values moved.

What moved is one cascade, not 60 findings: different WARM traces are selected,
so different analogue years are drawn, so every daily and evaluation statistic
downstream shifts. `gen.rlzN.n_distinct_sources` rises in all four cases
(3796→3999, 3901→4014, ...) — slightly more diverse resampling.

Evaluation MAEs move in both directions. Some improve (`mae_spell_Wet`
0.641→0.412, `mae_days_Dry` 0.569→0.433), some worsen (`mae_mean_precip`
0.253→0.309, `mae_sd_precip` 0.864→1.189). The worsening ones are expected and
are the point: 60% of the old candidate pool had its mean set to *exactly* the
observed mean, so whichever traces were selected scored near-perfectly on mean.
That score was measuring the artifact. With `n_realizations = 2` these metrics
are dominated by which two traces were drawn, so the direction of any single one
carries little signal.

Direct evidence on the baseline driver (20-year area-average annual precip,
3,000 candidates):

| | pinned mean | sd(trace means) | mean filter | sd filter |
|---|---|---|---|---|
| old | 1805 / 3000 | 0.059 | 99.5% | 72.4% |
| new | 0 / 3000 | 0.093 | 98.4% | 71.2% |

Re-record with `Rscript tools/record_baseline.R --force` once accepted.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C9 and § C10 — the measured
  numbers (1048/2000 traces mean-pinned; every clamped trace at exactly
  `target * sqrt(40/39)`, +1.274 percent).
- `.compute_gws_batch()` (`R/warm_filtering.R:1060`) already applies
  `* (n1/(n1-1))` and calls it "sample variance" — the convention to match.
