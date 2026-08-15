---
title: Make the WARM filter's observed benchmark deterministic
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
> **What** — filter_warm_pool() scores candidates against a single random contiguous window of the observations, and returns full-length traces that were scored on a window.
> **Why** — Acceptance thresholds move with the seed rather than only the draw, and a trace admitted on a 40-year slice ships as a 100-year trace whose full-length statistics were never tested.
> **Effort** — Small code change, but it moves numeric output and the right statistic (median over windows? full record?) is a judgment call, not a lookup.

## Progress

- [x] Quantified the defect before choosing a fix. **Observed side:** across the windows of the 21-year record truncated to 15, the observed sd spans **19.1%** and the 80th percentile 7.6% — against tolerances of 3%. **Simulated side, which is worse:** with 60-year traces scored on 20-year windows, **no trace passes in more than 52% of its 40 possible windows** (median 14%), and for **186 of 200** at least one window says pass and another says fail. The accept decision was mostly a property of the draw, not of the trace.
- [x] Found the sharper defect underneath: each candidate got its *own independent* random window, so candidates were never compared on the same footing at all. That is worse than the criterion merely moving with the seed.
- [x] Replaced both draws with a deterministic tail window — the last `n_use` values of each series, the same position for every candidate. Kept harmonisation rather than removing it: the spectral criterion genuinely needs like-for-like, since a 60-year spectrum resolves low-frequency power a 20-year one cannot, and comparing them unmatched would penalise length rather than shape.
- [x] Hardened `compute_tailmass_metrics()`, which divided **both** sides by `n_use * scale_obs` where `n_use = length(obs_use)`. That silently required equal lengths — a column sum over `n_sim` rows normalised by the observed length. Now each side is normalised by its own length. Numerically identical while the caller harmonises, so no output moves; it makes the function correct for unequal lengths rather than merely unused with them.
- [x] **No baseline movement.** Verified the change is live rather than inert: the calendar_year scenario does exercise obs windowing (21 obs, 20 sim, 2 windows) and now reports `obs_window = 2-21`. The old random draw agreed with that for only 1 of 5 sample seeds, so the recorded baseline happened to sit on a seed where the two coincide. The seed dependence is gone at zero delta.
- [ ] Documenting the harmonisation in `filter_warm_pool()`'s `@details` — deferred with the open question below, since what to document depends on that answer.

## Verification

- `WEATHERGENR_BASELINE=1 devtools::test()` — 960 pass, 0 fail, **0 skip**.
- `record_baseline.R --dry-run` — "No change: the current output matches the
  stored baseline." No re-record needed.
- `check_only()` — 0 errors, 0 warnings, 0 notes.
- `tools/lint.R --changed` — 9 findings, all pre-existing on `t2608061641`.
- Three new tests: the observed window is seed-independent, every candidate now
  shares one simulated window, and tail mass is normalised per side.

## Left open — see `t2608151507`

Harmonisation is only needed by the **spectral** criterion. Mean, sd and tail
mass are length-robust estimators (tail mass now genuinely so, after the
denominator fix), so they could be computed on the full record and the full
trace instead of on a 20-year slice of each. That would stop a 60-year trace
being judged on a third of its length for three of the four criteria, and use
strictly more data. It moves numeric output, so it is a separate landing.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C3.
- Note this interacts with `t2608151253a`: the mean-filter inertness measured
  there holds cleanly only on the unwindowed path, which is the default since
  `n_years` defaults to the observed year count (`R/generator.R:330`).
