---
title: Decide whether the WARM peak criterion should be a similarity filter
type: todo-item
status: backlog
effort: 2
area: warm
origin: review-2026-08-15
queue: 3
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — The peak criterion was a similarity filter implemented as a hypothesis test. Now that the invalid significance basis is gone it admits every candidate, so decide whether to leave it vacuous or re-express it as an explicit top-N spectral-peak match with no cone of influence and no chi-square.
> **Why** — It previously rejected half the pool. That filtering power is now gone, and peak_match_frac_min = 1.0 still reads as maximally strict while testing nothing.
> **Effort** — Large in decision, small in diff. The code change is a different peak-selection rule; the question is whether an ensemble should be required to reproduce spectral features that a short record cannot establish are real.

## Progress

- [ ] Decide the intent first. Is the goal *inference* — only ask traces to reproduce structure the record supports — or *similarity*, matching the observed spectrum whether or not it is distinguishable from red noise? Both are defensible for a stress-test ensemble; they are different features and the current defaults do not say which is meant.
- [ ] If similarity: replace the significance test in `identify_significant_peaks()` (or add a sibling) with a top-N-by-power selection over a period band, dropping the cone of influence and the chi-square entirely. Do not reintroduce a threshold and call it significance.
- [ ] If inference: leave the criterion vacuous on short records and reword `peak_match_frac_min` in `filter_warm_bounds_defaults()`, which still reads "require all significant peaks found to match" while nothing is found. Consider dropping it from `relax_order`, since a filter at a 100 percent pass rate is never selected by `which.min(rates)` anyway.
- [ ] Either way, check the interaction with `spectral_cor_min = 0.60`, which is now the whole of the `wavelet` criterion. If the peak half stays vacuous, that threshold may want revisiting to carry the load.
- [ ] Baseline gate before and after if the selection rule changes.

## Refs

- Closed sibling `t2608151353` removed the invalid significance basis and
  measured the cost: the criterion admitted 149 of 300 candidates before and
  300 of 300 after.
- On the packaged 20-year driver only 4 of 14 scales are testable inside the
  cone of influence; the longest testable period is 3.47 yr, while the peak the
  criterion used to enforce sat at 5.84 yr.
- `dev/drafts/package-review-2026-08-15.md` for the surrounding review.
