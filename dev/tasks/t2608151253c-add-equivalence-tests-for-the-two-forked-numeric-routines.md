---
title: Add equivalence tests for the two forked numeric routines
type: todo-item
status: backlog
effort: 1
area: tests
origin: review-2026-08-15
queue: 3
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — Assert .compute_gws_batch() matches analyze_wavelet_spectrum() and .knn_draw_one_rank() matches knn_sample() on small fixtures.
> **Why** — Both forks are hand-copied for speed and each IS the scientific criterion in its path; a change to the canonical routine leaves the fork behind with the suite still green.
> **Effort** — Small and mechanical. The one unknown is whether the forks currently agree — if either has already drifted, this item uncovers a live numeric bug rather than just fencing one off.

## Progress

- [ ] Write the GWS test: run `.compute_gws_batch()` (`R/warm_filtering.R:1012`) and `analyze_wavelet_spectrum()` on the same short series, assert the global wavelet spectra match to tolerance.
- [ ] Write the KNN test: same candidates, target, `k` and seed through `.knn_draw_one_rank()` (`R/resample.R:~355`) and `knn_sample()`, assert the same index is drawn.
- [ ] If either fails, stop and treat it as a numeric bug — baseline gate, then fix — rather than adjusting the test to pass.
- [ ] Add a comment at both fork sites pointing at the test, so the next person to edit the canonical routine sees what holds them together.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § A3.
- The forks assert their own equivalence in comments only:
  "CWT grid (identical to analyze_wavelet_spectrum internals)" and
  "mirrors analyze_wavelet_spectrum exactly".
