---
title: Stop the WARM pool filter deciding from plotting-only spectra
type: todo-item
status: backlog
effort: 1
area: warm
origin: review-2026-08-15
queue: 3
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — filter_warm_pool() feeds identify_significant_peaks() the unmasked global wavelet spectrum, which analyze_wavelet_spectrum()'s own roxygen labels plotting only, while every inference decision elsewhere uses the COI-masked curve.
> **Why** — The two halves of the package disagree about whether the record has significant low-frequency structure: on the packaged fixture the run log reports no significant periodicities while the pool filter simultaneously requires candidates to reproduce a 5.84-year peak.
> **Effort** — Small in code (which curve is passed), but it moves numeric output, so the cost is the baseline gate and reviewing the delta. The judgement is whether COI masking is right for a 20-year record, where it discards much of the usable spectrum.

## Progress

- [ ] Decide which curve the pool filter should use. The COI-masked `gws`/`gws_signif` are what `analyze_wavelet_spectrum()` calls the inference curves and what `has_significance` is built from (`R/wavelet_cwt.R:695-699`); the unmasked pair is roxygen-documented "plotting only" (`:383`, `:386`) yet consumed at `R/warm_filtering.R:1285-1286`.
- [ ] Weigh the cost of masking on short records: the cone of influence removes the long-period end, which on 20 years is most of the low-frequency band the filter cares about. Masking may be correct and still leave the `wavelet` criterion with nothing to test — check before switching.
- [ ] Whichever is chosen, make the run log agree with it. Today `generate_weather()` logs "No significant low-frequency periodicities detected" (masked) while the filter requires a 5.84-year peak match (unmasked), which reads as a contradiction to anyone watching the console.
- [ ] Baseline gate before and after; re-record only once the delta is reviewed.
- [ ] If the unmasked curves stay in use for a decision, correct their "plotting only" roxygen.

## Refs

- `dev/drafts/package-review-2026-08-15.md` — surfaced while closing
  `t2608151253b`, not part of the original review.
- Measured on the packaged fixture (20-year water-year driver), significant
  observed peaks via the **unmasked** spectrum: 1 at `warm_signif` 0.80 and
  0.90 (period 5.84 yr), 0 at 0.95. Via the **masked** spectrum
  `has_significance` is FALSE at all three.
- Related: `t2608151253b` documented `warm_signif`'s live effects and named this
  as the reason the pool filter's threshold bites while the log says otherwise.
