---
title: Stop the WARM pool filter deciding from plotting-only spectra
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
> **What** — filter_warm_pool() feeds identify_significant_peaks() the unmasked global wavelet spectrum, which analyze_wavelet_spectrum()'s own roxygen labels plotting only, while every inference decision elsewhere uses the COI-masked curve.
> **Why** — The two halves of the package disagree about whether the record has significant low-frequency structure: on the packaged fixture the run log reports no significant periodicities while the pool filter simultaneously requires candidates to reproduce a 5.84-year peak.
> **Effort** — Small in code (which curve is passed), but it moves numeric output, so the cost is the baseline gate and reviewing the delta. The judgement is whether COI masking is right for a 20-year record, where it discards much of the usable spectrum.

## Progress

- [x] **The curve choice was not the defect.** Measured all three combinations and they gave the *same* answer — 1 peak at 5.84 yr — because the threshold is passed through `fill_nearest()` before use. Switching curves alone changes nothing. The real defect is that the significance threshold is fabricated at scales where no test exists.
- [x] Found three separate fabrication mechanisms, all of which had to go for the masked NAs to survive: `gws_regrid()`'s `rule = 2` flat extrapolation (`:1299`), `fill_nearest()` on the threshold (`:1300`), and a second `fill_nearest()` inside `identify_significant_peaks()` (`:788`). The last is in an exported function, so this is an API-behaviour change, not an internal one.
- [x] Split the two roles rather than picking one curve. **Peak power stays unmasked** — `compute_peak_match_metrics()` compares it against simulated spectra by log-ratio, and those are unmasked by construction in `.compute_gws_batch()`, so it must stay like-for-like. **The threshold becomes the masked inference curve**, unfilled, regridded with `rule = 1` so it stays `NA` outside its support. Comparison quantity and inference quantity are different things.
- [x] Measured the cost of masking: on the 20-year driver only **4 of 14 scales** have enough degrees of freedom inside the cone of influence, and the longest testable period is **3.47 yr**. The 5.84-yr peak was never testable; its threshold had been carried over from 3.47 yr.
- [x] Both baseline scenarios go from 1 significant peak to 0 (water_year 4/14 testable, calendar_year 5/14).
- [x] Log now states when the criterion has nothing to test, and reports the testable-scale count. Moved it out of `compute_spectral_metrics()` — which takes no `verbose` — into `filter_warm_pool()`, which does, after a first attempt raised six R-level warnings across the suite.
- [x] Corrected the roxygen on both unmasked fields. `gws_unmasked` is not "plotting only": it is the right curve for comparing two spectra computed the same way, which is what this filter does. `gws_signif_unmasked` genuinely is plotting-only, and now says why — it floors the effective sample size at 1 so a plot has no gaps, which makes it permissive exactly where it matters.
- [ ] **Owner decision pending** — baseline delta reviewed but *not* re-recorded. See below.

## Baseline delta — awaiting acceptance

`Rscript tools/record_baseline.R --dry-run`: **15 of 106 keys differ** — all in
`water_year`, all in realization 2. `calendar_year` is untouched, and
`water_year` realization 1 is untouched. Every `config.*` and structural key
holds.

That is a much smaller delta than it might have been, and the reason is worth
stating: removing a criterion only changes the outcome where it was rejecting
the trace that would otherwise have been selected. It was not, in three of the
four realization slots.

## What was actually removed, and the question it raises

The peak criterion was doing real work. On a 300-trace pool at the baseline
settings it admitted **149 of 300 (49.7 percent)**; it now admits **300 of 300**.
Half the pool was being rejected for failing to reproduce a 5.84-year spectral
peak whose significance could not be established from a 20-year record.

What survives is the other half of the `wavelet` criterion,
`spectral_cor >= 0.60` — a correlation between log spectra, which is a genuine
similarity test and needs no significance machinery. So the filter is not
toothless; it lost the half that rested on an invalid basis.

The underlying issue is that **the peak criterion is a similarity filter
implemented as a hypothesis test**. Two coherent futures, and they are different
features:

1. *Accept as landed.* The ensemble is no longer required to reproduce a
   spectral feature that cannot be shown to be real. Pool roughly doubles;
   `peak_match_frac_min = 1.0` becomes vacuous and its documented "require all
   significant peaks to match" should be reworded.
2. *Restore the filtering power on an honest basis.* Match the top-N observed
   spectral peaks by power, with no cone of influence and no chi-square —
   an explicit similarity criterion rather than a significance one.

Landing (2) inside a task titled "stop deciding from plotting-only spectra"
would be delivering a different feature than the title describes, so it is left
as a separate decision.

Re-record with `Rscript tools/record_baseline.R --force` once accepted.

## Refs

- `dev/drafts/package-review-2026-08-15.md` — surfaced while closing
  `t2608151253b`, not part of the original review.
- Measured on the packaged fixture (20-year water-year driver), significant
  observed peaks via the **unmasked** spectrum: 1 at `warm_signif` 0.80 and
  0.90 (period 5.84 yr), 0 at 0.95. Via the **masked** spectrum
  `has_significance` is FALSE at all three.
- Related: `t2608151253b` documented `warm_signif`'s live effects and named this
  as the reason the pool filter's threshold bites while the log says otherwise.
