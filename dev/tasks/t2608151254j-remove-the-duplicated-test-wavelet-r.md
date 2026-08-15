---
title: Remove the duplicated test-wavelet.R
type: todo-item
status: backlog
effort: 1
area: tests
origin: review-2026-08-15
queue: 1
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — All 8 test_that() titles in test-wavelet.R are duplicated verbatim among test-warm.R's 24, and its header still cites the deleted R/wavelet.R.
> **Why** — Leftover from the wavelet_cwt.R / wavelet_warm.R split; duplicate tests cost suite time and make a real coverage gap hard to see.
> **Effort** — Small, but diff the two files before deleting rather than trusting the matching titles — same name does not guarantee same body.

## Progress

- [x] Diffed all 8 shared blocks. **Six byte-identical, two not — and in both cases `test-wavelet.R` held the stronger version.** Its `analyze_wavelet_spectrum` structure test asserts eight more return fields and their lengths (`gws_unmasked`, `gws_signif`, `gws_signif_unmasked`, `has_significance`, `signif_periods`, `gws_n_coi`, `neff_unmasked`, `comps_names`); its `simulate_warm` determinism test carries a `suppressWarnings()` and a comment explaining that the warnings surfaced only once `test-warm.R` stopped leaking a mocked ARMA fitter.
- [x] **Did not delete `test-wavelet.R`** — the review had the two files backwards, and following this note as written would have removed the better copy. `test-warm.R` is the one with the stale `R/wavelet.R` header; `test-wavelet.R` cites the current `R/wavelet_cwt.R` and `R/wavelet_warm.R`, and `git log` confirms it was last touched by `d067c8c` ("stop test-warm.R leaking mocks") — i.e. it is the newer file.
- [x] Split along the module layout instead of deleting either file. Removed the seven CWT and plotting blocks from `test-warm.R`, which `test-wavelet.R` already covers, and moved the stronger `simulate_warm` determinism test *into* `test-warm.R`, which owns that function. Headers rewritten on both; no `test_that()` title now appears in both files.
- [x] Suite count: 977 → 949, a drop of 28. Every removed block was byte-identical to, or weaker than, one that survives, so the drop is duplicate coverage only. Verified directly: the stronger structure assertions are still in `test-wavelet.R`, the stronger `simulate_warm` test is now in `test-warm.R`, and neither file mentions the other's functions.

## Verification

- `devtools::test(filter = "warm")` 138 pass, `filter = "wavelet"` 86 pass.
- `WEATHERGENR_BASELINE=1 devtools::test()` — 949 pass, 0 fail, **0 skip**. No
  numeric movement; the parallel baseline scenario ran rather than skipping.
- No shared `test_that()` titles remain between the two files.

## Note for the review draft

`dev/drafts/package-review-2026-08-15.md` § B8 says "`test-wavelet.R` is a
strict subset of `test-warm.R`" and attributes the stale `R/wavelet.R` header
to it. Both halves are wrong: it is a superset in two blocks, and the stale
header is in `test-warm.R`. The claim was made from matching titles without
diffing the bodies — which is exactly what this note's Effort line warned
against, and the warning was the only reason the coverage survived.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § B8.
