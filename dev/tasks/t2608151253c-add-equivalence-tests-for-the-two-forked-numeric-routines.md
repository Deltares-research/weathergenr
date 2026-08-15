---
title: Add equivalence tests for the two forked numeric routines
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
> **What** — Assert .compute_gws_batch() matches analyze_wavelet_spectrum() and .knn_draw_one_rank() matches knn_sample() on small fixtures.
> **Why** — Both forks are hand-copied for speed and each IS the scientific criterion in its path; a change to the canonical routine leaves the fork behind with the suite still green.
> **Effort** — Small and mechanical. The one unknown is whether the forks currently agree — if either has already drifted, this item uncovers a live numeric bug rather than just fencing one off.

## Progress

- [x] GWS test written in `tests/testthat/test-warm_filtering.R`. `.compute_gws_batch()` reproduces `analyze_wavelet_spectrum()$gws_unmasked` to **2e-16 (detrend off) and 4.4e-15 (detrend on)** — machine precision, no drift. Added a second test for column independence and chunk invariance, since the batching is where a per-trace answer could silently pick up a neighbour's.
- [x] KNN test written in `tests/testthat/test-resample.R` — but **not** as the note framed it. The step said "assert the same index is drawn"; that assertion would have been false. The two functions deliberately differ in the draw: `knn_sample()` uses `sample.int(prob = )`, the fork inverts the cumulative rank probabilities against a uniform drawn once for the whole simulation. Same seed, different rank. What is genuinely shared — and what the tests now assert — is the weighted squared distance, the neighbour selection and ordering, and the rank probabilities, on **both** the partial-sort (`k < 0.2 * nc`) and full-order branches. Ordering is compared by driving the fork across rank bins and by ranking `knn_sample()`'s 40,000 draws by frequency.
- [x] Neither fork had drifted, so this fenced them off rather than uncovering a bug.
- [x] Comments added at both fork sites naming the test that holds them together and what to re-check when the canonical routine changes.
- [x] **Refactor required to make the KNN fork testable at all:** `.knn_draw_one_rank()` was a closure defined inside `resample_weather_dates()`, unreachable from a test. Confirmed it closes over nothing (arguments and base functions only) and promoted it to a file-level internal next to `knn_sample()`, so the two now sit adjacent. Behaviour-preserving.
- [x] Verified the promotion on the worker path, which is the part that could have broken silently: PSOCK workers load the *installed* package, so `devtools::install()` first, then the full suite with the baseline gate enabled — **977 pass, 0 fail, 0 skip**. Zero skips is the evidence: the gate's parallel scenario ran rather than skipping, so the moved function resolves in a worker and the baseline is unchanged.

## Verification

- `devtools::test(filter = "resample")` 72 pass, `filter = "warm_filtering"` 53 pass.
- `devtools::install()` then `WEATHERGENR_BASELINE=1 devtools::test()` — 977 pass,
  0 fail, **0 skip**. No baseline movement; no re-record needed.
- `Rscript tools/lint.R --changed` — 13 findings, all pre-existing and already
  logged on `t2608061641` (the `.Random.seed` bindings and the `.log()` glue
  false positives). The promoted function contributes none.
- No `man/` or `NAMESPACE` change: the new internal is `@noRd`.
- No `NEWS.md` entry: tests and an internal refactor are not user-visible.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § A3.
- The forks assert their own equivalence in comments only:
  "CWT grid (identical to analyze_wavelet_spectrum internals)" and
  "mirrors analyze_wavelet_spectrum exactly".
