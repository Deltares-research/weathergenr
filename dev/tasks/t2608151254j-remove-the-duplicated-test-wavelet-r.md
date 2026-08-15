---
title: Remove the duplicated test-wavelet.R
type: todo-item
status: backlog
effort: 1
area: tests
origin: review-2026-08-15
queue: 13
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — All 8 test_that() titles in test-wavelet.R are duplicated verbatim among test-warm.R's 24, and its header still cites the deleted R/wavelet.R.
> **Why** — Leftover from the wavelet_cwt.R / wavelet_warm.R split; duplicate tests cost suite time and make a real coverage gap hard to see.
> **Effort** — Small, but diff the two files before deleting rather than trusting the matching titles — same name does not guarantee same body.

## Progress

- [ ] Diff the 8 shared `test_that()` blocks in `tests/testthat/test-wavelet.R` against their namesakes in `test-warm.R`; keep anything that differs.
- [ ] Delete `test-wavelet.R` once it is confirmed redundant.
- [ ] Decide whether the surviving file should be split to match the module layout — `AGENTS.md` says `test-<module>.R`, and `test-warm.R` now covers both `wavelet_cwt.R` and `wavelet_warm.R` while `test-wavelet_modwt.R` sits alongside. The stale `R/wavelet.R` header comment is the fossil of that split.
- [ ] Run `devtools::test()` and confirm the count drops by exactly the removed tests.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § B8.
