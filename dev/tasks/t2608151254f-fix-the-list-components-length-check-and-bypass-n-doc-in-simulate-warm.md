---
title: Fix the list-components length check and bypass_n doc in simulate_warm()
type: todo-item
status: backlog
effort: 1
area: warm
origin: review-2026-08-15
queue: 2
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — Line 297 tests any(lens != n) where the matrix and data.frame branches correctly test n_fit, so list-form components are rejected whenever n differs from the observed length. Separately the roxygen says bypass_n defaults to 25 while the formal is 15L.
> **Why** — Verified: matrix input returns 25x3 where identical list input errors. The existing list test uses n == n_fit, which is why it has gone unnoticed.
> **Effort** — Small and self-contained: a one-token fix plus a doc correction. No numeric-output risk, since `generate_weather()` passes a matrix and never reaches the broken branch.

## Progress

- [ ] Change `any(lens != n)` to `any(lens != n_fit)` at `R/wavelet_warm.R:297`, matching the matrix (`:275`) and data.frame (`:284`) branches.
- [ ] Extend the existing list-input test to use `n != n_fit` — it currently passes only because it uses `n == n_fit`, which is why this survived.
- [ ] Correct the `bypass_n` roxygen (`:55`), which says the default is 25 while the formal (`:103`) is `15L`.
- [ ] Run `devtools::test(filter = "warm")`, then `devtools::document()` for the roxygen change.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C8 — verified against the working
  tree: 40 observed years with `n = 25` gives `25 x 3` for matrix input and
  "All component vectors must have length 'n'." for the identical list.
