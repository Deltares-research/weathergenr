---
title: Resolve the order-dependent Cholesky whitening in simulate_warm()
type: todo-item
status: backlog
effort: 2
area: warm
origin: review-2026-08-15
queue: 4
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — chol() whitens sequentially, so the fitted ARMAs are no longer scale-specific and the result depends on component ordering. Either switch to a symmetric square root or document the limitation; the header comment at lines 17-20 now contradicts the Cholesky path.
> **Why** — Scale-specific persistence is the premise of WARM, and the ordering that decides the mixing is a convention rather than a physical fact.
> **Effort** — Large because it is a method question, not a bug. The algebra as written is correct; what is in doubt is whether whitening across scales is the right move at all in a method whose premise is scale separation.

## Progress

- [ ] Confirm the framing before changing anything: `chol()` gives a triangular factor, so Z1 is D1 scaled while ZJ mixes every component, and the mixing depends on component order (D1..DJ, S). The re-correlation itself is correct — do not "fix" the algebra.
- [ ] Evaluate a symmetric square root (eigendecomposition, `Sigma^{-1/2}`) as a drop-in for `chol()` at `R/wavelet_warm.R:354-359`; same cost, no ordering dependence.
- [ ] Compare the two on the ensemble's global wavelet spectrum against the observed — if the difference is negligible, the honest outcome is to keep `chol()` and document the choice.
- [ ] Fix the stale header comment (`:17-20`), which still says independent modelling does not preserve cross-component covariance "unless you implement coupling externally" — the Cholesky path is exactly that coupling.
- [ ] Baseline gate before and after if the factorisation changes.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C2.
- The code already half-acknowledges the issue at `:407-412` ("no longer smooth
  regardless of its name").
