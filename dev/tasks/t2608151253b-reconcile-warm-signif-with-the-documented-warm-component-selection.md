---
title: Reconcile warm_signif with the documented WARM component selection
type: todo-item
status: backlog
effort: 2
area: warm
origin: review-2026-08-15
queue: 1
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — warm_signif is documented as retaining low-frequency components but only nudges n_levels; cwt_to_modwt_map is computed, returned and never consumed, so every MODWT band is ARMA-modelled regardless of significance.
> **Why** — The documented behaviour does not exist and the method differs from the Steinschneider and Brown 2013 framing the README cites, with nothing disclosing it.
> **Effort** — Large, and the size is in the decision rather than the diff: gating on significance changes the generator's numeric output and needs a defence against the source paper, whereas correcting the docs is an afternoon.

## Progress

- [x] Read both papers from the owner's Zotero library (main text, PDF extraction). **Kwon et al. (2007) §2, eq. (1)-(2):** the series is `sum_k R_kt + e_t` — K signal components *plus* a residual stochastic process — and the model is `sum_k AR(R_kt, p_k) + AR(e_t, p)`, i.e. an AR fit to each retained signal **and** one to the residual. Selection: Morlet CWT, global wavelet power tested against a red-noise (AR(1)) background at the one-sided 95% limit; scales exceeding it "are then selected as candidates for reconstruction". Orders by AIC, coefficients by ML. **Steinschneider & Brown (2013) §2.1** adopts this verbatim (their eq. 1-2) and their own application retained a single signal (H = 1, significant periods 1-4 yr). The band-selection detail sits in their supporting information, not the main text; Kwon is the citable source for the criterion. Nothing is discarded in either: the non-retained part *is* the residual term.
- [x] Decided: **do not gate; correct the documentation.** Two reasons, both measured below. (1) This package decomposes with MODWT MRA, which is exactly additive — reconstruction error 9.5e-13 on the fixture — so there is no residual term for a non-significant band to be relegated to. Dropping `D1` would discard 22.6% of the variance outright, where the original would have kept it as `e_t`. (2) The CWT test never fires on records this short.
- [x] Measured on the packaged fixture (20-year water-year driver, `analyze_wavelet_additive`): 2 components, `D1` (period 2.83 yr, 22.6% of variance) and `S1` (period 4 yr, 59.2%). `has_significance` is **FALSE at `warm_signif` 0.80, 0.90 and 0.95 alike**, and `cwt_to_modwt_map` is empty at all three — so the `n_levels` branch never fires either. A naive gate would keep only `S1`. The fitted structure on this record is therefore ARMA(`D1`) + ARMA(`S1`), the same *shape* as the original with one retained signal — but that is a consequence of how few levels 20 points support, not a correspondence between a MODWT band and `e_t`: a longer record gives `D1`, `D2`, `S2` and no band plays that role.
- [x] Variance identity checked rather than assumed: fractions sum to 0.818, and the shortfall is exactly the cross-covariance (`2 * cross_cov / total = 0.1818`), with `variance_identity_error` at -4.8e-13. MRA components are additive but not orthogonal.
- [x] No baseline gate needed — documentation only, no numeric change. Full suite 954 pass / 0 fail.
- [x] Decided **against** splitting `warm_signif` into a separate pool-filter argument. It has one genuinely live effect (the `filter_warm_pool` peak threshold) and one that is inert on realistic records (`n_levels`), so a second argument would add API surface for little gain — and `blueearth_cst` pins this package by Git tag, making every new argument a contract. Both effects are now documented instead. Overturn if a user needs the two thresholds to differ.
- [x] Rewrote the `warm_signif` roxygen and added a WARM section to `generate_weather()`'s Details giving the relation to Kwon (2007) and Steinschneider & Brown (2013), with both references. Marked `cwt_to_modwt_map` diagnostic-only in `analyze_wavelet_additive()` and corrected its `signif` parameter.
- [x] Surfaced a **separate** defect while measuring, opened as `t2608151353`: the pool filter reads the *unmasked* global wavelet spectrum, which the API documents as plotting-only, while `has_significance` uses the COI-masked curve. That is why the run log says "no significant periodicities" while the filter requires a 5.84-year peak match. Kept out of this item — it moves numeric output and is a different decision.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C1 — what is verified (the dead
  `cwt_to_modwt_map`) versus what is not (the practical consequence).
- `R/wavelet_modwt.R:527-546` (n_levels selection) and `:557-577` (the map that
  is computed and never consumed).
