---
title: Reconcile warm_signif with the documented WARM component selection
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
> **What** — warm_signif is documented as retaining low-frequency components but only nudges n_levels; cwt_to_modwt_map is computed, returned and never consumed, so every MODWT band is ARMA-modelled regardless of significance.
> **Why** — The documented behaviour does not exist and the method differs from the Steinschneider and Brown 2013 framing the README cites, with nothing disclosing it.
> **Effort** — Large, and the size is in the decision rather than the diff: gating on significance changes the generator's numeric output and needs a defence against the source paper, whereas correcting the docs is an afternoon.

## Progress

- [ ] Read Steinschneider & Brown (2013) and Kwon et al. (2007) on how the significant bands are selected and what happens to the remainder — do not reconstruct this from memory.
- [ ] Decide: gate component selection on `cwt_to_modwt_map`, or keep modelling all bands and correct the documentation to say the significance test is diagnostic. Both are defensible; only silence is not.
- [ ] Measure before choosing — simulate with and without the non-significant bands and compare annual variance and the global wavelet spectrum. I did not test this; the contribution may be small.
- [ ] If gating: baseline gate before and after, since it moves numeric output.
- [ ] Either way, split `warm_signif`'s third role — it is also passed as `filter_warm_pool`'s `wavelet_args$signif_level` (`R/generator.R:439`). One knob should not mean three things.
- [ ] Update the `warm_signif` roxygen (`R/generator.R:70`), which currently claims it retains components.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C1 — what is verified (the dead
  `cwt_to_modwt_map`) versus what is not (the practical consequence).
- `R/wavelet_modwt.R:527-546` (n_levels selection) and `:557-577` (the map that
  is computed and never consumed).
