---
title: set.seed() discipline is held by convention, not construction
type: watch-item
area: warm
origin: review-2026-08-15
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — .set_seed_fixed_kind() is called only from resample.R; six other set.seed() sites use the bare form. None currently runs on a PSOCK worker -- .compute_gws_batch()'s parLapplyLB body is fully deterministic -- so this is not a live bug.
> **Why** — AGENTS.md marks the rule IMPORTANT because a bare set.seed() under clusterSetRNGStream makes one seed mean two streams. Nothing enforces it today.
> **Trigger** — Any function containing a bare set.seed() gains a worker call site -- most likely simulate_warm() gaining parallel = TRUE.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § B1, § B2 — the six bare sites and
  the four competing save/restore idioms.
- Bare sites: `R/generator.R:258`, `R/wavelet_warm.R:215` and `:400`,
  `R/warm_filtering.R:203`, `R/wavelet_cwt.R:586`,
  `R/evaluate_generator.R:164`, `R/climate_perturbations.R:193`.
