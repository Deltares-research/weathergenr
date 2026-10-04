---
title: set.seed() discipline is held by convention, not construction
type: todo-item
status: backlog
effort: 1
area: warm
origin: review-2026-08-15
queue:
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What changes** — .set_seed_fixed_kind() is called only from resample.R; six other set.seed() sites use the bare form. None currently runs on a PSOCK worker -- .compute_gws_batch()'s parLapplyLB body is fully deterministic -- so this is not a live bug.
> **Why it matters** — AGENTS.md marks the rule IMPORTANT because a bare set.seed() under clusterSetRNGStream makes one seed mean two streams. Nothing enforces it today.
> **What it takes** — Start when: Any function containing a bare set.seed() gains a worker call site -- most likely simulate_warm() gaining parallel = TRUE.

## Progress

- [ ] Reassess once the start condition in Effort holds

## Refs

- `dev/drafts/package-review-2026-08-15.md` § B1, § B2 — the six bare sites and
  the four competing save/restore idioms.
- Bare sites: `R/generator.R:258`, `R/wavelet_warm.R:215` and `:400`,
  `R/warm_filtering.R:203`, `R/wavelet_cwt.R:586`,
  `R/evaluate_generator.R:164`, `R/climate_perturbations.R:193`.
