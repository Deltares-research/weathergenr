---
title: Document or reshape the Markov spell-length factors
type: todo-item
status: backlog
effort: 1
area: resample
origin: review-2026-08-15
queue: 9
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — dry_spell_factor and wet_spell_factor divide transition probabilities and renormalise, so they also shift the unconditional wet-day fraction, and the extreme-state row p20/p21/p22 is never adjusted.
> **Why** — dry_spell_factor = 1.5 does not mean 50 percent longer dry spells and the docs give no way to calibrate; the asymmetry on the extreme state is undocumented.
> **Effort** — Small if the answer is to document the transformation honestly; larger if the arguments become target spell lengths, which means solving for the multiplier and moves numeric output.

## Progress

- [ ] Decide: document what the factors actually do, or re-express them as target spell lengths (or a target mean spell multiplier) and solve for the transition adjustment.
- [ ] Either way, state the wet-day-fraction coupling — dividing p01/p02 and renormalising raises p00, which lengthens dry spells *and* lowers the unconditional wet-day probability. The two knobs are not orthogonal.
- [ ] Decide what should happen to the extreme-state row. `p20`/`p21`/`p22` is currently untouched by either factor (`R/resample.R:1073-1085`), so persistence out of the extreme state ignores both. Deliberate or oversight?
- [ ] Update the `dry_spell_factor`/`wet_spell_factor` roxygen in `generate_weather()` (`R/generator.R:88-91`).
- [ ] Baseline gate if the transformation changes; docs-only needs no gate.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C6.
