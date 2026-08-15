---
title: Give PET a diurnal-range response to warming
type: todo-item
status: backlog
effort: 1
area: perturbations
origin: review-2026-08-15
queue: 3
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — temp_delta is added to temp, temp_min and temp_max alike, so temp_max minus temp_min is invariant and Hargreaves responds to warming only through the (T + 17.8) term. Add a temp_range knob or document the choice.
> **Why** — PET drives the drought limb of the response surface; +4 degC currently raises PET by only about 12 percent, systematically damping the signal a bottom-up assessment is built to explore.
> **Effort** — Small to add the argument; the real question is what default is defensible, since observed DTR trends are regionally mixed and a wrong default is worse than none.

## Progress

- [ ] Decide between adding a `temp_range_factor` (or separate `temp_min_delta`/`temp_max_delta`) and documenting the fixed-DTR choice as a stated limitation.
- [ ] If adding a knob, default it to 1 (current behaviour) so existing runs are unchanged, and follow the length-12-or-matrix convention the other factors use.
- [ ] Check `AGENTS.md` on whether PET belongs in the Hargreaves family at all for this use — a radiation- or humidity-aware method would respond to warming properly, but adds dependencies the package deliberately avoids. Consult `~/workspace/brain` on PET method choice under climate change before deciding.
- [ ] Document the sensitivity either way in `apply_climate_perturbations()`'s `@details`: at T = 15 degC, +4 degC currently gives about +12 percent PET.
- [ ] `NEWS.md` entry — this is user-visible whichever way it goes.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C4.
- `R/climate_perturbations.R:465-467` (the three identical deltas);
  `R/pet.R:105` (the Hargreaves form — the implementation itself is correct
  FAO-56 and is not what needs changing).
