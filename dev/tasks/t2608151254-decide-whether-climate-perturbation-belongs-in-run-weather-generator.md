---
title: Decide whether climate perturbation belongs in run_weather_generator()
type: todo-item
status: backlog
effort: 2
area: api
origin: review-2026-08-15
queue: 3
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — apply_climate_perturbations() is unreachable from either entry point, and the interfaces do not meet: it needs temp_min/temp_max per cell, while generate_weather() validates only precip/temp and returns dates rather than values.
> **Why** — The README calls the three components coupled; component 3 is hand-assembled by every user. Either wire it up or stop claiming the coupling.
> **Effort** — Large if wired up, since `generate_weather()` would have to return values rather than dates and carry temp_min/temp_max; trivial if the answer is to correct the README instead. Decide that first.

*Deferred 2026-08-15 at the owner's request — moved to the back of the queue,
not blocked. Nothing external is waiting on it; the scope question below is a
design decision the owner wants to take separately.*

## Progress

- [ ] Decide the scope: full pipeline integration, a documented recipe function, or a README correction. The middle option — a thin `apply_perturbations_to_realization()` that does the date-to-values join — may be the honest amount of coupling.
- [ ] If integrating: settle how `temp_min`/`temp_max` reach stage 3. `generate_weather()` validates only `precip`/`temp` (`R/generator.R:208`) and returns resampled dates, so the join is real work, not a wiring change.
- [ ] Check what `vignettes/climate_perturbation.qmd` hand-assembles today — that sequence is the de facto spec for whatever replaces it.
- [ ] Whatever is chosen, make `README.md`'s "three coupled components" claim true, and add a `NEWS.md` entry if anything user-visible changes.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § A1.
- Interface gap: `R/climate_perturbations.R:150` (required columns) versus
  `R/generator.R:208` (validated vars) and `:576` (return shape).
