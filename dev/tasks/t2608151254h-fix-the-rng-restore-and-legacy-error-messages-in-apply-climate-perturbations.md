---
title: Fix the RNG restore and legacy error messages in apply_climate_perturbations()
type: todo-item
status: backlog
effort: 1
area: perturbations
origin: review-2026-08-15
queue: 13
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — It calls set.seed() with no restore, the only one of six seed sites that does not, and its validation messages name climate.data, sim.dates and change.factor.* which are not in the signature. The latitude check accepts lat or y but the message names only y.
> **Why** — A user-facing entry point silently hijacks the caller's RNG stream, and an error tells them to fix an argument that does not exist.
> **Effort** — Small; the only judgment is whether anyone downstream greps for the legacy message strings before they change.

## Progress

- [ ] Add save/restore around `set.seed(as.integer(seed))` at `R/climate_perturbations.R:193`. Copy the `adjust_precipitation_qm()` idiom (`R/quantile_mapping.R:260-272`) — it is the only one of the six that also handles the case where `.Random.seed` did not exist.
- [ ] Replace the legacy argument names in the validation messages (`:123-130`): `climate.data` to `data`, `sim.dates` to `date`, `change.factor.precip.mean` to `precip_mean_factor`, `change.factor.precip.variance` to `precip_var_factor`, `change.factor.temp.mean` to `temp_delta`.
- [ ] Fix the latitude message (`:145-148`) to name both accepted columns, `lat` or `y`.
- [ ] Grep `blueearth_cst` for the legacy strings before changing them, in case anything matches on message text.
- [ ] Add a test asserting the caller's `.Random.seed` is unchanged after a seeded call.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § B1, § B7.
- The four competing RNG idioms are tabulated in § B1; worth settling on one
  while this file is open.
