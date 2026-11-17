---
title: Journal figures for the perturbed climate
type: todo-item
sequence: 4
former_id: t2611161814
status: active
branch: feat/optional-temperature-extrema
area: plots
queue:
doc: dev/tasks/t004/brief.md
created: 2026-11-16
updated: 2026-11-16
---

> [!note] Overview
> **Problem** — The three perturbation-preview figures compare targets with outcomes for one cell of one realization, so they verify the code but do not show what the perturbed climate looks like, and none is usable in a journal article.
>
> **Outcome** — One dev script renders six ensemble-based figures of perturbed temperature and precipitation as PNG and PDF, and a short findings file says how much year-to-year precipitation variability the perturbation removes.
>
> **Size and unknowns** — Medium: one script, six figures, no package behaviour change. Unknown: whether forcing each month-year to the target mean flattens annual precipitation variability enough to need a separate fix.

## Progress

Next: owner reviews the figures and `findings.md` and decides on `enforce_target_mean`.

- [x] Measure the interannual-variability loss under `enforce_target_mean` and record it in `findings.md` —
  median CV ratio 0.012 (range 0.001–0.028) with enforcement on, 1.000 with it off; defect confirmed
- [x] Build `dev/scripts/perturbation-figures.R` with full and smoke configurations and a cached ensemble —
  `smoke` argument, `cache/` keyed on config and `SUMMARY_VERSION`, `run-info.txt` provenance
- [x] Render figures 1–3: mechanics, seasonal cycle, daily distribution — PNG, PDF and CSV in
  `output/perturbation-figures/`
- [x] Render figures 4–6: transient trajectories, preserved structure, response surfaces — plus
  `fig4b-transient-enforce-off`, written because phase 0 confirmed the loss
- [x] Inspect every PNG at final size, run the numeric sanity checks and lint — all seven PNGs inspected;
  checks in `findings.md` (fig 3 M1/M3 and fig 6 Rx1day fail, traced to enforcement); no lints
- [x] Document the script in `dev/README.md` and write `captions.md` — usage paragraph added; captions
  generated with the figures
- [ ] Owner reviews the figures and decides on the `enforce_target_mean` behaviour

## Work location

Current: primary checkout, `feat/optional-temperature-extrema`, uncommitted. The script builds on
`R/climate_perturbations_plots.R`, which exists only as untracked owner work in this checkout.

## Scope and acceptance

In scope: one dev script, six figures, a variance diagnostic. Out of scope: any change under `R/`,
including a fix to `enforce_target_mean`, and replacing `plot_climate_perturbations()`.

Done when: the script renders all six figures in PNG and PDF from the repository root, the sanity checks
in the brief are reported, and the owner has accepted the figures.

## Refs

- `dev/tasks/t004/brief.md` — the implementation brief: experiment design, per-figure
  specification, constraints and checks.
- `output/perturbation-preview/review.html` — the earlier figures this task replaces as article material.
