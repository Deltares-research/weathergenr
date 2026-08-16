---
title: Standardize figure export geometry and theming
type: todo-item
status: backlog
effort: 2
area: plots
queue: 1
doc: dev/drafts/plot-standardization-2026-08-16.md
created: 2026-08-16
updated: 2026-08-16
---

> [!note] Overview
> **What** — The 17 PNGs a run writes come from six independently-written `ggsave()`
> calls, so they ship in three different ggplot2 themes and five different page
> geometries; `annual_precip.png` gets no theme at all and renders in ggplot2's grey
> default. Replace this with one house theme and one rule that derives figure size
> from the panel grid, and route every write through a single function.
> **Why** — These are the figures anyone judges a run by, and they currently look like
> output from three different packages. Two are also wrong rather than merely
> inconsistent: `precip_conditional_correlations.png` squeezes three rows of panels
> into a two-row canvas, and the documented `plot_dpi`/`plot_device` settings silently
> never reach five of the files.
> **Effort** — Large in reach rather than difficulty: one new file plus edits to six
> existing ones, all mechanical once the two helpers exist. The real unknown is
> aesthetic, not technical — every figure changes dimensions, and whether the new
> proportions look right can only be settled by regenerating the set and looking at it.

## Progress

- [ ] Write `R/plot_export.R` — `theme_weathergenr()`, `.PANEL_SIZES`, `.PLOT_GEOM`,
      `.figure_size()`, `.panel_dims()`, `.export_figure()`
- [ ] Land `tests/testthat/test-plot_export.R` **before** the migration; the
      `.figure_size(2,2,"square") == c(8, 8.5)` assertion is what proves the square
      family is behaviour-preserving
- [ ] Drop `theme` from `plot_config` in `R/evaluate_generator.R`, so a computational
      file stops holding a ggplot2 object
- [ ] Migrate `R/evaluate_generator_plots.R` — nine builders onto `theme_weathergenr()`
      and `.export_figure()`, fix `.create_annual_precip_plot()`, modernize the legend
- [ ] Migrate `R/warm_filtering_plots.R` and `R/wavelet_plots.R` onto the house theme
      and the shared geom constants
- [ ] Convert the four `R/generator.R` writes; add `plot_dpi`/`plot_device` to
      `generate_weather()` and forward from `run_weather_generator()`
- [ ] Update the four affected test files, `devtools::document()`, `NEWS.md`, then
      regenerate the 17 PNGs and look at them

## Why this is not just tidying

Three of the items above fix defects, not style:

1. **`.create_annual_precip_plot()`** applies no theme, bypasses the shared exporter,
   and hardcodes `dpi = 300` — so `plot_dpi` and `plot_device` do not reach it.
2. **`.export_multipanel_plot()` is blind to `facet_grid`.** It reads
   `p$facet$params$ncol/nrow`, which `facet_grid` does not populate, so
   `precip_conditional_correlations.png` is always written as if it were 2×2.
3. **The export test has been blind since it was written.** The helper calls
   `do.call(ggplot2::ggsave, args)`; the namespace qualifier bypasses
   `local_mocked_bindings(.package = "weathergenr")`, so the mock at
   `tests/testthat/test-evaluate_generator.R:889` has only ever intercepted the one
   bare `ggsave` in `.create_annual_precip_plot()`. That is why defect 1 survived.
   **The new writer must call bare `ggsave`** or the tests stay blind.

## Decisions already settled

Locked with the owner before design, so implementation does not need to reopen them:

- Sizing derives from the panel grid, not fixed presets.
- House theme is `theme_light()`-based at `base_size = 12` — matching what the
  evaluation pipeline already produces, so this change moves geometry only.
- **No new dependency.** Keep ggplot2's default PNG device; no `ragg`, no `cowplot`.
  But pass `units`, `dpi` and `bg` explicitly at every write.
- Figure sizes stay internal — no new width/height arguments. `pl_width`/`pl_height`
  at `R/generator.R:303-304` are local constants, not public API; there is nothing to
  deprecate.
- `plot_dpi`/`plot_device` are extended to `generate_weather()` (additive, defaults
  preserve behaviour) so all 17 files honour them.
- Geom constants and `rel()` font sizes are in scope. Vignettes are not.

## Open calls for whoever picks this up

- The 0.6 in header allowance is the only intentional change to the square family;
  setting it to 0 freezes geometry exactly.
- `precip_conditional_correlations.png` becomes much taller once its rows are honoured.
  That is the honest rendering, but the page shape may want an explicit `size`.
- Whether `theme_weathergenr()` is exported interacts with [[t2608151254k]], which is
  trimming the export surface. Recommend exporting — `evaluate_weather_generator()`
  returns ggplot objects users re-render — but coordinate rather than deciding alone.

## Verification

`devtools::test()` → `devtools::document()` → `check_only()` → `tools/lint.R --changed`,
then a **manual visual pass over the regenerated PNGs**, which is the check that
actually matters here and the one a green suite cannot substitute for.

The end-to-end baseline gate **does not apply**: `tests/testthat/helper-baseline.R`
runs every scenario with `save_plots = FALSE` and excludes plots from the fingerprint
as ggplot2-version-sensitive. Say so in the commit message so the omission reads as a
decision rather than a skip.

## Refs

- `dev/drafts/plot-standardization-2026-08-16.md` — full design: the six call sites,
  verified `ggplot_build()$layout$layout` probe results for all six facet layouts, the
  panel-size table with its reasoning, the complete geom-constant mapping, and the
  test plan.
- Related open item [[t2608151254k]] (trim the export surface).
