---
title: Standardize figure export geometry and theming
type: todo-item
status: active
effort: 2
area: plots
queue: 1
branch: refactor/plot-export-standardization
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

- [x] Wrote `R/plot_export.R` — `theme_weathergenr()` (exported), `.PANEL_SIZES`,
      `.PLOT_GEOM`, `.figure_size()`, `.panel_dims()`, `.export_figure()`
- [x] Landed `tests/testthat/test-plot_export.R` before the migration — 49 assertions,
      including `.figure_size(2,2,"square") == c(8, 8.5)`, which is what proves the
      square family reproduces the old `ncol*4 x nrow*4+0.5` arithmetic
- [x] Dropped `theme` from `plot_config`; `R/evaluate_generator.R` now holds no ggplot2
      object, and every builder calls `theme_weathergenr()` itself
- [x] Migrated `R/evaluate_generator_plots.R` — nine builders, `.create_annual_precip_plot()`
      fixed (theme + shared exporter + honours `plot_dpi`), legend modernized to
      `legend.position = "inside"`
- [x] Migrated `R/warm_filtering_plots.R` and `R/wavelet_plots.R`; the lone `theme_bw()`
      panel is gone
- [x] Converted the four `R/generator.R` writes; `generate_weather()` gained
      `plot_dpi`/`plot_device`, forwarded from `run_weather_generator()`
- [x] Tests, `document()`, `NEWS.md`; regenerated all 17 PNGs from the shipped Ntoum
      NetCDF fixture and reviewed them. `devtools::test()` 1213 pass / 0 fail,
      `check_only()` 0/0/0, lint clean
- [ ] Merge to `master` — held for owner review of the regenerated figures

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
- **Settled 2026-08-16:** `theme_weathergenr()` **is exported**. Owner's call, taken
  with [[t2608151254k]] in view — that item trims 46 exports to 41 on an unmerged
  branch, and this adds one back. `evaluate_weather_generator()` returns ggplot objects
  users re-render, and without an exported theme they cannot match the package's look.
- **Settled 2026-08-16:** the header allowance stays at **0.6 in** rather than 0, so
  titled figures grow 0.6 in taller and today's cramped title block gets room.

## Found during implementation

Two things the brief did not anticipate, both settled in the diff:

1. **`.panel_dims()` sizes to the panels a figure actually draws, not the grid it
   requested.** `facet_wrap(~variable, ncol = 2, nrow = 2)` with fewer than three
   variables renders a smaller grid, and `ggplot_build()` reports that. The old
   exporter read the requested 2x2 and wrote a full-size canvas holding one panel and
   three empty slots. With the usual four variables both agree exactly; with fewer,
   the new behaviour is better. Asserted in `test-evaluate_generator.R` so it is a
   decision rather than an accident.

2. **`alpha_bundle = 0.5` was added as a role after looking at the output.** The brief
   folded realization-line bundles into `alpha_faint = 0.2`. That is right for a
   ribbon or a boxplot interior but leaves `monthly_cycle`'s hue-scaled lines almost
   invisible — the comparison the figure exists for stopped working. An area fill and
   a line bundle are genuinely different roles; the draft anticipated this exact
   escape hatch.

**Also worth knowing:** lintr's `object_usage_linter` resolves symbols against the
**installed** package, not the working tree, so every new internal object or function
reports "no visible binding" until `devtools::install()` runs. Three symbols were
flagged here purely for that reason. This is the same class of trap as the
`devtools::load_all()`-does-not-reach-PSOCK-workers note in AGENTS.md, and cost a
wrong diagnosis before it was spotted.

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
