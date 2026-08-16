# Figure export standardization — design

Backing draft for `t2608161143`. Written 2026-08-16 from a full inventory of every
PNG-writing path in the package plus direct probes of the installed ggplot2 4.0.3 and
patchwork 1.3.2.

## 1. What is inconsistent today

`weathergenr` writes 17 PNGs per run from **6 `ggsave()` call sites**, each written
independently.

| # | Call site | Writes | Size | dpi |
|---|---|---|---|---|
| 1 | `R/generator.R:468` | `obs_power_spectra.png` | 8 × 5 literal | omitted |
| 2–4 | `R/generator.R:521-526` | `warm_annual_precip/stats/wavelet.png` | 8 × 5 via `pl_width`/`pl_height` (`:303-304`) | omitted |
| 5 | `.export_multipanel_plot()`, `R/evaluate_generator_plots.R:161-196` | 9 fixed names + `annual_pattern_<var>.png` | `ncol*4` × `nrow*4+0.5` | `plot_config$dpi` |
| 6 | `.create_annual_precip_plot()`, `:826-832` | `annual_precip.png` | 8 × 4 literal | **300 hardcoded** |

Observed symptoms, confirmed against a real run in `C:\TEMP\ntoum\10`:

- **Three base themes ship in one run.** `annual_precip.png` renders in ggplot2's
  default `theme_grey()` — it applies no theme at all. The WARM diagnostics use
  `theme_light()`, except `warm_annual_stats.png` which uses `theme_bw()`. The ten
  evaluation figures use `theme_bw(base_size = 12)`.
- **Five geometries**: 8×8.5, 8×5, 8×4.5, 8×4, 12×8.5 in.
- **`dpi`, `units`, `bg` omitted at four of six sites.**
- **`plot_dpi` / `plot_device` reach only 12 of 17 files.** They are threaded into
  `evaluate_weather_generator()` only; `generate_weather()` never receives them
  (`R/generator.R:888-889`).
- **Absolute font literals** — title 14, subtitle 10, legend 12 — that do not track
  `base_size`.

### Defects, not just inconsistency

1. `.create_annual_precip_plot()` applies no theme **and** bypasses the shared exporter
   **and** hardcodes `dpi = 300`. Three problems in one function.
2. `.export_multipanel_plot()` reads `p$facet$params$ncol/nrow`. `facet_grid` does not
   populate those fields, so `precip_conditional_correlations.png` is always written as
   if it were 2×2 — its three regime rows are squeezed into a two-row canvas.
3. `evaluate_generator_plots.R:672` uses `legend.position = c(1, 0)`, the
   pre-ggplot2-3.5.0 form. **Verified it does not currently warn** on ggplot2 4.0.3 —
   neither `ggplot_build()` nor `ggsave()` emits a deprecation — so this is hygiene,
   not a live bug, and no `DESCRIPTION` version pin is required for it.
4. **The export test has been blind since it was written.** `.export_multipanel_plot()`
   calls `do.call(ggplot2::ggsave, args)`. The namespace qualifier bypasses
   `local_mocked_bindings(.package = "weathergenr")`, so the mock at
   `tests/testthat/test-evaluate_generator.R:889` has only ever intercepted the *one*
   bare `ggsave` in `.create_annual_precip_plot()`. `expect_true(calls >= 1L)` passed
   on a single call while eleven others went unobserved — which is exactly why defect 1
   survived. **The new writer must call bare `ggsave`.**

## 2. Locked decisions

Settled with the owner before design:

| Decision | Choice |
|---|---|
| Multi-panel sizing | Derived from the panel grid, not fixed presets |
| House theme | `theme_light()`-based `theme_weathergenr(base_size = 12)` |
| PNG device | ggplot2's default. **No new dependency** — no `ragg`, no `cowplot`. Pass `units`, `dpi`, `bg` explicitly |
| Figure sizes | Stay internal — no new width/height arguments |
| `plot_dpi` / `plot_device` | Extended to `generate_weather()` so all 17 files honour them |
| Geom constants | In scope; font literals become `rel()` |
| Scope | Package-written PNGs only. Vignettes out of scope |

`pl_width` / `pl_height` are **not** public arguments — they are local constants under
a comment reading "no user-level change needed for now". Nothing to deprecate; figure
sizes are already internal, which is the intended end state.

## 3. Where the shared code lives

**New file `R/plot_export.R`.** AGENTS.md binds *module* plotting code to its module's
`*_plots.R` sibling. This code has no module — it is consumed by
`evaluate_generator_plots.R`, `warm_filtering_plots.R`, `wavelet_plots.R` and the
writer for `generator.R`. Putting it in any one sibling would make the other two depend
on an unrelated module's file.

It also *strengthens* the second half of that rule. Afterwards `generator.R` contains
zero `ggsave` and zero ggplot2 construction, and `evaluate_generator.R` no longer holds
a ggplot2 theme object — both are computational files that currently violate it.

Contents: `theme_weathergenr()`, `.PLOT_GEOM`, `.PANEL_SIZES`, `.figure_size()`,
`.panel_dims()`, `.export_figure()`. Bare ggplot2 calls throughout (`import(ggplot2)`
is already in NAMESPACE and bare is the package majority). No `aes()` appears in this
file, so the `.data` pronoun rule does not bite.

## 4. Panel-grid detection — verified, not inferred

An earlier design pass concluded that detection had to go through `p$facet$params`,
with `ggplot2::wrap_dims()` for the `nrow`-only case and a level-counting helper for
`facet_grid`, on the grounds that `p$facet$compute_layout()` is not callable. That
premise is true but the conclusion does not follow: `ggplot_build()` exposes the
computed layout directly. **Probed on ggplot2 4.0.3:**

| Layout | `ggplot_build(p)$layout$layout` |
|---|---|
| `facet_wrap(~g, ncol=2, nrow=2)` | ncol 2, nrow 2 |
| `facet_wrap(~g, nrow=1)`, 4 levels | ncol 4, nrow 1 |
| `facet_wrap(~g)` bare, 4 levels | ncol 2, nrow 2 |
| `facet_grid(r ~ v)`, 3 regimes × 2 vars | ncol 2, nrow 3 |
| unfaceted | ncol 1, nrow 1 |
| patchwork 2-panel | ncol 1, nrow 1 (no error) |

Every case is correct, including the two that motivated the `params` route:
`warm_annual_stats`'s `facet_wrap(~par, nrow = 1)` and the `facet_grid` defect. So:

```r
#' Infer a plot's panel grid as c(ncol, nrow)
#' @keywords internal
.panel_dims <- function(p) {
  # A patchwork governs its own split via plot_layout(); the outer canvas does
  # not need to know how many plots are inside it. An explicit branch because a
  # patchwork is also class "ggplot", and because patchwork exposes no stable
  # accessor for its grid -- attr(p, "patches") is NULL, p@patches errors, and
  # p[["patches"]] raises "Patchworks can only be indexed with numeric indices".
  if (inherits(p, "patchwork")) return(c(ncol = 1L, nrow = 1L))

  # One expression for facet_wrap, facet_grid and single-panel alike. Reading
  # p$facet$params$ncol/nrow instead is what made facet_grid fall through to a
  # hardcoded 2x2.
  tryCatch({
    lay <- ggplot_build(p)$layout$layout
    c(ncol = max(lay$COL), nrow = max(lay$ROW))
  }, error = function(e) c(ncol = 1L, nrow = 1L))
}
```

No `params` poking, no `wrap_dims()`, no level-counting helper, and no
`ggplot2 (>= 3.5.0)` pin needed for detection. `ggplot_build()` runs once per saved
figure — negligible against the render — and the `tryCatch` degrades a malformed plot
to a single panel rather than erroring out of the save.

## 5. Panel sizes

Panel size is a property of what is *inside* a panel, so it does not vary with grid
shape. An N×M grid is `M * panel_w` by `N * panel_h` plus a fixed margin allowance.

| family | panel W × H (in) | margin H | used by |
|---|---|---|---|
| `square` | 4.0 × 4.0 | 0.5 | faceted obs-vs-sim scatter, boxplot and monthly-cycle grids |
| `wide` | 8.0 × 4.0 | 0.5 | single-panel time series and spectra |
| `narrow` | 2.0 × 4.0 | 0.5 | one-row strips of single-category panels |

Reasoning anchored to what exists, not to taste:

- **`square` reproduces `width <- ncol*4; height <- nrow*4 + 0.5` exactly.** That is
  the current implicit panel size and it is the right one: these panels plot `Observed`
  against `Simulated` with a `geom_abline()` 1:1 line, which is only readable at
  aspect 1. Seven evaluation figures keep today's geometry bit-for-bit, so any size
  change that *does* appear is one someone chose.
- **`wide` 8.0 × 4.0** is the existing `annual_precip` literal and the existing
  `generator.R` width. All three current wide figures already approximate 2:1.
- **`narrow` 2.0 × 4.0** exists for `warm_annual_stats`: `facet_wrap(~par, nrow = 1)`,
  four panels each holding a *single* jittered column at `x = "a"` with `axis.text.x`
  blanked. Under `square` it would render 16 in wide for four columns of dots. At 2.0
  in per panel the four total 8.0 in — its width today.

**Header allowance.** The current `+0.5` is unconditional, but `show_title = TRUE` is
the pipeline default and adds a title *and* subtitle, three of which are two lines
(`crossgrid`, `intergrid`, `monthly_cycle`). `.figure_size()` takes `header = FALSE`
and adds **0.6 in** when a title was actually attached. This is the one place the new
geometry deliberately departs from today's arithmetic.

```r
.PANEL_SIZES <- list(
  square = list(panel_w = 4.0, panel_h = 4.0, margin_w = 0.0, margin_h = 0.5),
  wide   = list(panel_w = 8.0, panel_h = 4.0, margin_w = 0.0, margin_h = 0.5),
  narrow = list(panel_w = 2.0, panel_h = 4.0, margin_w = 0.0, margin_h = 0.5)
)
.HEADER_ALLOWANCE_IN <- 0.6

.figure_size <- function(nrow = 1L, ncol = 1L,
                         family = c("square", "wide", "narrow"),
                         header = FALSE, max_in = 30) {
  family <- match.arg(family)
  g <- .PANEL_SIZES[[family]]
  nrow <- max(1L, as.integer(nrow)); ncol <- max(1L, as.integer(ncol))
  width  <- ncol * g$panel_w + g$margin_w
  height <- nrow * g$panel_h + g$margin_h +
    if (isTRUE(header)) .HEADER_ALLOWANCE_IN else 0
  # ggsave() errors above 50 in; cap well below so a pathological panel count
  # degrades to a squashed figure rather than a failed run.
  c(width = min(width, max_in), height = min(height, max_in))
}
```

## 6. House theme

```r
#' @export
theme_weathergenr <- function(base_size = 12, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      plot.title    = element_text(size = rel(1.2)),   # 14.4pt at base 12 (was flat 14)
      plot.subtitle = element_text(size = rel(0.85))   # 10.2pt at base 12 (was flat 10)
    )
}
```

`base_size = 12` is load-bearing: `theme_light()`'s native 11 would silently shrink
every label relative to today's `theme_bw(base_size = 12)`. Changing geometry and type
size in one pass makes the visual diff uninterpretable.

Deliberately **not** in the house theme: `legend.text`. Exactly one exported figure has
a visible legend (`annual_pattern_<var>`; `monthly_cycle` and the wavelet field both
call `guides(... = "none")`), so legend styling stays local to that plot as
`element_text(size = rel(1))` — reproducing today's flat 12 at `base_size = 12` while
now tracking `base_size`.

### `plot_config` loses its `theme` field entirely

```r
plot_config <- list(
  subtitle = "Value range and median from all simulations shown against observed",
  dpi      = plot_dpi,
  device   = plot_device,
  alpha    = .PLOT_GEOM$alpha_summary,
  colors   = stats::setNames(c("blue3", "gray40"), c("Observed", "Simulated"))
)
```

Each builder calls `theme_weathergenr()` itself. This resolves the placement tension
without moving code — `evaluate_generator.R` no longer constructs a ggplot2 object —
and `plot_config` becomes pure data.

It also answers the wavelet / warm-filtering half: those modules never needed
`plot_config` plumbed to them. They call `theme_weathergenr()` directly, replacing
`theme_light()` at `wavelet_plots.R:137, 164, 237` and the mixed
`theme_light()` / `theme_bw()` / `ggplot2::theme_light()` at
`warm_filtering_plots.R:131, 153, 179`. **That is what unifies the three themes: one
function, every plot calls it** — rather than a config object plumbed to three modules.

## 7. Geom constants

```r
# One entry per role, not per call site, so the observed reference series has
# one weight in every figure.
.PLOT_GEOM <- list(
  lwd_observed  = 1.0,   # observed reference series
  lwd_simulated = 0.5,   # one simulated member
  lwd_summary   = 1.5,   # stat_summary linerange spanning the ensemble
  lwd_reference = 0.8,   # abline / hline / significance threshold
  pt_summary    = 2,     # ensemble median marker
  pt_member     = 3,     # single-member marker
  alpha_faint   = 0.2,   # ribbons, member bundles, boxplot fill
  alpha_member  = 0.3,
  alpha_summary = 0.4    # today's plot_config$alpha
)
```

Full mapping of current literals:

| current | where | role | new |
|---|---|---|---|
| `linewidth = 1.5` ×7 | `evaluate_generator_plots.R` lineranges | ensemble range | `lwd_summary` (unchanged) |
| `linewidth = 1.25` ×2 | `:753` monthly_cycle, `:814` annual_precip | observed | `lwd_observed` (**1.25 → 1.0**) |
| `linewidth = 0.9` | `warm_filtering_plots.R:135` | observed | `lwd_observed` (**0.9 → 1.0**) |
| `linewidth = 0.7` | `warm_filtering_plots.R:133` | member | `lwd_simulated` (**0.7 → 0.5**) |
| `linewidth = 0.8` ×5 | `warm_filtering_plots.R:154,188,194,200` | reference | `lwd_reference` (unchanged) |
| `size = 2` ×7 | median points | summary marker | `pt_summary` (unchanged) |
| `size = 3` ×2 | `:665` mean diamond, `warm_filtering_plots.R:155` jitter | member marker | `pt_member` (unchanged) |
| `alpha = 0.5` | `:805` annual_precip bundle | member bundle | `alpha_faint` (**0.5 → 0.2**) |
| `alpha = 0.8` | `:743` monthly_cycle bundle | member bundle | `alpha_faint` (**0.8 → 0.2**) |

The two alpha moves are the visible ones — both are bundles of N realization lines
where 0.8 saturates to solid at 20 realizations. If the owner prefers to freeze
appearance, add an `alpha_bundle = 0.5` role rather than special-casing a call site.

## 8. The writer

```r
#' @keywords internal
.export_figure <- function(p, filename, output_dir,
                           save_plots = TRUE,
                           show_title = FALSE, title = NULL, subtitle = NULL,
                           family = c("square", "wide", "narrow"),
                           dims = NULL,   # c(ncol, nrow), overrides detection
                           size = NULL,   # c(width, height) in, overrides everything
                           dpi = 300, device = NULL,
                           plot_config = NULL) {
  family <- match.arg(family)
  header <- isTRUE(show_title) && !is.null(title)
  if (header) p <- p + labs(title = title, subtitle = subtitle)
  if (!isTRUE(save_plots) || is.null(output_dir)) return(invisible(p))

  if (!is.null(plot_config$dpi))    dpi    <- plot_config$dpi
  if (!is.null(plot_config$device)) device <- plot_config$device

  if (is.null(size)) {
    if (is.null(dims)) dims <- .panel_dims(p)
    size <- .figure_size(nrow = dims[["nrow"]], ncol = dims[["ncol"]],
                         family = family, header = header)
  }

  args <- list(filename = file.path(output_dir, filename), plot = p,
               width = unname(size[1L]), height = unname(size[2L]),
               units = "in", dpi = dpi, bg = "white")
  # Only pass `device` when set: ggsave() infers it from the extension
  # otherwise, and an explicit NULL is not the same as omitting it.
  if (!is.null(device)) args$device <- device

  # Bare `ggsave` (via import(ggplot2)), NOT ggplot2::ggsave -- the qualified
  # form bypasses local_mocked_bindings(.package = "weathergenr"), which is how
  # the export test came to be observing only one of the six call sites.
  do.call(ggsave, args)
  invisible(p)
}
```

## 9. Call-site changes

| Site | Change |
|---|---|
| `generator.R:468` | `.export_figure(..., size = c(9.0, 4.5))` — patchwork exposes no stable grid accessor, and the 3:1.2 layout means `panel_w * ncol` would overstate it. 9.0 splits as ≈6.4 in field + ≈2.6 in companion |
| `generator.R:521-526` | `family = "wide"` / `"narrow"` / `"wide"`; delete `pl_width`/`pl_height` at `:303-304` |
| `.export_multipanel_plot()` | Replaced by `.export_figure()`; the nine builders swap `plot_config$theme` → `theme_weathergenr()` (kept as the **first** `+` term, before any local `theme()` override — the existing `:658` then `:671` order is already correct) and rename the exporter call |
| `.create_annual_precip_plot()` | Route through the exporter, apply the theme, delete its own `show_title` block (the exporter owns it, which is what makes the header allowance correct) |
| `:672`, `:675` | `legend.position = "inside"` + `legend.position.inside` + `legend.justification.inside`; `legend.text = rel(1)` |

`generate_weather()` gains `plot_dpi = 300`, `plot_device = NULL`, NULL-coalesced
alongside `relax_priority`, forwarded from `run_weather_generator()` mirroring the
existing `:888` pattern.

### Resulting sizes (`show_title = TRUE`, the default)

| output | grid | old | new |
|---|---|---|---|
| daily_mean, daily_sd, crossgrid_correlations, annual_pattern_*, monthly_cycle | 2×2 | 8 × 8.5 | 8 × 9.1 |
| spell_length, wetdry_days_count | 2×1 | 8 × 4.5 | 8 × 5.1 |
| intergrid_correlations | 3×2 | 12 × 8.5 | 12 × 9.1 |
| precip_conditional_correlations | 3 rows × n pairs | 8 × 8.5 (**wrong**) | derived |
| annual_precip | 1×1 wide | 8 × 4 | 8 × 5.1 |
| warm_annual_precip, warm_annual_wavelet | 1×1 wide | 8 × 5 | 8 × 4.5 |
| warm_annual_stats | 4×1 narrow | 8 × 5 | 8 × 4.5 |
| obs_power_spectra | patchwork | 8 × 5 | 9 × 4.5 |

## 10. Tests

New `tests/testthat/test-plot_export.R`, plus edits to four existing files.

**`.figure_size()`** — `.figure_size(2,2,"square") == c(8, 8.5)`, `(1,2,"square") ==
c(8, 4.5)`, `(2,3,"square") == c(12, 8.5)`. These three pin the "square reproduces
today's arithmetic" claim, which is the intended-no-change assertion. Plus:
`header = TRUE` adds exactly 0.6 to height and nothing to width; panel size is
grid-shape invariant; `(40,40)` caps at `max_in`; `match.arg` errors on an unknown
family.

**`.panel_dims()`** — one fixture per layout. The `facet_grid(regime ~ variable)` case
with 3 regimes and 2 pairs returning `ncol 2, nrow 3` is the regression test for the
reported defect; today's code returns 2×2 there. Also the `facet_wrap(~par, nrow = 1)`
case, the unfaceted case (including `p$data` being a waiver), and patchwork returning
the fallback without erroring.

**`.export_figure()`** — mock **bare** `ggsave` with `.package = "weathergenr"`,
capturing `list(...)` keyed by `basename(filename)`. Assert on every captured call:
`units == "in"`, `bg == "white"`, non-NULL `dpi`. Assert `device` is **absent from the
names** when NULL (`expect_false("device" %in% names(args))`) and present otherwise.
Assert `size` overrides `family`, `dims` overrides detection, `plot_config$dpi` beats
an explicit `dpi`, and `save_plots = FALSE` writes nothing while still returning `p`.

**`test-evaluate_generator.R`** —
- Delete the `theme = ggplot2::theme_bw(base_size = 12)` fixtures at `:733` and `:882`.
  They go inert under the reshape, and an inert fixture that still passes is how this
  regresses.
- `expect_false(any(vapply(plot_config, inherits, logical(1), "theme")))` against the
  config the pipeline actually builds — binds the reshape.
- Loop over every plot from `create_all_diagnostic_plots()` asserting
  `expect_identical(p$theme$panel.border$colour, theme_weathergenr()$panel.border$colour)`.
  **Verified this discriminates all three current themes**: `theme_light` `#B3B3B3FF`,
  `theme_bw` `#333333FF`, `theme_grey` `NULL`. `panel.background$fill` does *not* work
  — `theme_light` and `theme_bw` both give `"white"`. This is the assertion that would
  have caught `annual_precip` rendering in `theme_grey`.
- Rewrite the mock test at `:877-907` from a **counter** to a **recorder**: exact
  filename set, per-file geometry, and one run with `plot_dpi = 150` asserting every
  captured `dpi` is 150 — which binds the threading `annual_precip` currently ignores.
  Note in the commit message *why* the count changes: the old mock was
  namespace-bypassed and only ever saw one site.

**`test-wavelet_plots.R`, `test-warm_filtering.R`** — the same house-theme assertion on
`plot_wavelet_global_spectrum()`, the `p_gws` panel of `plot_wavelet_power()`, and all
three panels of `plot_filter_diagnostics()`. Without these the two modules that most
needed unifying stay unbound — and `plot_filter_diagnostics`'s `stats` panel is exactly
the one that used `theme_bw()` while its siblings used `theme_light()`.

**Deliberately not done:** no `vdiffr` snapshots. `test-wavelet_plots.R:1-4` already
records the reasoning — new dependency, SVG drifts across three OSes and three R
versions — and it applies unchanged.

## 11. Verification

```r
testthat::test_file("tests/testthat/test-plot_export.R")
devtools::test()
devtools::document()
source("tools/build_site_tools.R"); check_only()   # exports + docs change
```

```bash
Rscript tools/lint.R --changed
```

Then **a manual visual pass over the 17 regenerated PNGs** — the check that actually
matters for a styling change, and the one a green suite cannot substitute for.
`inst/extcode/wegen_pipeline.R` is the fastest route; compare against `C:\TEMP\ntoum\10`.

**The end-to-end baseline gate does not apply.** `tests/testthat/helper-baseline.R`
runs every scenario with `save_plots = FALSE` and its documentation states plots are
excluded from the fingerprint as "ggplot2-version-sensitive". Nothing here touches
resampling, WARM, wavelets, calendar logic or evaluation statistics. State this in the
commit message so the omission reads as a decision rather than a skip.

## 12. NEWS.md surface

**New features** — `theme_weathergenr()`, the exported house theme now applied to every
figure. `generate_weather()` gains `plot_dpi` and `plot_device`, so all 17 PNGs honour
them rather than only the evaluation set.

**Bug fixes** — `annual_precip.png` rendered in ggplot2's default grey theme and
ignored `plot_dpi`/`plot_device`. `precip_conditional_correlations.png` was sized as a
2×2 grid regardless of its real panel count.

**Other** — figure sizes derived from the panel grid; sizes remain internal. Shared
line weights, point sizes and alphas. Legend position modernized.

## 13. Open judgement calls

1. **The 0.6 in header allowance** is the only intentional change to the square family.
   `.HEADER_ALLOWANCE_IN <- 0` freezes geometry exactly, at the cost of the crowding
   that motivated it.
2. **`precip_conditional_correlations.png` gets much taller** — 3 rows honoured
   properly. That is the honest rendering, but 4 in wide × 13 in tall for one variable
   pair is an awkward page shape. If it reads badly the fix is an explicit `size` at
   that call site, not a change to detection.
3. **Two alphas and three line weights change visibly.** Unification necessarily means
   someone's current appearance changes.
4. **Whether `theme_weathergenr()` is exported.** Recommend yes —
   `evaluate_weather_generator()` returns ggplot objects users re-render, and without an
   exported theme they cannot match the package's look. But it interacts with
   `t2608151254k`, which is trimming the export surface; coordinate.
5. **`.panel_dims()` reads a ggplot2 internal structure.** `ggplot_build()$layout$layout`
   is stable in practice but not a documented API. A rename would silently return the
   fallback and produce wrong-sized figures rather than erroring. Mitigated by the six
   detection tests, which fail loudly, and by the `dims`/`size` escape hatches.
6. **The patchwork branch cannot introspect at all** — three accessor styles verified
   failing on patchwork 1.3.2. `obs_power_spectra.png` carries a hardcoded
   `size = c(9.0, 4.5)`, a judgement about the 3:1.2 layout rather than a derivation.
   If `plot_layout(widths = ...)` changes, that number needs a manual update and no
   test will catch it.
7. **`.data` pronoun conversion is out of scope.** The rule governs new plot code; the
   new helpers contain no `aes()`, and converting existing mappings would balloon the
   diff and pull `R/globals.R` in for no gain.
