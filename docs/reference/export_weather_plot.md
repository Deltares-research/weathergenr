# Export a weather diagnostic figure

Saves a returned ggplot or patchwork using the package's automatic
geometry, or explicit dimensions in millimetres. No scientific data are
recomputed.

## Usage

``` r
export_weather_plot(
  plot,
  filename,
  width = NULL,
  height = NULL,
  dpi = 300,
  device = NULL,
  bg = "white"
)
```

## Arguments

- plot:

  A ggplot or patchwork object.

- filename:

  Output path, including extension. Its parent must exist.

- width, height:

  Positive dimensions in millimetres. Either may be omitted to use that
  dimension from automatic geometry. Package-produced plots retain their
  figure family and explicit canvas size; other plots use square panels.

- dpi:

  Positive raster resolution. Default 300.

- device:

  Graphics device function or name; NULL infers from extension.

- bg:

  Background color. Default white; use "transparent" for transparency.

## Value

The exported plot, invisibly, with effective export typography.

## Plot styles

Plotting entry points accept `plot_style = NULL` or a named list with:

- theme:

  A partial
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html).
  Applied after local settings, so unspecified layout choices survive.
  Explicit font sizes win over automatic export scaling. No global theme
  or options are changed.

- colors:

  A named character vector, with any of `observed`, `simulated`,
  `reference`. Observed covers measured series; simulated covers
  ensemble summaries, monochrome members and ribbons; reference covers
  identity/zero guides, significance thresholds, contours and
  cone-of-influence lines. Omitted roles retain each plot's original
  colors, which can differ between figures. Invisible layout points do
  not carry semantic colors.

- ensemble_palette:

  Color vector interpolated across realization IDs in the evaluation
  monthly-cycle plot. Other ensemble plots use the simulated color.
  Omitted: retain the existing hue scale.

- power_palette:

  Color vector interpolated across log2 wavelet power. Existing limits,
  squishing and white missing-value cells are preserved. Omitted: retain
  the existing viridis C scale.

These are presentation settings only. YAML parsing and any
consumer-specific translation belong to the caller. Unknown fields and
color roles are errors.

## Examples

``` r
p <- plot_wavelet_global_spectrum(1:4, rep(2, 4), c(1, 3, 2, 1))
export_weather_plot(p, tempfile(fileext = ".png"), width = 160, height = 100)
```
