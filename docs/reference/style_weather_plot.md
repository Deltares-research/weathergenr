# Apply the compact diagnostic presentation profile

Restyles returned package plots without regenerating weather or
statistics. This explicit opt-in leaves all plotting entry-point
defaults unchanged.

## Usage

``` r
style_weather_plot(plot, id, plot_style = NULL, panel_layout = NULL)
```

## Arguments

- plot:

  A returned weather diagnostic ggplot or patchwork.

- id:

  Stable returned plot ID (not a filename).

- plot_style:

  Partial theme and semantic colors, as documented in
  [`export_weather_plot`](https://deltares-research.github.io/weathergenr/reference/export_weather_plot.md).
  Explicit theme entries win over this profile.

- panel_layout:

  Optional named arguments for
  [`arrange_weather_panels`](https://deltares-research.github.io/weathergenr/reference/arrange_weather_panels.md),
  excluding `plots`. NULL uses its defaults for faceted diagnostics.
  Local wavelet and WARM summary layouts are preserved.

## Value

A styled ggplot or patchwork with automatic export geometry.

## Details

Observations and agreement guides use `observed` (blue by default).
Members and ribbons use `simulated` (gray); central ensemble summaries
use black. Significance thresholds and cone-of-influence guides use
`reference` (red); significant wavelet regions retain black contours.
These semantics apply only to this profile. Boxplots retain medians,
quartiles, whiskers and outliers; their diamonds are means. Comparison
circles are ensemble medians. No new ensemble summary is computed. Line
widths and point sizes are fixed in millimetres, fonts in points.
Existing explicit font sizes in `plot_style` override the profile. Use
ordinary ggplot additions for further adjustments. Unit labels use
Unicode; for vector export choose a Unicode-capable device such as
[`grDevices::cairo_pdf`](https://rdrr.io/r/grDevices/cairo.html).

## Examples

``` r
# Minimal synthetic annual diagnostic with two supplied series:
annual <- data.frame(year = 1:6, Simulated = c(4, 5, 3, 4, 6, 5),
  Observed = c(4.2, 5.1, 3.4, 4.3, 5.8, 5.2))
p <- ggplot2::ggplot() + theme_weathergenr() +
  ggplot2::geom_line(data = annual[c("year", "Simulated")],
    ggplot2::aes(year, Simulated)) +
  ggplot2::geom_line(data = annual[c("year", "Observed")],
    ggplot2::aes(year, Observed))
styled <- style_weather_plot(p, "annual_precip")
```
