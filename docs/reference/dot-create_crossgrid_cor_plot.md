# Create cross-grid correlation diagnostic plot

Observed-vs-simulated comparison of cross-grid correlations, faceted by
the first variable in each correlation pair. Uses dummy points to
enforce symmetric axes.

## Usage

``` r
.create_crossgrid_cor_plot(
  stats_crosscor,
  plot_config,
  show_title,
  save_plots,
  out_dir
)
```

## Arguments

- stats_crosscor:

  Data frame of cross-grid correlation summaries with columns
  `Observed`, `Simulated`, and `variable1`.

- plot_config:

  List of plotting configuration options (subtitle, alpha, colors, dpi,
  device). The theme is not among them: every builder calls
  [`theme_weathergenr`](https://deltares-research.github.io/weathergenr/reference/theme_weathergenr.md)
  directly.

- show_title:

  Logical; if `TRUE`, adds title/subtitle.

- save_plots:

  Logical; if `TRUE`, writes plot to `out_dir`.

- out_dir:

  Character; output directory for saved plots.

## Value

ggplot object (returned invisibly by `.export_figure()`).
