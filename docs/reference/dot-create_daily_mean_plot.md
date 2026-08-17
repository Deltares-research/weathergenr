# Create daily mean diagnostic plot

Faceted observed-vs-simulated comparison of daily mean values by
variable, using summary ranges and medians across grid cells/months.
Optionally saves the plot.

## Usage

``` r
.create_daily_mean_plot(
  daily_stats_season,
  plot_config,
  show_title,
  save_plots,
  out_dir
)
```

## Arguments

- daily_stats_season:

  Data frame of seasonal daily statistics including columns `stat`,
  `Observed`, `Simulated`, and `variable`.

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
