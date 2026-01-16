# Create daily standard deviation diagnostic plot

Faceted observed-vs-simulated comparison of daily standard deviations by
variable, using summary ranges and medians across grid cells/months.
Optionally saves the plot.

## Usage

``` r
.create_daily_sd_plot(
  daily_stats_season,
  plot_config,
  show_title,
  save_plots,
  output_path
)
```

## Arguments

- daily_stats_season:

  Data frame of seasonal daily statistics including columns `stat`,
  `Observed`, `Simulated`, and `variable`.

- plot_config:

  List of plotting configuration options (theme, alpha, subtitle, etc.).

- show_title:

  Logical; if `TRUE`, adds title/subtitle.

- save_plots:

  Logical; if `TRUE`, writes plot to `output_path`.

- output_path:

  Character; output directory for saved plots.

## Value

ggplot object (returned invisibly by
[`.export_multipanel_plot()`](https://deltares-research.github.io/weathergenr/reference/dot-export_multipanel_plot.md)).
