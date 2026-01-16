# Create monthly cycle diagnostic plot

Plots annual cycles by month for each variable: simulated cycles are
shown for all realizations, while the observed cycle is shown as a
single reference line.

## Usage

``` r
.create_monthly_cycle_plot(
  daily_stats_season,
  plot_config,
  show_title,
  save_plots,
  output_path
)
```

## Arguments

- daily_stats_season:

  Data frame of daily seasonal statistics including mean values.

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
