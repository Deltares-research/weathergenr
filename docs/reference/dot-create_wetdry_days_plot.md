# Create wet/dry day count diagnostic plot

Observed-vs-simulated comparison of average monthly wet and dry day
counts, faceted by statistic. Uses dummy points to enforce symmetric
axes per facet.

## Usage

``` r
.create_wetdry_days_plot(
  stats_wetdry,
  plot_config,
  show_title,
  save_plots,
  output_path
)
```

## Arguments

- stats_wetdry:

  Data frame of wet/dry diagnostics including day-count statistics.

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
