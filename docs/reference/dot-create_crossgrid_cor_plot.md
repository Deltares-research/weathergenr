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
  output_path
)
```

## Arguments

- stats_crosscor:

  Data frame of cross-grid correlation summaries with columns
  `Observed`, `Simulated`, and `variable1`.

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
