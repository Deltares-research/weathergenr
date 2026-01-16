# Create inter-variable correlation diagnostic plot

Observed-vs-simulated comparison of inter-variable correlations, faceted
by variable. Uses dummy points to enforce symmetric axes per facet.

## Usage

``` r
.create_intergrid_cor_plot(
  stats_intercor,
  plot_config,
  show_title,
  save_plots,
  output_path
)
```

## Arguments

- stats_intercor:

  Data frame of inter-variable correlation summaries with columns
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
