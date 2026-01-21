# Create conditional precipitation correlation diagnostic plot

Observed-vs-simulated within-grid correlations between precipitation and
other variables, stratified by precipitation regime (e.g., all/wet/dry).
Faceted by regime (rows) and variable pair (columns).

## Usage

``` r
.create_precip_cond_cor_plot(
  stats_precip_cor_cond,
  plot_config,
  show_title,
  save_plots,
  output_dir
)
```

## Arguments

- stats_precip_cor_cond:

  Data frame of conditional correlation summaries with columns
  `variable1`, `variable2`, `id1`, `id2`, `regime`, `Observed`, and
  `Simulated`.

- plot_config:

  List of plotting configuration options (theme, alpha, subtitle, etc.).

- show_title:

  Logical; if `TRUE`, adds title/subtitle.

- save_plots:

  Logical; if `TRUE`, writes plot to `output_dir`.

- output_dir:

  Character; output directory for saved plots.

## Value

ggplot object (returned invisibly by
[`.export_multipanel_plot()`](https://deltares-research.github.io/weathergenr/reference/dot-export_multipanel_plot.md)).
