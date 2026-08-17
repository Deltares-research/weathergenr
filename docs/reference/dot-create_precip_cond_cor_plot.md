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
  out_dir
)
```

## Arguments

- stats_precip_cor_cond:

  Data frame of conditional correlation summaries with columns
  `variable1`, `variable2`, `id1`, `id2`, `regime`, `Observed`, and
  `Simulated`.

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
