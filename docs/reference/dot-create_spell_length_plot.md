# Create wet/dry spell length diagnostic plot

Observed-vs-simulated comparison of average wet and dry spell lengths,
faceted by spell type/statistic. Uses dummy points to enforce symmetric
axes per facet.

## Usage

``` r
.create_spell_length_plot(
  stats_wetdry,
  plot_config,
  show_title,
  save_plots,
  out_dir
)
```

## Arguments

- stats_wetdry:

  Data frame of wet/dry diagnostics including spell statistics.

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
