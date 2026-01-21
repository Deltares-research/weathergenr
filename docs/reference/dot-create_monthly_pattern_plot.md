# Create monthly pattern plot for a single variable

Compares observed and simulated monthly distributions (boxplots) for a
given variable, faceted by statistic. Simulated distributions are shown
across realizations.

## Usage

``` r
.create_monthly_pattern_plot(
  stats_mon_aavg_sim,
  stats_mon_aavg_obs,
  variable,
  plot_config,
  show_title,
  save_plots,
  output_dir
)
```

## Arguments

- stats_mon_aavg_sim:

  Data frame of simulated monthly aggregated statistics.

- stats_mon_aavg_obs:

  Data frame of observed monthly aggregated statistics.

- variable:

  Character; variable name to plot (must exist in the inputs).

- plot_config:

  List of plotting configuration options (theme, alpha, colors,
  subtitle, etc.).

- show_title:

  Logical; if `TRUE`, adds title/subtitle.

- save_plots:

  Logical; if `TRUE`, writes plot to `output_dir`.

- output_dir:

  Character; output directory for saved plots.

## Value

ggplot object (returned invisibly by
[`.export_multipanel_plot()`](https://deltares-research.github.io/weathergenr/reference/dot-export_multipanel_plot.md)).
