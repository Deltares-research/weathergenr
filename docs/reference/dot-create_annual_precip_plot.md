# Create annual mean precipitation diagnostic plot

Compares annual mean precipitation time series: simulated realizations
are plotted as multiple lines/points and the observed series is overlaid
as a single reference series.

## Usage

``` r
.create_annual_precip_plot(
  stats_annual_aavg_sim,
  stats_annual_aavg_obs,
  plot_config,
  show_title,
  save_plots,
  output_path
)
```

## Arguments

- stats_annual_aavg_sim:

  Data frame of simulated annual aggregated statistics.

- stats_annual_aavg_obs:

  Data frame of observed annual aggregated statistics.

- plot_config:

  List of plotting configuration options (theme, alpha, subtitle, etc.).

- show_title:

  Logical; if `TRUE`, adds title to the plot.

- save_plots:

  Logical; if `TRUE`, writes plot to `output_path`.

- output_path:

  Character; output directory for saved plots.

## Value

ggplot object returned invisibly.
