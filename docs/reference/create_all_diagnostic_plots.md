# Create all diagnostic plots

Runs the full set of diagnostic plotting routines and returns a named
list of ggplot objects. Optionally saves plots to disk (delegated to the
individual plot exporters).

## Usage

``` r
create_all_diagnostic_plots(
  plot_data,
  plot_config,
  variables,
  show_title,
  save_plots,
  output_path
)
```

## Arguments

- plot_data:

  List of precomputed diagnostic datasets produced by the evaluation
  pipeline.

- plot_config:

  List of plotting configuration options (theme, alpha, colors,
  subtitle).

- variables:

  Character vector of variable names to loop over for monthly pattern
  plots.

- show_title:

  Logical; if `TRUE`, titles/subtitles are added to plots where
  supported.

- save_plots:

  Logical; if `TRUE`, plots are written to `output_path`.

- output_path:

  Character; output directory for saved plots.

## Value

Named list of ggplot objects for all diagnostics created.

## Details

This helper expects the precomputed plot data returned by the evaluation
pipeline. It does not validate plot input structure beyond basic use in
downstream plotting.

## Examples

``` r
if (FALSE) { # \dontrun{
  plot_data <- list()
  plot_config <- list(
    subtitle = "Example",
    alpha = 0.4,
    colors = c(Observed = "blue3", Simulated = "gray40"),
    theme = ggplot2::theme_bw()
  )
  plots <- create_all_diagnostic_plots(
    plot_data = plot_data,
    plot_config = plot_config,
    variables = c("precip", "temp"),
    show_title = FALSE,
    save_plots = FALSE,
    output_path = NULL
  )
} # }
```
