# Export a faceted (multi-panel) ggplot

Adds title/subtitle (optional) and saves the plot to disk (optional),
using facet layout to infer width/height. Returns the plot invisibly.

## Usage

``` r
.export_multipanel_plot(
  p,
  filename,
  show_title,
  save_plots,
  title = NULL,
  subtitle = NULL,
  output_dir
)
```

## Arguments

- p:

  ggplot object, typically faceted via `facet_wrap()` or `facet_grid()`.

- filename:

  Character; output filename (e.g., `"daily_mean.png"`).

- show_title:

  Logical; if `TRUE`, adds `title`/`subtitle` via `labs()`.

- save_plots:

  Logical; if `TRUE`, writes plot to `output_dir`.

- title:

  Character; plot title (only used when `show_title = TRUE`).

- subtitle:

  Character; plot subtitle (only used when `show_title = TRUE`).

- output_dir:

  Character; output directory for saved plots.

## Value

The ggplot object `p`, returned invisibly.
