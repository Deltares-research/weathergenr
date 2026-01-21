# Display fit summary table for all realizations

Formats and prints a compact summary of fit metrics across realizations.
This is primarily intended for interactive review of model performance.

## Usage

``` r
.print_fit_summary_table(fit_summary)
```

## Arguments

- fit_summary:

  Data frame returned by \`.summarize_realization_fit()\`.

## Value

Invisibly returns \`fit_summary\`.
