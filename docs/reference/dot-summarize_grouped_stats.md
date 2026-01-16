# Summarize grouped statistics

Summarize grouped statistics

## Usage

``` r
.summarize_grouped_stats(df, variables, group_vars, stat_fns)
```

## Arguments

- df:

  Data frame to summarize.

- variables:

  Character vector of variables to summarize.

- group_vars:

  Character vector of grouping columns.

- stat_fns:

  Named list of summary functions.

## Value

Long-format data frame with \`variable\`, \`stat\`, and \`value\`.
