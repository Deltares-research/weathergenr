# Fast column standard deviations

The \*sample\* standard deviation, matching \`stats::sd()\`. The Bessel
factor is not cosmetic here: the targets these values are compared
against and rescaled to are produced by \`stats::sd()\`, so omitting it
left every corrected column at \`target_sd \* sqrt(n / (n - 1))\` – a
one-sided +1.3 \`filter_warm_pool()\` sd tolerance of 3

\`m2 - m1^2\` rather than a centred sum: with annual totals the
cancellation costs about two significant digits out of ~16, which is far
from mattering at the tolerances involved.

## Usage

``` r
.fast_col_sd(x)
```

## Arguments

- x:

  Numeric matrix.

## Value

Numeric vector of column sample standard deviations. \`NA\` per column
when \`x\` has fewer than two rows, which no caller can correct anyway.
