# Package-wide internal logger

Unified internal logging helper for the package.

Features: - Single entry point for all logging - Glue interpolation
resolved in caller environment - Safe with logger::formatter_glue (no
raw \`\` passed downstream) - Supports log levels (info, warn, error) -
Silent unless verbose = TRUE - Falls back to message()/warning()/stop()
if logger is unavailable

## Usage

``` r
.log(msg, ..., level = c("info", "warn", "error"), verbose = TRUE, tag = NULL)
```

## Arguments

- msg:

  Character scalar. Log message template.

- ...:

  Values interpolated into \`msg\` via glue.

- level:

  Character scalar. One of "info", "warn", "error".

- verbose:

  Logical. If FALSE, suppress output.

- tag:

  Optional character scalar. Component tag (e.g. "WARM", "KNN").

## Value

Invisibly returns NULL.
