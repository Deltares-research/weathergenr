# Package-wide internal logger

Unified internal logging helper for the package.

Features: - Single entry point for all logging - Brace interpolation
resolved in caller environment (base R, no glue) - Supports log levels
(info, warn, error) - Silent unless verbose = TRUE - Lines are formatted
as `HH:MM:SS - tag - message`, matching the downstream \`blueearth_cst\`
console syntax

## Usage

``` r
.log(msg, level = c("info", "warn", "error"), verbose = TRUE, tag = NULL)
```

## Arguments

- msg:

  Character scalar. Log message template with `{variable}` syntax.

- level:

  Character scalar. One of "info", "warn", "error".

- verbose:

  Logical. If FALSE, suppress output.

- tag:

  Optional character scalar. Component tag (e.g. "WARM", "KNN"), emitted
  lower-cased.

## Value

Invisibly returns NULL.
