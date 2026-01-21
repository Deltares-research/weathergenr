# Package-wide internal logger

Unified internal logging helper for the package.

Features: - Single entry point for all logging - Brace interpolation
resolved in caller environment (base R, no glue) - Supports log levels
(info, warn, error) - Silent unless verbose = TRUE - Timestamps in ISO
format

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

  Optional character scalar. Component tag (e.g. "WARM", "KNN").

## Value

Invisibly returns NULL.
