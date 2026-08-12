# Tidy/NSE column names used inside dplyr and ggplot2 calls. Declared here so
# R CMD check does not flag them as undefined globals.
#
# Treat this as a shrinking fallback, not a growing manifest: each entry
# disables typo detection for that name across every file in the package.
# Prefer the `.data` pronoun in new code (see `R/warm_filtering_plots.R`),
# which needs no declaration at all.
utils::globalVariables(c(
  "Observed", "Simulated",
  "day", "month", "year", "wyear", "dateo",
  "mon", "rlz", "type",
  "variable", "variable1", "variable2",
  "id1", "id2",
  "idx", "series", "obs",
  "value", "minval", "maxval",
  "x", "y",
  "xmin", "xmax", "ymin", "ymax",
  ".min", ".max"
))
