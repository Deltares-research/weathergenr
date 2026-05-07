# Tidy/NSE column names used inside dplyr and ggplot2 calls. Declared here so
# R CMD check does not flag them as undefined globals.
utils::globalVariables(c(
  "Observed", "Simulated",
  "day", "month", "year", "wyear",
  "mon", "rlz", "type",
  "variable", "variable1", "variable2",
  "id1", "id2",
  "x"
))
