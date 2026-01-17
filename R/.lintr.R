

linters: linters_with_defaults(
  line_length_linter = line_length_linter(100),
  object_name_linter = object_name_linter(styles = "snake_case"),
  object_usage_linter = object_usage_linter(),
  cyclocomp_linter = cyclocomp_linter(25),
  undesirable_function_linter = undesirable_function_linter(
    c("attach", "assign", "<<-")
  )
)

exclusions: list(
  "R/zzz.R",              # startup hooks
  "inst/extdata/.*",      # raw data
  "man/.*",               # Rd files
  "vignettes/.*"          # Quarto / Rmd
)
