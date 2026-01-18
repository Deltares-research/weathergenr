#!/usr/bin/env Rscript
# tools/lint.R
#
# Minimal lintr helper for local use or CI.
#
# Usage:
#   Rscript tools/lint.R            # lint whole package
#   Rscript tools/lint.R R/file.R   # lint specific file(s)
#   Rscript tools/lint.R --changed  # lint changed R files (git)
#
# Exit codes:
#   0 = no lints
#   1 = lints found
#   2 = execution error

args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  if (!requireNamespace("lintr", quietly = TRUE)) {
    message("lintr not installed. Run: install.packages('lintr')")
    quit(status = 2)
  }
})

lint_paths <- function(paths) {
  lints <- unlist(lapply(paths, lintr::lint), recursive = FALSE)
  if (length(lints) == 0) {
    message("No lints found.")
    return(0L)
  }
  print(lints)
  1L
}

lint_changed <- function() {
  if (Sys.which("git") == "") {
    message("git not found on PATH.")
    return(2L)
  }

  files <- system(
    "git diff --name-only --diff-filter=ACMRTUXB HEAD",
    intern = TRUE
  )

  r_files <- files[grepl("\\.R$", files)]
  r_files <- r_files[file.exists(r_files)]

  if (length(r_files) == 0) {
    message("No changed R files to lint.")
    return(0L)
  }

  message("Linting changed files:\n- ", paste(r_files, collapse = "\n- "))
  lint_paths(r_files)
}

status <- tryCatch(
  {
    if (length(args) == 0) {
      lints <- lintr::lint_package()
      if (length(lints) == 0) {
        message("No lints found.")
        0L
      } else {
        print(lints)
        1L
      }
    } else if (identical(args, "--changed")) {
      lint_changed()
    } else {
      lint_paths(args)
    }
  },
  error = function(e) {
    message("Lint failed: ", conditionMessage(e))
    2L
  }
)

quit(status = status)
