
# ---- linting ----
lint_pkg <- function() lintr::lint_package()
lint_changed <- function() system("Rscript tools/lint.R --changed")
lint_file <- function(path) {
  if (!requireNamespace("lintr", quietly = TRUE)) {
    stop("lintr is not installed. Run: install.packages('lintr')", call. = FALSE)
  }
  lintr::lint(path)
}
# ---- testing ----
test_pkg <- function(...) testthat::test_local(...)
test_file <- function(path) testthat::test_file(path)

# ---- dev helpers ----
load_pkg <- function() pkgload::load_all(reset = TRUE)
document_pkg <- function() roxygen2::roxygenise()

# ---- checks ----
check_encoding <- function() {
  if (!isTRUE(l10n_info()$`UTF-8`)) {
    warning("Session is not UTF-8.", call. = FALSE)
  }
}

if (interactive()) {
  message("weathergenr project loaded. Helpers available.")
}

message("Loaded project .Rprofile from: ", normalizePath(".Rprofile", mustWork = FALSE))
