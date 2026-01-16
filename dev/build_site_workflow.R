# build_site_workflow.R
# Clean, reproducible documentation + vignette + pkgdown build workflow
# Intended usage:
#   source("tools/build_site_workflow.R")
#
# Assumptions:
# - You are running from the package root (where DESCRIPTION lives).
# - Quarto CLI is installed if you render .qmd vignettes.
# - Your Quarto vignettes are configured for R vignette building:
#     DESCRIPTION: VignetteBuilder: quarto
#     vignette YAML: %\VignetteEngine{quarto::html}

.local_stop <- function(...) stop(paste0(...), call. = FALSE)

.assert_pkg_root <- function() {
  if (!file.exists("DESCRIPTION")) .local_stop("Run this from the package root (DESCRIPTION not found).")
}

.safe_unlink <- function(path) {
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
}

.check_quarto <- function() {
  qt <- Sys.which("quarto")
  if (qt == "") {
    message("WARN: Quarto CLI not found on PATH. Quarto vignettes may not build.")
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

.check_vignette_config <- function() {
  desc <- readLines("DESCRIPTION", warn = FALSE)
  has_vb <- any(grepl("^VignetteBuilder:\\s*quarto\\s*$", desc))
  if (!has_vb) {
    message("WARN: DESCRIPTION missing 'VignetteBuilder: quarto'. Quarto vignettes may be ignored by build_vignettes().")
  }
  invisible(has_vb)
}

.check_git_clean <- function() {
  if (Sys.which("git") == "") {
    message("WARN: git not found; skipping repo cleanliness check.")
    return(invisible(NULL))
  }
  st <- system("git status --porcelain", intern = TRUE)
  if (length(st) > 0) {
    message("WARN: Working tree is not clean. Consider committing/stashing before publish.")
  }
  invisible(st)
}

# Optional: user-provided Quarto preview builder (define elsewhere if you want)
# your_quarto_vignette_script <- function() { ... }

build_pkgdown_site <- function(run_quarto_preview = FALSE,
                               run_checks = TRUE,
                               clean_site = TRUE,
                               clean_vignettes = TRUE,
                               open_site = FALSE) {

  .assert_pkg_root()

  if (run_checks) {
    .check_git_clean()
    .check_quarto()
    .check_vignette_config()
  }

  # 0) Restart session (best-effort; if rsession not installed, skip)
  if (requireNamespace("rsession", quietly = TRUE)) {
    rsession::restart()
  } else {
    message("INFO: Package 'rsession' not installed; skipping session restart. Restart R manually for strict cleanliness.")
  }

  # 1) Load tools
  if (!requireNamespace("devtools", quietly = TRUE)) .local_stop("Missing 'devtools'. Install it first.")
  if (!requireNamespace("pkgdown", quietly = TRUE)) .local_stop("Missing 'pkgdown'. Install it first.")

  # 2) Regenerate Rd + NAMESPACE
  devtools::document()

  # 3) Clean old site outputs
  if (clean_site) {
    .safe_unlink("docs")
    .safe_unlink("pkgdown")
  }

  # 4) Optional Quarto vignette preview (does not replace build_vignettes)
  if (isTRUE(run_quarto_preview)) {
    if (exists("your_quarto_vignette_script", mode = "function")) {
      your_quarto_vignette_script()
    } else {
      message("WARN: run_quarto_preview=TRUE but no 'your_quarto_vignette_script()' found; skipping.")
    }
  }

  # 5) Build vignettes in R CMD build-compatible way
  # Optionally clean inst/doc to avoid stale artifacts
  if (clean_vignettes) {
    .safe_unlink(file.path("inst", "doc"))
  }
  devtools::build_vignettes()

  # 6) Run package checks before publishing docs
  # Keep as a gate: documentation should not be published if check fails.
  if (run_checks) {
    devtools::check(
      args = c("--no-manual"),
      check_dir = tempdir()
    )
  }

  # 7) Install (ensures pkgdown documents the installed package)
  devtools::install(upgrade = "never", build_vignettes = TRUE)

  # 8) Build pkgdown site (clean = TRUE forces a full rebuild)
  pkgdown::build_site(clean = TRUE)

  # 9) Optional local preview
  if (isTRUE(open_site)) {
    pkgdown::preview_site()
  }

  invisible(TRUE)
}

# Convenience wrapper with sensible defaults for "ready to push"
publish_docs <- function() {
  build_pkgdown_site(
    run_quarto_preview = FALSE,
    run_checks = TRUE,
    clean_site = TRUE,
    clean_vignettes = TRUE,
    open_site = FALSE
  )
}
