# dev/build_site_workflow.R
# Clean, reproducible documentation + vignette + pkgdown build workflow
#
# Intended usage (from package root):
#   source("dev/build_site_workflow.R")
#   publish_docs()
#
# Notes:
# - Designed for Quarto (.qmd) vignettes with DESCRIPTION: VignetteBuilder: quarto
# - Works on Windows by pre-setting TMPDIR/TMP/TEMP before any tooling loads.

# ------------------------------------------------------------------
# FIX Windows Quarto / system2 TMPDIR issue
# MUST run before pkgdown/devtools/quarto/rmarkdown loads and spawns processes.
# ------------------------------------------------------------------
.set_safe_tmpdir <- function() {
  td <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
  Sys.setenv(
    TMPDIR = td,
    TMP    = td,
    TEMP   = td
  )
  invisible(td)
}
.set_safe_tmpdir()


.check_quarto_cli <- function(strict = FALSE) {
  qt <- Sys.which("quarto")
  if (qt == "") {
    msg <- "Quarto CLI not found on PATH."
    if (isTRUE(strict)) .local_stop(msg) else .local_message(msg, prefix = "WARN")
    return(invisible(FALSE))
  }


  # Version check - use Sys.setenv, NOT system2 env argument (Windows bug)
  v <- tryCatch({
    # Suppress the spurious Windows warning about env argument
    suppressWarnings(
      system2("quarto", "--version", stdout = TRUE, stderr = TRUE)
    )
  }, error = function(e) NULL)

  if (!is.null(v) && length(v) > 0L && !grepl("ERROR:", v[1])) {
    .local_message("Quarto CLI version: ", v[1])
  }

  invisible(TRUE)
}

# ------------------------------------------------------------------
# Local helpers
# ------------------------------------------------------------------
.local_stop <- function(...) stop(paste0(...), call. = FALSE)

.assert_pkg_root <- function() {
  if (!file.exists("DESCRIPTION")) {
    .local_stop("Run this from the package root (DESCRIPTION not found).")
  }
}

.safe_unlink <- function(path) {
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  if (file.exists(path)) unlink(path, force = TRUE)
  invisible(TRUE)
}

.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    .local_stop("Missing required package: '", pkg, "'. Install it first.")
  }
  invisible(TRUE)
}

.check_quarto <- function(strict = FALSE) {
  qt <- Sys.which("quarto")
  if (qt == "") {
    msg <- "Quarto CLI not found on PATH. Quarto vignettes will fail to build."
    if (isTRUE(strict)) .local_stop(msg) else message("WARN: ", msg)
    return(invisible(FALSE))
  }
  # best-effort: show version (do not pass env=... here)
  v <- tryCatch(system2("quarto", "-V", stdout = TRUE, stderr = TRUE), error = function(e) e)
  if (inherits(v, "error")) {
    msg <- paste0("Quarto CLI exists but version check failed: ", conditionMessage(v))
    if (isTRUE(strict)) .local_stop(msg) else message("WARN: ", msg)
  }
  invisible(TRUE)
}

.check_vignette_config <- function(strict = FALSE) {
  desc <- readLines("DESCRIPTION", warn = FALSE)

  has_vb <- any(grepl("^VignetteBuilder:\\s*quarto\\s*$", desc))
  if (!has_vb) {
    msg <- "DESCRIPTION missing 'VignetteBuilder: quarto'. Quarto vignettes may not build via devtools::build_vignettes()."
    if (isTRUE(strict)) .local_stop(msg) else message("WARN: ", msg)
  }

  invisible(has_vb)
}

.check_pkgdown_config <- function() {
  # Informational only: many repos intentionally use docs/ for GitHub Pages.
  if (!file.exists("_pkgdown.yml")) {
    message("WARN: _pkgdown.yml not found. pkgdown will use defaults.")
  }
  invisible(TRUE)
}

.check_git_clean <- function() {
  if (Sys.which("git") == "") {
    message("WARN: git not found; skipping repo cleanliness check.")
    return(invisible(NULL))
  }
  st <- system("git status --porcelain", intern = TRUE)
  if (length(st) > 0) {
    message("WARN: Working tree is not clean. Consider committing/stashing before publishing.")
  }
  invisible(st)
}

# Optional: define this in a separate file if you want fast local preview rendering
# without affecting R CMD build-compatible vignette build.
# your_quarto_vignette_script <- function() { ... }

# ------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------
build_pkgdown_site <- function(run_quarto_preview = FALSE,
                               run_checks = TRUE,
                               strict_quarto = TRUE,
                               clean_site = TRUE,
                               clean_vignettes = TRUE,
                               clean_check_artifacts = TRUE,
                               open_site = FALSE) {

  .assert_pkg_root()

  # Preflight (do before loading devtools/pkgdown so failure is earlier/clearer)
  .check_pkgdown_config()
  if (run_checks) .check_git_clean()
  .check_quarto(strict = strict_quarto)
  .check_vignette_config(strict = FALSE)

  # Tooling
  .require_pkg("devtools")
  .require_pkg("pkgdown")

  # 1) Regenerate Rd + NAMESPACE (roxygen)
  devtools::document()

  # 2) Clean old outputs that commonly create stale sites / stale vignettes
  if (isTRUE(clean_site)) {
    # pkgdown output (choose one or both depending on your repo history)
    .safe_unlink("docs")
    .safe_unlink("pkgdown")

    # pkgdown cache dirs (safe to delete; pkgdown will recreate)
    .safe_unlink(file.path("inst", "pkgdown"))
    .safe_unlink(file.path("pkgdown", "assets"))
  }

  if (isTRUE(clean_vignettes)) {
    # Built vignettes (HTML/PDF) end up here after build/install
    .safe_unlink(file.path("inst", "doc"))
  }

  if (isTRUE(clean_check_artifacts)) {
    # Common leftovers
    .safe_unlink(".Rcheck")
    .safe_unlink("..Rcheck")
    .safe_unlink(file.path("vignettes", ".quarto"))
    .safe_unlink(file.path("vignettes", "_cache"))
    .safe_unlink(file.path("vignettes", "_freeze"))
  }

  # 3) Optional local Quarto preview step (non-authoritative)
  if (isTRUE(run_quarto_preview)) {
    if (exists("your_quarto_vignette_script", mode = "function")) {
      your_quarto_vignette_script()
    } else {
      message("WARN: run_quarto_preview=TRUE but no 'your_quarto_vignette_script()' found; skipping preview.")
    }
  }

  # 4) Build vignettes in R CMD build-compatible way
  # This is the step that populates inst/doc.
  devtools::build_vignettes()

  # 5) Gate: run checks before publishing docs
  # If check fails, stop here. Do not publish broken docs.
  if (isTRUE(run_checks)) {
    devtools::check(
      args = c("--no-manual"),
      check_dir = tempdir()
    )
  }

  # 6) Install (pkgdown generally documents the installed version best)
  devtools::install(upgrade = "never", build_vignettes = TRUE)

  # 7) Build pkgdown site
  # NOTE: some pkgdown versions do NOT support build_site(clean = TRUE).
  # We already cleaned docs/inst/doc above, so this is a clean rebuild in practice.
  pkgdown::build_site()

  # 8) Optional local preview
  if (isTRUE(open_site)) {
    pkgdown::preview_site()
  }

  invisible(TRUE)
}

# ------------------------------------------------------------------
# Convenience wrapper: "ready to push"
# ------------------------------------------------------------------
publish_docs <- function() {
  build_pkgdown_site(
    run_quarto_preview = FALSE,
    run_checks = TRUE,
    strict_quarto = TRUE,
    clean_site = TRUE,
    clean_vignettes = TRUE,
    clean_check_artifacts = TRUE,
    open_site = FALSE
  )
}
