# tools/build_site_tools.R
# Clean, reproducible documentation + vignette + pkgdown build workflow
#
# Intended usage (from package root):
#   source("tools/build_site_tools.R")
#   publish_docs()
#
# For development iteration (faster, no checks):
#   quick_site()
#
# Notes:
# - Designed for Quarto (.qmd) vignettes with DESCRIPTION: VignetteBuilder: quarto
# - Works on Windows by pre-setting TMPDIR/TMP/TEMP before any tooling loads.

# ==============================================================================
# WINDOWS TMPDIR FIX
# MUST run before pkgdown/devtools/quarto/rmarkdown loads and spawns processes.
# ==============================================================================

.install_quarto_cli_shim <- function() {
  if (.Platform$OS.type != "windows") return(invisible(FALSE))

  real_quarto <- Sys.which("quarto")
  if (real_quarto == "") {
    warning("Quarto CLI not found on PATH; cannot install shim.", call. = FALSE)
    return(invisible(FALSE))
  }

  shim_dir <- file.path(tempdir(), "quarto_shim")
  dir.create(shim_dir, recursive = TRUE, showWarnings = FALSE)

  shim_path <- file.path(shim_dir, "quarto.cmd")

  # A minimal shim: remove any TMP*, TEMP* tokens that were incorrectly passed as args.
  shim <- c(
    "@echo off",
    "setlocal EnableExtensions EnableDelayedExpansion",
    "set \"REAL_QUARTO=" %+% real_quarto %+% "\"",
    "set \"ARGS=\"",
    ":loop",
    "if \"%~1\"==\"\" goto run",
    "set \"A=%~1\"",
    "set \"SKIP=0\"",
    "echo(!A!| findstr /I /B /C:\"TMPDIR=\" /C:\"TMP=\" /C:\"TEMP=\" >nul && set \"SKIP=1\"",
    "if \"!SKIP!\"==\"0\" set \"ARGS=!ARGS! \"\"!A!\"\"\"",
    "shift",
    "goto loop",
    ":run",
    "call \"%REAL_QUARTO%\" %ARGS%",
    "exit /b %ERRORLEVEL%"
  )

  # helper for concatenation
  `%+%` <- function(a, b) paste0(a, b)

  writeLines(shim, shim_path, useBytes = TRUE)

  # Prepend shim_dir to PATH so system2("quarto", ...) hits the shim.
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(shim_dir, old_path, sep = .Platform$path.sep))

  invisible(TRUE)
}

# Call this as early as possible (before any quarto/rmarkdown/pkgdown/devtools loads)
.install_quarto_cli_shim()


# .set_safe_tmpdir <- function() {
#   td <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
#   Sys.setenv(TMPDIR = td, TMP = td, TEMP = td)
#   invisible(td)
# }
# .set_safe_tmpdir()

# ==============================================================================
# LOCAL HELPERS
# ==============================================================================

.local_stop <- function(...) {
  stop(paste0("[ERROR] ", ...), call. = FALSE)
}

.local_message <- function(..., prefix = "INFO") {
  message(paste0("[", prefix, "] ", ...))
}

.local_success <- function(...) {
  .local_message(..., prefix = "OK")
}

.local_warn <- function(...) {
  .local_message(..., prefix = "WARN")
}

.local_step <- function(step_num, total, description) {
  .local_message(
    sprintf("Step %d/%d: %s", step_num, total, description),
    prefix = "----"
  )
}

.assert_pkg_root <- function() {
  if (!file.exists("DESCRIPTION")) {
    .local_stop("Run this from the package root (DESCRIPTION not found).")
  }
  invisible(TRUE)
}

.safe_unlink <- function(path) {
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  if (file.exists(path)) unlink(path, force = TRUE)
  invisible(TRUE)
}

.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    .local_stop(
      "Missing required package: '", pkg,
      "'. Install with install.packages('", pkg, "')"
    )
  }
  invisible(TRUE)
}

.get_pkg_name <- function() {
  desc <- read.dcf("DESCRIPTION", fields = "Package")
  as.character(desc[1, "Package"])
}

.get_pkg_version <- function() {
  desc <- read.dcf("DESCRIPTION", fields = "Version")
  as.character(desc[1, "Version"])
}

# ==============================================================================
# PREFLIGHT CHECKS
# ==============================================================================

.check_quarto_cli <- function(strict = FALSE) {
  qt <- Sys.which("quarto")
  if (qt == "") {
    msg <- "Quarto CLI not found on PATH. Quarto vignettes will fail to build."
    if (isTRUE(strict)) .local_stop(msg) else .local_warn(msg)
    return(invisible(FALSE))
  }

  v <- tryCatch({
    suppressWarnings(system2("quarto", "--version", stdout = TRUE, stderr = TRUE))
  }, error = function(e) NULL)

  if (!is.null(v) && length(v) > 0L && !grepl("ERROR:", v[1])) {
    .local_message("Quarto CLI version: ", v[1])
  }

  invisible(TRUE)
}

.check_quarto_r_pkg <- function(strict = FALSE) {
  if (!requireNamespace("quarto", quietly = TRUE)) {
    msg <- "R package 'quarto' is not installed. Quarto vignettes may fail to build."
    if (isTRUE(strict)) .local_stop(msg) else .local_warn(msg)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

.check_vignette_config <- function(strict = FALSE) {
  desc <- readLines("DESCRIPTION", warn = FALSE)

  has_vb <- any(grepl("^VignetteBuilder:", desc))
  if (!has_vb) {
    msg <- "DESCRIPTION missing 'VignetteBuilder'. Add 'VignetteBuilder: quarto' for .qmd vignettes."
    if (isTRUE(strict)) .local_stop(msg) else .local_warn(msg)
    return(invisible(FALSE))
  }

  has_quarto <- any(grepl("^VignetteBuilder:.*quarto", desc))
  if (!has_quarto) {
    .local_message(
      "VignetteBuilder is not 'quarto'. Using: ",
      gsub("^VignetteBuilder:\\s*", "", desc[grepl("^VignetteBuilder:", desc)])
    )
  }

  invisible(has_vb)
}

.check_pkgdown_config <- function(strict = FALSE) {
  if (!file.exists("_pkgdown.yml") && !file.exists("pkgdown/_pkgdown.yml")) {
    msg <- "_pkgdown.yml not found. pkgdown will use defaults (may miss some functions in reference)."
    if (isTRUE(strict)) .local_stop(msg) else .local_warn(msg)
    return(invisible(FALSE))
  }
  .local_success("pkgdown configuration found")
  invisible(TRUE)
}

.check_git_clean <- function(strict = FALSE) {
  if (Sys.which("git") == "") {
    .local_warn("git not found; skipping repo cleanliness check.")
    return(invisible(NULL))
  }

  st <- tryCatch(
    system("git status --porcelain", intern = TRUE, ignore.stderr = TRUE),
    error = function(e) character(0)
  )

  if (length(st) > 0) {
    msg <- sprintf(
      "Working tree has %d uncommitted change(s). Consider committing before publishing.",
      length(st)
    )
    if (isTRUE(strict)) .local_stop(msg) else .local_warn(msg)
  } else {
    .local_success("Git working tree is clean")
  }

  invisible(st)
}

.check_namespace_current <- function() {
  if (!file.exists("NAMESPACE")) {
    .local_warn("NAMESPACE file not found. Run devtools::document() first.")
    return(invisible(FALSE))
  }

  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  if (length(r_files) == 0) return(invisible(TRUE))

  ns_time <- file.mtime("NAMESPACE")
  r_times <- file.mtime(r_files)

  if (any(r_times > ns_time)) {
    .local_warn("Some R files are newer than NAMESPACE. Will regenerate with document().")
  }

  invisible(TRUE)
}

.check_license <- function() {
  has_license <- file.exists("LICENSE") || file.exists("LICENSE.md") || file.exists("LICENCE")
  if (!has_license) {
    .local_warn("No LICENSE file found. Consider adding one.")
  }
  invisible(has_license)
}

.check_readme <- function() {
  has_readme <- file.exists("README.md") || file.exists("README.Rmd")
  if (!has_readme) {
    .local_warn("No README found. Consider adding README.md or README.Rmd.")
  }
  invisible(has_readme)
}

# ==============================================================================
# BUILD STEPS
# ==============================================================================

#' Clean All Build Artifacts
#'
#' Comprehensive removal of intermediate and output files from previous
#' documentation, vignette, and package builds.
#'
#' @param clean_site Logical. Remove pkgdown site output?
#' @param clean_vignettes Logical. Remove vignette build artifacts?
#' @param clean_check Logical. Remove R CMD check artifacts?
#' @param clean_build Logical. Remove package build artifacts (tarballs/binaries)?
#' @param clean_temp Logical. Remove temp/cache files?
#' @param verbose Logical. Print progress messages?
.clean_build_artifacts <- function(clean_site = TRUE,
                                   clean_vignettes = TRUE,
                                   clean_check = TRUE,
                                   clean_build = TRUE,
                                   clean_temp = TRUE,
                                   verbose = TRUE) {

  removed_count <- 0L

  # PKGDOWN SITE ARTIFACTS
  if (isTRUE(clean_site)) {
    site_artifacts <- c(
      "docs",
      "pkgdown",
      file.path("inst", "pkgdown"),
      "reference",
      "_site",
      "_book"
    )

    for (path in site_artifacts) {
      if (dir.exists(path) || file.exists(path)) {
        .safe_unlink(path)
        removed_count <- removed_count + 1L
      }
    }

    if (verbose) .local_message("Cleaned pkgdown/site artifacts")
  }

  # VIGNETTE BUILD ARTIFACTS
  if (isTRUE(clean_vignettes)) {
    vignette_dirs <- c(
      file.path("inst", "doc"),
      file.path("vignettes", ".quarto"),
      file.path("vignettes", "_cache"),
      file.path("vignettes", "_freeze"),
      file.path("vignettes", "_site"),
      file.path("vignettes", "rsconnect"),
      "doc"
    )

    for (path in vignette_dirs) {
      if (dir.exists(path)) {
        .safe_unlink(path)
        removed_count <- removed_count + 1L
      }
    }

    if (dir.exists("vignettes")) {
      vignette_html <- list.files("vignettes", pattern = "\\.html$", full.names = TRUE)
      lapply(vignette_html, .safe_unlink)
      removed_count <- removed_count + length(vignette_html)

      files_dirs <- list.files("vignettes", pattern = "_files$", full.names = TRUE)
      files_dirs <- files_dirs[dir.exists(files_dirs)]
      lapply(files_dirs, .safe_unlink)
      removed_count <- removed_count + length(files_dirs)

      cache_dirs <- list.files("vignettes", pattern = "_cache$", full.names = TRUE)
      cache_dirs <- cache_dirs[dir.exists(cache_dirs)]
      lapply(cache_dirs, .safe_unlink)
      removed_count <- removed_count + length(cache_dirs)

      quarto_files <- list.files("vignettes", pattern = "\\.(tex|log|aux)$", full.names = TRUE)
      lapply(quarto_files, .safe_unlink)
      removed_count <- removed_count + length(quarto_files)
    }

    if (verbose) .local_message("Cleaned vignette build artifacts")
  }

  # R CMD CHECK ARTIFACTS
  if (isTRUE(clean_check)) {
    check_dirs <- list.files(".", pattern = "\\.Rcheck$", full.names = TRUE, all.files = TRUE)
    lapply(check_dirs, .safe_unlink)
    removed_count <- removed_count + length(check_dirs)

    if (verbose) .local_message("Cleaned R CMD check artifacts")
  }

  # PACKAGE BUILD ARTIFACTS
  if (isTRUE(clean_build)) {
    tarballs <- list.files(".", pattern = "\\.tar\\.gz$", full.names = TRUE)
    pkg_name <- tryCatch(.get_pkg_name(), error = function(e) NULL)
    if (!is.null(pkg_name)) {
      pkg_tarballs <- tarballs[grepl(paste0("^\\./", pkg_name, "_"), tarballs)]
      lapply(pkg_tarballs, .safe_unlink)
      removed_count <- removed_count + length(pkg_tarballs)
    }

    binaries <- list.files(".", pattern = "\\.(zip|tgz)$", full.names = TRUE)
    if (!is.null(pkg_name)) {
      pkg_binaries <- binaries[grepl(paste0("^\\./", pkg_name, "_"), binaries)]
      lapply(pkg_binaries, .safe_unlink)
      removed_count <- removed_count + length(pkg_binaries)
    }

    .safe_unlink("build")
    .safe_unlink("Meta")

    if (verbose) .local_message("Cleaned package build artifacts")
  }

  # TEMPORARY AND CACHE FILES
  if (isTRUE(clean_temp)) {
    temp_files <- c(
      ".Rhistory",
      ".RData",
      ".DS_Store",
      "Thumbs.db"
    )

    for (f in temp_files) {
      if (file.exists(f)) {
        .safe_unlink(f)
        removed_count <- removed_count + 1L
      }
    }

    root_quarto <- c(".quarto", "_cache", "_freeze")
    for (d in root_quarto) {
      if (dir.exists(d)) {
        .safe_unlink(d)
        removed_count <- removed_count + 1L
      }
    }

    covr_dirs <- list.files(".", pattern = "^covr", full.names = TRUE)
    lapply(covr_dirs, .safe_unlink)
    .safe_unlink("codecov")
    .safe_unlink("coverage.html")

    for (f in Sys.glob(file.path("src", "*.{o,so,dll}"))) {
      .safe_unlink(f)
    }

    if (verbose) .local_message("Cleaned temporary and cache files")
  }

  if (verbose) {
    .local_success(sprintf("Artifact cleanup complete (%d items processed)", removed_count))
  }

  invisible(removed_count)
}

#' Deep Clean All Artifacts
#'
#' Aggressive cleaning of ALL build artifacts. Use with caution.
#' This removes everything that can be regenerated.
#'
#' @export
deep_clean <- function() {
  .assert_pkg_root()

  .local_warn("Deep clean will remove ALL generated files.")
  .local_warn("This includes: docs/, inst/doc/, vignette HTML, .Rcheck/, package builds")

  if (interactive()) {
    response <- readline("Continue? (y/N): ")
    if (!tolower(response) %in% c("y", "yes")) {
      .local_message("Aborted.")
      return(invisible(FALSE))
    }
  }

  .clean_build_artifacts(
    clean_site = TRUE,
    clean_vignettes = TRUE,
    clean_check = TRUE,
    clean_build = TRUE,
    clean_temp = TRUE,
    verbose = TRUE
  )

  .local_success("Deep clean complete. Run publish_docs() to rebuild everything.")
  invisible(TRUE)
}

.render_readme <- function() {
  if (file.exists("README.Rmd")) {
    .local_message("Rendering README.Rmd -> README.md")
    .require_pkg("rmarkdown")

    tryCatch({
      rmarkdown::render("README.Rmd", output_format = "github_document", quiet = TRUE)
      .local_success("README rendered successfully")
    }, error = function(e) {
      .local_warn("README rendering failed: ", conditionMessage(e))
    })
  }
  invisible(TRUE)
}

.run_tests <- function(stop_on_failure = TRUE) {
  .require_pkg("testthat")
  .local_message("Running package tests...")

  ok <- tryCatch({
    testthat::test_local(reporter = "summary")
    TRUE
  }, error = function(e) {
    .local_warn("Tests failed: ", conditionMessage(e))
    FALSE
  })

  if (!ok && stop_on_failure) .local_stop("Package tests failed. Fix before publishing.")
  invisible(ok)
}

.run_spell_check <- function(stop_on_failure = FALSE) {
  if (!requireNamespace("spelling", quietly = TRUE)) {
    .local_message("spelling package not installed; skipping spell check")
    return(invisible(NULL))
  }

  .local_message("Running spell check...")

  errors <- tryCatch({
    spelling::spell_check_package()
  }, error = function(e) {
    .local_warn("Spell check failed: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(errors) && nrow(errors) > 0) {
    .local_warn(
      sprintf("Found %d potential spelling errors. Review with spelling::spell_check_package()", nrow(errors))
    )
    if (stop_on_failure) {
      print(errors)
      .local_stop("Spell check found errors")
    }
  } else {
    .local_success("No spelling errors found")
  }

  invisible(errors)
}

.copy_built_vignettes_to_inst_doc <- function() {
  if (!dir.exists("doc")) {
    .local_warn("No 'doc/' directory found after vignette build; nothing to copy.")
    return(invisible(FALSE))
  }

  if (!dir.exists(file.path("inst", "doc"))) {
    dir.create(file.path("inst", "doc"), recursive = TRUE)
  }

  html_files <- list.files("doc", pattern = "\\.html$", full.names = TRUE)
  if (length(html_files) > 0) {
    file.copy(html_files, file.path("inst", "doc"), overwrite = TRUE)
  }

  invisible(TRUE)
}

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

#' Build pkgdown Site with Full Validation
#'
#' Complete workflow for building package documentation and pkgdown site.
#'
#' @param run_tests Logical. Run package tests before building? Default TRUE.
#' @param run_checks Logical. Run R CMD check before publishing? Default TRUE.
#' @param run_spell_check Logical. Run spelling check? Default FALSE.
#' @param strict_quarto Logical. Fail if Quarto CLI / R package not found? Default TRUE.
#' @param clean_site Logical. Remove old pkgdown output? Default TRUE.
#' @param clean_vignettes Logical. Remove old vignette builds? Default TRUE.
#' @param clean_check_artifacts Logical. Remove .Rcheck directories? Default TRUE.
#' @param clean_build_artifacts Logical. Remove package tarballs/binaries? Default TRUE.
#' @param clean_temp_files Logical. Remove temp files (.RData, etc.)? Default TRUE.
#' @param render_readme Logical. Render README.Rmd if present? Default TRUE.
#' @param check_args Character vector. Additional arguments for R CMD check.
#' @param check_build_args Character vector. Build arguments for `R CMD build/check`.
#'   Default includes `--no-build-vignettes` because vignettes are built explicitly
#'   in this workflow before the check step.
#' @param open_site Logical. Open site in browser when done? Default FALSE.
#' @param verbose Logical. Print detailed progress? Default TRUE.
#'
#' @details
#' ## Windows/Quarto stability approach
#'
#' On some Windows setups, rebuilding Quarto vignettes inside `R CMD check` can be
#' brittle due to process spawning / cached paths. This workflow:
#' 1. Builds vignettes explicitly via `tools::buildVignettes()` before the check step.
#' 2. Runs `devtools::check()` with `build_args = "--no-build-vignettes"` so vignettes
#'    are not rebuilt during check.
#'
#' @return Invisible TRUE on success.
#' @export
build_pkgdown_site <- function(run_tests = TRUE,
                               run_checks = TRUE,
                               run_spell_check = FALSE,
                               strict_quarto = TRUE,
                               clean_site = TRUE,
                               clean_vignettes = TRUE,
                               clean_check_artifacts = TRUE,
                               clean_build_artifacts = TRUE,
                               clean_temp_files = TRUE,
                               render_readme = TRUE,
                               check_args = c("--no-manual", "--as-cran"),
                               check_build_args = c("--no-build-vignettes"),
                               open_site = FALSE,
                               verbose = TRUE) {

  start_time <- Sys.time()
  .assert_pkg_root()

  pkg_name <- .get_pkg_name()
  pkg_version <- .get_pkg_version()
  .local_message("Building documentation for ", pkg_name, " v", pkg_version)

  total_steps <- 6 + as.integer(run_tests) + as.integer(run_checks) +
    as.integer(run_spell_check) + as.integer(render_readme)
  current_step <- 0

  # PREFLIGHT CHECKS
  current_step <- current_step + 1
  .local_step(current_step, total_steps, "Preflight checks")

  .check_pkgdown_config(strict = FALSE)
  .check_git_clean(strict = FALSE)
  .check_quarto_cli(strict = strict_quarto)
  .check_quarto_r_pkg(strict = strict_quarto)
  .check_vignette_config(strict = FALSE)
  .check_namespace_current()
  .check_license()
  .check_readme()

  # LOAD TOOLING
  .require_pkg("devtools")
  .require_pkg("pkgdown")
  .require_pkg("roxygen2")

  # CLEAN OLD ARTIFACTS
  current_step <- current_step + 1
  .local_step(current_step, total_steps, "Cleaning build artifacts")

  .clean_build_artifacts(
    clean_site = clean_site,
    clean_vignettes = clean_vignettes,
    clean_check = clean_check_artifacts,
    clean_build = clean_build_artifacts,
    clean_temp = clean_temp_files,
    verbose = verbose
  )

  # RENDER README
  if (isTRUE(render_readme)) {
    current_step <- current_step + 1
    .local_step(current_step, total_steps, "Rendering README")
    .render_readme()
  }

  # GENERATE DOCUMENTATION
  current_step <- current_step + 1
  .local_step(current_step, total_steps, "Generating Rd documentation (roxygen2)")

  devtools::document()
  .local_success("Documentation generated")

  # RUN TESTS
  if (isTRUE(run_tests)) {
    current_step <- current_step + 1
    .local_step(current_step, total_steps, "Running package tests")
    .run_tests(stop_on_failure = TRUE)
  }

  # SPELL CHECK
  if (isTRUE(run_spell_check)) {
    current_step <- current_step + 1
    .local_step(current_step, total_steps, "Running spell check")
    .run_spell_check(stop_on_failure = FALSE)
  }

  # BUILD VIGNETTES
  current_step <- current_step + 1
  .local_step(current_step, total_steps, "Building vignettes")

  tryCatch({
    # Use tools::buildVignettes instead of devtools::build_vignettes for stability on Windows.
    tools::buildVignettes(dir = ".")
    .copy_built_vignettes_to_inst_doc()
    .local_success("Vignettes built successfully")
  }, error = function(e) {
    .local_stop("Vignette build failed: ", conditionMessage(e))
  })

  # R CMD CHECK
  if (isTRUE(run_checks)) {
    current_step <- current_step + 1
    .local_step(current_step, total_steps, "Running R CMD check")

    if (any(check_build_args == "--no-build-vignettes")) {
      .local_message("Using --no-build-vignettes (vignettes already built above)")
    }

    check_result <- tryCatch({
      devtools::check(
        args = check_args,
        build_args = check_build_args,
        check_dir = tempdir()
      )
    }, error = function(e) {
      .local_stop("R CMD check failed: ", conditionMessage(e))
    })

    if (!is.null(check_result)) {
      n_errors <- length(check_result$errors)
      n_warnings <- length(check_result$warnings)
      n_notes <- length(check_result$notes)

      .local_message(sprintf(
        "Check results: %d errors, %d warnings, %d notes",
        n_errors, n_warnings, n_notes
      ))

      if (n_errors > 0) {
        .local_stop("R CMD check produced errors. Fix before publishing.")
      }
    }
  }

  # INSTALL PACKAGE
  current_step <- current_step + 1
  .local_step(current_step, total_steps, "Installing package")

  devtools::install(upgrade = "never", build_vignettes = FALSE, quiet = TRUE)
  .local_success("Package installed")

  # BUILD PKGDOWN SITE
  current_step <- current_step + 1
  .local_step(current_step, total_steps, "Building pkgdown site")

  tryCatch({
    pkgdown::build_site()
    .local_success("pkgdown site built successfully")
  }, error = function(e) {
    .local_stop("pkgdown site build failed: ", conditionMessage(e))
  })

  # SUMMARY
  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
  .local_message("", prefix = "====")
  .local_success(sprintf(
    "Build complete for %s v%s in %.1f minutes",
    pkg_name, pkg_version, elapsed
  ))

  if (file.exists("docs/index.html")) {
    .local_message("Site output: ", normalizePath("docs", mustWork = FALSE))
  }

  # OPTIONAL: OPEN SITE
  if (isTRUE(open_site)) {
    pkgdown::preview_site()
  }

  invisible(TRUE)
}

# ==============================================================================
# CONVENIENCE WRAPPERS
# ==============================================================================

#' Publish Documentation (Full Workflow)
#'
#' Production-ready build with all checks enabled. Use before pushing to GitHub
#' or releasing.
#'
#' @details
#' Vignettes are built explicitly before `R CMD check`. The check step uses
#' `build_args = "--no-build-vignettes"` to avoid rebuilding vignettes during check.
#'
#' @export
publish_docs <- function() {
  build_pkgdown_site(
    run_tests = TRUE,
    run_checks = TRUE,
    run_spell_check = FALSE,
    strict_quarto = TRUE,
    clean_site = TRUE,
    clean_vignettes = TRUE,
    clean_check_artifacts = TRUE,
    clean_build_artifacts = TRUE,
    clean_temp_files = TRUE,
    render_readme = TRUE,
    check_args = c("--no-manual", "--as-cran"),
    check_build_args = c("--no-build-vignettes"),
    open_site = FALSE
  )
}

#' Quick Site Build (Development)
#'
#' Fast iteration build without tests or R CMD check. Use during active development
#' to preview documentation changes.
#'
#' @param open Logical. Open site in browser? Default TRUE.
#' @export
quick_site <- function(open = TRUE) {
  build_pkgdown_site(
    run_tests = FALSE,
    run_checks = FALSE,
    run_spell_check = FALSE,
    strict_quarto = FALSE,
    clean_site = FALSE,
    clean_vignettes = FALSE,
    clean_check_artifacts = FALSE,
    clean_build_artifacts = FALSE,
    clean_temp_files = FALSE,
    render_readme = FALSE,
    open_site = open
  )
}

#' Build Vignettes Only
#'
#' Rebuild vignettes without the full site build. Useful when iterating on vignette content.
#'
#' @param clean Logical. Clean old vignette builds first? Default TRUE.
#' @export
build_vignettes_only <- function(clean = TRUE) {
  .assert_pkg_root()

  if (isTRUE(clean)) {
    .clean_build_artifacts(clean_site = FALSE, clean_vignettes = TRUE, clean_check = FALSE)
  }

  .local_message("Building vignettes...")

  tools::buildVignettes(dir = ".")
  .copy_built_vignettes_to_inst_doc()

  .local_success("Vignettes built. Output in inst/doc/ (HTML)")

  invisible(TRUE)
}

#' Check Package Only
#'
#' Run R CMD check without building the pkgdown site.
#'
#' @param as_cran Logical. Use CRAN check settings? Default TRUE.
#' @param build_vignettes Logical. Build vignettes during check build phase? Default FALSE.
#'   If FALSE, passes `--no-build-vignettes` via `build_args`.
#' @param document Logical. Run `devtools::document()` before check? Default TRUE.
#'
#' @export
check_only <- function(as_cran = TRUE, build_vignettes = FALSE, document = TRUE) {
  .assert_pkg_root()
  .require_pkg("devtools")

  if (isTRUE(document)) devtools::document()

  args <- c("--no-manual")
  if (isTRUE(as_cran)) args <- c(args, "--as-cran")

  build_args <- if (isTRUE(build_vignettes)) NULL else "--no-build-vignettes"

  devtools::check(args = args, build_args = build_args)
}

# ==============================================================================
# PRINT USAGE ON SOURCE
# ==============================================================================

.local_message("Package build workflow loaded. Available commands:")
.local_message("  publish_docs()         - Full production build with all checks")
.local_message("  quick_site()           - Fast dev build (no checks)")
.local_message("  build_vignettes_only() - Rebuild vignettes only")
.local_message("  check_only()           - Run R CMD check only")
.local_message("  deep_clean()           - Remove ALL generated artifacts")
