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

.set_safe_tmpdir <- function() {
  td <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
  Sys.setenv(TMPDIR = td, TMP = td, TEMP = td)
  invisible(td)
}
.set_safe_tmpdir()

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
  .local_message(sprintf("Step %d/%d: %s", step_num, total, description), prefix = "----")
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
    .local_stop("Missing required package: '", pkg, "'. Install with install.packages('", pkg, "')")
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

  # Version check
  v <- tryCatch({
    suppressWarnings(
      system2("quarto", "--version", stdout = TRUE, stderr = TRUE)
    )
  }, error = function(e) NULL)

  if (!is.null(v) && length(v) > 0L && !grepl("ERROR:", v[1])) {
    .local_message("Quarto CLI version: ", v[1])
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

  # Check if quarto is the builder
  has_quarto <- any(grepl("^VignetteBuilder:.*quarto", desc))
  if (!has_quarto) {
    .local_message("VignetteBuilder is not 'quarto'. Using: ",
                   gsub("^VignetteBuilder:\\s*", "", desc[grepl("^VignetteBuilder:", desc)]))
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
    msg <- sprintf("Working tree has %d uncommitted change(s). Consider committing before publishing.", length(st))
    if (isTRUE(strict)) .local_stop(msg) else .local_warn(msg)
  } else {
    .local_success("Git working tree is clean")
  }

  invisible(st)
}

.check_namespace_current <- function() {
  # Check if NAMESPACE might be stale
  if (!file.exists("NAMESPACE")) {
    .local_warn("NAMESPACE file not found. Run devtools::document() first.")
    return(invisible(FALSE))
  }

  # Get modification times
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
#' Comprehensive removal of all intermediate and output files from previous
#' documentation, vignette, and package builds.
#'
#' @param clean_site Logical. Remove pkgdown site output?
#' @param clean_vignettes Logical. Remove vignette build artifacts?
#' @param clean_check Logical. Remove R CMD check artifacts?
#' @param clean_build Logical
.clean_build_artifacts <- function(clean_site = TRUE,
                                   clean_vignettes = TRUE,
                                   clean_check = TRUE,
                                   clean_build = TRUE,
                                   clean_temp = TRUE,
                                   verbose = TRUE) {

  removed_count <- 0L

  # --------------------------------------------------------------------------
  # PKGDOWN SITE ARTIFACTS
  # --------------------------------------------------------------------------
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

  # --------------------------------------------------------------------------
  # VIGNETTE BUILD ARTIFACTS
  # --------------------------------------------------------------------------
  if (isTRUE(clean_vignettes)) {
    # Fixed paths
    vignette_dirs <- c(
      file.path("inst", "doc"),
      file.path("vignettes", ".quarto"),
      file.path("vignettes", "_cache"),
      file.path("vignettes", "_freeze"),
      file.path("vignettes", "_site"),
      file.path("vignettes", "rsconnect"),
      file.path("doc")
    )

    for (path in vignette_dirs) {
      if (dir.exists(path)) {
        .safe_unlink(path)
        removed_count <- removed_count + 1L
      }
    }

    # Pattern-matched vignette artifacts
    if (dir.exists("vignettes")) {
      # HTML output files (not .qmd or .Rmd source)
      vignette_html <- list.files("vignettes", pattern = "\\.html$", full.names = TRUE)
      lapply(vignette_html, .safe_unlink)
      removed_count <- removed_count + length(vignette_html)

      # *_files directories (knitr/rmarkdown figure directories)
      files_dirs <- list.files("vignettes", pattern = "_files$", full.names = TRUE)
      files_dirs <- files_dirs[dir.exists(files_dirs)]
      lapply(files_dirs, .safe_unlink)
      removed_count <- removed_count + length(files_dirs)

      # *_cache directories (knitr cache)
      cache_dirs <- list.files("vignettes", pattern = "_cache$", full.names = TRUE)
      cache_dirs <- cache_dirs[dir.exists(cache_dirs)]
      lapply(cache_dirs, .safe_unlink)
      removed_count <- removed_count + length(cache_dirs)

      # Quarto intermediate files
      quarto_files <- list.files("vignettes", pattern = "\\.(tex|log|aux)$", full.names = TRUE)
      lapply(quarto_files, .safe_unlink)
      removed_count <- removed_count + length(quarto_files)
    }

    if (verbose) .local_message("Cleaned vignette build artifacts")
  }

  # --------------------------------------------------------------------------
  # R CMD CHECK ARTIFACTS
  # --------------------------------------------------------------------------
  if (isTRUE(clean_check)) {
    # .Rcheck directories (anywhere in package root)
    check_dirs <- list.files(".", pattern = "\\.Rcheck$", full.names = TRUE, all.files = TRUE)
    lapply(check_dirs, .safe_unlink)
    removed_count <- removed_count + length(check_dirs)

    # Common check output locations
    check_artifacts <- c(
      ".Rcheck",
      "..Rcheck"
    )

    for (path in check_artifacts) {
      if (dir.exists(path)) {
        .safe_unlink(path)
        removed_count <- removed_count + 1L
      }
    }

    if (verbose) .local_message("Cleaned R CMD check artifacts")
  }

  # --------------------------------------------------------------------------
  # PACKAGE BUILD ARTIFACTS
  # --------------------------------------------------------------------------
  if (isTRUE(clean_build)) {
    # Source package tarballs (*.tar.gz)
    tarballs <- list.files(".", pattern = "\\.tar\\.gz$", full.names = TRUE)
    # Only remove if they look like package builds (pkgname_version.tar.gz)
    pkg_name <- tryCatch(.get_pkg_name(), error = function(e) NULL)
    if (!is.null(pkg_name)) {
      pkg_tarballs <- tarballs[grepl(paste0("^\\./", pkg_name, "_"), tarballs)]
      lapply(pkg_tarballs, .safe_unlink)
      removed_count <- removed_count + length(pkg_tarballs)
    }

    # Binary packages (*.zip on Windows, *.tgz on macOS)
    binaries <- list.files(".", pattern = "\\.(zip|tgz)$", full.names = TRUE)
    if (!is.null(pkg_name)) {
      pkg_binaries <- binaries[grepl(paste0("^\\./", pkg_name, "_"), binaries)]
      lapply(pkg_binaries, .safe_unlink)
      removed_count <- removed_count + length(pkg_binaries)
    }

    # Build directory
    .safe_unlink("build")

    # Meta directory (sometimes created during build)
    .safe_unlink("Meta")

    if (verbose) .local_message("Cleaned package build artifacts")
  }

  # --------------------------------------------------------------------------
  # TEMPORARY AND CACHE FILES
  # --------------------------------------------------------------------------
  if (isTRUE(clean_temp)) {
    temp_files <- c(
      ".Rhistory",
      ".RData",
      ".Rproj.user",
      ".DS_Store",
      "Thumbs.db"
    )

    for (f in temp_files) {
      if (file.exists(f)) {
        .safe_unlink(f)
        removed_count <- removed_count + 1L
      }
    }

    # Root-level Quarto artifacts
    root_quarto <- c(".quarto", "_cache", "_freeze")
    for (d in root_quarto) {
      if (dir.exists(d)) {
        .safe_unlink(d)
        removed_count <- removed_count + 1L
      }
    }

    # Coverage artifacts
    covr_dirs <- list.files(".", pattern = "^covr", full.names = TRUE)
    lapply(covr_dirs, .safe_unlink)
    .safe_unlink("codecov")
    .safe_unlink("coverage.html")

    # Rcpp artifacts (if applicable)
    .safe_unlink(file.path("src", "*.o"))
    .safe_unlink(file.path("src", "*.so"))
    .safe_unlink(file.path("src", "*.dll"))

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

  # Prompt for confirmation in interactive mode
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

  result <- tryCatch({
    devtools::test(stop_on_failure = FALSE)
  }, error = function(e) {
    .local_warn("Test execution failed: ", conditionMessage(e))
    return(NULL)
  })

  if (is.null(result)) {
    if (stop_on_failure) .local_stop("Tests failed to execute")
    return(invisible(FALSE))
  }

  # Check for failures
  if (inherits(result, "testthat_results")) {
    n_fail <- sum(vapply(result, function(x) x$failed, integer(1)))
    n_skip <- sum(vapply(result, function(x) x$skipped, integer(1)))
    n_pass <- sum(vapply(result, function(x) x$passed, integer(1)))

    .local_message(sprintf("Tests: %d passed, %d failed, %d skipped", n_pass, n_fail, n_skip))

    if (n_fail > 0 && stop_on_failure) {
      .local_stop("Package tests failed. Fix tests before publishing.")
    }
  }

  invisible(TRUE)
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
    .local_warn(sprintf("Found %d potential spelling errors. Review with spelling::spell_check_package()", nrow(errors)))
    if (stop_on_failure) {
      print(errors)
      .local_stop("Spell check found errors")
    }
  } else {
    .local_success("No spelling errors found")
  }

  invisible(errors)
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
#' @param strict_quarto Logical. Fail if Quarto not found? Default TRUE.
#' @param clean_site Logical. Remove old pkgdown output? Default TRUE.
#' @param clean_vignettes Logical. Remove old vignette builds? Default TRUE.
#' @param clean_check_artifacts Logical. Remove .Rcheck directories? Default TRUE.
#' @param clean_build_artifacts Logical. Remove package tarballs/binaries? Default TRUE.
#' @param clean_temp_files Logical. Remove temp files (.RData, etc.)? Default TRUE.
#' @param render_readme Logical. Render README.Rmd if present? Default TRUE.
#' @param check_args Character vector. Additional arguments for R CMD check.
#'   Default includes `--ignore-vignettes` to work around a known Windows/Quarto
#'   bug where system2() incorrectly passes env arguments. Vignettes are built
#'   separately via devtools::build_vignettes() before the check runs.
#' @param open_site Logical. Open site in browser when done? Default FALSE.
#' @param verbose Logical. Print detailed progress? Default TRUE.
#'
#' @details
#' ## Windows Quarto Workaround
#'
#' On Windows, R CMD check with Quarto vignettes fails with:
#' ```
#' ERROR: Unknown command "TMPDIR=...". Did you mean command "create"?
#' ```
#'
#' This is a known bug in how system2() handles the `env` argument on Windows.
#' The workaround used here is:
#' 1. Build vignettes separately with `devtools::build_vignettes()` (before check)
#' 2. Run R CMD check with `--ignore-vignettes` to skip re-building them
#'
#' The vignettes are still validated and included in the package; they're just
#' not re-built during the check phase.
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
                               check_args = c("--no-manual", "--as-cran", "--ignore-vignettes"),
                               open_site = FALSE,
                               verbose = TRUE) {

  start_time <- Sys.time()
  .assert_pkg_root()

  pkg_name <- .get_pkg_name()
  pkg_version <- .get_pkg_version()
  .local_message("Building documentation for ", pkg_name, " v", pkg_version)

  # Count total steps for progress reporting
  total_steps <- 6 + run_tests + run_checks + run_spell_check + render_readme
  current_step <- 0

  # --------------------------------------------------------------------------
  # PREFLIGHT CHECKS
  # --------------------------------------------------------------------------
  current_step <- current_step + 1
  .local_step(current_step, total_steps, "Preflight checks")

  .check_pkgdown_config(strict = FALSE)
  .check_git_clean(strict = FALSE)
  .check_quarto_cli(strict = strict_quarto)
  .check_vignette_config(strict = FALSE)
  .check_namespace_current()
  .check_license()
  .check_readme()

  # --------------------------------------------------------------------------
  # LOAD TOOLING
  # --------------------------------------------------------------------------
  .require_pkg("devtools")
  .require_pkg("pkgdown")
  .require_pkg("roxygen2")

  # --------------------------------------------------------------------------
  # CLEAN OLD ARTIFACTS
  # --------------------------------------------------------------------------
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

  # --------------------------------------------------------------------------
  # RENDER README (if applicable)
  # --------------------------------------------------------------------------
  if (isTRUE(render_readme)) {
    current_step <- current_step + 1
    .local_step(current_step, total_steps, "Rendering README")
    .render_readme()
  }

  # --------------------------------------------------------------------------
  # GENERATE DOCUMENTATION (roxygen2)
  # --------------------------------------------------------------------------
  current_step <- current_step + 1
  .local_step(current_step, total_steps, "Generating Rd documentation (roxygen2)")

  devtools::document()
  .local_success("Documentation generated")

  # --------------------------------------------------------------------------
  # RUN TESTS
  # --------------------------------------------------------------------------
  if (isTRUE(run_tests)) {
    current_step <- current_step + 1
    .local_step(current_step, total_steps, "Running package tests")
    .run_tests(stop_on_failure = TRUE)
  }

  # --------------------------------------------------------------------------
  # SPELL CHECK
  # --------------------------------------------------------------------------
  if (isTRUE(run_spell_check)) {
    current_step <- current_step + 1
    .local_step(current_step, total_steps, "Running spell check")
    .run_spell_check(stop_on_failure = FALSE)
  }

  # --------------------------------------------------------------------------
  # BUILD VIGNETTES
  # --------------------------------------------------------------------------
  current_step <- current_step + 1
  .local_step(current_step, total_steps, "Building vignettes")

  tryCatch({
    devtools::build_vignettes()
    .local_success("Vignettes built successfully")
  }, error = function(e) {
    .local_stop("Vignette build failed: ", conditionMessage(e))
  })

  # --------------------------------------------------------------------------
  # R CMD CHECK
  # --------------------------------------------------------------------------
  if (isTRUE(run_checks)) {
    current_step <- current_step + 1
    .local_step(current_step, total_steps, "Running R CMD check")

    # Log if using vignette workaround
    if ("--ignore-vignettes" %in% check_args) {
      .local_message("Using --ignore-vignettes (vignettes already built above)")
      .local_message("This avoids the Windows/Quarto system2() env bug")
    }

    check_result <- tryCatch({
      devtools::check(args = check_args, check_dir = tempdir())
    }, error = function(e) {
      .local_stop("R CMD check failed: ", conditionMessage(e))
    })

    # Check for errors/warnings
    if (!is.null(check_result)) {
      n_errors <- length(check_result$errors)
      n_warnings <- length(check_result$warnings)
      n_notes <- length(check_result$notes)

      .local_message(sprintf("Check results: %d errors, %d warnings, %d notes",
                             n_errors, n_warnings, n_notes))

      if (n_errors > 0) {
        .local_stop("R CMD check produced errors. Fix before publishing.")
      }
    }
  }

  # --------------------------------------------------------------------------
  # INSTALL PACKAGE
  # --------------------------------------------------------------------------
  current_step <- current_step + 1
  .local_step(current_step, total_steps, "Installing package")

  devtools::install(upgrade = "never", build_vignettes = TRUE, quiet = TRUE)
  .local_success("Package installed")

  # --------------------------------------------------------------------------
  # BUILD PKGDOWN SITE
  # --------------------------------------------------------------------------
  current_step <- current_step + 1
  .local_step(current_step, total_steps, "Building pkgdown site")

  tryCatch({
    pkgdown::build_site()
    .local_success("pkgdown site built successfully")
  }, error = function(e) {
    .local_stop("pkgdown site build failed: ", conditionMessage(e))
  })

  # --------------------------------------------------------------------------
  # SUMMARY
  # --------------------------------------------------------------------------
  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
  .local_message("", prefix = "====")
  .local_success(sprintf("Build complete for %s v%s in %.1f minutes", pkg_name, pkg_version, elapsed))

  if (file.exists("docs/index.html")) {
    .local_message("Site output: ", normalizePath("docs", mustWork = FALSE))
  }

  # --------------------------------------------------------------------------
  # OPTIONAL: OPEN SITE
  # --------------------------------------------------------------------------
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
#' Production-ready build with all checks enabled.
#' Use before pushing to GitHub or releasing.
#'
#' @details
#' Vignettes are built separately before R CMD check runs, then check is run
#' with `--ignore-vignettes` to avoid the Windows/Quarto system2() bug.
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
    check_args = c("--no-manual", "--as-cran", "--ignore-vignettes"),
    open_site = FALSE
  )
}

#' Quick Site Build (Development)
#'
#' Fast iteration build without tests or R CMD check.
#' Use during active development to preview changes.
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
#' Rebuild vignettes without full site build.
#' Useful when iterating on vignette content.
#'
#' @param clean Logical. Clean old vignette builds first? Default TRUE.
#' @export
build_vignettes_only <- function(clean = TRUE) {
  .assert_pkg_root()
  .require_pkg("devtools")

  if (isTRUE(clean)) {
    .clean_build_artifacts(clean_site = FALSE, clean_vignettes = TRUE, clean_check = FALSE)
  }

  .local_message("Building vignettes...")
  devtools::build_vignettes()
  .local_success("Vignettes built. Output in inst/doc/")

  invisible(TRUE)
}

#' Check Package Only
#'
#' Run R CMD check without building site.
#'
#' @param as_cran Logical. Use CRAN check settings? Default TRUE.
#' @param ignore_vignettes Logical. Skip vignette re-building during check?
#'   Default TRUE to avoid Windows/Quarto system2() bug. Set to FALSE if you
#'   don't use Quarto vignettes or are on Linux/macOS.
#'
#' @details
#' By default, uses `--ignore-vignettes` to work around a known Windows bug
#' where Quarto vignettes fail during R CMD check due to system2() incorrectly
#' passing environment variables as command arguments.
#'
#' If you need to test vignette building during check, either:
#' - Set `ignore_vignettes = FALSE` (may fail on Windows with Quarto)
#' - Build vignettes separately with `build_vignettes_only()`
#'
#' @export
check_only <- function(as_cran = TRUE, ignore_vignettes = TRUE) {
  .assert_pkg_root()
  .require_pkg("devtools")

  devtools::document()

  args <- "--no-manual"
  if (isTRUE(as_cran)) args <- c(args, "--as-cran")
  if (isTRUE(ignore_vignettes)) {
    args <- c(args, "--ignore-vignettes")
    .local_message("Using --ignore-vignettes to avoid Windows/Quarto bug")
  }

  devtools::check(args = args)
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
