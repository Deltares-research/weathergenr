# tools/build_quarto_vignettes.R

build_quarto_vignettes <- function(vignettes_dir = "vignettes",
                                   output_dir = file.path("inst", "doc"),
                                   pattern = "\\.qmd$",
                                   quiet = FALSE,
                                   stop_on_error = TRUE,
                                   copy_assets = TRUE,
                                   clean_tmp = TRUE,
                                   verify_rcmd_build = FALSE) {
  if (!file.exists("DESCRIPTION")) {
    stop("Run from the package root (DESCRIPTION not found).", call. = FALSE)
  }

  quarto_bin <- Sys.which("quarto")
  if (!nzchar(quarto_bin)) {
    stop("Quarto CLI not found on PATH. Install Quarto or add it to PATH.", call. = FALSE)
  }

  if (!dir.exists(vignettes_dir)) {
    return(invisible(NULL))
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  qmd_files <- list.files(vignettes_dir, pattern = pattern, full.names = TRUE)
  qmd_files <- qmd_files[!grepl("[/\\\\]_", qmd_files)]
  if (length(qmd_files) == 0) {
    return(invisible(NULL))
  }

  if (isTRUE(verify_rcmd_build)) {
    if (!nzchar(Sys.which("R"))) {
      stop("R executable not found on PATH; cannot run R CMD build.", call. = FALSE)
    }
    pkg_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
    pkg_name <- basename(pkg_root)
    parent_dir <- normalizePath("..", winslash = "/", mustWork = TRUE)

    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(parent_dir)

    out <- system2("R", c("CMD", "build", "--no-manual", pkg_name), stdout = TRUE, stderr = TRUE)
    tarballs <- list.files(parent_dir, pattern = paste0("^", pkg_name, "_.*\\.tar\\.gz$"), full.names = TRUE)
    if (length(tarballs) == 0) {
      stop("verify_rcmd_build: R CMD build did not produce a tarball.\n", paste(out, collapse = "\n"), call. = FALSE)
    }
    setwd(old_wd)
  }

  failures <- character(0)

  copy_dir_recursive <- function(src_dir, dst_dir) {
    if (!dir.exists(dst_dir)) dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
    all <- list.files(src_dir, full.names = TRUE, recursive = TRUE, all.files = TRUE, include.dirs = TRUE)

    dirs <- all[dir.exists(all)]
    if (length(dirs) > 0) {
      rel_dirs <- substring(dirs, nchar(src_dir) + 2)
      for (d in rel_dirs) dir.create(file.path(dst_dir, d), recursive = TRUE, showWarnings = FALSE)
    }

    files <- all[file.exists(all) & !dir.exists(all)]
    if (length(files) > 0) {
      rel_files <- substring(files, nchar(src_dir) + 2)
      for (i in seq_along(files)) {
        dst <- file.path(dst_dir, rel_files[i])
        dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
        file.copy(files[i], dst, overwrite = TRUE)
      }
    }

    invisible(TRUE)
  }

  run_quarto_render <- function(out_tmp, qmd_basename, out_html, quiet) {
    # Write Quarto output to a log so we can show it on failure.
    log_file <- file.path(out_tmp, "quarto_render.log")

    args <- c("render", qmd_basename, "--to", "html", "--output", out_html)
    if (isTRUE(quiet)) args <- c(args, "--quiet")

    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(out_tmp)

    # Capture both stdout+stderr to a file; also return it.
    out <- system2(quarto_bin, args, stdout = TRUE, stderr = TRUE)

    # Persist output for debugging
    writeLines(out, con = log_file, useBytes = TRUE)

    list(status = 0L, out = out, log_file = log_file)
  }

  for (qmd in qmd_files) {
    out_tmp <- tempfile("quarto_vignette_")
    dir.create(out_tmp, recursive = TRUE, showWarnings = FALSE)
    if (isTRUE(clean_tmp)) on.exit(unlink(out_tmp, recursive = TRUE, force = TRUE), add = TRUE)

    qmd_tmp <- file.path(out_tmp, basename(qmd))
    if (!file.copy(qmd, qmd_tmp, overwrite = TRUE)) {
      msg <- paste0("Failed to stage vignette into temp dir: ", qmd)
      if (isTRUE(stop_on_error)) stop(msg, call. = FALSE)
      warning(msg, call. = FALSE)
      failures <- c(failures, qmd)
      next
    }

    out_name <- sub("\\.qmd$", ".html", basename(qmd_tmp))

    # Render (with real output capture)
    old_quiet <- quiet
    res <- tryCatch(
      run_quarto_render(out_tmp, basename(qmd_tmp), out_name, quiet = old_quiet),
      error = function(e) e
    )

    rendered <- file.path(out_tmp, out_name)

    if (inherits(res, "error") || !file.exists(rendered)) {
      # Re-run once without --quiet to get full error context (if user asked for quiet)
      if (isTRUE(old_quiet)) {
        res2 <- tryCatch(run_quarto_render(out_tmp, basename(qmd_tmp), out_name, quiet = FALSE),
                         error = function(e) e)
        if (!inherits(res2, "error")) res <- res2
      }

      log_path <- if (!inherits(res, "error")) res$log_file else file.path(out_tmp, "quarto_render.log")
      log_txt <- if (file.exists(log_path)) paste(readLines(log_path, warn = FALSE), collapse = "\n") else ""

      msg <- paste0(
        "Quarto render failed for: ", qmd, "\n",
        "Quarto did not produce expected HTML: ", rendered, "\n\n",
        "Quarto log:\n", log_txt
      )

      if (isTRUE(stop_on_error)) stop(msg, call. = FALSE)
      warning(msg, call. = FALSE)
      failures <- c(failures, qmd)
      next
    }

    # Copy HTML
    if (!file.copy(rendered, file.path(output_dir, out_name), overwrite = TRUE)) {
      msg <- paste0("Failed to copy rendered HTML for: ", qmd)
      if (isTRUE(stop_on_error)) stop(msg, call. = FALSE)
      warning(msg, call. = FALSE)
      failures <- c(failures, qmd)
      next
    }

    # Copy common Quarto dependency dir: <stem>_files/
    if (isTRUE(copy_assets)) {
      stem <- sub("\\.html$", "", out_name)
      dep_dir <- file.path(out_tmp, paste0(stem, "_files"))
      if (dir.exists(dep_dir)) {
        dst_dep <- file.path(output_dir, basename(dep_dir))
        if (dir.exists(dst_dep)) unlink(dst_dep, recursive = TRUE, force = TRUE)
        copy_dir_recursive(dep_dir, dst_dep)
      }
    }
  }

  if (length(failures) > 0) {
    return(invisible(list(success = FALSE, failures = failures)))
  }
  invisible(list(success = TRUE, failures = character(0)))
}



# Standard use
#source("tools/.build_quarto_vignettes.R")
# build_quarto_vignettes(quiet = FALSE)
# enforce "standard toolchain sanity" occasionally:
#build_quarto_vignettes(verify_rcmd_build = TRUE)
