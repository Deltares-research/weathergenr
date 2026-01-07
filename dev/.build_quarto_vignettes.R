# .build_quarto_vignettes.R

# tools/build_quarto_vignettes.R

build_quarto_vignettes <- function() {
  if (!nzchar(Sys.which("quarto"))) {
    stop("Quarto not found on PATH. Install Quarto or add it to PATH.", call. = FALSE)
  }
  if (!requireNamespace("quarto", quietly = TRUE)) {
    stop("R package 'quarto' not installed. Install it to render qmd vignettes.", call. = FALSE)
  }

  if (!dir.exists("inst/doc")) dir.create("inst/doc", recursive = TRUE)

  qmd_files <- list.files("vignettes", pattern = "\\.qmd$", full.names = TRUE)
  if (length(qmd_files) == 0) return(invisible(NULL))

  for (qmd in qmd_files) {
    quarto::quarto_render(input = qmd)

    html_file <- sub("\\.qmd$", ".html", basename(qmd))
    file.copy(
      from = file.path("vignettes", html_file),
      to   = file.path("inst/doc", html_file),
      overwrite = TRUE
    )
  }

  invisible(TRUE)
}

# source("tools/build_quarto_vignettes.R")
# build_quarto_vignettes()
