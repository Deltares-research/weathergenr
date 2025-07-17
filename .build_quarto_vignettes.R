# .build_quarto_vignettes.R

# Set TMPDIR
Sys.setenv(TMPDIR = tempdir())

# Ensure output directory exists
if (!dir.exists("inst/doc")) dir.create("inst/doc", recursive = TRUE)

# List all .qmd vignette source files
qmd_files <- list.files("vignettes", pattern = "\\.qmd$", full.names = TRUE)

# Render each .qmd in place and copy HTML output to inst/doc
for (qmd in qmd_files) {
  # Render the .qmd to HTML in the same directory
  quarto::quarto_render(input = qmd)

  # Construct HTML output filename
  html_file <- sub("\\.qmd$", ".html", basename(qmd))

  # Copy rendered HTML to inst/doc/
  file.copy(
    from = file.path("vignettes", html_file),
    to   = file.path("inst/doc", html_file),
    overwrite = TRUE
  )
}
