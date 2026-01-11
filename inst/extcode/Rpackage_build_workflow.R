
mypackage <- "weathergenr"

# Set your package directory
pkg_dir <- "."

# Restart session (works only in RStudio)
if ("rstudioapi" %in% rownames(installed.packages())) {
  rstudioapi::restartSession()
}


# Clean start
rs.restartR()

# Remove all objects
rm(list = ls())

# Document package
devtools::document(pkg_dir)

# Run package checks
devtools::check(pkg_dir)

# Run tests (if any)
devtools::test(pkg_dir)

# Build and install locally
devtools::build(pkg_dir)
devtools::install(pkg_dir)

# Optionally, spell check
devtools::spell_check(pkg_dir)

# Optionally, build README and site
if (requireNamespace("pkgdown", quietly = TRUE)) {
  devtools::build_readme(pkg_dir)
  pkgdown::build_site(pkg_dir)
}

# Try loading and running a function (replace 'yourpackagename' and 'your_function')
library(mypackage, character.only = TRUE)


# Before releasing
devtools::check_rhub()     # Checks on remote platforms, mimics CRAN servers
devtools::check_win_release()  # Checks on Windows

devtools::build_vignettes()
