# Load required devtools ecosystem
if (!requireNamespace("devtools")) install.packages("devtools")
if (!requireNamespace("roxygen2")) install.packages("roxygen2")
if (!requireNamespace("testthat")) install.packages("testthat")
if (!requireNamespace("usethis")) install.packages("usethis")
if (!requireNamespace("styler")) install.packages("styler")
if (!requireNamespace("spelling")) install.packages("spelling")
if (!requireNamespace("pkgdown")) install.packages("pkgdown")

library(devtools)
library(usethis)
library(styler)
library(spelling)
library(pkgdown)

# -------- CONFIGURATION ----------
version_bump <- "patch" # or "minor", "major"
# ---------------------------------

# 1. Bump version
if (interactive()) {
  usethis::use_version(version_bump)  # or "minor", "major"
}

# 2. Style all code
message("[STYLE] Formatting code")
styler::style_pkg()

# 3. Document functions and rebuild NAMESPACE
message("[DOCS] Rebuilding documentation")
devtools::document()

# 4. Run spell check on documentation and vignettes
message("[SPELL CHECK] Checking spelling")
spelling::spell_check_package()

# 5. Run unit tests
message("[TESTS] Running tests")
devtools::test()

# 6. Run full CMD check
message("[CHECK] Running devtools::check()")
devtools::check()

# 7. Build PDF manual (requires LaTeX)
message("[MANUAL] Building manual")
devtools::build_manual()

# 8. Build source tarball for submission
message("[BUILD] Building package tarball")
devtools::build()

# 9. Rebuild site (if pkgdown is used)
if (file.exists("pkgdown.yml")) {
  message("[PKGDOWN] Building documentation website")
  pkgdown::build_site()
}

# 10. Optional: Run GitHub Action checks (if configured)
if (file.exists(".github/workflows/R-CMD-check.yaml")) {
  message("[GITHUB ACTIONS] Reminder: CI will recheck on push.")
}

message("\n??? All done. Review the following before submitting:")
message("- Check CRAN submission policies: https://cran.r-project.org/web/packages/policies.html")
message("- Review R CMD Check output: no ERRORs, and Warnings/Notes are explainable.")
message("- Check that version and NEWS.md are updated.")
