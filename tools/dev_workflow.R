# ==============================================================================
# R PACKAGE DEVELOPMENT WORKFLOW
# ==============================================================================
#
# This script documents the complete workflow for developing, testing, and
# publishing the weathergenr R package.
#
# Prerequisites:
#   - R (>= 4.1)
#   - Quarto CLI installed and on PATH
#   - Required packages: devtools, pkgdown, roxygen2, testthat, quarto
#
# Usage:
#   1. Open this file in RStudio
#   2. Run sections as needed using Ctrl+Enter (line) or Ctrl+Shift+Enter (chunk)
#   3. Follow the workflow steps in order when preparing a release
#
# ==============================================================================


# Load the build workflow tools
source("tools/build_site_tools.R")

# ------------------------------------------------------------------------------

# 1. Full publish workflow (RECOMMENDED)
publish_docs()

# The publish_docs() function executes these steps in order:
#
#   Step 1/9: Preflight checks (git status, quarto, config files)
#   Step 2/9: Clean all build artifacts
#   Step 3/9: Render README.Rmd (if present)
#   Step 4/9: Generate Rd documentation (roxygen2)
#   Step 5/9: Run package tests (testthat)
#   Step 6/9: Build vignettes (Quarto)
#   Step 7/9: Run R CMD check (with --ignore-vignettes workaround)
#   Step 8/9: Install package locally
#   Step 9/9: Build pkgdown site to docs/
#
# Output:
#   - docs/          : pkgdown website (ready for GitHub Pages)
#   - inst/doc/      : Built vignettes
#   - man/           : Rd documentation files

# 4b. (Alternative) Quick site build for development preview (no checks)
# quick_site(open = TRUE)


# ==============================================================================
# STEP 2: PUSH TO GITHUB
# ==============================================================================
#
# After publish_docs() completes successfully, commit and push
#
# ------------------------------------------------------------------------------

# 2a. Stage all changes (run in terminal or use RStudio Git pane)
# system("git add -A")

# 2b. Commit with descriptive message
# system('git commit -m "Your descriptive commit message"')

# 2c. Push to remote
# system("git push origin main")

# ==============================================================================
# MAINTENANCE COMMANDS
# ==============================================================================

# Deep clean: Remove ALL generated artifacts
# Use when builds are broken or site looks wrong
# deep_clean()

# Check only: Run R CMD check without site build
# check_only(as_cran = TRUE, ignore_vignettes = TRUE)

# Build vignettes only: Faster iteration on vignette content
# build_vignettes_only(clean = TRUE)

# Quick site: Fast site rebuild for previewing (no tests/checks)
# quick_site(open = TRUE)


# ==============================================================================
# TROUBLESHOOTING
# ==============================================================================
#
# Problem: Stale documentation or broken site
# Solution:
#   deep_clean()
#   publish_docs()
#
# Problem: Tests failing
# Solution:
#   devtools::test()
#   # Fix failing tests, then re-run
#
# Problem: Vignette build errors
# Solution:
#   # Check .qmd syntax in vignettes/
#   quarto::quarto_render("vignettes/problem_vignette.qmd")
#   # Fix errors, then:
#   build_vignettes_only()
#
# Problem: R CMD check warnings/notes
# Solution:
#   # Review output from check_only()
#   # Address issues or document why they're acceptable
#
# Problem: Windows Quarto error "Unknown command TMPDIR=..."
# Solution:
#   # This is handled automatically by --ignore-vignettes flag
#   # Vignettes are built separately before R CMD check
#
# Problem: pkgdown site missing functions
# Solution:
#   # Check _pkgdown.yml reference section
#   # Ensure all exported functions are documented with @export
#
# ==============================================================================


# ==============================================================================
# QUICK REFERENCE
# ==============================================================================
#
# Daily Development:
#   devtools::document()
#   devtools::load_all()
#   devtools::test()
#
# Before Commit:
#   check_only()
#
# Before Push:
#   publish_docs()
#
# Maintenance:
#   deep_clean()
#   quick_site()
#   build_vignettes_only()
#
# ==============================================================================
