# weathergenr (development version)

## Bug fixes

* `evaluate_weather_generator()` no longer emits tidyselect deprecation warnings
  ("Use of `.data` in tidyselect expressions was deprecated in tidyselect
  1.2.0"). Thirty `.data$` references inside `dplyr::select()`,
  `dplyr::rename()`, `tidyr::pivot_longer()` and `tidyr::pivot_wider()` now pass
  column names as strings. Results are unchanged; the `.data` pronoun is still
  used in data-masking verbs, where it remains correct.

## Internal

* `DESCRIPTION` declares authorship with `Authors@R` instead of the deprecated
  `Author`/`Maintainer` pair.
* The abandoned Zarr I/O prototypes moved from `inst/experiments/` to
  `dev/scripts/zarr-prototype/`, so no part of the development tree sits under
  `inst/` any more.
* Added `testthat` coverage for `R/qm_diagnostics.R` and the wavelet plotting
  functions, including regression tests for the spell-length and dry-day
  diagnostics' `NA` handling.
* Added pkgdown and lint GitHub Actions workflows, plus a `.lintr` baseline
  that scopes linting to semantic findings.

# weathergenr 1.2.0

## Bug fixes

* Fixed warnings raised by the evaluation and diagnostic plotting functions.
* Fixed netCDF template propagation when writing gridded output.

## Internal

* `inst/experiments/` and `inst/extcode/` are no longer installed with the
  package. Both held development scratch material — including abandoned Zarr
  I/O prototypes — that was previously shipping into user libraries.
* Consolidated the two `utils::globalVariables()` declarations into
  `R/globals.R` and removed `R/zzz_globals.R`. Of the 42 names the removed
  file declared, 28 were unused. `colorRampPalette` was among them: declaring
  a function there also silences R CMD check's "no visible global function
  definition", which would have masked the missing `grDevices` import had a
  call been reintroduced.
* Dropped the `LazyData` field. The package ships no `data/` directory, so it
  produced an R CMD check note without effect.
* `AGENTS.md` is now the canonical agent-instruction file, with `CLAUDE.md`
  importing it.

Versions before 1.2.0 predate this changelog; see the commit history.
