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
