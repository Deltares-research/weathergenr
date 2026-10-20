# AGENTS.md

`AGENTS.md` is canonical; `CLAUDE.md` imports it with `@AGENTS.md`.

## Overview

`weathergenr` is an R package implementing a gridded, semiparametric, multivariate,
multisite stochastic weather generator (Steinschneider & Brown 2013) for bottom-up
climate vulnerability assessment. See `README.md` for the method framing.

## Repo Map

- `R/` — package source; change package behavior here.
- `man/`, `NAMESPACE` — roxygen2 output.
- `docs/` — generated pkgdown site.
- `tools/` — repo-only dev scripts, `.Rbuildignore`d and not shipped.
- `dev/` — development process: the work board (`tasks/`, `TODO.md`), the closed-work
  ledger (`LOG.md`), task-backing drafts, and `scripts/` for exploratory code that must
  not ship. Never shipped; see `dev/README.md`.
- `dev/scripts/zarr-prototype/` — abandoned Zarr I/O prototypes. Not package code, not
  tested, not on the maintenance path; do not extend unless the task says so.
- `vignettes/` — Quarto `.qmd` only, built by `quarto`, not `knitr`.

## Key Commands

For a fresh clone or an R/library mismatch, follow `dev/README.md` before running checks.

```r
devtools::load_all()                                   # iterate
testthat::test_file("tests/testthat/test-calendar.R")  # single file
devtools::test()                                       # full suite
devtools::document()                                   # after any roxygen change
source("tools/build_site_tools.R"); check_only()       # R CMD check
check_only(build_vignettes = TRUE)                     # release gate — matches CI
```

```powershell
Rscript tools/lint.R --changed   # lint changed R files
Rscript tools/lint.R             # lint whole package — what CI runs
```

`object_usage_linter` resolves symbols against the installed package. After adding
internal objects, run `devtools::install()` before interpreting binding warnings;
do not mask stale-install findings with code changes or additions to `R/globals.R`.

For the opt-in numeric baseline gate (see Workflow), use PowerShell:

```powershell
$env:WEATHERGENR_BASELINE = "1"
Rscript -e 'devtools::test(filter = "baseline-e2e")'
Rscript tools/record_baseline.R --dry-run             # inspect the delta
Rscript tools/record_baseline.R --force               # re-record, once accepted
```

IMPORTANT: never `source("tools/dev_workflow.R")`: its top level calls `publish_docs()`
and runs the full publishing pipeline. Source `tools/build_site_tools.R` for helpers.

## Conventions

- Exports come from roxygen `@export`; regenerate with `devtools::document()`.
- Dates run on a 365-day no-leap calendar (`R/calendar.R`). Assume Feb 29 is absent from
  any date vector reaching the generator; keep new date logic leap-free rather than
  special-casing leap days downstream.
- Realization iteration runs in a PSOCK cluster (`parallel::makeCluster` +
  `doParallel::registerDoParallel`, `R/generator.R`) — realizations, not grid
  cells; grid cells are iterated sequentially. Workers inherit no globals:
  reference package functions or pass objects explicitly, and thread the `seed`
  argument through so runs stay reproducible.
- IMPORTANT: seed with `.set_seed_fixed_kind()`, never a bare `set.seed()`, in
  any code that can run on a worker. `clusterSetRNGStream()` switches workers to
  L'Ecuyer-CMRG, and `set.seed()` applies the seed to whatever generator is
  active, so a bare call can give serial and parallel runs different streams.
- `devtools::load_all()` does not reach PSOCK workers; they load the installed
  package. Run `devtools::install()` before worker tests, or use `R CMD check`,
  which installs before testing. A skipped parallel baseline is not verification.
- Plotting lives in the `*_plots.R` sibling of its module
  (`warm_filtering_plots.R`, `wavelet_plots.R`, `evaluate_generator_plots.R`); keep
  ggplot2 code out of the computational files.
- New dev-only files at the repo root must be added to `.Rbuildignore`, or `R CMD check`
  raises a non-standard-file NOTE.
- `.agent-system/agent-manifest.yml` owns project agent selection. Provider directories
  (`.agents/`, `.claude/`, `.codex/`, `.github/agents/`) are links or generated
  projections; edit canonical sources, never these projections.

## Workflow

- Change behavior → update or add the matching `tests/testthat/test-<module>.R`.
- Change *user-visible* behavior → add a `NEWS.md` entry under the
  `# weathergenr (development version)` heading, in the same commit. User-visible means
  a user of the package could notice: exported functions, arguments, defaults, outputs,
  errors, or documented workflows. Internal-only work — refactors, tests, CI, tooling,
  roxygen fixes — gets no entry.
- For package changes, run the affected test file, then `devtools::test()`. Run
  `check_only()` when exports, package documentation, or dependencies change.
  Instruction-only edits need reference and diff checks, not package execution.
- Change anything that could move numeric output — resampling, WARM, wavelets,
  calendar logic, evaluation statistics, or a "behavior-preserving" refactor of any of
  them — run the end-to-end baseline gate before and after the change. It is opt-in
  and local; ordinary tests, R CMD check, and CI skip it. Keep scenario and fingerprint
  definitions in `tests/testthat/helper-baseline.R`, shared by the test and recorder.
  Investigate moved keys; re-record only after the delta is reviewed and intended.
- Versioning is manual under `.git-workflow.yml`; hold development versions steady
  until an explicitly requested release.
- `blueearth_cst` consumes Git tags. A local release commit does not deliver a release;
  follow `dev/workflows/release-and-publish.md` for publication and downstream handoff.

### Work board

- Admit before you track: open a `dev/` board item only when work spans sessions or
  needs visibility beyond the current one. Work that its diff and commit message fully
  explain gets no board item, no task ID, and no draft note.
- `dev/tasks/` is the live board; `dev/TODO.md` is generated. Use `todo-board` for
  board operations and `project-system` for lifecycle; see `dev/README.md` for setup.
- A `dev/drafts/` note that anything durable cites — a test, a module, a tracked
  document — is promoted at closure with its citations updated, never deleted.

## Hard Constraints

- IMPORTANT: do not hand-edit generated output — `man/`, `NAMESPACE`, `docs/`, or
  built vignette artifacts.
- Do not add to `Imports` in `DESCRIPTION` without justification; keep dependencies small.
- Do not touch `.github/workflows/`, the `DESCRIPTION` version, or release tags unless
  the task asks for it.

## References

- Before changing WARM, resampling, or perturbation behavior, read `README.md`.
- Before running or debugging a docs or site build, read `tools/build_site_tools.R`.
- Before adding or renaming a vignette, read the explicit articles list in `_pkgdown.yml`.
- When a change might behave differently across OSes or R versions, read
  `.github/workflows/R-CMD-check.yaml`.
- For development setup, agent activation, and board setup, see `dev/README.md`.
- Before releasing, read `dev/workflows/release-and-publish.md`; use the apply and
  validate commands defined in `.git-workflow.yml`.
