# AGENTS.md

> **Canonical.** Source-of-truth agent instructions for every runtime: Codex reads
> `AGENTS.md` directly; `CLAUDE.md` imports it verbatim (`@AGENTS.md`).

## Overview

`weathergenr` is an R package implementing a gridded, semiparametric, multivariate,
multisite stochastic weather generator (Steinschneider & Brown 2013) for bottom-up
climate vulnerability assessment. See `README.md` for the method framing.

## Repo Map

- `R/` — package source; the only place to change behavior.
- `man/`, `NAMESPACE` — roxygen2 output. Regenerate with `devtools::document()`.
- `docs/` — generated pkgdown site.
- `tools/` — repo-only dev scripts, `.Rbuildignore`d and not shipped.
- `dev/` — development process: the work board (`tasks/`, `TODO.md`), the closed-work
  ledger (`LOG.md`), task-backing drafts, and `scripts/` for exploratory code that must
  not ship. Never shipped; see `dev/README.md`.
- `dev/scripts/zarr-prototype/` — abandoned Zarr I/O prototypes. Not package code, not
  tested, not on the maintenance path; do not extend unless the task says so.
- `vignettes/` — Quarto `.qmd` only, built by `quarto`, not `knitr`.

## Key Commands

```r
devtools::load_all()                                   # iterate
testthat::test_file("tests/testthat/test-calendar.R")  # single file
devtools::test()                                       # full suite
devtools::document()                                   # after any roxygen change
source("tools/build_site_tools.R"); check_only()       # R CMD check
check_only(build_vignettes = TRUE)                     # release gate — matches CI
```

```bash
Rscript tools/lint.R --changed   # lint changed R files; exit 1 = lints found
Rscript tools/lint.R             # lint whole package — what CI runs
```

End-to-end numeric baseline — opt-in local gate, ~40 s (see *Workflow* below).
The gate reads the `WEATHERGENR_BASELINE` environment variable, so the run command
is shell-specific; the repo's default shell is PowerShell 7:

```powershell
$env:WEATHERGENR_BASELINE = "1"                       # PowerShell — NOT VAR=1 cmd
Rscript -e 'devtools::test(filter = "baseline-e2e")'
Rscript tools/record_baseline.R --dry-run             # inspect the delta
Rscript tools/record_baseline.R --force               # re-record, once accepted
```

From an R console (RStudio), no shell variable needed:

```r
Sys.setenv(WEATHERGENR_BASELINE = "1")
devtools::test(filter = "baseline-e2e")
```

IMPORTANT: never `source("tools/dev_workflow.R")`. It is a runnable notes file whose top
level calls `publish_docs()` — the full 9-step clean, render, check, install, and site
build. The functions it documents are defined in `tools/build_site_tools.R`
(`check_only`, `quick_site`, `build_vignettes_only`, `deep_clean`, `publish_docs`);
source that file instead.

## Conventions

- Exports come from roxygen `@export` plus `devtools::document()`; never hand-edit
  `NAMESPACE` or `man/*.Rd`.
- Dates run on a 365-day no-leap calendar (`calendar.R`). Assume Feb 29 is absent from
  any date vector reaching the generator; keep new date logic leap-free rather than
  special-casing leap days downstream.
- Grid-cell iteration runs in a PSOCK cluster (`parallel::makeCluster` +
  `doParallel::registerDoParallel`, `R/generator.R`). Workers inherit no globals:
  reference package functions or pass objects explicitly, and thread the `seed`
  argument through so runs stay reproducible.
- Plotting lives in the `*_plots.R` sibling of its module
  (`warm_filtering_plots.R`, `wavelet_plots.R`, `evaluate_generator_plots.R`); keep
  ggplot2 code out of the computational files.
- New dev-only files at the repo root must be added to `.Rbuildignore`, or `R CMD check`
  raises a non-standard-file NOTE.
- Vendored agent pins (`.claude/skills/`, `.agents/skills/`, `.claude/agents/`,
  `.codex/agents/`) are resync artifacts of `.claude/agent-manifest.yml` — never edit
  them in place; see the global `AGENTS.md` self-improvement rule.

## Workflow

- Change behavior → update or add the matching `tests/testthat/test-<module>.R`.
- Change *user-visible* behavior → add a `NEWS.md` entry under the
  `# weathergenr (development version)` heading, in the same commit. User-visible means
  a user of the package could notice: exported functions, arguments, defaults, outputs,
  errors, or documented workflows. Internal-only work — refactors, tests, CI, tooling,
  roxygen fixes — gets no entry. Writing entries as the work lands is what keeps a
  release cheap; do not reconstruct them from `git log` at release time.
- Before finishing: run the affected test file, then `devtools::test()`. Run
  `check_only()` when the change touches exports, documentation, or dependencies.
- Change anything that could move **numeric output** — resampling, WARM, wavelets,
  calendar logic, evaluation statistics, or a "behavior-preserving" refactor of any of
  them — run the end-to-end baseline gate **before and after** the change. It is opt-in
  and local: `devtools::test()`, `R CMD check`, and CI all skip it, so **a green CI does
  not mean the baseline was checked**. Definitions live in
  `tests/testthat/helper-baseline.R`, which `tools/record_baseline.R` sources — there is
  one definition of the scenarios and the fingerprint, never two. A failure names which
  keys moved; treat it as a question, not a chore, and re-record only once the delta is
  reviewed and intended. Recording blesses whatever the package currently produces.
- Report which commands were run and what they returned; never describe a check that
  did not run.
- The `DESCRIPTION` version is **not** bumped per commit. `.git-workflow.yml` sets
  `cadence: manual`, and `1.2.0.9000`-style development versions hold steady until an
  explicitly requested release. See that file for the release gate and apply command.
- IMPORTANT: a downstream consumer pins this package by **Git tag**. The `blueearth_cst`
  toolkit installs `tanerumit/weathergenr@v1.2.0` via `remotes::install_github()` from
  `dev/scripts/install_weathergenr.R`, driven by `pixi run install-rdeps`. Consequences:
  a release is not delivered until its tag is **pushed to `origin`** — with
  `auto_push: false` that is a deliberate, owner-run final step, so never report a
  release as complete at the commit. Tags are consumed artifacts: never move, delete, or
  re-point a published one. And user-visible change stranded on an unreleased `.9000`
  is invisible downstream, which is the reason to cut releases promptly rather than let
  development versions drift.

### Work board

- Admit before you track: open a `dev/` board item only when work spans sessions or
  needs visibility beyond the current one. Work that its diff and commit message fully
  explain gets no board item, no task ID, and no draft note.
- `dev/tasks/` is the live board — one note per open item. `dev/TODO.md` is generated
  by `brain todo render`; never hand-edit it. Close with `brain todo done <id>`, which
  appends the `dev/LOG.md` row and removes the note.
- A `dev/drafts/` note that anything durable cites — a test, a module, a tracked
  document — is promoted at closure with its citations updated, never deleted.
- The `todo-board` skill owns the board mechanic; `project-system` owns the
  surrounding lifecycle. Neither is restated here.

## Hard Constraints

- IMPORTANT: do not hand-edit generated output — `man/`, `NAMESPACE`, `docs/`, or
  built vignette artifacts.
- Do not add to `Imports` in `DESCRIPTION` without justification; the dependency set is
  deliberately small and CI builds on three OSes and three R versions.
- Do not touch `.github/workflows/`, the `DESCRIPTION` version, or release tags unless
  the task asks for it.

## References

- `README.md` — method framing and the three coupled components; read before changing
  WARM, resampling, or perturbation behavior.
- `tools/build_site_tools.R` — read before running or debugging any docs or site build.
- `_pkgdown.yml` — read when adding or renaming a vignette; its `articles` list is
  explicit and breaks the site build when stale.
- `.github/workflows/R-CMD-check.yaml` — read when a change might behave differently on
  macOS, Windows, or oldrel R.
