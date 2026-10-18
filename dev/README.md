# dev/

Development process and project record. Everything here supports building the
project; none of it ships. The board mechanic is the `todo-board` skill; the
surrounding lifecycle (handoffs, reviews, closure promotion) is `project-system`.

## Layout

| Path | Holds |
|---|---|
| `TODO.base` | The open-work board, rendered live by Obsidian (Obsidian projects only) |
| `TODO.md` | The open-work board as a generated read-only table (non-Obsidian projects only; `todoboard render` writes it) |
| `LOG.md` | The closed-work ledger -- one row per item closed off the board |
| `tasks/` | Board notes -- one `t*.md` per **open** item (`type: todo-item`); the source the board renders |
| `drafts/` | Ephemeral task-backing notes: designs, reviews, audits, plans, handoffs -- pruned on ship or by age |
| `decisions/` | Decision records (context, alternatives, consequences) |
| `reviews/` | Periodic and milestone review summaries |
| `workflows/` | Reusable multi-stage process definitions |
| `scripts/` | Runnable developer scripts -- build, lint, profile, and exploratory one-offs |

Disposable scratch goes in the project-root `.tmp/` (gitignored), not under `dev/`.
Generated results, figures, and model outputs go in the project-root `output/`
(gitignored), not `dev/`. Shard `reviews/` into `<year>/` subfolders only if a
flat folder ever grows unwieldy.

## Development setup

Use R >= 4.1, Git, and the Quarto CLI. On Windows, install the matching Rtools when building
dependencies from source. Use a shell with this project's intended R installation; an active pixi
environment from another project can silently select a different R and library.

From the repository root, inspect the interpreter and library before installing or testing:

```powershell
Get-Command Rscript | Select-Object Source
Rscript --vanilla -e 'cat(R.version.string, "\n"); cat("R home:", R.home(), "\n"); print(.libPaths()); cat("cwd:", getwd(), "\n")'
quarto --version
```

If the selected R belongs to another project, open a shell outside that environment or explicitly
prepend the intended R `bin` directory to the process PATH. Do not install into the other project's
library. For example, use your installed R version in place of `<version>`:

```powershell
$env:PATH = "${env:ProgramFiles}\R\R-<version>\bin;$env:PATH"
```

On Windows, inherited `C.UTF-8` locale settings can cause R startup warnings. For check commands in
that shell, use the portable C locale without changing persistent user settings:

```powershell
$env:LC_ALL = 'C'
$env:LC_CTYPE = 'C'
$env:LANG = 'C'
```

Install development dependencies into the selected R's writable user library. The R package `quarto`
and the Quarto CLI are separate requirements; `DESCRIPTION` declares the package dependencies.

```r
install.packages(c("devtools", "lintr", "pkgdown", "rcmdcheck"))
devtools::install_deps(dependencies = TRUE, upgrade = "never")
devtools::install(build_vignettes = FALSE, upgrade = FALSE)
stopifnot(requireNamespace("quarto", quietly = TRUE))
```

Installing the current package matters for PSOCK workers and lint symbol resolution. Run the test,
lint, and check commands in `AGENTS.md`; use `document = FALSE` with `check_only()` when validating
without regenerating documentation. Package dependencies follow `DESCRIPTION` and the CI R matrix;
this repository does not promise an exactly locked research environment.

## Agent activation

R package development does not require Brain. For the curated agent setup, install the knowledge and
infrastructure checkouts first and set persistent `BRAIN_ROOT` and `BRAIN_INFRA_ROOT` environment
variables to their respective checkout roots. Follow the infrastructure checkout's
`docs/machine-setup.md` for Python/Conda prerequisites and user-runtime installation. Do not commit
machine-specific roots here.

After those prerequisites, these PowerShell commands activate and inspect only this project, even
when the `brain` shell command is unavailable:

```powershell
if (!(Test-Path "$env:BRAIN_ROOT/wiki/index.md")) { throw 'Check BRAIN_ROOT' }
if (!(Test-Path "$env:BRAIN_INFRA_ROOT/brain.py")) { throw 'Check BRAIN_INFRA_ROOT' }
python "$env:BRAIN_INFRA_ROOT/brain.py" refresh --agent-system --project "$PWD"
python "$env:BRAIN_INFRA_ROOT/brain.py" status --agent-system --project "$PWD" --detail
```

The tracked `.agent-system/agent-manifest.yml` owns the project selection. `.agents/`, `.claude/`,
`.codex/`, and `.github/agents/` are ignored links or generated projections; edit their canonical
sources instead. Start a fresh agent session after activation to verify discovery. A passing status
check establishes configuration compatibility, not execution by every runtime.

For the board CLI without an installed `todoboard` command, use the activated skill:

```powershell
$env:PYTHONPATH = (Resolve-Path '.agents/skills/todo-board').Path
python -m todoboard --help
```

This repository retains its existing R layout without Copier adoption. There is no
`.copier-answers.yml`, so `copier update` is not a supported maintenance command.

## Board

`dev/tasks/` is the source of truth for open work -- one note per item, its
frontmatter the row and its body the working draft. The view is the only
per-project variable:

- **In Obsidian** -- open `TODO.base` for the Queue / Active / Backlog / Blocked
  views.
- **Outside Obsidian (agents, CLI)** -- read the notes in `tasks/` directly,
  sorted by their `queue:` field, or run `todoboard render` to regenerate the
  read-only `dev/TODO.md` table. That table is generated output; never hand-edit
  it.

Both views show the same lean columns: `ID` · `Item` · `Area` · `#` · `State`.

## Working rules

- **Admit before you track.** Add an item when work spans sessions, is
  substantial enough that its scope, progress, or decisions need tracking, or
  needs visibility beyond the current session. Small work fully explained by its
  diff and Git history creates no board item.
- **One note per item** at `tasks/t<YYMMDDHHMM>[letter]-<slug>.md` -- the local
  creation minute, with `a`, `b`, … only for a same-minute collision. IDs are
  never reused or renamed. `todoboard add "<title>"` scaffolds one.
- **Body opens with a `> [!note] Overview`** -- What / Why / Effort, written for
  a reader who has forgotten the context: gloss every name, **What** = the change
  only, **Why** = what it costs to leave it undone, **Effort** = the shape of the
  work and its main unknown.
- **`## Progress` follows the Overview** -- a flat `- [ ]` / `- [x]` checklist,
  one line per step, no nesting, about seven steps at most. A ticked step names
  what it produced; a `blocked` item states what must clear in one italic line
  under the heading. It is evidence for a reader, not a state source: an
  all-ticked list never flips `status:` and never closes the item. Everything
  after it is free-form working draft, `## Refs` by convention for pointers.
- **State = the `status` field** -- `backlog` (default) · `active` · `blocked`.
  Completeness is not a frontmatter field; it lives in `## Progress`.
  Items sort by `queue` (lowest first, unranked last).
- **Long-form backing material lives in `drafts/`**, pointed at by the note's
  `doc:` field -- not inlined into the task note.
- **Closing = `todoboard done <id>`** -- appends one row to `LOG.md`, `git rm`s
  the note and its eligible tracked `doc:` draft, and compacts the queue. There
  is no "done" state on the board; Git history carries the detail.
- **Handoffs are self-contained.** A note handed to another session or runtime
  states objective, state, decisions, location, validation, next action, and
  blockers.
- **Record exact validation** -- the commands run and their outcomes.
- **Log shipped features in the root changelog** (`CHANGELOG.md`, or `NEWS.md`
  for R packages) -- feature-level entries only, linking `decisions/` for the
  detail. It is not a `dev/` file; it lives at the project root.

Create optional folders only when first needed.
