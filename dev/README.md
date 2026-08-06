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
