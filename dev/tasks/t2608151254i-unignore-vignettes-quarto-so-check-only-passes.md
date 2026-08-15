---
title: Unignore vignettes/.quarto so check_only() passes
type: todo-item
status: backlog
effort: 1
area: ci
origin: review-2026-08-15
queue: 10
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — Add a .Rbuildignore entry for vignettes/.quarto so the 18 cached _freeze artifacts stop shipping into the built package.
> **Why** — R CMD check is 0 errors and 0 notes but 2 warnings, and check_only() errors on warnings, so the documented release gate is currently red on master.
> **Effort** — One line, but confirm it is the whole story: the second warning names a missing `inst/doc`, which may be normal for a Quarto `VignetteBuilder` or may need its own fix.

## Progress

- [ ] Add `^vignettes/\.quarto$` to `.Rbuildignore`.
- [ ] Re-run `source("tools/build_site_tools.R"); check_only()` and confirm it now exits 0.
- [ ] If the "Directory 'inst/doc' does not exist" warning survives, decide whether that is expected under `VignetteBuilder: quarto` or needs `check_only(build_vignettes = TRUE)` to be the gate instead.
- [ ] Check whether `vignettes/.quarto/` should also be in `.gitignore` — 18 `_freeze` artifacts are currently tracked.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § B11.
- Full check result: 0 errors, 0 notes, 2 warnings, 3m34s.
- `AGENTS.md` names `check_only(build_vignettes = TRUE)` as the release gate,
  and the wrapper errors on warnings — so this is currently red on master.
