---
title: Unignore vignettes/.quarto so check_only() passes
type: todo-item
status: backlog
effort: 1
area: ci
origin: review-2026-08-15
queue: 1
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — Add a .Rbuildignore entry for vignettes/.quarto so the 18 cached _freeze artifacts stop shipping into the built package.
> **Why** — R CMD check is 0 errors and 0 notes but 2 warnings, and check_only() errors on warnings, so the documented release gate is currently red on master.
> **Effort** — One line, but confirm it is the whole story: the second warning names a missing `inst/doc`, which may be normal for a Quarto `VignetteBuilder` or may need its own fix.

## Progress

- [x] **The premise was wrong.** `^vignettes/\.quarto$` was already in `.Rbuildignore`. There are *two* freeze directories — `vignettes/.quarto/_freeze` (ignored) and `vignettes/_freeze` (not) — and the warning named the second. Added `^vignettes/_freeze$` with a comment saying both exist, since the near-identical names are exactly what made this look already-fixed.
- [x] Left `vignettes/_freeze/` tracked in Git. It is Quarto's freeze cache and committing it is what stops the vignettes re-executing the generator on every build; `vignettes/.gitignore` already excludes `.quarto/`. The problem was only that it shipped in the tarball, not that it is versioned.
- [x] Resolved the `inst/doc` warning: it was never a defect. `check_only()` defaults to `build_vignettes = FALSE`, which passed `--no-build-vignettes` but left the vignette *checks* running — so every quick check warned about output it had just been told not to produce. The documented release gate, `check_only(build_vignettes = TRUE)`, renders them and is **0 errors / 0 warnings / 0 notes**.
- [x] Made the quick form honest too: it now also passes `--ignore-vignettes`, so it skips vignettes on both the build and the check side. `check_only()` is now **0/0/0** as well. Two unavoidable warnings on every quick check are worse than none — they train you to read past the warning count, which is how a real one slips through. That happened to me earlier in this sweep: a genuine new `withr` warning was only caught because I had recorded the expected count beforehand.

## Verification

- `check_only()` — 0 errors, 0 warnings, 0 notes (1m13s).
- `check_only(build_vignettes = TRUE)` — 0 errors, 0 warnings, 0 notes (1m48s),
  vignettes rebuilt OK. The release gate is green on master for the first time
  this sweep.
- No `NEWS.md` entry: `.Rbuildignore` and `tools/` are packaging and repo-only
  tooling, not user-visible behaviour.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § B11.
- Full check result: 0 errors, 0 notes, 2 warnings, 3m34s.
- `AGENTS.md` names `check_only(build_vignettes = TRUE)` as the release gate,
  and the wrapper errors on warnings — so this is currently red on master.
