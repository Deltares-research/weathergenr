---
title: warm_filtering.R and evaluate_generator.R are god modules
type: todo-item
status: backlog
effort: 1
area: architecture
origin: review-2026-08-15
queue:
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What changes** — 1903 and 1836 lines, 14 and 17 functions, together 28 percent of the package; evaluate_weather_generator() alone spans lines 102-360.
> **Why it matters** — The *_plots.R separation convention is held consistently but is not enough on its own at this size. Not worth a speculative refactor; worth knowing before the next substantial change lands in either file.
> **What it takes** — Nothing for now. Split them when a change to either file needs more than a small, local edit.

## Progress

- [ ] Check whether the condition in What it takes has happened

## Refs

- `dev/drafts/package-review-2026-08-15.md` § A2.
