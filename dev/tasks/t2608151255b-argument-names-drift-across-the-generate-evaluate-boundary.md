---
title: Argument names drift across the generate/evaluate boundary
type: watch-item
area: api
origin: review-2026-08-15
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — out_dir versus output_dir, warm_filter_bounds becoming filter_bounds, relax_priority becoming relax_order. run_weather_generator() absorbs the renames; anyone calling the pieces directly hits them.
> **Why** — Renaming is a breaking change and blueearth_cst pins by Git tag, so it cannot be done casually -- but the drift should not grow.
> **Trigger** — A major version is cut, or a new argument is added across the same boundary.

## Trigger fired 2026-08-17 — undecided

A major version **was** cut: `efc071d` bumps to 2.0.0 for the export trim
(`t2608151254k`). The renames were not made, and this stays a watch item only
because nobody has decided; the trigger is spent, not unmet.

The window is still open. 2.0.0 is committed but **untagged and unpushed**, and
`blueearth_cst` pins by tag, so nothing downstream has seen it. Renaming now
costs the same one break already being paid. Once the tag is pushed, the next
chance is 3.0.0.

Decide before tagging: rename and fold into 2.0.0, or accept the drift and
re-arm this for the next major.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § B6.
- Closed item `t2608151254k` (trim the export surface) — the break this would
  have shared a release window with.
