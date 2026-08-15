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

## Refs

- `dev/drafts/package-review-2026-08-15.md` § B6.
- Related open item `t2608151254k` (trim the export surface) — both are
  breaking-change decisions that want the same release window.
