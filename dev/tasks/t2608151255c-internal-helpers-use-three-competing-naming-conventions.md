---
title: Internal helpers use three competing naming conventions
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
> **What changes** — Dot-prefixed (.markov_month_probs, .date_parts), bare-unexported (fit_monthly_distributions, find_local_maxima, log_filter_iteration), and @keywords internal but exported (criteria_string_compact, get_result_index). All three appear within single files.
> **Why it matters** — Costs nothing today, but it means neither the prefix nor the export list reliably signals whether a function is internal.
> **What it takes** — Waits on choosing one naming convention and writing it into the Conventions section of AGENTS.md.

## Progress

- [ ] Check whether the condition in What it takes has happened

## Refs

- `dev/drafts/package-review-2026-08-15.md` § A5.
- Related open item `t2608151254k` (trim the export surface) — settling the
  convention would make the export list self-evident.
