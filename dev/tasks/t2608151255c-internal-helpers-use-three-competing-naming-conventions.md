---
title: Internal helpers use three competing naming conventions
type: watch-item
area: architecture
origin: review-2026-08-15
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — Dot-prefixed (.markov_month_probs, .date_parts), bare-unexported (fit_monthly_distributions, find_local_maxima, log_filter_iteration), and @keywords internal but exported (criteria_string_compact, get_result_index). All three appear within single files.
> **Why** — Costs nothing today, but it means neither the prefix nor the export list reliably signals whether a function is internal.
> **Trigger** — A convention is chosen and written into the AGENTS.md Conventions section.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § A5.
- Related open item `t2608151254k` (trim the export surface) — settling the
  convention would make the export list self-evident.
