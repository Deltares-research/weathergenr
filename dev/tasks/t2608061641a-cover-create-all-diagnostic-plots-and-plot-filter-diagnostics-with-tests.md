---
title: Cover create_all_diagnostic_plots and plot_filter_diagnostics with tests
type: todo-item
status: backlog
effort: 2
area: tests
queue:
created: 2026-08-06
updated: 2026-08-06
---

> [!note] Overview
> **What** — Add tests for the two exported plot functions left uncovered: create_all_diagnostic_plots (needs plot_data + plot_config + variables) and plot_filter_diagnostics (10 arguments).
> **Why** — Both are exported and untested; a fake fixture built only to reach expect_s3_class would assert nothing, so they need real fixtures.
> **Effort** — large

## Progress

- [ ] <first step>
