---
title: LOG
type: log
---

# LOG — closed work

Terminal ledger for items closed off the board — one line per closed item,
newest first. A browsable index *over* Git history, not a copy of it: full task
bodies stay recoverable with `git log --diff-filter=D -- dev/tasks/`.

Rows are appended by `todoboard done`, never hand-edited.

| Closed     | ID           | Item                                     | Area         |
| ---------- | ------------ | ---------------------------------------- | ------------ |
| 2026-08-15 | t2608151253c | Add equivalence tests for the two forked numeric routines | tests        |
| 2026-08-15 | t2608151253b | Reconcile warm_signif with the documented WARM component selection | warm         |
| 2026-08-15 | t2608151253a | Fix the mean pinning and sample-sd bias in .variance_match_matrix() | warm         |
| 2026-08-15 | t2608151253  | Fix the NetCDF noleap calendar round trip | io           |
