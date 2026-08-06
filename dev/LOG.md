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
