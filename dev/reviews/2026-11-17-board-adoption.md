# Persistent board adoption

Owner requested migration to the current persistent board, with open task IDs renumbered to the
sequential format. The primary is `master`; `.git-workflow.yml` now sets `trunk: master` and
`worktree_policy: concurrent`. Historical `dev/LOG.md` is unchanged.

## Renumbering

| New | Former ID | Task |
|---|---|---|
| t001 | t2608151255 | warm_filtering.R and evaluate_generator.R are god modules |
| t002 | t2608151255a | set.seed() discipline is held by convention, not construction |
| t003 | t2608151255c | Internal helpers use three competing naming conventions |
| t004 | t2611161814 | Journal figures for the perturbed climate |

Each note moved to `dev/tasks/<id>/<id>-<slug>.md` with `sequence:` and `former_id:` frontmatter.
References in t004's brief and findings, `dev/README.md` and the 2026-08-15 package review now use
the new IDs. t004's code is on `feat/optional-temperature-extrema` (`e4ae00d`, commit message uses
t004); its task record moved to `master`, where board edits belong.

## Reviewed checks

- Infrastructure's portable `test_adoption.py`, `test_persistent.py`, `test_projection.py` and
  `test_runtime.py`: **14 passed** in 90.6 seconds on 2026-11-17.
- Read-only `todoboard adopt` inventory on `master`: 4 live records, 32 historical identities,
  `next_sequence: 5`; every timestamp identity is recognized as legacy.

## Limits

Existing linked worktrees (`article-refresh`, `fresh-2026-11-05`, prunable `t2610241523c`) were not
reconciled; they predate adoption and keep their own checkouts. Fresh interactive-agent sandbox
behavior is not claimed.
