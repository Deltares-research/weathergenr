---
title: Trim the export surface
type: todo-item
status: backlog
effort: 2
area: api
origin: review-2026-08-15
queue: 4
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — NAMESPACE exports criteria_string_compact, generate_symmetric_dummy_points, get_result_index and match_transition_positions, which have no plausible external audience; several carry @keywords internal and are exported anyway.
> **Why** — blueearth_cst pins this package by Git tag, so every export is a contract that cannot be quietly withdrawn. Deciding now is cheaper than after another consumer appears.
> **Effort** — Large not in code but in coordination: unexporting is breaking, the downstream consumer pins by Git tag, and the release has to be sequenced so nothing is pulled from under it.

## Progress

- [ ] Grep `blueearth_cst` for every candidate before touching anything — an export nobody calls is a different decision from one that is called.
- [ ] Confirm the shortlist: `criteria_string_compact`, `generate_symmetric_dummy_points`, `get_result_index`, `match_transition_positions`. The wavelet utilities (`fill_nearest`, `gws_regrid`, `extract_signif_curve`) are general-purpose enough to keep.
- [ ] Reconcile the `@keywords internal`-but-exported cases: either drop `@export` or drop the keyword. Both is incoherent.
- [ ] Decide the release path — deprecate for one version or drop at the next minor. `AGENTS.md` notes releases are only delivered once the tag is pushed, so sequence accordingly.
- [ ] `devtools::document()`, then `check_only()` since this touches exports and documentation.
- [ ] `NEWS.md` entry — removing an export is user-visible.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § A4.
- Related watch-item `t2608151255c` (three competing internal-naming
  conventions) — settling that would make the export list self-evident.
