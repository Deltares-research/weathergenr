---
title: Fix the NetCDF noleap calendar round trip
type: todo-item
status: backlog
effort: 1
area: io
origin: review-2026-08-15
queue: 1
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — read_netcdf() ignores the calendar attribute and decodes proleptic-Gregorian, so a file write_netcdf() stamped calendar=noleap reads back one day short per leap year. Also give write_netcdf() a way to verify origin_date is the series start, and validate the calendar argument.
> **Why** — Verified round trip loses ~24 days over a century run, and the mismatch crosses the blueearth_cst boundary where CF-aware readers disagree with weathergenr's own.
> **Effort** — Small in code; the open question is which side is wrong. Making `read_netcdf()` CF-correct is a behaviour change for anyone already compensating for it.

## Progress

- [ ] Decide the direction: teach `.parse_time_to_date()` to honour a `noleap`/`365_day` calendar attribute (CF-correct), or drop the attribute from `write_netcdf()`. CF-correct is the right answer unless a downstream consumer depends on the current decode.
- [ ] Implement the decode branch in `.parse_time_to_date()` (`R/io_netcdf.R:116-150`) — build the date axis by counting 365-day years rather than adding days to a POSIXct origin.
- [ ] Give `write_netcdf()` the series dates (or an explicit `n_years`/start check) so it can verify `origin_date` is the series start; `vals = 0:(nt-1)` is only correct when it is (`:546-553`).
- [ ] Validate the `calendar` argument against a known set instead of accepting free text (`:409`, `:650`).
- [ ] Fix the roxygen example, which suggests `origin_date = "1970-01-01"` (`:391`) — against `0:(nt-1)` that silently relocates a 2020 run to 1970.
- [ ] Add a round-trip test: write N noleap days, read back, assert the date axis matches what was written.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C5 — the verified round trip
  (730 days written as 2020-01-01..2021-12-31 read back as ..2021-12-30).
- `blueearth_cst` consumes these files; see the tag-pin note in `AGENTS.md`.
