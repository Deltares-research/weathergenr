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

- [x] Decide the direction — **make the reader CF-correct**. The writer's encoding was already right: the data are 365-day, so `0:(nt-1)` plus `calendar=noleap` is correct CF, and xarray/Wflow already decoded it that way. Dropping the attribute would have made the files wrong for every external consumer. Confirmed no downstream dependency on the old decode: the packaged fixture is `proleptic_gregorian` and the existing round-trip test writes `proleptic_gregorian`, so neither exercised the broken path.
- [x] Implement the decode branch — added `.cf_calendar_kind()`, `.noleap_ymd_to_index()`, `.noleap_index_to_ymd()`, `.noleap_offset_to_date()`, `.noleap_date_to_offset()` and `.cf_time_to_date()` as file-level internals in `R/io_netcdf.R`. `.parse_time_to_date()` is now a thin wrapper that reads the `units` and `calendar` attributes and delegates.
- [x] Give `write_netcdf()` the series dates — added an optional `dates` argument. When supplied it validates length, `origin_date == dates[1]`, and strict monotonicity, and computes the time coordinate on the declared calendar via `.noleap_date_to_offset()` instead of assuming contiguity. Default `NULL` keeps the old behaviour, so existing callers are unaffected.
- [x] Validate the `calendar` argument — `write_netcdf()` now rejects anything outside the writable set. Also removed the `try(..., silent = TRUE)` around the calendar `ncatt_put`: losing that attribute silently is the whole bug.
- [x] Fix the roxygen example — `origin_date` is now documented as *the date of the first time step*, with the `"1970-01-01"` trap named explicitly.
- [x] Add a round-trip test — 5 new `test_that()` blocks in `tests/testthat/test-io_netcdf.R` covering calendar classification, index round-trip (including negative offsets and fractional truncation), decode divergence at the leap day, the unrepresentable-calendar errors, the noleap write/read round trip, and the `origin_date`/`dates` disagreement.
- [x] Removed a second, independent parse of the time axis. `leap_keep_idx` was rebuilt by re-calling `.parse_time_to_date()`, so the date vector and the Feb-29 row mask could have disagreed about the calendar while still matching on row count — nothing downstream would have caught it. Both now derive from one decode.

## Verification

- `devtools::test(filter = "io_netcdf")` — 91 pass, 0 fail. Includes the
  pre-existing `proleptic_gregorian` round-trip test, which is the evidence the
  standard branch is untouched.
- `devtools::test()` — 939 pass, 0 fail, 4 skips (the opt-in baseline).
- Baseline gate with `WEATHERGENR_BASELINE=1` — 8 pass, 0 fail. No numeric
  output moved, as expected given the fixture is Gregorian, but confirmed
  rather than assumed since this touches calendar logic.
- `Rscript tools/lint.R --changed` — no new findings; the single `nt` warning in
  `io_netcdf.R` predates this change and is logged on `t2608061641`.
- Direct reproduction of the review's case: 730 noleap days written from
  2020-01-01 now read back as 2020-01-01..2021-12-31 (was ..2021-12-30). A
  100-year encode/decode round trip is exact; the drift the old decode carried
  measures 25 days.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § C5 — the verified round trip
  (730 days written as 2020-01-01..2021-12-31 read back as ..2021-12-30).
- `blueearth_cst` consumes these files; see the tag-pin note in `AGENTS.md`.
