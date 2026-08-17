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

## Resolved 2026-08-17 in 2.0.0 — two renamed, one accepted

Trigger fired when `efc071d` cut 2.0.0 for the export trim. Acted on before
tagging, so the break is absorbed into a release nothing downstream had seen.

**Checking the consumer changed the answer.** `blueearth_cst`'s
`config/defaults/weathergen_config.yml` is an explicit name-mirroring contract —
"each section is a weathergenr FUNCTION name and each key is one of that
function's ARGUMENT names, so any key here can be checked against `?<section>`".
`warm_filter_bounds: []` is live in it and `out_dir` is injected into the same
section by rule 3.10. So the three pairs are not symmetric, and the
low-blast-radius direction is to move the *inner* functions onto the *outer*
names rather than the reverse the note implied.

Renamed (`84d7ebd`), old names deprecated but working:

- `evaluate_weather_generator(output_dir=)` and
  `create_all_diagnostic_plots(output_dir=)` → `out_dir`. The real defect:
  one concept, two names, two exported functions. `out_dir` already won 2-1,
  and the consumer never calls the evaluator directly. 88 references.
- `generate_weather(relax_priority=)` → `relax_order`, matching
  `filter_warm_pool()` which it forwards to. Also `config$relax_order`.
  Pure synonym drift, neither name live downstream; `relax_order` chosen
  because the argument is a sequence, not a ranking.

**Accepted, not fixed:** `generate_weather(warm_filter_bounds=)` versus
`filter_warm_pool(filter_bounds=)`. The `warm_` prefix is doing real work —
`warm_var`, `warm_signif`, `warm_pool_size`, `warm_filter_bounds` are a family,
and breaking it to match one niche function costs more than the mismatch does.
It is also the single name live in the downstream config. Recorded here so a
future reader finds a decision rather than an oversight.

Mechanism is a hand-rolled `.resolve_renamed_arg()` in `utils_internal.R`, not
`lifecycle` — one renamed argument does not justify an Import. It resolves on
suppliedness rather than value, and errors when both names are supplied.

## Found while doing this — the downstream note is stale

`blueearth_cst`'s config says `relax_priority` is deliberately absent because
"that wrapper forwards every generate_weather argument EXCEPT this one — so a
key here would read as a setting and reach nothing. Restore it if upstream
forwards it." `R/generator.R` **does** forward it (and now forwards
`relax_order`). The restore condition is met, so a disabled stress-test knob is
available downstream and the config comment is wrong. Not actioned here — that
is a `blueearth_cst` change, and it will want the new spelling.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § B6.
- Closed item `t2608151254k` (trim the export surface) — the break this shared
  its release window with.
