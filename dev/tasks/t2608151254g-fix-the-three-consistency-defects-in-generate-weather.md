---
title: Fix the three consistency defects in generate_weather()
type: todo-item
status: backlog
effort: 1
area: api
origin: review-2026-08-15
queue: 9
created: 2026-08-15
updated: 2026-08-15
---

> [!note] Overview
> **What** — Seed arithmetic can overflow to NA, save_plots = FALSE still builds the filter plots, and one .log() call drops verbose.
> **Why** — The overflow fails hard under a specific seed deep into a long run, and the other two make documented switches quietly not mean what they say.
> **Effort** — Small and mechanical, three independent one-liners in one file. The only care needed is that changing how seeds are derived shifts the RNG stream and therefore numeric output.

## Progress

- [ ] Fix the seed overflow: `warm_seed`/`daily_seed` come from `sample.int(.Machine$integer.max)` (`R/generator.R:261-262`) and are then incremented — `warm_seed + 1L` (`:434`), `daily_seed + n` (`:508`, `:536`) — which can reach `NA`. Same pattern at `R/wavelet_warm.R:400` (`seed + k * 1000L`). Draw from a smaller range, or take `%%` before adding.
- [ ] Baseline gate around the seed change — it shifts the RNG stream even though it is not a method change.
- [ ] Pass `make_plots = save_plots` to `filter_warm_pool()` (`:437`) so `save_plots = FALSE` stops building three panels it will not write.
- [ ] Add the missing `verbose = verbose` to the `.log()` call at `:457`; its sibling at `:401` has it.
- [ ] Consider dropping the dead self-recursion guard at `:753` — `identical(generate_weather, run_weather_generator)` cannot fire in a namespaced package.

## Refs

- `dev/drafts/package-review-2026-08-15.md` § B3, § B4, § B5.
