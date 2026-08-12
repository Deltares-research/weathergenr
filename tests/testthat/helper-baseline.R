# ==============================================================================
# End-to-end baseline: shared configuration and fingerprinting
# ==============================================================================
#
# This helper is the SINGLE definition of what the baseline covers and how it is
# fingerprinted. It is loaded automatically by testthat and is also sourced by
# tools/record_baseline.R, so the recorder and the checker can never drift apart.
#
# Depends only on the package's existing Imports (rlang) and base R. No new
# dependency is introduced.
#
# See tests/testthat/test-baseline-e2e.R for the gate and
# tools/record_baseline.R for how to re-record.

# ------------------------------------------------------------------------------
# Scenarios
# ------------------------------------------------------------------------------

#' Baseline scenarios
#'
#' Two scenarios, deliberately chosen: the water-year and calendar-year paths
#' diverge in `compute_water_year()`, `.filter_grid_years()` and
#' `estimate_monthly_markov_probs()`, so a single scenario would leave half of
#' each of those branches unguarded.
#'
#' `n_years = 20` is bounded below by the wavelet minimum: `filter_warm_pool()`
#' analyses each *simulated* series, whose length is `n_years`, and
#' `analyze_wavelet_spectrum()` rejects any series shorter than 16
#' (`wavelet_cwt.R:437`). The observed record is ~30 years here and is not the
#' binding constraint. Do not lower `n_years` below 16.
#'
#' `warm_pool_size` is small because WARM is under 1% of runtime -- shrinking it
#' does not make the baseline meaningfully faster and only weakens the pool.
#' Runtime scales with `n_years` x `n_realizations`, not with the pool.
#' @noRd
baseline_scenarios <- function() {
  base <- list(
    vars             = c("precip", "temp"),
    n_years          = 20L,
    start_year       = 2020L,
    n_realizations   = 2L,
    warm_var         = "precip",
    warm_signif      = 0.80,
    warm_pool_size   = 300L,
    annual_knn_n     = 100L,
    wet_q            = 0.20,
    extreme_q        = 0.80,
    seed             = 100L,
    eval_max_grids   = 5L,
    eval_seed        = 42L
  )

  list(
    water_year    = utils::modifyList(base, list(year_start_month = 10L)),
    calendar_year = utils::modifyList(base, list(year_start_month = 1L))
  )
}

#' Path to the fixture NetCDF used by every baseline scenario
#' @noRd
baseline_fixture_path <- function() {
  system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
}

#' Path to the stored baseline artifact
#' @noRd
baseline_path <- function() {
  testthat::test_path("baseline", "e2e-baseline.csv")
}

# ------------------------------------------------------------------------------
# Running a scenario
# ------------------------------------------------------------------------------

#' Run one baseline scenario end to end
#'
#' Returns the raw objects rather than a fingerprint, so the recorder and the
#' checker fingerprint identical inputs.
#'
#' `save_plots = FALSE` throughout: plot rendering is roughly a quarter of the
#' pipeline's runtime and contributes nothing to the numeric fingerprint.
#' @noRd
baseline_run_scenario <- function(cfg, ncdata, out_dir = tempfile("wgr-baseline-"),
                                  parallel = FALSE, n_cores = NULL) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # markov_next_state() latches a global option on an invalid state; isolate it
  # so a preceding test cannot change this run's warning behaviour.
  old_opt <- options(weathergenr.markov_invalid_state_warned = NULL)
  on.exit(options(old_opt), add = TRUE)

  gen <- generate_weather(
    obs_data         = ncdata$data,
    obs_grid         = ncdata$grid,
    obs_dates        = ncdata$date,
    vars             = cfg$vars,
    n_years          = cfg$n_years,
    start_year       = cfg$start_year,
    year_start_month = cfg$year_start_month,
    n_realizations   = cfg$n_realizations,
    warm_var         = cfg$warm_var,
    warm_signif      = cfg$warm_signif,
    warm_pool_size   = cfg$warm_pool_size,
    annual_knn_n     = cfg$annual_knn_n,
    wet_q            = cfg$wet_q,
    extreme_q        = cfg$extreme_q,
    out_dir          = out_dir,
    seed             = cfg$seed,
    parallel         = parallel,
    n_cores          = n_cores,
    verbose          = FALSE,
    save_plots       = FALSE
  )

  eval_data <- prepare_evaluation_data(
    gen_output       = gen,
    obs_data         = ncdata$data,
    obs_dates        = ncdata$date,
    grid_ids         = ncdata$grid$id,
    variables        = cfg$vars,
    year_start_month = cfg$year_start_month,
    verbose          = FALSE
  )

  assessment <- evaluate_weather_generator(
    daily_sim        = eval_data$sim_data,
    daily_obs        = eval_data$obs_data,
    vars             = cfg$vars,
    variable_labels  = NULL,
    n_realizations   = cfg$n_realizations,
    year_start_month = cfg$year_start_month,
    eval_max_grids   = cfg$eval_max_grids,
    wet_q            = cfg$wet_q,
    extreme_q        = cfg$extreme_q,
    output_dir       = out_dir,
    save_plots       = FALSE,
    seed             = cfg$eval_seed
  )

  list(gen = gen, fit_summary = attr(assessment, "fit_summary"))
}

# ------------------------------------------------------------------------------
# Fingerprinting
# ------------------------------------------------------------------------------

# Tier-2 statistics are rounded before they are recorded, so a difference in the
# last bits of a double does not read as a regression. Eight significant digits
# is far tighter than any real change to the resampling logic would produce.
BASELINE_SIGNIF <- 8L

#' Format a numeric value for stable text storage
#' @noRd
.baseline_num <- function(x) {
  if (!length(x) || all(is.na(x))) return(NA_character_)
  formatC(signif(as.numeric(x), BASELINE_SIGNIF),
          format = "e", digits = BASELINE_SIGNIF - 1L)
}

#' Fingerprint one scenario's results
#'
#' Two tiers:
#'   * `hash_*` keys are exact -- the tripwire. Any change to the resampled
#'     dates or the fit metrics moves them.
#'   * every other key is a rounded summary statistic whose only job is to say
#'     *how* the output moved once the tripwire fires.
#'
#' `evaluate_weather_generator()` returns ggplot objects with `fit_summary` as an
#' attribute, and its metadata embeds `Sys.Date()`. Only `fit_summary` is
#' fingerprinted: the plots are ggplot2-version-sensitive and the date is not
#' reproducible.
#' @noRd
baseline_fingerprint <- function(scenario, res, cfg = NULL) {

  rows <- list()
  add <- function(key, value) {
    rows[[length(rows) + 1L]] <<- data.frame(
      scenario = scenario, key = key, value = as.character(value),
      stringsAsFactors = FALSE
    )
  }

  # --- scenario configuration ------------------------------------------------
  # Recorded AND compared. A fingerprint is meaningless without the settings
  # that produced it, so a change to baseline_scenarios() must fail the gate
  # rather than silently compare a different experiment against old numbers.
  if (!is.null(cfg)) {
    for (nm in names(cfg)) {
      add(paste0("config.", nm), paste(as.character(cfg[[nm]]), collapse = ","))
    }
  }

  # --- generate_weather() ----------------------------------------------------
  resampled <- res$gen$resampled
  # Compare the underlying integer date codes, not the printed representation:
  # Date formatting is locale- and timezone-sensitive, the integers are not.
  resampled_int <- lapply(resampled, function(x) as.integer(unclass(x)))
  names(resampled_int) <- names(resampled)

  add("gen.hash_resampled", rlang::hash(resampled_int))
  add("gen.hash_sim_dates", rlang::hash(as.integer(unclass(res$gen$dates))))
  add("gen.n_days", length(res$gen$dates))
  add("gen.n_realizations", length(resampled_int))

  for (nm in names(resampled_int)) {
    v <- resampled_int[[nm]]
    add(paste0("gen.", nm, ".n_distinct_sources"), length(unique(v)))
    add(paste0("gen.", nm, ".mean_source_date"), .baseline_num(mean(v)))
    add(paste0("gen.", nm, ".sd_source_date"), .baseline_num(stats::sd(v)))
    add(paste0("gen.", nm, ".n_na"), sum(is.na(v)))
  }

  # --- evaluate_weather_generator() -----------------------------------------
  fit <- res$fit_summary

  if (is.null(fit) || !nrow(fit)) {
    add("eval.fit_summary_present", "FALSE")
    return(do.call(rbind, rows))
  }

  fit <- as.data.frame(fit, stringsAsFactors = FALSE)
  fit <- fit[order(fit$rlz), , drop = FALSE]

  # Hash the rounded frame, not the raw one, so the exact tier carries the same
  # precision contract as the summary tier.
  fit_rounded <- fit
  num_cols <- names(fit)[vapply(fit, is.numeric, logical(1))]
  for (nm in num_cols) fit_rounded[[nm]] <- signif(fit_rounded[[nm]], BASELINE_SIGNIF)

  add("eval.hash_fit_summary", rlang::hash(fit_rounded))
  add("eval.fit_n_rows", nrow(fit))
  add("eval.fit_columns", paste(names(fit), collapse = "|"))

  for (i in seq_len(nrow(fit))) {
    rlz <- fit$rlz[i]
    for (nm in num_cols) {
      if (identical(nm, "rlz")) next
      add(paste0("eval.rlz", rlz, ".", nm), .baseline_num(fit[[nm]][i]))
    }
  }

  do.call(rbind, rows)
}

# ------------------------------------------------------------------------------
# Artifact I/O
# ------------------------------------------------------------------------------

#' Do PSOCK workers run the same build as the master?
#'
#' `devtools::load_all()` does not propagate to PSOCK workers: they resolve the
#' package from the installed library instead. A parallel test run under
#' `devtools::test()` therefore exercises whatever is installed, not the working
#' tree, and would report a pass or a failure that says nothing about the code
#' being edited. `R CMD check` installs first, so there the two agree.
#'
#' Compares the deparsed body of `resample_weather_dates()` on a worker against
#' the master's. Returns TRUE only when they match.
#' @noRd
baseline_workers_match_master <- function(n_cores = 2L) {

  master_body <- paste(deparse(body(resample_weather_dates)), collapse = "\n")

  cl <- parallel::makeCluster(n_cores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  worker_body <- try(
    parallel::clusterEvalQ(cl, {
      fn <- try(get("resample_weather_dates",
                    envir = asNamespace("weathergenr")), silent = TRUE)
      if (inherits(fn, "try-error")) NA_character_
      else paste(deparse(body(fn)), collapse = "\n")
    })[[1]],
    silent = TRUE
  )

  !inherits(worker_body, "try-error") &&
    !is.na(worker_body) &&
    identical(worker_body, master_body)
}

#' Read the stored baseline, or NULL when it does not exist
#' @noRd
baseline_read <- function(path = baseline_path()) {
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, colClasses = "character", stringsAsFactors = FALSE)
}

#' Split a baseline frame into provenance rows and comparable fingerprint rows
#'
#' Provenance is recorded for attribution and is deliberately NOT compared:
#' a baseline recorded on a different commit is still a valid baseline.
#' @noRd
baseline_split <- function(df) {
  is_prov <- df$scenario == "provenance"
  list(provenance = df[is_prov, , drop = FALSE],
       fingerprint = df[!is_prov, , drop = FALSE])
}
