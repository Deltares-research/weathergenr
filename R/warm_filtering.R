#' Filter and sample WARM realizations using distributional, tail, and spectral criteria
#'
#' @description
#' Filters an ensemble of annual WARM realizations against an observed annual
#' reference series. Filtering is based on three families of criteria:
#' distributional moments, tail mass behavior, and spectral similarity from the
#' global wavelet spectrum.
#'
#' The function computes pass or fail vectors for each filter family and builds
#' a candidate pool. If the pool is smaller than the requested sample size, the
#' function relaxes thresholds iteratively until enough candidates are found or
#' the maximum number of relaxation iterations is reached.
#'
#' Relaxation is adaptive. At each iteration, the currently most restrictive
#' active filter is relaxed, defined as the filter with the lowest pass rate
#' among the active filters. This is constrained by the relaxation order given
#' in \code{relax_order}, which sets which filters are eligible to be relaxed
#' and how wavelet relaxation is parameterized.
#'
#' If the pool is still smaller than \code{n_select} after relaxation, a
#' deterministic fallback selects the \code{n_select} realizations with the
#' smallest absolute relative mean difference.
#'
#' @param obs_series Numeric vector. Observed annual series used as the
#'   reference target.
#' @param sim_series Numeric matrix. Simulated annual realizations with years
#'   in rows and realizations in columns.
#' @param n_select Integer scalar. Number of realizations to return in the
#'   \code{selected} element.
#' @param seed Integer scalar or NULL. If provided, sets the random seed used
#'   for selecting alignment windows when \code{obs_series} and \code{sim_series}
#'   have different lengths, and for sampling within the final candidate pool.
#' @param relax_order Character vector. Relaxation priority ordering for the
#'   filter families. Must contain exactly: \code{"mean"}, \code{"sd"},
#'   \code{"tail_low"}, \code{"tail_high"}, \code{"wavelet"}.
#' @param filter_bounds Named list. Overrides for filtering thresholds and
#'   relaxation controls. Keys must be snake case and match those returned by
#'   \code{filter_warm_bounds_defaults()}.
#' @param wavelet_args Named list. Parameters passed to
#'   \code{analyze_wavelet_spectrum()} for observed and simulated series.
#'   Expected entries include \code{signif_level}, \code{noise_type},
#'   \code{period_lower_limit}, and \code{detrend}.
#' @param modwt_n_levels Integer or NULL. Number of MODWT levels used in WARM.
#'   This value is used only for diagnostics. If NULL, a value is estimated from
#'   the series length.
#' @param make_plots Logical scalar. If TRUE, compute diagnostic plots for the
#'   selected realizations and return them in \code{plots}.
#' @param parallel Logical scalar. If TRUE, allows spectral metrics to be
#'   computed in parallel inside \code{compute_spectral_metrics()}.
#' @param n_cores Integer scalar or NULL. Number of worker processes to use for
#'   parallel execution. If NULL, an internal default is used.
#' @param cache_gws Logical scalar. If TRUE, cache simulated global wavelet
#'   spectra for diagnostics; this can increase memory use.
#' @param verbose Logical scalar. If TRUE, logs setup information, per iteration
#'   pass rates, and relaxation actions.
#'
#' @return Named list with the following elements:
#' \describe{
#'   \item{pool}{Numeric matrix. Candidate pool of realizations that passed the
#'   final set of filters. Subset of columns from \code{sim_series}.}
#'   \item{selected}{Numeric matrix. The \code{n_select} realizations chosen from
#'   the pool. Subset of columns from \code{sim_series}.}
#'   \item{summary}{Data frame. Pass counts and pass rates for each filter
#'   family, plus the selection mode used.}
#'   \item{diagnostics}{List. Window metadata, relaxation log, spectral
#'   diagnostics, spectral metrics, and the final bounds used. When
#'   \code{cache_gws = TRUE}, includes \code{gws_cache}.}
#'   \item{plots}{NULL or a named list of ggplot objects when
#'   \code{make_plots = TRUE}.}
#' }
#'
#' @importFrom utils modifyList
#' @importFrom stats acf median runif setNames cor var quantile
#' @export
filter_warm_pool <- function(
    obs_series = NULL,
    sim_series = NULL,
    n_select = 5,
    seed = NULL,
    relax_order = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
    filter_bounds = list(),
    wavelet_args = list(
      signif_level = 0.80,
      noise_type = "red",
      period_lower_limit = 2,
      detrend = TRUE
    ),
    modwt_n_levels = NULL,
    make_plots = FALSE,
    parallel = FALSE,
    n_cores = NULL,
    cache_gws = FALSE,
    verbose = FALSE
) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (is.null(obs_series) || !is.numeric(obs_series) || !is.vector(obs_series)) {
    stop("'obs_series' must be a numeric vector.", call. = FALSE)
  }
  if (is.null(sim_series) || !is.matrix(sim_series) || !is.numeric(sim_series)) {
    stop("'sim_series' must be a numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(n_select) || length(n_select) != 1L || !is.finite(n_select) || n_select < 1) {
    stop("'n_select' must be a positive integer.", call. = FALSE)
  }
  n_select <- as.integer(n_select)

  make_plots <- isTRUE(make_plots)
  verbose <- isTRUE(verbose)

  parallel <- isTRUE(parallel)
  cache_gws <- isTRUE(cache_gws)
  if (!is.null(n_cores)) {
    if (!is.numeric(n_cores) || length(n_cores) != 1L || !is.finite(n_cores) || n_cores < 1) {
      stop("'n_cores' must be a positive integer or NULL.", call. = FALSE)
    }
    n_cores <- as.integer(n_cores)
  }

  n_obs0 <- length(obs_series)
  n_sim0 <- nrow(sim_series)
  n_rlz  <- ncol(sim_series)

  if (n_rlz < n_select) {
    warning("n_select exceeds number of series; returning at most ncol(sim_series).", call. = FALSE)
    n_select <- n_rlz
  }

  # ---------------------------------------------------------------------------
  # Relaxation order
  # ---------------------------------------------------------------------------
  allowed <- c("mean", "sd", "tail_low", "tail_high", "wavelet")
  if (is.null(relax_order) || !is.character(relax_order)) {
    stop("'relax_order' must be a character vector.", call. = FALSE)
  }
  relax_order <- as.character(relax_order)

  if (anyNA(relax_order) || any(!nzchar(relax_order))) {
    stop("'relax_order' must not contain NA/empty strings.", call. = FALSE)
  }
  if (any(duplicated(relax_order))) {
    stop("'relax_order' must not contain duplicates.", call. = FALSE)
  }
  if (!all(relax_order %in% allowed) || !setequal(relax_order, allowed)) {
    stop(
      "'relax_order' must contain exactly these filters once each: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }

  RELAX_ORDER <- relax_order

  # ---------------------------------------------------------------------------
  # Bounds defaults (snake_case keys)
  # ---------------------------------------------------------------------------
  b_list <- modifyList(filter_warm_bounds_defaults(), filter_bounds)
  b <- list2env(b_list, parent = environment())
  b$relax_max_iter <- as.integer(b$relax_max_iter)

  # Validate tail parameters (snake_case)
  b$tail_low_p   <- as.numeric(b$tail_low_p)
  b$tail_high_p  <- as.numeric(b$tail_high_p)
  b$tail_tol_log <- as.numeric(b$tail_tol_log)
  b$tail_eps     <- as.numeric(b$tail_eps)

  if (!is.finite(b$tail_low_p) || b$tail_low_p <= 0 || b$tail_low_p >= 0.5) {
    stop("filter_bounds$tail_low_p must be in (0, 0.5).", call. = FALSE)
  }
  if (!is.finite(b$tail_high_p) || b$tail_high_p <= 0.5 || b$tail_high_p >= 1) {
    stop("filter_bounds$tail_high_p must be in (0.5, 1).", call. = FALSE)
  }
  if (!is.finite(b$tail_tol_log) || b$tail_tol_log <= 0) {
    stop("filter_bounds$tail_tol_log must be a positive finite number.", call. = FALSE)
  }
  if (!is.finite(b$tail_eps) || b$tail_eps <= 0) {
    stop("filter_bounds$tail_eps must be a positive finite number.", call. = FALSE)
  }

  if (verbose) {
    log_filtering_start(
      n_obs = n_obs0,
      n_sim = n_sim0,
      n_realizations = n_rlz,
      sample_target = n_select,
      relax_priority = RELAX_ORDER
    )
  }

  # ---------------------------------------------------------------------------
  # RNG management
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- .Random.seed
      has_seed <- TRUE
    } else {
      has_seed <- FALSE
    }
    on.exit({ if (has_seed) .Random.seed <<- old_seed }, add = TRUE)
    set.seed(seed)
  }

  # ---------------------------------------------------------------------------
  # Length harmonisation
  # ---------------------------------------------------------------------------
  n_use <- min(n_obs0, n_sim0)

  if (n_obs0 > n_use) {
    max_start_obs <- n_obs0 - n_use + 1L
    start_obs <- sample.int(max_start_obs, size = 1L)
    end_obs <- start_obs + n_use - 1L
    obs_use <- obs_series[start_obs:end_obs]
    obs_window <- c(start = start_obs, end = end_obs)
  } else {
    obs_use <- obs_series
    obs_window <- NULL
  }

  if (n_sim0 > n_use) {
    max_start_sim <- n_sim0 - n_use + 1L
    starts <- sample.int(max_start_sim, size = n_rlz, replace = TRUE)
    ends <- starts + n_use - 1L
    window_index <- lapply(seq_len(n_rlz), function(j) c(start = starts[j], end = ends[j]))

    sim_use <- vapply(
      seq_len(n_rlz),
      FUN = function(j) sim_series[starts[j]:ends[j], j],
      FUN.VALUE = numeric(n_use)
    )
    colnames(sim_use) <- colnames(sim_series)
  } else {
    window_index <- NULL
    sim_use <- sim_series
  }

  if (!is.matrix(sim_use)) sim_use <- as.matrix(sim_use)
  if (nrow(sim_use) != n_use) stop("Internal error: sim_use row mismatch.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # Base statistics
  # ---------------------------------------------------------------------------
  .log("Computing statistics (mean, sd, tail mass)", verbose = verbose, tag = "FILTER")

  obs_mean <- mean(obs_use)
  obs_sd   <- stats::sd(obs_use)

  sim_means <- colMeans(sim_use)
  sim_sds   <- apply(sim_use, 2, stats::sd)

  eps0 <- 1e-10
  rel_diff_vec <- function(sim, obs) {
    if (!is.finite(obs) || abs(obs) < eps0) rep(0, length(sim)) else (sim - obs) / obs
  }
  mean_rel_diff <- rel_diff_vec(sim_means, obs_mean)
  sd_rel_diff   <- rel_diff_vec(sim_sds,   obs_sd)

  # ---------------------------------------------------------------------------
  # Tail metrics (snake_case args)
  # ---------------------------------------------------------------------------
  tail_stats <- compute_tailmass_metrics(
    obs_use = obs_use,
    sim_series_stats = sim_use,
    tail_low_p = b$tail_low_p,
    tail_high_p = b$tail_high_p,
    tail_eps = b$tail_eps
  )

  .recompute_tailmass <- function() {
    tail_stats <<- compute_tailmass_metrics(
      obs_use = obs_use,
      sim_series_stats = sim_use,
      tail_low_p = b$tail_low_p,
      tail_high_p = b$tail_high_p,
      tail_eps = b$tail_eps
    )
  }

  # ---------------------------------------------------------------------------
  # Spectral matching metrics
  # ---------------------------------------------------------------------------
  .log("Computing spectral matching metrics", verbose = verbose, tag = "FILTER")

  spectral_results <- compute_spectral_metrics(
    obs_use = obs_use,
    sim_series_stats = sim_use,
    wavelet_pars = wavelet_args,
    modwt_n_levels = modwt_n_levels,
    n_sig_peaks_max = b$n_sig_peaks_max,
    peak_period_tol = b$peak_period_tol,
    peak_mag_tol_log = b$peak_mag_tol_log,
    eps = b$spectral_eps,
    parallel = parallel,
    n_cores = n_cores,
    cache_gws = cache_gws
  )

  wavelet_active <- spectral_results$active
  spectral_metrics <- spectral_results$metrics
  spectral_diag <- spectral_results$diagnostics

  period <- spectral_results$period
  gws_obs <- spectral_results$gws_obs
  gws_signif <- spectral_results$gws_signif
  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  .pool_from_pass <- function(pass_list) {
    ok <- rep(TRUE, n_rlz)
    for (nm in names(pass_list)) {
      v <- pass_list[[nm]]
      if (!is.null(v)) ok <- ok & as.logical(v)
    }
    which(ok)
  }

  .relax_one <- function(filter_name) {
    relax_bounds_one_filter(
      filter_name = filter_name,
      bounds_env = b,
      wavelet_active_env = environment(),
      recompute_tailmass_fn = .recompute_tailmass
    )
  }

  .pass_vectors_current <- function() {
    pass_mean <- abs(mean_rel_diff) <= b$mean
    pass_sd   <- abs(sd_rel_diff)   <= b$sd

    pass_tail_low  <- is.finite(tail_stats$logdiff_low)  & (tail_stats$logdiff_low  <= b$tail_tol_log)
    pass_tail_high <- is.finite(tail_stats$logdiff_high) & (tail_stats$logdiff_high <= b$tail_tol_log)

    pass_wavelet <- rep(TRUE, n_rlz)
    if (isTRUE(wavelet_active)) {

      pass_cor <- is.finite(spectral_metrics$spectral_cor) &
        (spectral_metrics$spectral_cor >= b$spectral_cor_min)

      pass_peak <- rep(TRUE, n_rlz)
      n_sig <- 0L
      if (!is.null(spectral_diag$n_sig_peaks_found)) n_sig <- as.integer(spectral_diag$n_sig_peaks_found)

      if (n_sig > 0L) {
        pass_peak <- is.finite(spectral_metrics$peak_match_frac) &
          (spectral_metrics$peak_match_frac >= b$peak_match_frac_min)
      }

      pass_wavelet <- pass_cor & pass_peak
    }

    list(
      mean      = pass_mean,
      sd        = pass_sd,
      tail_low  = pass_tail_low,
      tail_high = pass_tail_high,
      wavelet   = pass_wavelet
    )
  }

  passes <- .pass_vectors_current()
  pool_idx <- .pool_from_pass(passes)

  if (verbose) {
    log_filter_iteration(
      iter = 0L,
      passes = passes,
      pool = pool_idx,
      n_total = n_rlz,
      target = n_select,
      bounds = b,
      tail_metrics = tail_stats,
      wavelet_active = wavelet_active,
      spectral_diag = spectral_diag,
      note = "Initial evaluation."
    )
  }

  # ---------------------------------------------------------------------------
  # Relaxation loop
  # ---------------------------------------------------------------------------
  relax_log <- character(0)
  selection_mode <- "tiered"

  iter <- 0L
  while (length(pool_idx) < n_select && iter < b$relax_max_iter) {
    iter <- iter + 1L

    active <- names(passes)[!vapply(passes, is.null, logical(1))]
    active <- intersect(RELAX_ORDER, active)
    if (length(active) == 0) break

    rates <- vapply(active, function(nm) mean(passes[[nm]]), numeric(1))
    filter_to_relax <- active[which.min(rates)]

    before_pool <- length(pool_idx)

    rel <- .relax_one(filter_to_relax)
    if (!isTRUE(rel$changed)) break

    passes <- .pass_vectors_current()
    pool_idx <- .pool_from_pass(passes)

    note <- sprintf(
      "Relaxed '%s': %s | pool %s -> %s",
      filter_to_relax,
      rel$msg,
      format(before_pool, big.mark = ","),
      format(length(pool_idx), big.mark = ",")
    )
    relax_log <- c(relax_log, sprintf("Iteration %d: %s", iter, note))

    if (verbose) {
      log_filter_iteration(
        iter = iter,
        passes = passes,
        pool = pool_idx,
        n_total = n_rlz,
        target = n_select,
        bounds = b,
        tail_metrics = tail_stats,
        wavelet_active = wavelet_active,
        spectral_diag = spectral_diag,
        note = note
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Fallback if needed
  # ---------------------------------------------------------------------------
  if (length(pool_idx) < n_select) {
    selection_mode <- "closest_mean_fallback"
    pool_idx <- order(abs(mean_rel_diff))[seq_len(n_select)]
    if (verbose) {
      .log(sprintf(
        "Fallback activated: closest_mean to guarantee pool size = %s (from %s total).",
        format(n_select, big.mark = ","),
        format(n_rlz, big.mark = ",")
      ))
    }
  }

  # ---------------------------------------------------------------------------
  # Select from pool
  # ---------------------------------------------------------------------------
  .log("Selecting {format(n_select, big.mark = ',')} from pool of {format(length(pool_idx), big.mark = ',')}",
       tag = "FILTER", verbose = verbose)

  # Prefer best peak matching first, then overall spectral shape.
  # If no significant peaks exist, fall back to spectral correlation only.
  idx_pool <- pool_idx

  peak_score <- spectral_metrics$peak_match_frac[idx_pool]
  cor_score  <- spectral_metrics$spectral_cor[idx_pool]

  has_sig_peaks <- !is.null(spectral_diag$n_sig_peaks_found) &&
    is.finite(spectral_diag$n_sig_peaks_found) &&
    as.integer(spectral_diag$n_sig_peaks_found) > 0L

  if (has_sig_peaks) {

    # Order: highest peak_match_frac, then highest spectral_cor
    ord <- order(
      -ifelse(is.finite(peak_score), peak_score, -Inf),
      -ifelse(is.finite(cor_score),  cor_score,  -Inf)
    )

  } else {

    # No significant peaks to enforce: rank by overall spectral shape only
    ord <- order(-ifelse(is.finite(cor_score), cor_score, -Inf))
  }

  idx_select <- idx_pool[ord][seq_len(n_select)]
  selection_mode <- if (selection_mode == "tiered") "ranked_wavelet_priority" else selection_mode


  # ---------------------------------------------------------------------------
  # Summary table
  # ---------------------------------------------------------------------------
  report_filters <- c("mean", "sd", "tail_low", "tail_high", "wavelet")
  n_passed <- vapply(report_filters, function(nm) {
    if (!is.null(passes[[nm]])) sum(passes[[nm]]) else NA_integer_
  }, integer(1))
  pct_passed <- vapply(report_filters, function(nm) {
    if (!is.null(passes[[nm]])) mean(passes[[nm]]) * 100 else NA_real_
  }, numeric(1))

  summary <- data.frame(
    filter = report_filters,
    n_passed = n_passed,
    pct_passed = pct_passed,
    selection_mode = selection_mode,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    log_final_summary(
      pool_size = length(pool_idx),
      n_total = n_rlz,
      n_sampled = length(idx_select),
      relaxation_level = selection_mode
    )
  }

  # ---------------------------------------------------------------------------
  # Plots
  # ---------------------------------------------------------------------------
  plots_out <- NULL
  if (make_plots) {

    .log("Creating diagnostic plots for the selected realizations",
         tag = "FILTER", verbose = verbose)

    if (length(pool_idx) < 1L) {
      warning("make_plots=TRUE but final pool is empty; plots=NULL.", call. = FALSE)
    } else {

      plots_out <- plot_filter_diagnostics(
        obs_series    = obs_use,
        sim_series    = sim_use,
        pool          = idx_select,
        rel_diff_mean = mean_rel_diff,
        rel_diff_sd   = sd_rel_diff,
        tail_metrics  = tail_stats,
        power_period  = period,
        power_obs     = gws_obs,
        power_signif  = gws_signif,
        wavelet_pars  = wavelet_args,
        wavelet_q     = b$plot_wavelet_q
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  diagnostics <- list(
    n_obs_original = n_obs0,
    n_sim_original = n_sim0,
    n_use = n_use,
    n_periods = length(period),
    obs_window = obs_window,
    window_index = window_index,
    pool_idx = pool_idx,
    selected_idx = idx_select,
    relax_log = relax_log,
    spectral_diag = spectral_diag,
    spectral_metrics = spectral_metrics,
    final_bounds = as.list(b)
  )

  if (cache_gws) {
    diagnostics$gws_cache <- spectral_results$gws_cache
  }

  list(
    pool = sim_series[, pool_idx, drop = FALSE],
    selected = sim_series[, idx_select, drop = FALSE],
    summary = summary,
    diagnostics = diagnostics,
    plots = plots_out
  )
}


#' WARM filtering default bounds
#'
#' @description
#' Returns internal default bounds for \code{filter_warm_pool()}.
#'
#' The defaults are designed to be moderately selective for annual records on the
#' order of 50 to 100 years. Typical usage is to override only a small subset of
#' entries via \code{filter_bounds = list(...)} while keeping the remaining
#' defaults unchanged.
#'
#' @details
#' The returned list contains thresholds and controls used by the filtering and
#' relaxation logic.
#'
#' Distributional tolerances:
#' \describe{
#'   \item{mean}{Maximum absolute relative difference in the mean between a
#'   simulated realization and the observed series. Default 0.03.}
#'   \item{sd}{Maximum absolute relative difference in standard deviation between
#'   a simulated realization and the observed series. Default 0.03.}
#' }
#'
#' Tail mass behavior:
#' \describe{
#'   \item{tail_low_p}{Lower quantile used to define the low tail threshold in the
#'   observed series. Default 0.20.}
#'   \item{tail_high_p}{Upper quantile used to define the high tail threshold in
#'   the observed series. Default 0.80.}
#'   \item{tail_tol_log}{Maximum absolute log difference between simulated and
#'   observed tail mass metrics. Default log(1.03).}
#'   \item{tail_eps}{Positive constant used for numerical stability in log
#'   transforms. Default 1e-5.}
#' }
#'
#' Spectral similarity:
#' \describe{
#'   \item{spectral_cor_min}{Minimum correlation between log transformed observed
#'   and simulated global wavelet spectra. Default 0.60.}
#'   \item{spectral_eps}{Positive constant used for numerical stability in
#'   spectral log transforms. Default 1e-10.}
#' }
#'
#' Peak matching:
#' \describe{
#'   \item{n_sig_peaks_max}{Maximum number of significant observed spectral peaks
#'   to enforce. Default 2.}
#'   \item{peak_period_tol}{Tolerance for matching peak periods in log2 period
#'   space. Default 0.50.}
#'   \item{peak_mag_tol_log}{Tolerance for matching peak magnitudes as the
#'   absolute log ratio between simulated and observed peak power. Default
#'   log(1.5).}
#'   \item{peak_match_frac_min}{Minimum fraction of significant observed peaks
#'   that must be matched. Default 1.0.}
#' }
#'
#' Plot controls:
#' \describe{
#'   \item{plot_wavelet_q}{Two probabilities used to summarize simulated spectra
#'   in diagnostic plots. Default c(0.50, 0.95).}
#' }
#'
#' Relaxation controls:
#' \describe{
#'   \item{relax_mult}{Multiplicative factor applied when relaxing some bounds.
#'   Default 1.25.}
#'   \item{relax_mean_max}{Maximum allowed mean tolerance during relaxation.
#'   Default 0.25.}
#'   \item{relax_sd_max}{Maximum allowed sd tolerance during relaxation.
#'   Default 0.25.}
#'   \item{relax_tail_tol_log_max}{Maximum tail tolerance during relaxation.
#'   Default log(2.0).}
#'   \item{relax_tail_p_step}{Step size for relaxing tail quantiles. Default
#'   0.02.}
#'   \item{relax_tail_p_low_max}{Maximum lower tail quantile during relaxation.
#'   Default 0.40.}
#'   \item{relax_tail_p_high_min}{Minimum upper tail quantile during relaxation.
#'   Default 0.40.}
#'   \item{relax_spectral_cor_step}{Decrement applied to spectral_cor_min during
#'   wavelet relaxation. Default 0.05.}
#'   \item{relax_spectral_cor_min}{Minimum spectral_cor_min allowed during
#'   relaxation. Default 0.30.}
#'   \item{relax_peak_match_frac_step}{Decrement applied to peak_match_frac_min
#'   during wavelet relaxation. Default 0.10.}
#'   \item{relax_peak_match_frac_min}{Minimum peak_match_frac_min allowed during
#'   relaxation. Default 0.00.}
#'   \item{relax_max_iter}{Maximum number of relaxation iterations. Default 20.}
#' }
#'
#' @return Named list of default bounds using snake case keys.
#' @keywords internal
#' @export
filter_warm_bounds_defaults <- function() {
  list(
    # --- distributional tolerances (relative diff) ---
    mean = 0.03,
    sd   = 0.03,

    # --- tail behaviour ---
    tail_low_p   = 0.20,
    tail_high_p  = 0.80,
    tail_tol_log = log(1.03),
    tail_eps     = 1e-5,

    # --- spectral matching (overall shape) ---
    spectral_cor_min = 0.60,
    spectral_eps     = 1e-10,

    # --- peak matching (significant observed peaks only) ---
    n_sig_peaks_max      = 2L,
    peak_period_tol      = 0.50,     # tolerance in log2(period) (octaves)
    peak_mag_tol_log     = log(1.5), # abs(log(sim/obs)) <= log(1.5) => within 50%
    peak_match_frac_min  = 1.0,      # require all significant peaks found to match (set <1 to relax)

    # --- plotting diagnostics ---
    plot_wavelet_q = c(0.50, 0.95),

    # --- relaxation controls ---
    relax_mult     = 1.25,
    relax_mean_max = 0.25,
    relax_sd_max   = 0.25,

    relax_tail_tol_log_max = log(2.0),
    relax_tail_p_step      = 0.02,
    relax_tail_p_low_max   = 0.40,
    relax_tail_p_high_min  = 0.40,

    # Spectral relaxation
    relax_spectral_cor_step     = 0.05,
    relax_spectral_cor_min      = 0.30,
    relax_peak_match_frac_step  = 0.10,
    relax_peak_match_frac_min   = 0.00,
    relax_max_iter = 20L
  )
}



# ==============================================================================
# Spectral Matching Metrics
# ==============================================================================

#' Find local maxima indices in a numeric vector
#'
#' @description
#' Scans a numeric vector and returns indices of local maxima. A local maximum is
#' defined by comparing each element to its immediate neighbors.
#'
#' @param x Numeric vector. Values to scan for local maxima.
#' @param strict Logical scalar. If TRUE, requires strict inequality on both
#'   sides. If FALSE, allows ties but still requires at least one strict
#'   inequality so flat plateaus do not produce multiple peaks.
#'
#' @return Integer vector of peak indices. Returns an empty integer vector if
#'   fewer than three values are available.
#' @keywords internal
find_local_maxima <- function(x, strict = TRUE) {
  n <- length(x)
  if (n < 3L) return(integer(0))

  peaks <- integer(0)
  for (i in 2:(n - 1)) {
    if (strict) {
      if (x[i] > x[i - 1] && x[i] > x[i + 1]) {
        peaks <- c(peaks, i)
      }
    } else {
      if (x[i] >= x[i - 1] && x[i] >= x[i + 1] && (x[i] > x[i - 1] || x[i] > x[i + 1])) {
        peaks <- c(peaks, i)
      }
    }
  }
  peaks
}


#' Identify significant peaks in an observed global wavelet spectrum
#'
#' @description
#' Detects local maxima in an observed global wavelet spectrum and keeps only
#' those peaks whose power exceeds the corresponding significance curve value.
#'
#' Peaks are ranked by signal to noise ratio defined as power divided by
#' significance, then by power. Up to \code{n_max} peaks are returned.
#'
#' If the significance curve is missing or is not aligned to \code{gws}, the
#' function returns an empty result, treating the significance information as
#' unavailable.
#'
#' @param gws Numeric vector. Observed global wavelet spectrum values.
#' @param gws_signif Numeric vector or NULL. Significance curve aligned to
#'   \code{gws}. Must have the same length as \code{gws}.
#' @param period Numeric vector. Period values associated with \code{gws}. Must
#'   have the same length as \code{gws}.
#' @param n_max Integer scalar. Maximum number of significant peaks to return.
#'
#' @return Data frame with one row per selected peak and columns:
#'   \code{idx} index in the spectrum,
#'   \code{period} the period at the peak,
#'   \code{power} the peak power,
#'   \code{signif} the significance curve value at the peak,
#'   \code{snr} signal to noise ratio.
#'   Returns a data frame with zero rows if no significant peaks are found.
#' @keywords internal
identify_significant_peaks <- function(gws, gws_signif, period, n_max = 3L) {

  if (length(gws) < 3L || length(period) != length(gws)) {
    return(data.frame(idx = integer(0), period = numeric(0),
                      power = numeric(0), signif = numeric(0), snr = numeric(0)))
  }

  gws_clean <- fill_nearest(as.numeric(gws))

  # If signif is missing or wrong length: treat as "no significance available"
  if (is.null(gws_signif) || !is.numeric(gws_signif) || length(gws_signif) != length(gws_clean)) {
    return(data.frame(idx = integer(0), period = numeric(0),
                      power = numeric(0), signif = numeric(0), snr = numeric(0)))
  }

  signif_clean <- fill_nearest(as.numeric(gws_signif))

  peak_idx <- find_local_maxima(gws_clean, strict = FALSE)
  if (length(peak_idx) == 0L) {
    return(data.frame(idx = integer(0), period = numeric(0),
                      power = numeric(0), signif = numeric(0), snr = numeric(0)))
  }

  # Significant peaks only
  keep <- is.finite(gws_clean[peak_idx]) & is.finite(signif_clean[peak_idx]) &
    (gws_clean[peak_idx] > signif_clean[peak_idx])

  peak_idx <- peak_idx[keep]
  if (length(peak_idx) == 0L) {
    return(data.frame(idx = integer(0), period = numeric(0),
                      power = numeric(0), signif = numeric(0), snr = numeric(0)))
  }

  snr <- gws_clean[peak_idx] / pmax(signif_clean[peak_idx], 1e-12)

  # Rank by SNR then by power
  ord <- order(snr, gws_clean[peak_idx], decreasing = TRUE)
  peak_idx <- peak_idx[ord]
  snr <- snr[ord]

  n_keep <- min(as.integer(n_max), length(peak_idx))
  peak_idx <- peak_idx[seq_len(n_keep)]
  snr <- snr[seq_len(n_keep)]

  data.frame(
    idx = peak_idx,
    period = period[peak_idx],
    power = gws_clean[peak_idx],
    signif = signif_clean[peak_idx],
    snr = snr
  )
}


#' Compute significant peak match metrics for one simulated spectrum
#'
#' @description
#' Compares a simulated global wavelet spectrum to a set of significant observed
#' peaks. For each observed peak, a match is declared if the simulated spectrum
#' contains a local maximum within a period tolerance and the matched power is
#' within a magnitude tolerance.
#'
#' Period tolerance is measured in log2 period space. Magnitude tolerance is
#' measured as the absolute log ratio between simulated and observed power.
#'
#' @param gws_sim Numeric vector. Simulated global wavelet spectrum aligned to
#'   \code{period}.
#' @param period Numeric vector. Period grid for \code{gws_sim}.
#' @param obs_peaks Data frame. Significant observed peaks returned by
#'   \code{identify_significant_peaks()}.
#' @param period_tol Numeric scalar. Maximum allowed absolute difference in
#'   log2 period between an observed peak period and a simulated candidate period.
#' @param mag_tol_log Numeric scalar. Maximum allowed absolute log ratio between
#'   simulated and observed peak power.
#' @param eps Numeric scalar. Positive constant used for numerical stability.
#'
#' @return Named list with:
#'   \code{peak_match_frac} fraction of observed peaks that were matched,
#'   \code{peak_mag_mean_abs_log_ratio} mean absolute log ratio for matched peaks,
#'   or NA if no peaks were matched.
#' @keywords internal
compute_peak_match_metrics <- function(gws_sim, period, obs_peaks,
                                       period_tol = 0.5,
                                       mag_tol_log = log(1.5),
                                       eps = 1e-10) {

  if (is.null(obs_peaks) || nrow(obs_peaks) == 0L) {
    return(list(peak_match_frac = 1.0, peak_mag_mean_abs_log_ratio = NA_real_))
  }

  gws_sim_clean <- fill_nearest(as.numeric(gws_sim))

  n_matched <- 0L
  abs_log_ratios <- numeric(0)

  for (i in seq_len(nrow(obs_peaks))) {
    obs_period <- obs_peaks$period[i]
    obs_power  <- obs_peaks$power[i]

    # candidate indices within period tolerance (log2 space)
    log2_diff <- abs(log2(period) - log2(obs_period))
    within_tol <- which(is.finite(log2_diff) & (log2_diff <= period_tol))

    if (length(within_tol) == 0L) next

    sim_power <- max(gws_sim_clean[within_tol], na.rm = TRUE)
    if (!is.finite(sim_power) || !is.finite(obs_power)) next

    abs_log_ratio <- abs(log((sim_power + eps) / (obs_power + eps)))

    # Match requires magnitude agreement within tolerance
    if (abs_log_ratio <= mag_tol_log) {
      n_matched <- n_matched + 1L
      abs_log_ratios <- c(abs_log_ratios, abs_log_ratio)
    }
  }

  frac <- n_matched / nrow(obs_peaks)
  mean_abs <- if (length(abs_log_ratios) > 0L) mean(abs_log_ratios) else NA_real_

  list(
    peak_match_frac = frac,
    peak_mag_mean_abs_log_ratio = mean_abs
  )
}



#' Compute spectral match metrics for one realization
#'
#' @description
#' Computes two spectral similarity metrics between an observed and a simulated
#' global wavelet spectrum:
#' 1) correlation between log transformed spectra, and
#' 2) significant peak match metrics relative to significant observed peaks.
#'
#' @param gws_obs Numeric vector. Observed global wavelet spectrum aligned to
#'   \code{period}.
#' @param gws_sim Numeric vector. Simulated global wavelet spectrum aligned to
#'   \code{period}.
#' @param period Numeric vector. Period grid for both spectra.
#' @param obs_peaks Data frame. Significant observed peaks returned by
#'   \code{identify_significant_peaks()}.
#' @param peak_period_tol Numeric scalar. Period matching tolerance in log2
#'   period space.
#' @param peak_mag_tol_log Numeric scalar. Magnitude matching tolerance as an
#'   absolute log ratio.
#' @param eps Numeric scalar. Positive constant used to avoid log of zero.
#'
#' @return Named list with:
#'   \code{spectral_cor} correlation between log transformed spectra,
#'   \code{peak_match_frac} fraction of significant observed peaks matched,
#'   \code{peak_mag_mean_abs_log_ratio} mean absolute log ratio for matched peaks.
#' @keywords internal
compute_spectral_match_single <- function(gws_obs, gws_sim, period,
                                          obs_peaks,
                                          peak_period_tol = 0.5,
                                          peak_mag_tol_log = log(1.5),
                                          eps = 1e-10) {

  gws_obs <- pmax(as.numeric(gws_obs), eps)
  gws_sim <- pmax(as.numeric(gws_sim), eps)

  # Overall spectral shape match (log-space correlation)
  log_obs <- log(gws_obs)
  log_sim <- log(gws_sim)
  ok <- is.finite(log_obs) & is.finite(log_sim)

  spectral_cor <- if (sum(ok) < 3L) NA_real_ else stats::cor(log_sim[ok], log_obs[ok])

  # Significant-peak matching (period + magnitude)
  pk <- compute_peak_match_metrics(
    gws_sim = gws_sim,
    period = period,
    obs_peaks = obs_peaks,
    period_tol = peak_period_tol,
    mag_tol_log = peak_mag_tol_log,
    eps = eps
  )

  list(
    spectral_cor = spectral_cor,
    peak_match_frac = pk$peak_match_frac,
    peak_mag_mean_abs_log_ratio = pk$peak_mag_mean_abs_log_ratio
  )
}



#' Compute spectral metrics for an ensemble of realizations
#'
#' @description
#' Computes wavelet based spectral similarity metrics between an observed annual
#' series and each realization in a simulated ensemble.
#'
#' The observed series is analyzed once to obtain a period grid, an observed
#' global wavelet spectrum, a significance curve, and a set of significant peaks.
#' Each simulated realization is then analyzed on the same period grid and is
#' scored using:
#' 1) correlation between log transformed spectra, and
#' 2) significant peak match fraction and magnitude agreement.
#'
#' If \code{parallel = TRUE}, the per realization computation can run in parallel.
#' Parallel execution uses a PSOCK cluster created by the base parallel package,
#' which is typically the most portable option across operating systems.
#'
#' @param obs_use Numeric vector. Observed annual series after any window
#'   alignment performed by the caller.
#' @param sim_series_stats Numeric matrix. Simulated ensemble after any window
#'   alignment performed by the caller. Rows are years and columns are
#'   realizations.
#' @param wavelet_pars Named list. Parameters passed to
#'   \code{analyze_wavelet_spectrum()}. Expected entries include \code{signif_level},
#'   \code{noise_type}, \code{period_lower_limit}, and \code{detrend}.
#' @param modwt_n_levels Integer or NULL. Number of MODWT levels used in WARM.
#'   Used only for diagnostics. If NULL, a value is estimated from series length.
#' @param n_sig_peaks_max Integer scalar. Maximum number of significant observed
#'   peaks to enforce when computing peak matching metrics.
#' @param peak_period_tol Numeric scalar. Period matching tolerance in log2
#'   period space.
#' @param peak_mag_tol_log Numeric scalar. Magnitude matching tolerance as an
#'   absolute log ratio.
#' @param eps Numeric scalar. Positive constant used for numerical stability in
#'   log transforms and ratios.
#' @param parallel Logical scalar. If TRUE, compute per realization spectra and
#'   metrics in parallel.
#' @param n_cores Integer scalar or NULL. Number of worker processes to use when
#'   \code{parallel = TRUE}. If NULL, an internal default is used.
#' @param cache_gws Logical scalar. If TRUE, store simulated global wavelet
#'   spectra in \code{gws_cache}. When FALSE, \code{gws_cache} is returned as
#'   NULL and simulated spectra are not retained.
#'
#' @return Named list with:
#' \describe{
#'   \item{active}{Logical scalar. TRUE when wavelet metrics were computed.}
#'   \item{period}{Numeric vector. Period grid used for all spectra.}
#'   \item{gws_obs}{Numeric vector. Observed global wavelet spectrum on \code{period}.}
#'   \item{gws_signif}{Numeric vector or NULL. Significance curve on \code{period}.}
#'   \item{gws_cache}{Numeric matrix or NULL. Simulated spectra on \code{period}
#'   for each realization when \code{cache_gws = TRUE}; otherwise NULL.}
#'   \item{metrics}{List with numeric vectors \code{spectral_cor},
#'   \code{peak_match_frac}, and \code{peak_mag_mean_abs_log_ratio}.}
#'   \item{diagnostics}{List. Summary information including significant peaks
#'   and observed spectrum summary statistics.}
#' }
#'
#' @keywords internal
#' @export
compute_spectral_metrics <- function(obs_use, sim_series_stats, wavelet_pars,
                                     modwt_n_levels = NULL,
                                     n_sig_peaks_max = 2L,
                                     peak_period_tol = 0.50,
                                     peak_mag_tol_log = log(1.5),
                                     eps = 1e-10,
                                     parallel = FALSE,
                                     n_cores = NULL,
                                     cache_gws = FALSE) {


  parallel <- isTRUE(parallel)
  cache_gws <- isTRUE(cache_gws)

  if (!is.null(n_cores)) {
    if (!is.numeric(n_cores) || length(n_cores) != 1L || !is.finite(n_cores) || n_cores < 1) {
      stop("'n_cores' must be a positive integer or NULL.", call. = FALSE)
    }
    n_cores <- as.integer(n_cores)
  }

  n_use <- length(obs_use)
  n_realizations <- ncol(sim_series_stats)

  if (is.null(modwt_n_levels) || !is.finite(modwt_n_levels) || modwt_n_levels < 1) {
    modwt_n_levels <- max(2L, floor(log2(n_use)) - 1L)
  }
  modwt_n_levels <- as.integer(modwt_n_levels)

  # -------------------------------------------------------------------------
  # Observed wavelet analysis (sequential; only once)
  # -------------------------------------------------------------------------
  wv_obs <- analyze_wavelet_spectrum(
    obs_use,
    signif = wavelet_pars$signif_level,
    noise  = wavelet_pars$noise_type,
    min_period = wavelet_pars$period_lower_limit,
    detrend = isTRUE(wavelet_pars$detrend),
    mode = "fast"
  )

  if (is.null(wv_obs$period) || !is.numeric(wv_obs$period)) {
    stop("analyze_wavelet_spectrum(obs) must return numeric $period.", call. = FALSE)
  }

  period <- as.numeric(wv_obs$period)

  # Observed GWS (use unmasked when available)
  gws_obs <- if (!is.null(wv_obs$gws_unmasked) && is.numeric(wv_obs$gws_unmasked)) {
    as.numeric(wv_obs$gws_unmasked)
  } else {
    as.numeric(wv_obs$gws)
  }
  gws_obs <- gws_regrid(wv_obs, period, use_unmasked = TRUE)
  gws_obs <- fill_nearest(gws_obs)

  # Significance curve aligned to 'period' (single regrid only)
  gws_signif <- wv_obs$gws_signif_unmasked
  if (is.null(gws_signif)) gws_signif <- wv_obs$gws_signif

  if (!is.null(gws_signif) && is.numeric(gws_signif)) {
    wv_sig <- list(period = wv_obs$period, gws = as.numeric(gws_signif))
    gws_signif <- gws_regrid(wv_sig, target_period = period, use_unmasked = FALSE)
    gws_signif <- fill_nearest(as.numeric(gws_signif))
  } else {
    gws_signif <- NULL
  }

  # Significant observed peaks only
  obs_peaks_sig <- identify_significant_peaks(
    gws = gws_obs,
    gws_signif = gws_signif,
    period = period,
    n_max = n_sig_peaks_max
  )

  # -------------------------------------------------------------------------
  # Worker function: analyze one simulation column
  # -------------------------------------------------------------------------
  one_sim <- function(j, sim_series_stats, wavelet_pars, period, gws_obs, obs_peaks_sig,
                      peak_period_tol, peak_mag_tol_log, eps, cache_gws) {

    wv_sim <- analyze_wavelet_spectrum(
      sim_series_stats[, j],
      signif = wavelet_pars$signif_level,
      noise  = wavelet_pars$noise_type,
      min_period = wavelet_pars$period_lower_limit,
      detrend = isTRUE(wavelet_pars$detrend),
      mode = "fast"
    )

    gws_sim <- gws_regrid(wv_sim, period, use_unmasked = TRUE)
    gws_sim <- fill_nearest(gws_sim)

    m <- compute_spectral_match_single(
      gws_obs = gws_obs,
      gws_sim = gws_sim,
      period = period,
      obs_peaks = obs_peaks_sig,
      peak_period_tol = peak_period_tol,
      peak_mag_tol_log = peak_mag_tol_log,
      eps = eps
    )

    out <- list(
      spectral_cor = m$spectral_cor,
      peak_match_frac = m$peak_match_frac,
      peak_mag_mean_abs_log_ratio = m$peak_mag_mean_abs_log_ratio
    )
    if (isTRUE(cache_gws)) out$gws_sim <- gws_sim
    out
  }

  # -------------------------------------------------------------------------
  # Run loop: sequential or parallel
  # -------------------------------------------------------------------------
  gws_cache <- NULL
  if (cache_gws) {
    gws_cache <- matrix(NA_real_, nrow = length(period), ncol = n_realizations)
  }
  spectral_cor <- rep(NA_real_, n_realizations)
  peak_match_frac <- rep(NA_real_, n_realizations)
  peak_mag_mean_abs_log_ratio <- rep(NA_real_, n_realizations)

  if (!parallel || n_realizations < 2L) {

    for (j in seq_len(n_realizations)) {
        res <- one_sim(
          j = j,
          sim_series_stats = sim_series_stats,
          wavelet_pars = wavelet_pars,
          period = period,
          gws_obs = gws_obs,
          obs_peaks_sig = obs_peaks_sig,
          peak_period_tol = peak_period_tol,
          peak_mag_tol_log = peak_mag_tol_log,
          eps = eps,
          cache_gws = cache_gws
        )

      if (cache_gws) gws_cache[, j] <- res$gws_sim
      spectral_cor[j] <- res$spectral_cor
      peak_match_frac[j] <- res$peak_match_frac
      peak_mag_mean_abs_log_ratio[j] <- res$peak_mag_mean_abs_log_ratio
    }

  } else {

    # robust core selection
    max_cores <- parallel::detectCores(logical = TRUE)
    if (!is.finite(max_cores) || max_cores < 1L) max_cores <- 1L

    use_cores <- if (is.null(n_cores)) max_cores else min(n_cores, max_cores)
    use_cores <- max(1L, as.integer(use_cores))

    if (use_cores == 1L) {
      for (j in seq_len(n_realizations)) {
        res <- one_sim(
          j = j,
          sim_series_stats = sim_series_stats,
          wavelet_pars = wavelet_pars,
          period = period,
          gws_obs = gws_obs,
          obs_peaks_sig = obs_peaks_sig,
          peak_period_tol = peak_period_tol,
          peak_mag_tol_log = peak_mag_tol_log,
          eps = eps,
          cache_gws = cache_gws
        )

        if (cache_gws) gws_cache[, j] <- res$gws_sim
        spectral_cor[j] <- res$spectral_cor
        peak_match_frac[j] <- res$peak_match_frac
        peak_mag_mean_abs_log_ratio[j] <- res$peak_mag_mean_abs_log_ratio
      }

    } else {

      cl <- parallel::makeCluster(use_cores)
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

      # Ensure worker sessions can see required functions/objects.
      # Use explicit var list; this avoids "object not found" on Windows PSOCK.
      parallel::clusterExport(
        cl,
        varlist = c(
          "analyze_wavelet_spectrum",
          "gws_regrid",
          "fill_nearest",
          "compute_spectral_match_single",
          "compute_peak_match_metrics",
          "identify_significant_peaks",
          "find_local_maxima",
          "one_sim"
        ),
        envir = environment()
      )

      # Load required namespaces on workers (safe even if already loaded)
      parallel::clusterEvalQ(cl, {
        NULL
      })

      idx <- as.integer(seq_len(n_realizations))

      res_list <- parallel::parLapply(
        cl,
        X = idx,
        fun = function(j, sim_series_stats, wavelet_pars, period, gws_obs, obs_peaks_sig,
                       peak_period_tol, peak_mag_tol_log, eps, cache_gws) {
          one_sim(
            j = j,
            sim_series_stats = sim_series_stats,
            wavelet_pars = wavelet_pars,
            period = period,
            gws_obs = gws_obs,
            obs_peaks_sig = obs_peaks_sig,
            peak_period_tol = peak_period_tol,
            peak_mag_tol_log = peak_mag_tol_log,
            eps = eps,
            cache_gws = cache_gws
          )
        },
        sim_series_stats = sim_series_stats,
        wavelet_pars = wavelet_pars,
        period = period,
        gws_obs = gws_obs,
        obs_peaks_sig = obs_peaks_sig,
        peak_period_tol = peak_period_tol,
        peak_mag_tol_log = peak_mag_tol_log,
        eps = eps,
        cache_gws = cache_gws
      )

      for (j in seq_len(n_realizations)) {
        res <- res_list[[j]]
        if (cache_gws) gws_cache[, j] <- res$gws_sim
        spectral_cor[j] <- res$spectral_cor
        peak_match_frac[j] <- res$peak_match_frac
        peak_mag_mean_abs_log_ratio[j] <- res$peak_mag_mean_abs_log_ratio
      }
    }
  }

  spectral_diag <- list(
    modwt_n_levels = modwt_n_levels,
    n_periods = length(period),
    obs_peaks_sig = obs_peaks_sig,
    n_sig_peaks_found = nrow(obs_peaks_sig),
    gws_obs_summary = c(
      min = min(gws_obs, na.rm = TRUE),
      median = stats::median(gws_obs, na.rm = TRUE),
      max = max(gws_obs, na.rm = TRUE)
    )
  )

  list(
    active = TRUE,
    period = period,
    gws_obs = gws_obs,
    gws_signif = gws_signif,
    gws_cache = gws_cache,
    metrics = list(
      spectral_cor = spectral_cor,
      peak_match_frac = peak_match_frac,
      peak_mag_mean_abs_log_ratio = peak_mag_mean_abs_log_ratio
    ),
    diagnostics = spectral_diag
  )
}



# ==============================================================================
# Tail-mass Metrics
# ==============================================================================

#' Compute tail mass metrics for filtering
#'
#' @description
#' Computes tail mass metrics for an observed series and a simulated ensemble.
#' Tail mass metrics quantify how much probability mass lies beyond observed
#' quantile thresholds, expressed as normalized deficit or excess mass.
#'
#' The method uses robust scale estimation for normalization. The scale is chosen
#' as IQR when available, otherwise MAD, otherwise standard deviation, otherwise 1.
#' Tail masses are normalized by series length and scale to make values more
#' comparable across datasets.
#'
#' The returned log difference vectors are used by \code{filter_warm_pool()} to
#' decide whether a realization passes tail filters.
#'
#' @param obs_use Numeric vector. Observed annual values after any window
#'   alignment performed by the caller.
#' @param sim_series_stats Numeric matrix. Simulated values aligned to
#'   \code{obs_use}. Rows are years and columns are realizations.
#' @param tail_low_p Numeric scalar. Lower tail quantile probability used to
#'   compute the low threshold on the observed series.
#' @param tail_high_p Numeric scalar. Upper tail quantile probability used to
#'   compute the high threshold on the observed series.
#' @param tail_eps Numeric scalar. Positive constant added inside log transforms
#'   to avoid log of zero.
#'
#' @return Named list with tail thresholds, scale, observed and simulated tail
#'   mass metrics, and log difference vectors:
#' \describe{
#'   \item{thr_low}{Lower threshold from the observed series.}
#'   \item{thr_high}{Upper threshold from the observed series.}
#'   \item{scale_obs}{Robust scale used for normalization.}
#'   \item{M_obs_low}{Observed low tail deficit mass.}
#'   \item{M_obs_high}{Observed high tail excess mass.}
#'   \item{M_sim_low}{Vector of simulated low tail deficit masses.}
#'   \item{M_sim_high}{Vector of simulated high tail excess masses.}
#'   \item{logdiff_low}{Vector of absolute log differences for low tail mass.}
#'   \item{logdiff_high}{Vector of absolute log differences for high tail mass.}
#' }
#'
#' @keywords internal
#' @export
compute_tailmass_metrics <- function(obs_use, sim_series_stats,
                                     tail_low_p, tail_high_p, tail_eps) {

  n_use <- length(obs_use)
  n_realizations <- ncol(sim_series_stats)

  # Robust scale with fallback hierarchy: IQR -> MAD -> SD -> 1
  scale_obs <- stats::IQR(obs_use, na.rm = TRUE, type = 7)
  if (!is.finite(scale_obs) || scale_obs <= 0) {
    scale_obs <- stats::mad(obs_use, constant = 1, na.rm = TRUE)
  }
  if (!is.finite(scale_obs) || scale_obs <= 0) {
    scale_obs <- stats::sd(obs_use, na.rm = TRUE)
  }
  if (!is.finite(scale_obs) || scale_obs <= 0) scale_obs <- 1.0

  thr_low <- stats::quantile(obs_use, probs = tail_low_p, names = FALSE, type = 7, na.rm = TRUE)
  thr_high <- stats::quantile(obs_use, probs = tail_high_p, names = FALSE, type = 7, na.rm = TRUE)

  if (!is.finite(thr_low) || !is.finite(thr_high) || !is.finite(scale_obs)) {
    return(list(
      thr_low = NA_real_, thr_high = NA_real_, scale_obs = NA_real_,
      M_obs_low = NA_real_, M_obs_high = NA_real_,
      M_sim_low = rep(NA_real_, n_realizations),
      M_sim_high = rep(NA_real_, n_realizations),
      logdiff_low = rep(Inf, n_realizations),
      logdiff_high = rep(Inf, n_realizations)
    ))
  }

  denom <- max(n_use * scale_obs, 1e-12)

  # Lower-tail deficit mass
  X_low <- thr_low - sim_series_stats
  S_sim_low <- colSums(pmax(X_low, 0), na.rm = TRUE)
  S_obs_low <- sum(pmax(thr_low - obs_use, 0), na.rm = TRUE)

  # Upper-tail excess mass
  X_high <- sim_series_stats - thr_high
  S_sim_high <- colSums(pmax(X_high, 0), na.rm = TRUE)
  S_obs_high <- sum(pmax(obs_use - thr_high, 0), na.rm = TRUE)

  M_obs_low <- S_obs_low / denom
  M_obs_high <- S_obs_high / denom
  M_sim_low <- S_sim_low / denom
  M_sim_high <- S_sim_high / denom

  logdiff_low <- abs(log(M_sim_low + tail_eps) - log(M_obs_low + tail_eps))
  logdiff_high <- abs(log(M_sim_high + tail_eps) - log(M_obs_high + tail_eps))

  list(
    thr_low = thr_low,
    thr_high = thr_high,
    scale_obs = scale_obs,
    M_obs_low = M_obs_low,
    M_obs_high = M_obs_high,
    M_sim_low = M_sim_low,
    M_sim_high = M_sim_high,
    logdiff_low = logdiff_low,
    logdiff_high = logdiff_high
  )
}


# ==============================================================================
# Relaxation Logic
# ==============================================================================

#' Relax bounds for one filter family
#'
#' @description
#' Applies one relaxation step for a single filter family and updates the bounds
#' environment in place. This function is called by \code{filter_warm_pool()}.
#'
#' Relaxation behavior by filter family:
#' \describe{
#'   \item{mean}{Increases the mean tolerance up to \code{relax_mean_max}.}
#'   \item{sd}{Increases the sd tolerance up to \code{relax_sd_max}.}
#'   \item{tail_low}{First increases tail tolerance, then increases
#'   \code{tail_low_p} in steps, recomputing tail metrics after changes.}
#'   \item{tail_high}{First increases tail tolerance, then decreases
#'   \code{tail_high_p} in steps, recomputing tail metrics after changes.}
#'   \item{wavelet}{Relaxes spectral correlation threshold, then relaxes peak
#'   match fraction, then disables peak matching, then disables wavelet filtering.}
#' }
#'
#' @param filter_name Character scalar. Filter family to relax. Must be one of
#'   \code{"mean"}, \code{"sd"}, \code{"tail_low"}, \code{"tail_high"},
#'   \code{"wavelet"}.
#' @param bounds_env Environment. Environment containing the current bounds using
#'   snake case keys. Updated in place.
#' @param wavelet_active_env Environment. Environment containing a logical
#'   \code{wavelet_active} flag. May be updated in place.
#' @param recompute_tailmass_fn Function. Callback used to recompute tail mass
#'   metrics when tail quantile thresholds change.
#'
#' @return Named list with:
#'   \code{changed} logical scalar indicating whether a bound was changed,
#'   \code{msg} character message describing the applied change.
#' @keywords internal
#' @export
relax_bounds_one_filter <- function(filter_name, bounds_env, wavelet_active_env,
                                    recompute_tailmass_fn) {

  msg <- NULL
  changed <- FALSE
  b <- bounds_env

  if (filter_name == "mean") {
    old <- b$mean
    b$mean <- min(b$mean * b$relax_mult, b$relax_mean_max)
    changed <- (b$mean > old + 1e-15)
    msg <- sprintf("mean tol %.4f -> %.4f", old, b$mean)
  }

  if (filter_name == "sd") {
    old <- b$sd
    b$sd <- min(b$sd * b$relax_mult, b$relax_sd_max)
    changed <- (b$sd > old + 1e-15)
    msg <- sprintf("sd tol %.4f -> %.4f", old, b$sd)
  }

  if (filter_name == "tail_low") {
    if (b$tail_tol_log < b$relax_tail_tol_log_max - 1e-15) {
      old <- b$tail_tol_log
      b$tail_tol_log <- min(b$tail_tol_log * b$relax_mult, b$relax_tail_tol_log_max)
      changed <- TRUE
      msg <- sprintf("tail_tol_log %.4f -> %.4f", old, b$tail_tol_log)
    } else if (b$tail_low_p < b$relax_tail_p_low_max - 1e-15) {
      old <- b$tail_low_p
      b$tail_low_p <- min(b$tail_low_p + b$relax_tail_p_step, b$relax_tail_p_low_max)
      recompute_tailmass_fn()
      changed <- TRUE
      msg <- sprintf("tail_low_p %.2f -> %.2f (recompute)", old, b$tail_low_p)
    } else {
      msg <- "tail_low at limit"
    }
  }

  if (filter_name == "tail_high") {
    if (b$tail_tol_log < b$relax_tail_tol_log_max - 1e-15) {
      old <- b$tail_tol_log
      b$tail_tol_log <- min(b$tail_tol_log * b$relax_mult, b$relax_tail_tol_log_max)
      changed <- TRUE
      msg <- sprintf("tail_tol_log %.4f -> %.4f", old, b$tail_tol_log)
    } else if (b$tail_high_p > b$relax_tail_p_high_min + 1e-15) {
      old <- b$tail_high_p
      b$tail_high_p <- max(b$tail_high_p - b$relax_tail_p_step, b$relax_tail_p_high_min)
      recompute_tailmass_fn()
      changed <- TRUE
      msg <- sprintf("tail_high_p %.2f -> %.2f (recompute)", old, b$tail_high_p)
    } else {
      msg <- "tail_high at limit"
    }
  }

  if (filter_name == "wavelet") {
    wavelet_active <- get("wavelet_active", envir = wavelet_active_env)

    if (!isTRUE(wavelet_active)) {
      msg <- "wavelet already inactive"

    } else if (b$spectral_cor_min > b$relax_spectral_cor_min + 1e-15) {
      # 1) Relax spectral correlation threshold
      old <- b$spectral_cor_min
      b$spectral_cor_min <- max(
        b$spectral_cor_min - b$relax_spectral_cor_step,
        b$relax_spectral_cor_min
      )
      changed <- TRUE
      msg <- sprintf("spectral_cor_min %.2f -> %.2f", old, b$spectral_cor_min)

    } else if (isTRUE(b$peak_match_enabled) &&
               b$peak_match_frac_min > b$relax_peak_match_frac_min + 1e-15) {
      # 2) Relax peak match fraction (only while enabled)
      old <- b$peak_match_frac_min
      b$peak_match_frac_min <- max(
        b$peak_match_frac_min - b$relax_peak_match_frac_step,
        b$relax_peak_match_frac_min
      )
      changed <- TRUE
      msg <- sprintf("peak_match_frac_min %.2f -> %.2f", old, b$peak_match_frac_min)

    } else if (isTRUE(b$peak_match_enabled)) {
      # 3) Disable peak matching
      b$peak_match_enabled <- FALSE
      changed <- TRUE
      msg <- "peak_match_enabled TRUE -> FALSE"

    } else {
      # 4) Disable wavelet filter entirely
      assign("wavelet_active", FALSE, envir = wavelet_active_env)
      changed <- TRUE
      msg <- "wavelet filter disabled"
    }
  }


  list(changed = changed, msg = msg)
}


# ==============================================================================
# Logging Helpers
# ==============================================================================

#' Create compact criteria text for logging
#'
#' @description
#' Builds a compact one line string describing the current criteria for a given
#' filter family. Used by the logging helpers in this script.
#'
#' @param filter_name Character scalar. Filter family name.
#' @param bounds List or environment. Bounds values used to build the criteria
#'   string.
#' @param tail_metrics List. Tail metrics produced by \code{compute_tailmass_metrics()}.
#' @param wavelet_active Logical scalar. TRUE if wavelet filtering is active.
#' @param spectral_diag List. Spectral diagnostics returned by
#'   \code{compute_spectral_metrics()}.
#'
#' @return Character scalar. A compact criteria string.
#' @keywords internal
#' @export
criteria_string_compact <- function(filter_name, bounds, tail_metrics,
                                    wavelet_active, spectral_diag) {

  if (filter_name == "mean") return(sprintf("tol = %.4f", bounds$mean))
  if (filter_name == "sd")   return(sprintf("tol = %.4f", bounds$sd))

  if (filter_name == "tail_low") {
    return(sprintf("p=%.2f, log.tol=%.4f", bounds$tail_low_p, bounds$tail_tol_log))
  }

  if (filter_name == "tail_high") {
    return(sprintf("p=%.2f, log.tol=%.4f", bounds$tail_high_p, bounds$tail_tol_log))
  }

  if (filter_name == "wavelet") {
    if (!isTRUE(wavelet_active)) return("inactive")

    peak_str <- if (isTRUE(bounds$peak_match_enabled)) {
      sprintf(", pk>=%.0f%%", bounds$peak_match_frac_min * 100)
    } else {
      ""
    }

    if (filter_name == "wavelet") {
      if (!isTRUE(wavelet_active)) return("inactive")

      peak_str <- if (isTRUE(bounds$peak_match_enabled)) {
        sprintf(", pk>=%.0f%%", bounds$peak_match_frac_min * 100)
      } else {
        ""
      }

      return(sprintf("cor>=%.2f%s", bounds$spectral_cor_min, peak_str))
    }


  }

  "NA"
}


#' Log filtering setup information
#'
#' @description
#' Writes a header block describing the filtering configuration. This includes
#' series lengths, the number of candidate realizations, the selection target,
#' and the relaxation priority.
#'
#' @param n_obs Integer scalar. Length of the observed series in years.
#' @param n_sim Integer scalar. Length of the simulated series in years.
#' @param n_realizations Integer scalar. Number of candidate realizations.
#' @param sample_target Integer scalar. Number of realizations requested.
#' @param relax_priority Character vector. Relaxation priority ordering.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
log_filtering_start <- function(n_obs, n_sim, n_realizations, sample_target, relax_priority) {

  .log("[FILTER] ===========================================================================")
  .log("[FILTER] FILTERING SETUP")
  .log("[FILTER] ===========================================================================")
  .log("[FILTER] Observed series: {format(n_obs, big.mark = ',')} years")
  .log("[FILTER] Simulated series: {format(n_sim, big.mark = ',')} years x {format(n_realizations, big.mark = ',')} candidate realizations")
  .log("[FILTER] Target: select {format(sample_target, big.mark = ',')} realizations from pool")
  .log("[FILTER] Relaxation priority: {paste(relax_priority, collapse = ' > ')}")
  .log("[FILTER] Filters relax left-to-right: {relax_priority[1]} relaxes FIRST, {relax_priority[length(relax_priority)]} relaxes LAST")
  .log("[FILTER] ===========================================================================")

  invisible(NULL)
}


#' Log one filtering iteration summary
#'
#' @description
#' Writes a compact per iteration summary showing pass counts, pass rates, and
#' current criteria for each active filter family, followed by the current pool
#' size relative to the target.
#'
#' @param iter Integer scalar. Iteration number. Use 0 for the initial evaluation.
#' @param passes Named list. Logical vectors indicating pass or fail for each
#'   filter family.
#' @param pool Integer vector. Indices of realizations currently in the pool.
#' @param n_total Integer scalar. Total number of realizations evaluated.
#' @param target Integer scalar. Target pool size.
#' @param bounds Environment. Current bounds values.
#' @param tail_metrics List. Tail metrics produced by \code{compute_tailmass_metrics()}.
#' @param wavelet_active Logical scalar. TRUE if wavelet filtering is active.
#' @param spectral_diag List. Spectral diagnostics returned by
#'   \code{compute_spectral_metrics()}.
#' @param note Character scalar or NULL. Optional message describing the action
#'   taken in this iteration.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
log_filter_iteration <- function(iter, passes, pool, n_total, target, bounds,
                                 tail_metrics, wavelet_active, spectral_diag,
                                 note = NULL) {

  active <- names(passes)[!vapply(passes, is.null, logical(1))]
  pool_size <- length(pool)
  pool_pct <- if (n_total > 0) 100 * pool_size / n_total else 0

  if (iter == 0L) {
    .log("[FILTER] ---------------------------------------------------------------------------")
    .log("[FILTER] ITERATION {format(iter, big.mark = ',')} | Initial evaluation")
  } else {
    .log("[FILTER] ---------------------------------------------------------------------------")
    if (!is.null(note)) {
      .log("[FILTER] ITERATION {format(iter, big.mark = ',')} | {note}")
    } else {
      .log("[FILTER] ITERATION {format(iter, big.mark = ',')}")
    }
  }
  .log("[FILTER] ---------------------------------------------------------------------------")

  filter_order <- c("mean", "sd", "tail_low", "tail_high", "wavelet")
  show_filters <- intersect(filter_order, active)

  if (length(show_filters) > 0) {
    .log("[FILTER] {sprintf('%-12s %10s %8s  %-40s', 'Filter', 'Passed', 'Rate', 'Criteria')}")
    .log("[FILTER] ---------------------------------------------------------------------------")

    for (nm in show_filters) {
      n_pass <- sum(passes[[nm]])
      rate <- 100 * mean(passes[[nm]])
      crit <- criteria_string_compact(nm, bounds, tail_metrics, wavelet_active, spectral_diag)
      .log("[FILTER] {sprintf('%-12s %10s %7.1f%%  %-40s', nm, format(n_pass, big.mark = ','), rate, crit)}")
    }
    .log("[FILTER] ---------------------------------------------------------------------------")
  }

  status_icon <- if (pool_size >= target) "[OK]" else "[>>]"
  status_txt  <- if (pool_size >= target) "TARGET REACHED" else "Need more candidates"

  .log("[FILTER] {status_icon} Pool: {format(pool_size, big.mark = ',')} / {format(n_total, big.mark = ',')} ({sprintf('%.1f', pool_pct)}%) | Need: {format(target, big.mark = ',')} | Status: {status_txt}")

  invisible(NULL)
}


#' Log final filtering summary
#'
#' @description
#' Writes a footer block summarizing the final pool size, number sampled, and the
#' selection mode that was used.
#'
#' @param pool_size Integer scalar. Final candidate pool size.
#' @param n_total Integer scalar. Total number of candidate realizations.
#' @param n_sampled Integer scalar. Number of realizations sampled.
#' @param relaxation_level Character scalar. Label describing the selection mode.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
log_final_summary <- function(pool_size, n_total, n_sampled, relaxation_level) {

  pct <- if (n_total > 0) 100 * pool_size / n_total else NA_real_

  .log("[FILTER] ===========================================================================")
  .log("[FILTER] FILTERING COMPLETE")
  .log("[FILTER] ===========================================================================")
  .log("[FILTER] Final pool: {format(pool_size, big.mark = ',')} / {format(n_total, big.mark = ',')} ({sprintf('%.1f', pct)}%)")
  .log("[FILTER] Sampled: {format(n_sampled, big.mark = ',')} realizations")
  .log("[FILTER] Selection mode: {relaxation_level}")
  .log("[FILTER] ===========================================================================")

  invisible(NULL)
}
