#' Filter and sample WARM realizations using distributional, tail, and spectral criteria
#'
#' @description
#' Filters an ensemble of WARM-generated annual realizations against an observed
#' annual series using three criterion families:
#' \itemize{
#'   \item \strong{Distributional}: relative differences in mean and standard deviation.
#'   \item \strong{Tail behaviour}: lower/upper tail mass relative to observed quantile thresholds.
#'   \item \strong{Spectral matching}: direct comparison of simulated vs observed global
#'     wavelet spectra, focusing on reproducing low-frequency structure and dominant peaks.
#' }
#'
#' Spectral matching criteria include:
#' \itemize{
#'   \item \strong{Spectral correlation}: correlation between log-transformed GWS of
#'     simulated and observed series.
#'   \item \strong{Low-frequency fraction}: relative difference in variance fraction
#'     at periods beyond a threshold tied to MODWT decomposition levels.
#'   \item \strong{Peak matching}: verification that dominant peaks in observed spectrum
#'     are reproduced in simulated spectrum within tolerance.
#' }
#'
#' If fewer than \code{n_select} realizations pass, the function relaxes criteria
#' iteratively (up to \code{filter_bounds$relax_max_iter}) by loosening the currently
#' most restrictive active filter (lowest pass rate). If still insufficient, a
#' deterministic fallback returns the \code{n_select} realizations with the smallest
#' absolute relative mean difference.
#'
#' @param obs_series Numeric vector. Observed annual series used as the reference.
#' @param sim_series Numeric matrix. Simulated annual realizations with years in rows
#'   and realizations in columns.
#' @param n_select Integer scalar. Number of realizations to return in \code{selected}.
#' @param seed Optional integer scalar. Random seed used for window selection (if lengths differ)
#'   and for sampling from the final candidate pool.
#' @param relax_order Character vector. Relaxation priority ordering for criteria.
#'   Must contain each of \code{c("mean","sd","tail_low","tail_high","wavelet")} exactly once.
#' @param filter_bounds Named list. Filtering thresholds and relaxation controls. Any entry
#'   overrides internal defaults. Uses snake_case keys.
#' @param wavelet_args Named list passed to \code{\link{analyze_wavelet_spectrum}} for
#'   observed and simulated series (e.g., \code{signif_level}, \code{noise_type},
#'   \code{period_lower_limit}, \code{detrend}).
#' @param modwt_n_levels Integer or NULL. Number of MODWT levels used in WARM simulation.
#'   Used to determine the low-frequency threshold period as 2^J. If NULL, estimated from
#'   series length.
#' @param make_plots Logical scalar. If \code{TRUE}, returns diagnostic plots in \code{plots}.
#' @param verbose Logical scalar. If \code{TRUE}, logs per-iteration pass rates and relaxation steps.
#'
#' @return A list with:
#' \describe{
#'   \item{pool}{Numeric matrix. Final candidate pool (subset of columns from \code{sim_series}).}
#'   \item{selected}{Numeric matrix. \code{n_select} realizations selected from the final pool.}
#'   \item{summary}{Data frame summarising pass counts/rates and selection mode.}
#'   \item{diagnostics}{List with window metadata, indices, relaxation log, spectral metrics,
#'     and final bounds.}
#'   \item{plots}{NULL or a named list of ggplot objects when \code{make_plots = TRUE}.}
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
    n_top_peaks = b$n_top_peaks,
    peak_prominence_frac = b$peak_prominence_frac,
    peak_period_tol = b$peak_period_tol,
    eps = b$spectral_eps
  )

  wavelet_active <- spectral_results$active
  spectral_metrics <- spectral_results$metrics
  spectral_diag <- spectral_results$diagnostics

  period <- spectral_results$period
  gws_obs <- spectral_results$gws_obs
  gws_signif <- spectral_results$gws_signif
  gws_cache <- spectral_results$gws_cache

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

      # Primary criterion 1: spectral correlation
      pass_cor <- is.finite(spectral_metrics$spectral_cor) &
        (spectral_metrics$spectral_cor >= b$spectral_cor_min)

      # Primary criterion 2: low-frequency fraction difference
      pass_lf <- is.finite(spectral_metrics$lf_frac_diff) &
        (spectral_metrics$lf_frac_diff <= b$lf_frac_tol)

      # Criterion 3: peak matching
      if (isTRUE(b$peak_match_enabled) && !is.null(spectral_metrics$peak_match_frac)) {
        pass_peak <- is.finite(spectral_metrics$peak_match_frac) &
          (spectral_metrics$peak_match_frac >= b$peak_match_frac_min)
      } else {
        pass_peak <- rep(TRUE, n_rlz)
      }

      pass_wavelet <- pass_cor & pass_lf & pass_peak
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

  idx_select <- if (length(pool_idx) == n_select) pool_idx else sample(pool_idx, size = n_select, replace = FALSE)

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
        gws_cache     = gws_cache,
        wavelet_q     = b$plot_wavelet_q
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  list(
    pool = sim_series[, pool_idx, drop = FALSE],
    selected = sim_series[, idx_select, drop = FALSE],
    summary = summary,
    diagnostics = list(
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
    ),
    plots = plots_out
  )
}


# ==============================================================================
# Default Bounds
# ==============================================================================

#' WARM filtering default bounds
#'
#' @description
#' Internal defaults for filter_warm_pool() bounds. Users should usually override
#' only a few entries via filter_bounds = list(...).
#'
#' @return Named list of defaults (snake_case keys).
#' @keywords internal
#' @export
filter_warm_bounds_defaults <- function() {
  list(
    # --- distributional tolerances (relative diff) ---
    mean = 0.03,
    sd   = 0.03,

    # --- tail behaviour (quantile-defined tails + log-distance tol) ---
    tail_low_p   = 0.20,
    tail_high_p  = 0.80,
    tail_tol_log = log(1.03),
    tail_eps     = 1e-5,

    # --- spectral matching (primary criteria) ---
    spectral_cor_min = 0.60,  # Min correlation of log-GWS
    lf_frac_tol      = 0.25, # Max relative diff in LF variance fraction
    spectral_eps     = 1e-10,

    # --- peak matching ---
    peak_match_enabled   = TRUE,
    peak_match_frac_min  = 1.0,  # At least 50% of observed peaks must match
    n_top_peaks          = 2L,    # Number of prominent peaks to detect
    peak_prominence_frac = 0.10,  # Min prominence as fraction of max power
    peak_period_tol      = 0.50,   # Period tolerance in log2 scale (octaves)

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
    relax_lf_frac_tol_step      = 0.10,
    relax_lf_frac_tol_max       = 0.60,
    relax_peak_match_frac_step  = 0.10,
    relax_peak_match_frac_min   = 0.00,

    relax_max_iter = 20L
  )
}


# ==============================================================================
# Spectral Matching Metrics
# ==============================================================================

#' Find local maxima (peaks) in a numeric vector
#'
#' @description
#' Identifies indices of local maxima in a vector.
#'
#' @param x Numeric vector.
#' @param strict Logical. If TRUE, requires strict inequality on both sides.
#'
#' @return Integer vector of peak indices.
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


#' Identify prominent peaks in observed GWS
#'
#' @description
#' Finds the top N peaks in the observed global wavelet spectrum that exceed
#' a prominence threshold (relative to the mean power).
#'
#' @param gws Numeric vector. Global wavelet spectrum.
#' @param period Numeric vector. Corresponding periods.
#' @param n_top Integer. Maximum number of peaks to return.
#' @param prominence_frac Numeric. Minimum prominence as fraction of max power.
#'
#' @return Data frame with columns: idx, period, power, prominence.
#' @keywords internal
identify_prominent_peaks <- function(gws, period, n_top = 3, prominence_frac = 0.10) {

  if (length(gws) < 3 || all(!is.finite(gws))) {
    return(data.frame(idx = integer(0), period = numeric(0),
                      power = numeric(0), prominence = numeric(0)))
  }

  gws_clean <- fill_nearest(gws)
  peak_idx <- find_local_maxima(gws_clean, strict = FALSE)

  if (length(peak_idx) == 0) {
    # No local maxima found - use global maximum
    max_idx <- which.max(gws_clean)
    return(data.frame(idx = max_idx, period = period[max_idx],
                      power = gws_clean[max_idx], prominence = 1.0))
  }

  # Compute prominence for each peak
  mean_power <- mean(gws_clean, na.rm = TRUE)
  max_power <- max(gws_clean, na.rm = TRUE)
  min_threshold <- prominence_frac * max_power

  prominences <- numeric(length(peak_idx))
  for (i in seq_along(peak_idx)) {
    pi <- peak_idx[i]
    peak_val <- gws_clean[pi]

    # Find lowest point between this peak and neighbors
    left_min <- if (pi > 1) min(gws_clean[1:(pi - 1)], na.rm = TRUE) else peak_val
    right_min <- if (pi < length(gws_clean)) min(gws_clean[(pi + 1):length(gws_clean)], na.rm = TRUE) else peak_val

    # Prominence is height above the higher of the two valleys
    base <- max(left_min, right_min)
    prominences[i] <- peak_val - base
  }

  # Filter by prominence and sort by power
  keep <- prominences >= min_threshold & gws_clean[peak_idx] >= min_threshold
  if (!any(keep)) {
    # Keep at least the highest peak
    keep[which.max(gws_clean[peak_idx])] <- TRUE
  }

  peak_idx <- peak_idx[keep]
  prominences <- prominences[keep]

  # Sort by power (descending) and take top N
  ord <- order(gws_clean[peak_idx], decreasing = TRUE)
  peak_idx <- peak_idx[ord]
  prominences <- prominences[ord]

  n_keep <- min(n_top, length(peak_idx))
  peak_idx <- peak_idx[seq_len(n_keep)]
  prominences <- prominences[seq_len(n_keep)]

  data.frame(
    idx = peak_idx,
    period = period[peak_idx],
    power = gws_clean[peak_idx],
    prominence = prominences
  )
}


#' Check if simulated spectrum matches observed peaks
#'
#' @description
#' For each observed peak, checks if the simulated spectrum has elevated power
#' at or near that period.
#'
#' @param gws_sim Numeric vector. Simulated GWS.
#' @param gws_obs Numeric vector. Observed GWS.
#' @param period Numeric vector. Periods.
#' @param obs_peaks Data frame. Observed peak information from identify_prominent_peaks.
#' @param period_tol Numeric. Relative tolerance for period matching (log2 scale).
#'
#' @return Numeric. Fraction of observed peaks matched (0 to 1).
#' @keywords internal
compute_peak_match_fraction <- function(gws_sim, gws_obs, period, obs_peaks, period_tol = 0.5) {

  if (nrow(obs_peaks) == 0) return(1.0)

  gws_sim_clean <- fill_nearest(gws_sim)
  n_matched <- 0L

  for (i in seq_len(nrow(obs_peaks))) {
    obs_period <- obs_peaks$period[i]
    obs_power <- obs_peaks$power[i]

    # Find periods within tolerance (in log2 scale)
    log2_diff <- abs(log2(period) - log2(obs_period))
    within_tol <- which(log2_diff <= period_tol)

    if (length(within_tol) > 0) {
      # Check if simulated power in this region is reasonably elevated
      sim_power_region <- max(gws_sim_clean[within_tol], na.rm = TRUE)

      # Match if simulated power is at least 50% of observed peak power
      # OR if it's a local maximum in the region
      if (sim_power_region >= 0.5 * obs_power) {
        n_matched <- n_matched + 1L
      } else {
        # Check if there's a local maximum in the region
        region_max_idx <- within_tol[which.max(gws_sim_clean[within_tol])]
        neighbors <- c(region_max_idx - 1L, region_max_idx + 1L)
        neighbors <- neighbors[neighbors >= 1 & neighbors <= length(gws_sim_clean)]

        if (length(neighbors) > 0 && all(gws_sim_clean[region_max_idx] >= gws_sim_clean[neighbors])) {
          n_matched <- n_matched + 1L
        }
      }
    }
  }

  n_matched / nrow(obs_peaks)
}


#' Compute spectral matching metrics for a single realization
#'
#' @description
#' Computes direct spectral matching metrics between observed and simulated
#' global wavelet spectra, focusing on reproducing observed low-frequency
#' structure and dominant peaks.
#'
#' @param gws_obs Numeric vector. Observed global wavelet spectrum.
#' @param gws_sim Numeric vector. Simulated global wavelet spectrum (same length).
#' @param period Numeric vector. Periods corresponding to GWS values.
#' @param lf_period_threshold Numeric. Period threshold defining low-frequency.
#' @param obs_peaks Data frame. Observed peaks from identify_prominent_peaks.
#' @param peak_period_tol Numeric. Period tolerance for peak matching.
#' @param eps Numeric. Small constant for log stability.
#'
#' @return Named list with spectral matching metrics.
#' @keywords internal
compute_spectral_match_single <- function(gws_obs, gws_sim, period,
                                          lf_period_threshold, obs_peaks,
                                          peak_period_tol = 0.5, eps = 1e-10) {

  gws_obs <- pmax(as.numeric(gws_obs), eps)
  gws_sim <- pmax(as.numeric(gws_sim), eps)

  # 1. Spectral shape correlation (log-space)
  log_obs <- log(gws_obs)
  log_sim <- log(gws_sim)

  ok <- is.finite(log_obs) & is.finite(log_sim)
  if (sum(ok) < 3L) {
    spectral_cor <- NA_real_
  } else {
    spectral_cor <- stats::cor(log_sim[ok], log_obs[ok])
  }

  # 2. Low-frequency variance fraction matching
  lf_idx <- which(period >= lf_period_threshold)
  total_obs <- sum(gws_obs, na.rm = TRUE)
  total_sim <- sum(gws_sim, na.rm = TRUE)

  if (length(lf_idx) > 0 && total_obs > eps && total_sim > eps) {
    lf_frac_obs <- sum(gws_obs[lf_idx], na.rm = TRUE) / total_obs
    lf_frac_sim <- sum(gws_sim[lf_idx], na.rm = TRUE) / total_sim
    lf_frac_diff <- abs(lf_frac_sim - lf_frac_obs) / max(lf_frac_obs, eps)
  } else {
    lf_frac_obs <- NA_real_
    lf_frac_sim <- NA_real_
    lf_frac_diff <- NA_real_
  }

  # 3. Peak matching
  peak_match_frac <- compute_peak_match_fraction(
    gws_sim = gws_sim,
    gws_obs = gws_obs,
    period = period,
    obs_peaks = obs_peaks,
    period_tol = peak_period_tol
  )

  list(
    spectral_cor = spectral_cor,
    lf_frac_obs = lf_frac_obs,
    lf_frac_sim = lf_frac_sim,
    lf_frac_diff = lf_frac_diff,
    peak_match_frac = peak_match_frac
  )
}


#' Compute spectral metrics for all realizations
#'
#' @description
#' Performs wavelet analysis on observed series and all simulated realizations,
#' then computes spectral matching metrics for filtering.
#'
#' @param obs_use Numeric vector. Observed values.
#' @param sim_series_stats Numeric matrix. Simulated values (n_use x n_realizations).
#' @param wavelet_pars List of wavelet parameters (signif_level, noise_type, etc.).
#' @param modwt_n_levels Integer or NULL. MODWT levels from WARM. If NULL, estimated.
#' @param n_top_peaks Integer. Number of top peaks to match.
#' @param peak_prominence_frac Numeric. Minimum prominence for peak detection.
#' @param peak_period_tol Numeric. Period tolerance for peak matching.
#' @param eps Numeric. Small constant for stability.
#'
#' @return List with spectral filter diagnostics, metrics, and cached spectra.
#' @keywords internal
#' @export
compute_spectral_metrics <- function(obs_use, sim_series_stats, wavelet_pars,
                                     modwt_n_levels = NULL,
                                     n_top_peaks = 3L,
                                     peak_prominence_frac = 0.10,
                                     peak_period_tol = 0.50,
                                     eps = 1e-10) {

  n_use <- length(obs_use)
  n_realizations <- ncol(sim_series_stats)

  # Estimate MODWT levels if not provided
  if (is.null(modwt_n_levels) || !is.finite(modwt_n_levels) || modwt_n_levels < 1) {
    modwt_n_levels <- max(2L, floor(log2(n_use)) - 1L)
  }
  modwt_n_levels <- as.integer(modwt_n_levels)

  # Low-frequency threshold: periods >= 2^J (smooth component threshold)
  lf_period_threshold <- 2^modwt_n_levels

  # Observed wavelet analysis
  wv_obs <- analyze_wavelet_spectrum(
    obs_use,
    signif = wavelet_pars$signif_level,
    noise = wavelet_pars$noise_type,
    min_period = wavelet_pars$period_lower_limit,
    detrend = isTRUE(wavelet_pars$detrend),
    mode = "fast"
  )

  if (is.null(wv_obs$period) || !is.numeric(wv_obs$period)) {
    stop("analyze_wavelet_spectrum(obs) must return numeric $period.", call. = FALSE)
  }

  period <- as.numeric(wv_obs$period)

  # Use unmasked GWS for matching
  gws_obs <- if (!is.null(wv_obs$gws_unmasked) && is.numeric(wv_obs$gws_unmasked)) {
    as.numeric(wv_obs$gws_unmasked)
  } else {
    as.numeric(wv_obs$gws)
  }

  gws_signif <- wv_obs$gws_signif_unmasked
  if (is.null(gws_signif)) gws_signif <- wv_obs$gws_signif

  # Ensure matching lengths
  if (length(gws_obs) != length(period)) {
    gws_obs <- gws_regrid(wv_obs, period, use_unmasked = TRUE)
  }

  gws_obs <- fill_nearest(gws_obs)

  # Identify prominent peaks in observed spectrum
  obs_peaks <- identify_prominent_peaks(
    gws = gws_obs,
    period = period,
    n_top = n_top_peaks,
    prominence_frac = peak_prominence_frac
  )

  # Cache simulated GWS and compute metrics
  gws_cache <- matrix(NA_real_, nrow = length(period), ncol = n_realizations)

  spectral_cor    <- rep(NA_real_, n_realizations)
  lf_frac_obs_vec <- rep(NA_real_, n_realizations)
  lf_frac_sim_vec <- rep(NA_real_, n_realizations)
  lf_frac_diff    <- rep(NA_real_, n_realizations)
  peak_match_frac <- rep(NA_real_, n_realizations)

  for (j in seq_len(n_realizations)) {
    wv_sim <- analyze_wavelet_spectrum(
      sim_series_stats[, j],
      signif = wavelet_pars$signif_level,
      noise = wavelet_pars$noise_type,
      min_period = wavelet_pars$period_lower_limit,
      detrend = isTRUE(wavelet_pars$detrend),
      mode = "fast"
    )

    gws_sim <- gws_regrid(wv_sim, period, use_unmasked = TRUE)
    gws_sim <- fill_nearest(gws_sim)
    gws_cache[, j] <- gws_sim

    m <- compute_spectral_match_single(
      gws_obs = gws_obs,
      gws_sim = gws_sim,
      period = period,
      lf_period_threshold = lf_period_threshold,
      obs_peaks = obs_peaks,
      peak_period_tol = peak_period_tol,
      eps = eps
    )

    spectral_cor[j]    <- m$spectral_cor
    lf_frac_obs_vec[j] <- m$lf_frac_obs
    lf_frac_sim_vec[j] <- m$lf_frac_sim
    lf_frac_diff[j]    <- m$lf_frac_diff
    peak_match_frac[j] <- m$peak_match_frac
  }

  # Diagnostics
  spectral_diag <- list(
    modwt_n_levels = modwt_n_levels,
    lf_period_threshold = lf_period_threshold,
    n_periods = length(period),
    obs_peaks = obs_peaks,
    gws_obs_summary = c(min = min(gws_obs, na.rm = TRUE),
                        median = stats::median(gws_obs, na.rm = TRUE),
                        max = max(gws_obs, na.rm = TRUE))
  )

  list(
    active = TRUE,
    period = period,
    gws_obs = gws_obs,
    gws_signif = gws_signif,
    gws_cache = gws_cache,
    metrics = list(
      spectral_cor = spectral_cor,
      lf_frac_obs = lf_frac_obs_vec,
      lf_frac_sim = lf_frac_sim_vec,
      lf_frac_diff = lf_frac_diff,
      peak_match_frac = peak_match_frac
    ),
    diagnostics = spectral_diag
  )
}


# ==============================================================================
# Tail-mass Metrics
# ==============================================================================

#' Compute tail-mass metrics for filtering
#'
#' @description
#' Computes tail-mass metrics based on lower and upper tail quantile thresholds.
#' Uses robust scale estimation (IQR -> MAD -> SD fallback) and normalizes
#' tail deficit/excess masses by series length and scale.
#'
#' @param obs_use Numeric vector of observed values.
#' @param sim_series_stats Numeric matrix of simulated values (n_use x n_realizations).
#' @param tail_low_p Lower tail quantile probability (e.g., 0.10).
#' @param tail_high_p Upper tail quantile probability (e.g., 0.90).
#' @param tail_eps Epsilon for log transform to avoid log(0).
#'
#' @return List with tail thresholds, scale, masses, and log-distance metrics.
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

#' Relax bounds for one filter
#'
#' @description
#' Applies one relaxation step to a single filter. Updates bounds in place.
#'
#' @param filter_name Character. Name of filter to relax.
#' @param bounds_env Environment containing bounds (snake_case keys).
#' @param wavelet_active_env Environment containing wavelet_active flag.
#' @param recompute_tailmass_fn Function to recompute tail mass when thresholds change.
#'
#' @return List with changed (logical) and msg (character).
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
      # First relax: spectral correlation threshold
      old <- b$spectral_cor_min
      b$spectral_cor_min <- max(b$spectral_cor_min - b$relax_spectral_cor_step,
                                b$relax_spectral_cor_min)
      changed <- TRUE
      msg <- sprintf("spectral_cor_min %.2f -> %.2f", old, b$spectral_cor_min)

    } else if (b$lf_frac_tol < b$relax_lf_frac_tol_max - 1e-15) {
      # Second relax: LF fraction tolerance
      old <- b$lf_frac_tol
      b$lf_frac_tol <- min(b$lf_frac_tol + b$relax_lf_frac_tol_step,
                           b$relax_lf_frac_tol_max)
      changed <- TRUE
      msg <- sprintf("lf_frac_tol %.2f -> %.2f", old, b$lf_frac_tol)

    } else if (isTRUE(b$peak_match_enabled) &&
               b$peak_match_frac_min > b$relax_peak_match_frac_min + 1e-15) {
      # Third relax: peak match fraction
      old <- b$peak_match_frac_min
      b$peak_match_frac_min <- max(b$peak_match_frac_min - b$relax_peak_match_frac_step,
                                   b$relax_peak_match_frac_min)
      changed <- TRUE
      msg <- sprintf("peak_match_frac_min %.2f -> %.2f", old, b$peak_match_frac_min)

    } else if (isTRUE(b$peak_match_enabled)) {
      # Fourth relax: disable peak matching
      b$peak_match_enabled <- FALSE
      changed <- TRUE
      msg <- "peak_match_enabled TRUE -> FALSE"

    } else {
      # Final: disable wavelet filter entirely
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

#' Compact criteria string for a filter
#'
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

    return(sprintf("cor>=%.2f, lf<=%.2f%s",
                   bounds$spectral_cor_min, bounds$lf_frac_tol, peak_str))
  }

  "NA"
}


#' Log initial setup information
#'
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


#' Log filter iteration details
#'
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
