#' Filter and sample WARM realizations using distributional, tail, and wavelet criteria
#'
#' @description
#' Filters an ensemble of WARM-generated annual realizations against an observed
#' annual series using three criterion families:
#' \itemize{
#'   \item \strong{Distributional}: relative differences in mean and standard deviation.
#'   \item \strong{Tail behaviour}: lower/upper tail mass relative to observed quantile thresholds.
#'   \item \strong{Spectral}: observed-relevant global wavelet spectrum (GWS) filtering.
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
#' @param pad_periods Logical scalar. If \code{TRUE}, expands the observed significant-period
#'   band by one index on each side when checking simulated presence in the observed-relevant band.
#' @param relax_order Character vector. Relaxation priority ordering for criteria.
#'   Must contain each of \code{c("mean","sd","tail_low","tail_high","wavelet")} exactly once.
#' @param filter_bounds Named list. Filtering thresholds and relaxation controls. Any entry
#'   overrides internal defaults. Uses snake_case keys (e.g. \code{tail_low_p}, not \code{tail.low.p}).
#' @param wavelet_args Named list passed to \code{\link{analyze_wavelet_spectrum}} for
#'   observed and simulated series (e.g., \code{signif_level}, \code{noise_type}, \code{period_lower_limit}, \code{detrend}).
#' @param make_plots Logical scalar. If \code{TRUE}, returns diagnostic plots in \code{plots}.
#' @param verbose Logical scalar. If \code{TRUE}, logs per-iteration pass rates and relaxation steps.
#'
#' @return A list with:
#' \describe{
#'   \item{pool}{Numeric matrix. Final candidate pool (subset of columns from \code{sim_series}).}
#'   \item{selected}{Numeric matrix. \code{n_select} realizations selected from the final pool.}
#'   \item{summary}{Data frame summarising pass counts/rates and selection mode.}
#'   \item{diagnostics}{List with window metadata, indices, relaxation log, and final bounds.}
#'   \item{plots}{NULL or a named list of ggplot objects when \code{make_plots = TRUE}.}
#' }
#'
#' @importFrom utils modifyList
#' @importFrom stats acf median runif setNames
#' @export
filter_warm_pool <- function(
    obs_series = NULL,
    sim_series = NULL,
    n_select = 5,
    seed = NULL,
    pad_periods = TRUE,
    relax_order = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
    filter_bounds = list(),
    wavelet_args = list(
      signif_level = 0.80,
      noise_type = "red",
      period_lower_limit = 2,
      detrend = TRUE
    ),
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

  if (!is.logical(pad_periods) || length(pad_periods) != 1L) {
    stop("'pad_periods' must be TRUE/FALSE.", call. = FALSE)
  }
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
  if (verbose) {
    log_step(sprintf(
      "Computing mean, sd, tail mass for observed and %s simulated series",
      format(n_rlz, big.mark = ",")
    ))
  }

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
  # Wavelet metrics + caching (snake_case inside diagnostics)
  # ---------------------------------------------------------------------------
  if (verbose) {
    log_step(sprintf(
      "Computing wavelet spectra for the observed and %s simulated series",
      format(n_rlz, big.mark = ",")
    ))
  }

  wavelet_results <- compute_wavelet_metrics(
    obs_use = obs_use,
    sim_series_stats = sim_use,
    wavelet_pars = wavelet_args,
    padding = pad_periods,
    min_bg = b$wavelet_min_bg
  )

  wavelet_active <- wavelet_results$active
  wavelet_diag <- wavelet_results$diagnostics

  period <- wavelet_results$power_period
  obs_power <- wavelet_results$power_obs
  signif_unmasked <- wavelet_results$power_signif_unmasked

  p_sim_reg <- wavelet_results$p_sim_reg
  p_sim_bg <- wavelet_results$p_sim_bg
  presence_rpad <- wavelet_results$presence_rpad

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
      if (!is.null(wavelet_diag$regions) && length(wavelet_diag$regions) > 0 &&
          !is.null(p_sim_reg) && !is.null(p_sim_bg) &&
          !is.null(wavelet_diag$p_obs_reg) && length(wavelet_diag$p_obs_reg) == ncol(p_sim_reg)) {

        p_obs_reg <- as.numeric(wavelet_diag$p_obs_reg)
        p_obs_bg  <- max(sum(fill_nearest(as.numeric(wavelet_diag$power_obs))[wavelet_diag$bg_idx], na.rm = TRUE), b$wavelet_min_bg)

        c_obs_reg <- p_obs_reg / p_obs_bg
        C_sim_reg <- p_sim_reg / pmax(p_sim_bg, b$wavelet_min_bg)

        reg_power_ok <- abs(p_sim_reg - rep(p_obs_reg, each = n_rlz)) <= (b$wavelet_region_tol * rep(p_obs_reg, each = n_rlz))
        reg_power_ok <- matrix(reg_power_ok, nrow = n_rlz)

        reg_contrast_ok <- abs(C_sim_reg - matrix(rep(c_obs_reg, each = n_rlz), nrow = n_rlz)) <=
          (b$wavelet_contrast_tol * matrix(rep(c_obs_reg, each = n_rlz), nrow = n_rlz))

        reg_ok <- reg_power_ok & reg_contrast_ok
        frac_ok <- rowMeans(reg_ok)

        pass_wavelet <- frac_ok >= b$sig_frac

        if (isTRUE(b$wavelet_require_presence) && !is.null(presence_rpad)) {
          pass_wavelet <- pass_wavelet & as.logical(presence_rpad)
        }
      }
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
      wavelet_pars = wavelet_args,
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
        wavelet_pars = wavelet_args,
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
  if (verbose) {
    log_step(
      "Selecting series",
      sprintf(
        "selecting %s from pool of %s",
        format(n_select, big.mark = ","),
        format(length(pool_idx), big.mark = ",")
      )
    )
  }

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
  # plots
  # ---------------------------------------------------------------------------
  plots_out <- NULL
  if (make_plots) {
    if (verbose) {
      log_step("Creating diagnostic plots", "time-series, statistics, and wavelet spectra")
    }

    if (length(pool_idx) < 1L) {
      warning("make_plots=TRUE but final pool is empty; plots=NULL.", call. = FALSE)
    } else {

      plots_out <- plot_filter_diagnostics(
        obs_series   = obs_use,
        sim_series   = sim_use,
        pool         = idx_select,
        rel_diff_mean = mean_rel_diff,
        rel_diff_sd   = sd_rel_diff,
        tail_metrics  = tail_stats,
        power_period  = period,
        power_obs     = obs_power,
        power_signif  = signif_unmasked,
        gws_cache     = wavelet_results$gws_cache_unmasked,
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
      final_bounds = as.list(b)
    ),
    plots = plots_out
  )
}


# ==============================================================================
# Helper Functions (snake_case keys + args)
# ==============================================================================

#' WARM filtering default bounds
#'
#' @description
#' Internal defaults for filter_warm_pool() bounds. Users should usually override
#' only a few entries via filter_bounds = list(...).
#'
#' @return Named list of defaults (snake_case keys).
#' @keywords internal
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

    # --- wavelet (observed-relevant regions) ---
    sig_frac                = 0.60,
    wavelet_region_tol      = 0.50,
    wavelet_contrast_tol    = 0.30,
    wavelet_min_bg          = 1e-12,
    wavelet_require_presence = TRUE,
    wavelet_presence_frac   = NULL,

    # --- plotting diagnostics ---
    plot_wavelet_q = c(0.50, 0.95),

    # --- relaxation controls ---
    relax_mult = 1.25,
    relax_mean_max = 0.25,
    relax_sd_max   = 0.25,

    relax_tail_tol_log_max = log(2.0),
    relax_tail_p_step      = 0.02,
    relax_tail_p_low_max   = 0.40,
    relax_tail_p_high_min  = 0.40,

    relax_wavelet_sig_frac_step     = 0.05,
    relax_wavelet_sig_frac_min      = 0.30,
    relax_wavelet_region_tol_step   = 0.10,
    relax_wavelet_region_tol_max    = 1.00,
    relax_wavelet_contrast_tol_step = 0.10,
    relax_wavelet_contrast_tol_max  = 1.00,

    relax_max_iter = 20L
  )
}

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

#' Compute wavelet metrics for all realizations
#'
#' @description
#' performs wavelet analysis on observed series and all simulated realizations.
#' Identifies significant periods, computes regional power and contrast metrics,
#' and caches both masked (for filtering) and unmasked (for plotting) GWS.
#'
#' @param obs_use Numeric vector of observed values.
#' @param sim_series_stats Numeric matrix of simulated values.
#' @param wavelet_pars List of wavelet parameters (signif_level, noise_type, etc.).
#' @param padding Logical for period padding.
#' @param min_bg Minimum background power threshold.
#'
#' @return List with wavelet filter diagnostics and cached spectra.
#' @keywords internal
#' @export
compute_wavelet_metrics <- function(obs_use, sim_series_stats, wavelet_pars,
                                    padding, min_bg) {

  n_realizations <- ncol(sim_series_stats)

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
  if (is.null(wv_obs$gws) || !is.numeric(wv_obs$gws)) {
    stop("analyze_wavelet_spectrum(obs) must return numeric $gws.", call. = FALSE)
  }

  power_period <- as.numeric(wv_obs$period)

  # Use unmasked observed GWS when available
  power_obs <- if (!is.null(wv_obs$gws_unmasked) && is.numeric(wv_obs$gws_unmasked)) {
    as.numeric(wv_obs$gws_unmasked)
  } else {
    as.numeric(wv_obs$gws)
  }

  power_signif <- wv_obs$gws_signif
  power_signif_unmasked <- wv_obs$gws_signif_unmasked

  if (length(power_obs) != length(power_period)) {
    power_obs <- gws_regrid(wv_obs, power_period, use_unmasked = TRUE)
  }
  if (length(power_signif) != length(power_period)) {
    stop("Observed significance curve length does not match observed period.", call. = FALSE)
  }

  signif_grid <- fill_nearest(as.numeric(power_signif))
  gws_obs_grid <- fill_nearest(as.numeric(power_obs))

  R_obs <- gws_obs_grid / signif_grid
  periods_sig_core <- which(is.finite(R_obs) & (R_obs > 1))

  split_regions <- function(idx) {
    idx <- sort(unique(as.integer(idx)))
    if (length(idx) == 0) return(list())
    brk <- c(TRUE, diff(idx) != 1L)
    split(idx, cumsum(brk))
  }

  wavelet_active <- TRUE
  wavelet_diag <- list(
    power_period = power_period,
    power_obs = power_obs,
    power_signif = power_signif,
    periods_sig_core = periods_sig_core,
    regions = list(),
    p_obs_reg = numeric(0),
    c_obs_reg = numeric(0)
  )

  p_sim_reg <- NULL
  p_sim_bg  <- NULL
  presence_rpad <- NULL

  # ALWAYS cache
  gws_cache <- matrix(NA_real_, nrow = length(power_period), ncol = n_realizations)
  gws_cache_unmasked <- matrix(NA_real_, nrow = length(power_period), ncol = n_realizations)

  for (j in seq_len(n_realizations)) {
    wv <- analyze_wavelet_spectrum(
      sim_series_stats[, j],
      signif = wavelet_pars$signif_level,
      noise = wavelet_pars$noise_type,
      min_period = wavelet_pars$period_lower_limit,
      detrend = isTRUE(wavelet_pars$detrend),
      mode = "fast"
    )

    gws_cache[, j] <- gws_regrid(wv, power_period, use_unmasked = FALSE)
    gws_cache_unmasked[, j] <- gws_regrid(wv, power_period, use_unmasked = TRUE)
  }

  if (length(periods_sig_core) == 0) {
    wavelet_active <- FALSE
  } else {
    regions <- split_regions(periods_sig_core)
    n_regions <- length(regions)

    core_union <- sort(unique(unlist(regions, use.names = FALSE)))
    bg_idx <- setdiff(seq_along(power_period), core_union)

    p_obs_reg <- vapply(regions, function(ii) sum(gws_obs_grid[ii], na.rm = TRUE), FUN.VALUE = numeric(1))
    p_obs_bg <- max(sum(gws_obs_grid[bg_idx], na.rm = TRUE), min_bg)
    c_obs_reg <- p_obs_reg / p_obs_bg

    periods_sig_pad <- periods_sig_core
    if (isTRUE(padding)) {
      n_periods <- length(power_period)
      periods_sig_pad <- unique(sort(c(
        pmax(periods_sig_core - 1L, 1L),
        periods_sig_core,
        pmin(periods_sig_core + 1L, n_periods)
      )))
    }

    p_sim_reg <- matrix(NA_real_, nrow = n_realizations, ncol = n_regions)
    p_sim_bg  <- rep(NA_real_, n_realizations)
    presence_rpad <- rep(FALSE, n_realizations)

    for (j in seq_len(n_realizations)) {
      gj_masked <- gws_cache[, j]

      bg <- max(sum(gj_masked[bg_idx], na.rm = TRUE), min_bg)
      p_sim_bg[j] <- bg

      for (r in seq_len(n_regions)) {
        p_sim_reg[j, r] <- sum(gj_masked[regions[[r]]], na.rm = TRUE)
      }

      rpad <- gj_masked[periods_sig_pad] / signif_grid[periods_sig_pad]
      presence_rpad[j] <- any(is.finite(rpad) & (rpad > 1))
    }

    wavelet_diag$regions <- regions
    wavelet_diag$p_obs_reg <- p_obs_reg
    wavelet_diag$c_obs_reg <- c_obs_reg
    wavelet_diag$bg_idx <- bg_idx
    wavelet_diag$periods_sig_pad <- periods_sig_pad
  }

  list(
    active = wavelet_active,
    diagnostics = wavelet_diag,
    power_period = power_period,
    power_obs = power_obs,
    power_signif = power_signif,
    power_signif_unmasked = power_signif_unmasked,
    p_sim_reg = p_sim_reg,
    p_sim_bg = p_sim_bg,
    presence_rpad = presence_rpad,
    gws_cache = gws_cache,
    gws_cache_unmasked = gws_cache_unmasked
  )
}

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
    } else if (b$sig_frac > b$relax_wavelet_sig_frac_min + 1e-15) {
      old <- b$sig_frac
      b$sig_frac <- max(b$sig_frac - b$relax_wavelet_sig_frac_step, b$relax_wavelet_sig_frac_min)
      changed <- TRUE
      msg <- sprintf("sig_frac %.2f -> %.2f", old, b$sig_frac)
    } else if (b$wavelet_region_tol < b$relax_wavelet_region_tol_max - 1e-15) {
      old <- b$wavelet_region_tol
      b$wavelet_region_tol <- min(b$wavelet_region_tol + b$relax_wavelet_region_tol_step,
                                  b$relax_wavelet_region_tol_max)
      changed <- TRUE
      msg <- sprintf("wavelet_region_tol %.2f -> %.2f", old, b$wavelet_region_tol)
    } else if (b$wavelet_contrast_tol < b$relax_wavelet_contrast_tol_max - 1e-15) {
      old <- b$wavelet_contrast_tol
      b$wavelet_contrast_tol <- min(b$wavelet_contrast_tol + b$relax_wavelet_contrast_tol_step,
                                    b$relax_wavelet_contrast_tol_max)
      changed <- TRUE
      msg <- sprintf("wavelet_contrast_tol %.2f -> %.2f", old, b$wavelet_contrast_tol)
    } else if (isTRUE(b$wavelet_require_presence)) {
      b$wavelet_require_presence <- FALSE
      changed <- TRUE
      msg <- "wavelet_require_presence TRUE -> FALSE"
    } else {
      assign("wavelet_active", FALSE, envir = wavelet_active_env)
      changed <- TRUE
      msg <- "wavelet disabled"
    }
  }

  list(changed = changed, msg = msg)
}

# ==============================================================================
# Logging helpers: only changed criteria_string_compact() to snake_case bounds keys
# ==============================================================================

#' Compact criteria string for a filter
#'
#' @keywords internal
#' @export
criteria_string_compact <- function(filter_name, bounds, tail_metrics,
                                    wavelet_active, wavelet_pars) {

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
    return(sprintf("sig_frac >= %.2f", bounds$sig_frac))
  }

  "NA"
}


# ==============================================================================
# Logging helpers for filter_warm_pool()
# ==============================================================================

#' Log initial setup information
#'
#' @description
#' Internal helper. Displays general information at the start of filtering.
#'
#' @param n_obs Integer scalar. Number of observations in observed series.
#' @param n_sim Integer scalar. Number of years in simulated series used after windowing.
#' @param n_realizations Integer scalar. Number of candidate realizations.
#' @param sample_target Integer scalar. Target number to select.
#' @param relax_priority Character vector. Relaxation priority vector.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
log_filtering_start <- function(n_obs, n_sim, n_realizations, sample_target, relax_priority) {

  logger::log_info("[FILTERING] ===========================================================================")
  logger::log_info("[FILTERING] FILTERING SETUP")
  logger::log_info("[FILTERING] ===========================================================================")
  logger::log_info("[FILTERING] Observed series: {format(n_obs, big.mark = ',')} years")
  logger::log_info("[FILTERING] Simulated series: {format(n_sim, big.mark = ',')} years x {format(n_realizations, big.mark = ',')} candidate realizations")
  logger::log_info("[FILTERING] Target: select {format(sample_target, big.mark = ',')} realizations from pool")
  logger::log_info("[FILTERING] Relaxation priority: {paste(relax_priority, collapse = ' > ')}")
  logger::log_info("[FILTERING] Filters relax left-to-right: {relax_priority[1]} relaxes FIRST, {relax_priority[length(relax_priority)]} relaxes LAST")
  logger::log_info("[FILTERING] ===========================================================================")

  invisible(NULL)
}

#' Log major step progress
#'
#' @description
#' Internal helper. Displays progress message for major computational steps.
#'
#' @param step_name Character scalar. Name of the step.
#' @param details Optional character scalar. Details string.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
log_step <- function(step_name, details = NULL) {
  if (!is.null(details)) {
    logger::log_info("[FILTERING] {step_name} | {details}")
  } else {
    logger::log_info("[FILTERING] {step_name}")
  }
  invisible(NULL)
}

#' Log filter iteration details
#'
#' @description
#' Internal helper. Prints iteration diagnostics in a compact table-like format.
#'
#' @param iter Integer scalar. Iteration number.
#' @param passes Named list of logical vectors. Per-filter pass vectors.
#' @param pool Integer vector. Pool indices passing all active filters.
#' @param n_total Integer scalar. Total number of realizations.
#' @param target Integer scalar. Target pool size.
#' @param bounds Environment or list of bounds.
#' @param tail_metrics List. Tail metrics used for criteria display.
#' @param wavelet_active Logical scalar.
#' @param wavelet_pars List. Wavelet parameter list.
#' @param note Optional character scalar.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
log_filter_iteration <- function(iter, passes, pool, n_total, target, bounds,
                                 tail_metrics, wavelet_active, wavelet_pars,
                                 note = NULL) {

  active <- names(passes)[!vapply(passes, is.null, logical(1))]
  pool_size <- length(pool)
  pool_pct <- if (n_total > 0) 100 * pool_size / n_total else 0

  if (iter == 0L) {
    logger::log_info("[FILTERING] ---------------------------------------------------------------------------")
    logger::log_info("[FILTERING] ITERATION {format(iter, big.mark = ',')} | Initial evaluation")
  } else {
    logger::log_info("[FILTERING] ---------------------------------------------------------------------------")
    if (!is.null(note)) {
      logger::log_info("[FILTERING] ITERATION {format(iter, big.mark = ',')} | {note}")
    } else {
      logger::log_info("[FILTERING] ITERATION {format(iter, big.mark = ',')}")
    }
  }
  logger::log_info("[FILTERING] ---------------------------------------------------------------------------")

  filter_order <- c("mean", "sd", "tail_low", "tail_high", "wavelet")
  show_filters <- intersect(filter_order, active)

  if (length(show_filters) > 0) {
    logger::log_info("[FILTERING] {sprintf('%-12s %10s %8s  %-30s', 'Filter', 'Passed', 'Rate', 'Criteria')}")
    logger::log_info("[FILTERING] ---------------------------------------------------------------------------")

    for (nm in show_filters) {
      n_pass <- sum(passes[[nm]])
      rate <- 100 * mean(passes[[nm]])
      crit <- criteria_string_compact(nm, bounds, tail_metrics, wavelet_active, wavelet_pars)
      logger::log_info("[FILTERING] {sprintf('%-12s %10s %7.1f%%  %-30s', nm, format(n_pass, big.mark = ','), rate, crit)}")
    }
    logger::log_info("[FILTERING] ---------------------------------------------------------------------------")
  }

  status_icon <- if (pool_size >= target) "[OK]" else "[>>]"
  status_txt  <- if (pool_size >= target) "TARGET REACHED" else "Need more candidates"

  logger::log_info("[FILTERING] {status_icon} Pool: {format(pool_size, big.mark = ',')} / {format(n_total, big.mark = ',')} ({sprintf('%.1f', pool_pct)}%) | Need: {format(target, big.mark = ',')} | Status: {status_txt}")

  invisible(NULL)
}

#' Log final filtering summary
#'
#' @description
#' Internal helper. Displays final filtering results.
#'
#' @param pool_size Integer scalar. Final pool size.
#' @param n_total Integer scalar. Total realizations.
#' @param n_sampled Integer scalar. Number sampled.
#' @param relaxation_level Character scalar. Selection mode / relaxation level reached.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
log_final_summary <- function(pool_size, n_total, n_sampled, relaxation_level) {

  pct <- if (n_total > 0) 100 * pool_size / n_total else NA_real_

  logger::log_info("[FILTERING] ===========================================================================")
  logger::log_info("[FILTERING] FILTERING COMPLETE")
  logger::log_info("[FILTERING] ===========================================================================")
  logger::log_info("[FILTERING] Final pool: {format(pool_size, big.mark = ',')} / {format(n_total, big.mark = ',')} ({sprintf('%.1f', pct)}%)")
  logger::log_info("[FILTERING] Sampled: {format(n_sampled, big.mark = ',')} realizations")
  logger::log_info("[FILTERING] Selection mode: {relaxation_level}")
  logger::log_info("[FILTERING] ===========================================================================")

  invisible(NULL)
}
