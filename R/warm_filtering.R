#' Filter and sample WARM simulations using distributional, tail, and wavelet criteria
#'
#' @description
#' Filters a set of WARM realizations (simulated annual series) against an observed
#' annual series using three tiers of checks:
#' \itemize{
#'   \item \strong{Distributional}: relative differences in mean and standard deviation,
#'   \item \strong{Tail behavior}: lower/upper tail mass relative to observed quantile thresholds,
#'   \item \strong{Spectral}: an observed-relevant global wavelet spectrum (GWS) filter.
#' }
#'
#' If fewer than \code{sample.num} realizations pass, the function relaxes criteria
#' iteratively (up to \code{bounds$relax.max.iter}) by loosening the currently
#' most restrictive active filter. If still insufficient, a deterministic fallback
#' returns the \code{sample.num} realizations with the smallest absolute relative
#' mean difference.
#'
#' @param series.obs Numeric vector. Observed annual series (no missing values recommended).
#' @param series.sim Numeric matrix. Simulated annual series with years in rows and
#'   realizations in columns.
#' @param sample.num Integer scalar. Number of realizations to return in \code{sampled}.
#'   If greater than \code{ncol(series.sim)}, it is reduced to \code{ncol(series.sim)} with a warning.
#' @param seed Optional integer. Random seed used for window selection (if lengths differ)
#'   and for sampling from the final candidate pool.
#' @param padding Logical scalar. If \code{TRUE}, expands the observed significant-period
#'   set by one index on each side when assessing presence of observed-relevant signal
#'   in simulated spectra.
#' @param relax.order Character vector. Relaxation priority ordering for criteria.
#'   Must contain each of \code{c("mean","sd","tail_low","tail_high","wavelet")} exactly once.
#'   When multiple filters have identical pass rates, this order is used to break ties.
#' @param bounds Named list. Filtering thresholds and relaxation controls. Any entry
#'   provided overrides the defaults. Common entries include:
#'   \itemize{
#'     \item \code{mean}, \code{sd}: relative-difference tolerances for mean and sd.
#'     \item \code{tail.low.p}, \code{tail.high.p}: quantiles defining lower/upper tail thresholds.
#'     \item \code{tail.tol.log}: tolerance on log-distance between simulated and observed tail mass.
#'     \item \code{tail.eps}: epsilon added to tail mass before log transform for stability.
#'     \item \code{sig.frac}: minimum fraction of observed-relevant regions required to pass the wavelet filter.
#'     \item \code{wavelet.region.tol}, \code{wavelet.contrast.tol}: tolerances for regional power and contrast ratios.
#'     \item \code{wavelet.require_presence}, \code{wavelet.presence.frac}: controls for requiring signal presence.
#'     \item \code{relax.*}: relaxation step sizes and caps (including \code{relax.max.iter}).
#'   }
#' @param wavelet.pars Named list passed to \code{\link{wavelet_spectral_analysis}} for
#'   observed and simulated series. Typical entries include \code{signif.level},
#'   \code{noise.type}, \code{period.lower.limit}, and \code{detrend}.
#' @param make.plots Logical scalar. If \code{TRUE}, returns diagnostic plots (time series overlay,
#'   relative-difference summaries, and wavelet GWS diagnostics) in \code{plots}.
#' @param verbose Logical scalar. If \code{TRUE}, logs per-iteration pass rates and relaxation steps
#'   using \code{logger::.log_info()}.
#'
#' @details
#' \strong{Length harmonization}:
#' The evaluation is performed over \code{n_use = min(length(series.obs), nrow(series.sim))}.
#' If the observed series is longer than \code{n_use}, one contiguous window of length \code{n_use}
#' is sampled and used for all comparisons. If simulated realizations are longer than \code{n_use},
#' each realization is evaluated on its own sampled contiguous window of length \code{n_use}.
#' Returned series in \code{subsetted} and \code{sampled} always use the original rows from
#' \code{series.sim} (no trimming in the outputs).
#'
#' \strong{Tail-mass metrics}:
#' Lower and upper thresholds are defined by the observed quantiles at
#' \code{bounds$tail.low.p} and \code{bounds$tail.high.p}. Tail behavior is summarized as
#' normalized deficit/excess mass relative to those thresholds, and compared using a
#' stabilized log-distance with \code{bounds$tail.eps}.
#'
#' \strong{Wavelet filter (observed-relevant GWS)}:
#' The observed GWS is tested against its significance curve to identify significant periods,
#' which are grouped into contiguous regions. For each realization, regional integrated power
#' and region-to-background contrast are compared to the observed within the specified
#' tolerances, requiring at least \code{bounds$sig.frac} of regions to pass. Optional
#' presence checks ensure simulated power exceeds a minimum fraction of observed power in
#' the observed-relevant band.
#'
#' @return A list with:
#' \describe{
#'   \item{subsetted}{Numeric matrix. All realizations that remain in the final candidate pool
#'     (columns subset of \code{series.sim}).}
#'   \item{sampled}{Numeric matrix. \code{sample.num} realizations sampled from the final pool
#'     (columns subset of \code{series.sim}).}
#'   \item{filter_summary}{Data frame summarizing pass counts and pass rates by criterion, and the
#'     final selection mode (tiered relaxation vs fallback).}
#'   \item{diagnostics}{List containing window metadata, pool/sample indices, relaxation log, and
#'     final parameter values used for filtering.}
#'   \item{plots}{NULL or a named list of ggplot objects when \code{make.plots = TRUE}.}
#' }
#'
#' @seealso \code{\link{wavelet_spectral_analysis}}
#'
#' @examples
#' \dontrun{
#' out <- filter_warm_simulations(
#'   series.obs = obs,
#'   series.sim = sim_mat,
#'   sample.num = 5,
#'   seed = 123,
#'   make.plots = TRUE,
#'   verbose = TRUE
#' )
#' out$filter_summary
#' out$plots$wavelet_gws
#' }
#' @importFrom utils modifyList
#' @importFrom stats acf median runif setNames
#' @export
filter_warm_simulations <- function(series.obs = NULL,
                                    series.sim = NULL,
                                    sample.num = 5,
                                    seed = NULL,
                                    padding = TRUE,
                                    relax.order = c("wavelet", "sd", "tail_low", "tail_high", "mean"),
                                    bounds = list(),
                                    wavelet.pars = list(
                                      signif.level = 0.80,
                                      noise.type = "red",
                                      period.lower.limit = 2,
                                      detrend = TRUE),
                                    make.plots = FALSE,
                                    verbose = FALSE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (is.null(series.obs) || !is.numeric(series.obs) || !is.vector(series.obs)) {
    stop("'series.obs' must be a numeric vector.", call. = FALSE)
  }
  if (is.null(series.sim) || !is.matrix(series.sim) || !is.numeric(series.sim)) {
    stop("'series.sim' must be a numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(sample.num) || length(sample.num) != 1L || !is.finite(sample.num) || sample.num < 1) {
    stop("'sample.num' must be a positive integer.", call. = FALSE)
  }
  sample.num <- as.integer(sample.num)

  if (!is.logical(padding) || length(padding) != 1L) stop("'padding' must be TRUE/FALSE.", call. = FALSE)
  make.plots <- isTRUE(make.plots)
  verbose <- isTRUE(verbose)

  n_years_obs0 <- length(series.obs)
  n_years_sim0 <- nrow(series.sim)
  n_realizations <- ncol(series.sim)

  if (n_realizations < sample.num) {
    warning("sample.num exceeds number of series; returning at most ncol(series.sim).", call. = FALSE)
    sample.num <- n_realizations
  }

  # ---------------------------------------------------------------------------
  # Relaxation order (NEW: user-configurable)
  # ---------------------------------------------------------------------------
  allowed <- c("mean", "sd", "tail_low", "tail_high", "wavelet")
  if (is.null(relax.order) || !is.character(relax.order)) {
    stop("'relax.order' must be a character vector.", call. = FALSE)
  }
  relax.order <- as.character(relax.order)

  if (anyNA(relax.order) || any(!nzchar(relax.order))) {
    stop("'relax.order' must not contain NA/empty strings.", call. = FALSE)
  }
  if (any(duplicated(relax.order))) {
    stop("'relax.order' must not contain duplicates.", call. = FALSE)
  }
  if (!all(relax.order %in% allowed)) {
    bad <- setdiff(relax.order, allowed)
    stop("Invalid entries in 'relax.order': ", paste(bad, collapse = ", "),
         ". Allowed: ", paste(allowed, collapse = ", "), call. = FALSE)
  }
  if (!setequal(relax.order, allowed)) {
    missing <- setdiff(allowed, relax.order)
    extra   <- setdiff(relax.order, allowed)
    msg <- character(0)
    if (length(missing) > 0) msg <- c(msg, paste0("missing: ", paste(missing, collapse = ", ")))
    if (length(extra) > 0)   msg <- c(msg, paste0("extra: ", paste(extra, collapse = ", ")))
    stop("'relax.order' must contain exactly these filters once each: ",
         paste(allowed, collapse = ", "), ". Problem: ", paste(msg, collapse = " | "),
         call. = FALSE)
  }

  RELAX_ORDER <- relax.order

  # ---------------------------------------------------------------------------
  # Bounds defaults
  # ---------------------------------------------------------------------------
  b_list <- modifyList(list(
    mean = 0.05,
    sd   = 0.05,

    tail.low.p  = 0.20,
    tail.high.p = 0.80,
    tail.tol.log = log(1.05),
    tail.eps = 1e-5,

    sig.frac = 0.60,
    wavelet.region.tol = 0.50,
    wavelet.contrast.tol = 0.30,
    wavelet.min_bg = 1e-12,
    wavelet.require_presence = TRUE,
    wavelet.presence.frac = NULL,

    plot.wavelet_q = c(0.50, 0.95),

    relax.mult = 1.25,
    relax.mean.max = 0.25,
    relax.sd.max = 0.25,

    relax.tail.tol.log.max = log(2.0),
    relax.tail.p.step = 0.02,
    relax.tail.p.low.max  = 0.40,
    relax.tail.p.high.min = 0.40,

    relax.wavelet.sig.frac.step = 0.05,
    relax.wavelet.sig.frac.min = 0.30,
    relax.wavelet.region.tol.step = 0.10,
    relax.wavelet.region.tol.max = 1.00,
    relax.wavelet.contrast.tol.step = 0.10,
    relax.wavelet.contrast.tol.max = 1.00,

    relax.max.iter = 20L
  ), bounds)

  b <- list2env(b_list, parent = environment())
  b$relax.max.iter <- as.integer(b$relax.max.iter)

  # Validate tail parameters
  b$tail.low.p  <- as.numeric(b$tail.low.p)
  b$tail.high.p <- as.numeric(b$tail.high.p)
  b$tail.tol.log <- as.numeric(b$tail.tol.log)
  b$tail.eps <- as.numeric(b$tail.eps)

  if (!is.finite(b$tail.low.p) || b$tail.low.p <= 0 || b$tail.low.p >= 0.5) {
    stop("bounds$tail.low.p must be in (0, 0.5).", call. = FALSE)
  }
  if (!is.finite(b$tail.high.p) || b$tail.high.p <= 0.5 || b$tail.high.p >= 1) {
    stop("bounds$tail.high.p must be in (0.5, 1).", call. = FALSE)
  }
  if (!is.finite(b$tail.tol.log) || b$tail.tol.log <= 0) {
    stop("bounds$tail.tol.log must be a positive finite number.", call. = FALSE)
  }
  if (!is.finite(b$tail.eps) || b$tail.eps <= 0) {
    stop("bounds$tail.eps must be a positive finite number.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Display initial setup information
  # ---------------------------------------------------------------------------
  if (verbose) {
    log_filtering_setup(
      n_obs = n_years_obs0,
      n_sim = n_years_sim0,
      n_realizations = n_realizations,
      sample_num = sample.num,
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
    set.seed(seed)                    # Only set.seed() needed
  }

  # ---------------------------------------------------------------------------
  # Length harmonization
  # ---------------------------------------------------------------------------
  n_use <- min(n_years_obs0, n_years_sim0)

  if (n_years_obs0 > n_use) {
    max_start_obs <- n_years_obs0 - n_use + 1L
    start_obs <- sample.int(max_start_obs, size = 1L)
    end_obs <- start_obs + n_use - 1L
    obs.use <- series.obs[start_obs:end_obs]
    obs_window <- c(start = start_obs, end = end_obs)
  } else {
    obs.use <- series.obs
    obs_window <- NULL
  }

  if (n_years_sim0 > n_use) {
    max_start_sim <- n_years_sim0 - n_use + 1L
    starts <- sample.int(max_start_sim, size = n_realizations, replace = TRUE)
    ends <- starts + n_use - 1L
    window_index <- lapply(seq_len(n_realizations), function(j) c(start = starts[j], end = ends[j]))

    series_sim_for_stats <- vapply(
      seq_len(n_realizations),
      FUN = function(j) series.sim[starts[j]:ends[j], j],
      FUN.VALUE = numeric(n_use)
    )
    colnames(series_sim_for_stats) <- colnames(series.sim)
  } else {
    window_index <- NULL
    series_sim_for_stats <- series.sim
  }

  if (!is.matrix(series_sim_for_stats)) series_sim_for_stats <- as.matrix(series_sim_for_stats)
  if (nrow(series_sim_for_stats) != n_use) stop("Internal error: series_sim_for_stats row mismatch.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # Base statistics
  # ---------------------------------------------------------------------------

  if (verbose) {
    log_step("Computing ", sprintf("mean, sd, tail mass for observed and %d simulated series", n_realizations))
  }

  obs_mean <- mean(obs.use)
  obs_sd   <- stats::sd(obs.use)

  sim_means <- colMeans(series_sim_for_stats)
  sim_sds   <- apply(series_sim_for_stats, 2, stats::sd)

  eps0 <- 1e-10
  rel_diff_vec <- function(sim, obs) {
    if (!is.finite(obs) || abs(obs) < eps0) rep(0, length(sim)) else (sim - obs) / obs
  }
  rel.diff.mean <- rel_diff_vec(sim_means, obs_mean)
  rel.diff.sd   <- rel_diff_vec(sim_sds,   obs_sd)

  # ---------------------------------------------------------------------------
  # IMPROVEMENT 1: Use helper function for tail-mass metrics
  # ---------------------------------------------------------------------------
  tail_metrics <- compute_tailmass_metrics(
    obs.use = obs.use,
    series_sim_for_stats = series_sim_for_stats,
    tail.low.p = b$tail.low.p,
    tail.high.p = b$tail.high.p,
    tail.eps = b$tail.eps
  )

  # For relaxation, we need to recompute tail metrics when tail.low.p/tail.high.p change
  .recompute_tailmass <- function() {
    tail_metrics <<- compute_tailmass_metrics(
      obs.use = obs.use,
      series_sim_for_stats = series_sim_for_stats,
      tail.low.p = b$tail.low.p,
      tail.high.p = b$tail.high.p,
      tail.eps = b$tail.eps
    )
  }

  # ---------------------------------------------------------------------------
  # IMPROVEMENT 2: Use helper function for wavelet metrics + ALWAYS cache
  # ---------------------------------------------------------------------------
  if (verbose) {
    log_step("Computing wavelet spectra for ", sprintf("observed and %d simulated series", n_realizations))
  }

  wavelet_results <- compute_wavelet_metrics(
    obs.use = obs.use,
    series_sim_for_stats = series_sim_for_stats,
    wavelet.pars = wavelet.pars,
    padding = padding,
    min_bg = b$wavelet.min_bg
  )

  wavelet_active <- wavelet_results$active
  wavelet_diag <- wavelet_results$diagnostics
  power.period <- wavelet_results$power.period
  power.obs <- wavelet_results$power.obs
  power.signif <- wavelet_results$power.signif
  power.signif_unmasked <- wavelet_results$power.signif_unmasked
  P_sim_reg <- wavelet_results$P_sim_reg
  P_sim_bg <- wavelet_results$P_sim_bg
  presence_rpad <- wavelet_results$presence_rpad
  gws_cache_mat <- wavelet_results$gws_cache

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  .pool_from_pass <- function(pass_list) {
    ok <- rep(TRUE, n_realizations)
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

  # ---------------------------------------------------------------------------
  # Initial evaluation
  # ---------------------------------------------------------------------------

  # --- build pass vectors from currently-available diagnostics -----------------

  .pass_vectors_current <- function() {
    # mean / sd
    pass_mean <- abs(rel.diff.mean) <= b$mean
    pass_sd   <- abs(rel.diff.sd)   <= b$sd

    # tails: use the log-difference metrics already returned by compute_tailmass_metrics()
    pass_tail_low  <- is.finite(tail_metrics$logdiff_low)  & (tail_metrics$logdiff_low  <= b$tail.tol.log)
    pass_tail_high <- is.finite(tail_metrics$logdiff_high) & (tail_metrics$logdiff_high <= b$tail.tol.log)

    # wavelet
    pass_wavelet <- rep(TRUE, n_realizations)
    if (isTRUE(wavelet_active)) {
      # If your compute_wavelet_metrics() computed region metrics, use them.
      # Otherwise (no significance) wavelet_active should be FALSE already.
      if (!is.null(wavelet_diag$regions) && length(wavelet_diag$regions) > 0 &&
          !is.null(P_sim_reg) && !is.null(P_sim_bg) &&
          !is.null(wavelet_diag$P_obs_reg) && length(wavelet_diag$P_obs_reg) == ncol(P_sim_reg)) {

        P_obs_reg <- as.numeric(wavelet_diag$P_obs_reg)
        P_obs_bg  <- max(sum(fill_nearest(as.numeric(wavelet_diag$power.obs))[wavelet_diag$bg_idx], na.rm = TRUE), b$wavelet.min_bg)

        # observed contrast per region
        C_obs_reg <- P_obs_reg / P_obs_bg

        # sim contrast per region
        C_sim_reg <- P_sim_reg / pmax(P_sim_bg, b$wavelet.min_bg)

        # region-wise pass: within tol for power AND contrast
        # (use multiplicative bands; symmetric)
        reg_power_ok <- abs(P_sim_reg - rep(P_obs_reg, each = n_realizations)) <= (b$wavelet.region.tol * rep(P_obs_reg, each = n_realizations))
        reg_power_ok <- matrix(reg_power_ok, nrow = n_realizations)

        reg_contrast_ok <- abs(C_sim_reg - matrix(rep(C_obs_reg, each = n_realizations), nrow = n_realizations)) <=
          (b$wavelet.contrast.tol * matrix(rep(C_obs_reg, each = n_realizations), nrow = n_realizations))

        reg_ok <- reg_power_ok & reg_contrast_ok

        frac_ok <- rowMeans(reg_ok)
        pass_wavelet <- frac_ok >= b$sig.frac

        # optional presence requirement
        if (isTRUE(b$wavelet.require_presence) && !is.null(presence_rpad)) {
          pass_wavelet <- pass_wavelet & as.logical(presence_rpad)
        }
      } else {
        # if wavelet_active TRUE but no region metrics exist, fail-safe: don't block
        pass_wavelet <- rep(TRUE, n_realizations)
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
  pool <- .pool_from_pass(passes)




  if (verbose) {
    log_filter_iteration(
      iter = 0L,
      passes = passes,
      pool = pool,
      n_total = n_realizations,
      target = sample.num,
      bounds = b,
      tail_metrics = tail_metrics,
      wavelet_active = wavelet_active,
      wavelet_pars = wavelet.pars,
      note = "Initial evaluation."
    )
  }

  # ---------------------------------------------------------------------------
  # SIMPLIFIED relaxation loop (deterministic ordering, no tie-breaking)
  # ---------------------------------------------------------------------------
  #if (verbose) {
  #  log_step("Applying filter criteria",
  #           sprintf("evaluating %d realizations against tolerance bounds", n_realizations))
  #}

  relax_log <- character(0)
  relaxation_level <- "tiered"

  iter <- 0L
  while (length(pool) < sample.num && iter < b$relax.max.iter) {
    iter <- iter + 1L

    # Get active filters in priority order
    active <- names(passes)[!vapply(passes, is.null, logical(1))]
    active <- intersect(RELAX_ORDER, active)

    if (length(active) == 0) break

    # Pick filter with lowest pass rate (ties broken by RELAX_ORDER)
    rates <- vapply(active, function(nm) mean(passes[[nm]]), numeric(1))
    filter_to_relax <- active[which.min(rates)]

    before_pool <- length(pool)

    # Relax that filter
    rel <- .relax_one(filter_to_relax)
    if (!isTRUE(rel$changed)) break  # All relaxation exhausted

    # Recompute passes
    passes <- .pass_vectors_current()
    pool <- .pool_from_pass(passes)

    note <- sprintf("Relaxed '%s': %s | pool %d -> %d",
                    filter_to_relax, rel$msg, before_pool, length(pool))
    relax_log <- c(relax_log, sprintf("Iteration %d: %s", iter, note))

    if (verbose) {
      log_filter_iteration(
        iter = iter,
        passes = passes,
        pool = pool,
        n_total = n_realizations,
        target = sample.num,
        bounds = b,
        tail_metrics = tail_metrics,
        wavelet_active = wavelet_active,
        wavelet_pars = wavelet.pars,
        note = note
      )
    }
  }

  # Fallback if needed
  if (length(pool) < sample.num) {
    relaxation_level <- "closest_mean_fallback"
    pool <- order(abs(rel.diff.mean))[seq_len(sample.num)]
    if (verbose) {
      .log_info(sprintf(
        "Fallback activated: closest_mean to guarantee pool size = %d (from %d total).",
        sample.num, n_realizations
      ))
    }
  }

  # Sample from pool
  if (verbose) {
    log_step("Sampling series",
             sprintf("selecting %d from pool of %d", sample.num, length(pool)))
  }

  idx_sampled <- if (length(pool) == sample.num) pool else sample(pool, size = sample.num, replace = FALSE)

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

  filter_summary <- data.frame(
    filter = report_filters,
    n_passed = n_passed,
    pct_passed = pct_passed,
    relaxation_level = relaxation_level,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    log_final_summary(
      pool_size = length(pool),
      n_total = n_realizations,
      n_sampled = length(idx_sampled),
      relaxation_level = relaxation_level
    )
  }

  # ---------------------------------------------------------------------------
  # Helper function for plotting
  # ---------------------------------------------------------------------------

  plots_out <- NULL
  if (make.plots) {
    if (verbose) {
      log_step("Creating diagnostic plots", "time-series, statistics, and wavelet spectra")
    }

    if (length(pool) < 1L) {
      warning("make.plots=TRUE but final pool is empty; plots=NULL.", call. = FALSE)

      } else {

      plots_out <- plot_filter_diagnostics(
        obs.use = obs.use,
        series_sim_for_stats = series_sim_for_stats,
        pool = idx_sampled,
        rel.diff.mean = rel.diff.mean,
        rel.diff.sd = rel.diff.sd,
        tail_metrics = tail_metrics,
        power.period = power.period,
        power.obs = power.obs,
        power.signif = power.signif_unmasked,
        gws_cache_mat = wavelet_results$gws_cache_unmasked,
        wavelet_q = b$plot.wavelet_q
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  list(
    subsetted = series.sim[, pool, drop = FALSE],
    sampled = series.sim[, idx_sampled, drop = FALSE],
    filter_summary = filter_summary,
    diagnostics = list(
      n_years_obs_original = n_years_obs0,
      n_years_sim_original = n_years_sim0,
      n_use = n_use,
      n_periods = length(power.period),
      obs_window = obs_window,
      window_index = window_index,
      idx_pool = pool,
      idx_sampled = idx_sampled,
      relaxation_log = relax_log,
      final_params = as.list(b)
    ),
    plots = plots_out
  )
}




# ==============================================================================
# Helper Functions for filter_warm_simulations()
# ==============================================================================

#' WARM filtering default bounds
#'
#' @description
#' Internal defaults for filter_warm_simulations() bounds. Users should usually
#' override only a few entries via generate_weather_series(warm.bounds = list(...)).
#'
#' @return Named list of defaults.
#' @keywords internal
filter_warm_bounds_defaults <- function() {
  list(
    # --- distributional tolerances (relative diff) ---
    mean = 0.05,
    sd   = 0.05,

    # --- tail behaviour (quantile-defined tails + log-distance tol) ---
    tail.low.p   = 0.20,
    tail.high.p  = 0.80,
    tail.tol.log = log(1.05),
    tail.eps     = 1e-5,

    # --- wavelet (observed-relevant regions) ---
    sig.frac               = 0.60,
    wavelet.region.tol     = 0.50,
    wavelet.contrast.tol   = 0.30,
    wavelet.min_bg         = 1e-12,
    wavelet.require_presence = TRUE,
    wavelet.presence.frac  = NULL,

    # --- plotting diagnostics ---
    plot.wavelet_q = c(0.50, 0.95),

    # --- relaxation controls ---
    relax.mult = 1.25,
    relax.mean.max = 0.25,
    relax.sd.max   = 0.25,

    relax.tail.tol.log.max = log(2.0),
    relax.tail.p.step      = 0.02,
    relax.tail.p.low.max   = 0.40,
    relax.tail.p.high.min  = 0.40,

    relax.wavelet.sig.frac.step     = 0.05,
    relax.wavelet.sig.frac.min      = 0.30,
    relax.wavelet.region.tol.step   = 0.10,
    relax.wavelet.region.tol.max    = 1.00,
    relax.wavelet.contrast.tol.step = 0.10,
    relax.wavelet.contrast.tol.max  = 1.00,

    relax.max.iter = 20L
  )
}

#' Compute tail-mass metrics for filtering
#'
#' @description
#' Computes tail-mass metrics based on lower and upper tail quantile thresholds.
#' Uses robust scale estimation (IQR -> MAD -> SD fallback) and normalizes
#' tail deficit/excess masses by series length and scale.
#'
#' @param obs.use Numeric vector of observed values
#' @param series_sim_for_stats Numeric matrix of simulated values (n_use x n_realizations)
#' @param tail.low.p Lower tail quantile probability (e.g., 0.10)
#' @param tail.high.p Upper tail quantile probability (e.g., 0.90)
#' @param tail.eps Epsilon for log transform to avoid log(0)
#'
#' @return List with elements:
#' \itemize{
#'   \item thr_low, thr_high: Observed quantile thresholds
#'   \item scale_obs: Robust scale estimate
#'   \item M_obs_low, M_obs_high: Observed tail masses
#'   \item M_sim_low, M_sim_high: Simulated tail masses (vectors)
#'   \item logdiff_low, logdiff_high: Log-distance metrics (vectors)
#' }
#'
#' @keywords internal
#' @export
compute_tailmass_metrics <- function(obs.use, series_sim_for_stats,
                                     tail.low.p, tail.high.p, tail.eps) {

  n_use <- length(obs.use)
  n_realizations <- ncol(series_sim_for_stats)

  # Robust scale with fallback hierarchy: IQR -> MAD -> SD -> 1
  scale_obs <- stats::IQR(obs.use, na.rm = TRUE, type = 7)
  if (!is.finite(scale_obs) || scale_obs <= 0) {
    scale_obs <- stats::mad(obs.use, constant = 1, na.rm = TRUE)
  }
  if (!is.finite(scale_obs) || scale_obs <= 0) {
    scale_obs <- stats::sd(obs.use, na.rm = TRUE)
  }
  if (!is.finite(scale_obs) || scale_obs <= 0) scale_obs <- 1.0

  # Compute thresholds
  thr_low <- stats::quantile(obs.use, probs = tail.low.p, names = FALSE, type = 7, na.rm = TRUE)
  thr_high <- stats::quantile(obs.use, probs = tail.high.p, names = FALSE, type = 7, na.rm = TRUE)

  # Handle invalid cases
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

  denom <- n_use * scale_obs
  denom <- max(denom, 1e-12)

  # Lower-tail deficit mass (how much obs is below threshold vs sim)
  X_low <- thr_low - series_sim_for_stats
  S_sim_low <- colSums(pmax(X_low, 0), na.rm = TRUE)
  S_obs_low <- sum(pmax(thr_low - obs.use, 0), na.rm = TRUE)

  # Upper-tail excess mass (how much obs is above threshold vs sim)
  X_high <- series_sim_for_stats - thr_high
  S_sim_high <- colSums(pmax(X_high, 0), na.rm = TRUE)
  S_obs_high <- sum(pmax(obs.use - thr_high, 0), na.rm = TRUE)

  # Normalize by length and scale
  M_obs_low <- S_obs_low / denom
  M_obs_high <- S_obs_high / denom
  M_sim_low <- S_sim_low / denom
  M_sim_high <- S_sim_high / denom

  # Log differences for scale-invariant comparison
  logdiff_low <- abs(log(M_sim_low + tail.eps) - log(M_obs_low + tail.eps))
  logdiff_high <- abs(log(M_sim_high + tail.eps) - log(M_obs_high + tail.eps))

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
#' Performs wavelet analysis on observed series and all simulated realizations.
#' Identifies significant periods, computes regional power and contrast metrics,
#' and ALWAYS caches both masked (for filtering) and unmasked (for plotting) GWS.
#'
#' @param obs.use Numeric vector of observed values
#' @param series_sim_for_stats Numeric matrix of simulated values
#' @param wavelet.pars List of wavelet parameters (signif.level, noise.type, etc.)
#' @param padding Logical for period padding
#' @param min_bg Minimum background power threshold
#'
#' @return List with elements:
#' \itemize{
#'   \item active: Logical, whether wavelet filter is active
#'   \item diagnostics: Detailed wavelet diagnostics
#'   \item power.period, power.obs, power.signif: Period and power vectors
#'   \item P_sim_reg: Regional power matrix (n_realizations x n_regions)
#'   \item P_sim_bg: Background power vector (n_realizations)
#'   \item presence_rpad: Presence indicators (n_realizations)
#'   \item gws_cache: ALWAYS cached masked GWS matrix (n_periods x n_realizations) for filtering
#'   \item gws_cache_unmasked: ALWAYS cached unmasked GWS matrix (n_periods x n_realizations) for plotting
#' }
#'
#' @keywords internal
#' @export
compute_wavelet_metrics <- function(obs.use, series_sim_for_stats, wavelet.pars,
                                    padding, min_bg) {

  n_realizations <- ncol(series_sim_for_stats)

  # Observed wavelet analysis
  wv_obs <- wavelet_spectral_analysis(
    obs.use,
    signif.level = wavelet.pars$signif.level,
    noise.type = wavelet.pars$noise.type,
    period.lower.limit = wavelet.pars$period.lower.limit,
    detrend = isTRUE(wavelet.pars$detrend),
    mode = "fast"
  )

  if (is.null(wv_obs$gws_period) || !is.numeric(wv_obs$gws_period)) {
    stop("wavelet_spectral_analysis(obs) must return numeric $gws_period.", call. = FALSE)
  }
  if (is.null(wv_obs$gws) || !is.numeric(wv_obs$gws)) {
    stop("wavelet_spectral_analysis(obs) must return numeric $gws.", call. = FALSE)
  }

  power.period <- as.numeric(wv_obs$gws_period)

  # Use unmasked observed GWS when available
  if (!is.null(wv_obs$gws_unmasked) && is.numeric(wv_obs$gws_unmasked)) {
    power.obs <- as.numeric(wv_obs$gws_unmasked)
  } else {
    power.obs <- as.numeric(wv_obs$gws)
  }

  # Extract significance curves
  power.signif <- wv_obs$gws_signif
  power.signif_unmasked <- wv_obs$gws_signif_unmasked

  # Regrid if needed
  if (length(power.obs) != length(power.period)) {
    power.obs <- gws_regrid(wv_obs, power.period, use_unmasked = TRUE)
  }
  if (length(power.signif) != length(power.period)) {
    stop("Observed significance curve length does not match observed gws_period.", call. = FALSE)
  }

  signif_grid <- fill_nearest(as.numeric(power.signif))
  gws_obs_grid <- fill_nearest(as.numeric(power.obs))

  # Identify significant periods
  R_obs <- gws_obs_grid / signif_grid
  periods_sig_core <- which(is.finite(R_obs) & (R_obs > 1))

  # Split into contiguous regions
  split_regions <- function(idx) {
    idx <- sort(unique(as.integer(idx)))
    if (length(idx) == 0) return(list())
    brk <- c(TRUE, diff(idx) != 1L)
    split(idx, cumsum(brk))
  }

  wavelet_active <- TRUE
  wavelet_diag <- list(
    power.period = power.period,
    power.obs = power.obs,
    power.signif = power.signif,
    periods_sig_core = periods_sig_core,
    regions = list(),
    P_obs_reg = numeric(0),
    C_obs_reg = numeric(0)
  )

  P_sim_reg <- NULL
  P_sim_bg  <- NULL
  presence_rpad <- NULL

  # -------------------------------------------------------------------------
  # FIX: ALWAYS allocate and compute caches (independent of significance)
  # -------------------------------------------------------------------------
  gws_cache_mat <- matrix(NA_real_, nrow = length(power.period), ncol = n_realizations)
  gws_cache_unmasked_mat <- matrix(NA_real_, nrow = length(power.period), ncol = n_realizations)

  for (j in seq_len(n_realizations)) {
    wv <- wavelet_spectral_analysis(
      series_sim_for_stats[, j],
      signif.level = wavelet.pars$signif.level,
      noise.type = wavelet.pars$noise.type,
      period.lower.limit = wavelet.pars$period.lower.limit,
      detrend = isTRUE(wavelet.pars$detrend),
      mode = "fast"
    )

    gws_cache_mat[, j] <- gws_regrid(wv, power.period, use_unmasked = FALSE)
    gws_cache_unmasked_mat[, j] <- gws_regrid(wv, power.period, use_unmasked = TRUE)
  }

  # -------------------------------------------------------------------------
  # The rest: only compute region-based metrics if observed has significant signal
  # -------------------------------------------------------------------------
  if (length(periods_sig_core) == 0) {
    wavelet_active <- FALSE
  } else {
    regions <- split_regions(periods_sig_core)
    n_regions <- length(regions)

    core_union <- sort(unique(unlist(regions, use.names = FALSE)))
    bg_idx <- setdiff(seq_along(power.period), core_union)

    # Observed regional power and contrast
    P_obs_reg <- vapply(regions, function(ii) sum(gws_obs_grid[ii], na.rm = TRUE), FUN.VALUE = numeric(1))
    P_obs_bg <- sum(gws_obs_grid[bg_idx], na.rm = TRUE)
    P_obs_bg <- max(P_obs_bg, min_bg)
    C_obs_reg <- P_obs_reg / P_obs_bg

    # Padded periods for presence check
    periods_sig_pad <- periods_sig_core
    if (isTRUE(padding)) {
      n_periods <- length(power.period)
      periods_sig_pad <- unique(sort(c(
        pmax(periods_sig_core - 1L, 1L),
        periods_sig_core,
        pmin(periods_sig_core + 1L, n_periods)
      )))
    }

    # Initialize storage
    P_sim_reg <- matrix(NA_real_, nrow = n_realizations, ncol = n_regions)
    P_sim_bg  <- rep(NA_real_, n_realizations)
    presence_rpad <- rep(FALSE, n_realizations)

    # Use the already-computed masked cache for metrics (no second wavelet call)
    for (j in seq_len(n_realizations)) {
      gj_masked <- gws_cache_mat[, j]

      bg <- sum(gj_masked[bg_idx], na.rm = TRUE)
      bg <- max(bg, min_bg)
      P_sim_bg[j] <- bg

      for (r in seq_len(n_regions)) {
        P_sim_reg[j, r] <- sum(gj_masked[regions[[r]]], na.rm = TRUE)
      }

      rpad <- gj_masked[periods_sig_pad] / signif_grid[periods_sig_pad]
      presence_rpad[j] <- any(is.finite(rpad) & (rpad > 1))
    }

    wavelet_diag$regions <- regions
    wavelet_diag$P_obs_reg <- P_obs_reg
    wavelet_diag$C_obs_reg <- C_obs_reg
    wavelet_diag$bg_idx <- bg_idx
    wavelet_diag$periods_sig_pad <- periods_sig_pad
  }

  list(
    active = wavelet_active,
    diagnostics = wavelet_diag,
    power.period = power.period,
    power.obs = power.obs,
    power.signif = power.signif,
    power.signif_unmasked = power.signif_unmasked,
    P_sim_reg = P_sim_reg,
    P_sim_bg = P_sim_bg,
    presence_rpad = presence_rpad,
    gws_cache = gws_cache_mat,
    gws_cache_unmasked = gws_cache_unmasked_mat
  )
}


#' Relax bounds for one filter
#'
#' @description
#' Applies one relaxation step to a single filter. Updates bounds in place
#' and returns status.
#'
#' @param filter_name Character. Name of filter to relax
#' @param bounds_env Environment containing bounds
#' @param wavelet_active_env Environment containing wavelet_active flag
#' @param recompute_tailmass_fn Function to recompute tail mass when thresholds change
#'
#' @return List with changed (logical) and msg (character)
#'
#' @keywords internal
#' @export
relax_bounds_one_filter <- function(filter_name, bounds_env, wavelet_active_env,
                                    recompute_tailmass_fn) {

  msg <- NULL
  changed <- FALSE
  b <- bounds_env

  if (filter_name == "mean") {
    old <- b$mean
    b$mean <- min(b$mean * b$relax.mult, b$relax.mean.max)
    changed <- (b$mean > old + 1e-15)
    msg <- sprintf("mean tol %.4f -> %.4f", old, b$mean)
  }

  if (filter_name == "sd") {
    old <- b$sd
    b$sd <- min(b$sd * b$relax.mult, b$relax.sd.max)
    changed <- (b$sd > old + 1e-15)
    msg <- sprintf("sd tol %.4f -> %.4f", old, b$sd)
  }

  if (filter_name == "tail_low") {
    # First try to relax tolerance
    if (b$tail.tol.log < b$relax.tail.tol.log.max - 1e-15) {
      old <- b$tail.tol.log
      b$tail.tol.log <- min(b$tail.tol.log * b$relax.mult, b$relax.tail.tol.log.max)
      changed <- TRUE
      msg <- sprintf("tail.tol.log %.4f -> %.4f", old, b$tail.tol.log)
    }
    # Then try to relax threshold
    else if (b$tail.low.p < b$relax.tail.p.low.max - 1e-15) {
      old <- b$tail.low.p
      b$tail.low.p <- min(b$tail.low.p + b$relax.tail.p.step, b$relax.tail.p.low.max)
      recompute_tailmass_fn()
      changed <- TRUE
      msg <- sprintf("tail.low.p %.2f -> %.2f (recompute)", old, b$tail.low.p)
    } else {
      msg <- "tail_low at limit"
    }
  }

  if (filter_name == "tail_high") {
    # First try to relax tolerance
    if (b$tail.tol.log < b$relax.tail.tol.log.max - 1e-15) {
      old <- b$tail.tol.log
      b$tail.tol.log <- min(b$tail.tol.log * b$relax.mult, b$relax.tail.tol.log.max)
      changed <- TRUE
      msg <- sprintf("tail.tol.log %.4f -> %.4f", old, b$tail.tol.log)
    }
    # Then try to relax threshold
    else if (b$tail.high.p > b$relax.tail.p.high.min + 1e-15) {
      old <- b$tail.high.p
      b$tail.high.p <- max(b$tail.high.p - b$relax.tail.p.step, b$relax.tail.p.high.min)
      recompute_tailmass_fn()
      changed <- TRUE
      msg <- sprintf("tail.high.p %.2f -> %.2f (recompute)", old, b$tail.high.p)
    } else {
      msg <- "tail_high at limit"
    }
  }

  if (filter_name == "wavelet") {
    wavelet_active <- get("wavelet_active", envir = wavelet_active_env)

    if (!isTRUE(wavelet_active)) {
      msg <- "wavelet already inactive"
    }
    # Relax sig.frac
    else if (b$sig.frac > b$relax.wavelet.sig.frac.min + 1e-15) {
      old <- b$sig.frac
      b$sig.frac <- max(b$sig.frac - b$relax.wavelet.sig.frac.step, b$relax.wavelet.sig.frac.min)
      changed <- TRUE
      msg <- sprintf("sig.frac %.2f -> %.2f", old, b$sig.frac)
    }
    # Relax region tolerance
    else if (b$wavelet.region.tol < b$relax.wavelet.region.tol.max - 1e-15) {
      old <- b$wavelet.region.tol
      b$wavelet.region.tol <- min(b$wavelet.region.tol + b$relax.wavelet.region.tol.step,
                                  b$relax.wavelet.region.tol.max)
      changed <- TRUE
      msg <- sprintf("wavelet.region.tol %.2f -> %.2f", old, b$wavelet.region.tol)
    }
    # Relax contrast tolerance
    else if (b$wavelet.contrast.tol < b$relax.wavelet.contrast.tol.max - 1e-15) {
      old <- b$wavelet.contrast.tol
      b$wavelet.contrast.tol <- min(b$wavelet.contrast.tol + b$relax.wavelet.contrast.tol.step,
                                    b$relax.wavelet.contrast.tol.max)
      changed <- TRUE
      msg <- sprintf("wavelet.contrast.tol %.2f -> %.2f", old, b$wavelet.contrast.tol)
    }
    # Disable presence requirement
    else if (isTRUE(b$wavelet.require_presence)) {
      b$wavelet.require_presence <- FALSE
      changed <- TRUE
      msg <- "wavelet.require_presence TRUE -> FALSE"
    }
    # Finally disable wavelet entirely
    else {
      assign("wavelet_active", FALSE, envir = wavelet_active_env)
      changed <- TRUE
      msg <- "wavelet disabled"
    }
  }

  list(changed = changed, msg = msg)
}


# ==============================================================================
# Logging Functions for filter_warm_simulations()
# ==============================================================================

#' Log initial setup information
#'
#' @description
#' Displays general information at the start of filtering.
#'
#' @param n_obs Number of observations in observed series
#' @param n_sim Number of years in simulated series
#' @param n_realizations Number of realizations
#' @param sample_num Target number to sample
#' @param relax_priority Relaxation priority vector
#'
#' @keywords internal
#' @export
log_filtering_setup <- function(n_obs, n_sim, n_realizations, sample_num, relax_priority) {
  .log_info(strrep("=", 75))
  .log_info("FILTERING SETUP")
  .log_info(strrep("=", 75))
  .log_info(sprintf("Observed series: %d years", n_obs))
  .log_info(sprintf("Simulated series: %d years x %d realizations", n_sim, n_realizations))
  .log_info(sprintf("Target: Sample %d realizations from pool", sample_num))
  .log_info(sprintf("Relaxation priority: %s", paste(relax_priority, collapse = " > ")))
  .log_info(sprintf("  Filters relax left-to-right: %s relaxes FIRST, %s relaxes LAST",
                   relax_priority[1], relax_priority[length(relax_priority)]))
  .log_info(strrep("=", 75))
  .log_info("")
}

#' Log major step progress
#'
#' @description
#' Displays progress message for major computational steps.
#'
#' @param step_name Name of the step
#' @param details Optional details string
#'
#' @keywords internal
#' @export
log_step <- function(step_name, details = NULL) {
  msg <- if (!is.null(details)) {
    sprintf(" %s - %s", step_name, details)
  } else {
    sprintf(" %s", step_name)
  }
  .log_info(msg)
}

#' Log filter iteration details with table format
#'
#' @description
#' Prints iteration diagnostics in table format for all iterations.
#'
#' @param iter Iteration number
#' @param passes List of pass vectors
#' @param pool Vector of pool indices
#' @param n_total Total number of realizations
#' @param target Target pool size
#' @param bounds Bounds environment or list
#' @param tail_metrics Tail metrics list
#' @param wavelet_active Logical
#' @param wavelet_pars Wavelet parameters list
#' @param note Optional note string
#'
#' @keywords internal
#' @export
log_filter_iteration <- function(iter, passes, pool, n_total, target, bounds,
                                 tail_metrics, wavelet_active, wavelet_pars,
                                 note = NULL) {

  active <- names(passes)[!vapply(passes, is.null, logical(1))]
  pool_size <- length(pool)
  pool_pct <- if (n_total > 0) 100 * pool_size / n_total else 0

  # Determine status
  if (pool_size >= target) {
    status <- "TARGET REACHED"
    status_icon <- "[OK]"
  } else {
    status <- "Need more candidates"
    status_icon <- "[>>]"
  }

  # Iteration header
  if (iter == 0L) {
    .log_info(strrep("-", 75))
    .log_info(sprintf("ITERATION %d - Initial Evaluation", iter))
  } else {
    .log_info("")
    .log_info(strrep("-", 75))
    .log_info(sprintf("ITERATION %d - %s", iter, note))
  }
  .log_info(strrep("-", 75))

  # Show filter table
  filter_order <- c("mean", "sd", "tail_low", "tail_high", "wavelet")
  show_filters <- intersect(filter_order, active)

  if (length(show_filters) > 0) {
    # Table header
    .log_info(sprintf("%-12s %10s %8s  %-30s", "Filter", "Passed", "Rate", "Criteria"))
    .log_info(strrep("-", 75))

    # Table rows
    for (nm in show_filters) {
      if (!is.null(passes[[nm]])) {
        n_pass <- sum(passes[[nm]])
        rate <- sprintf("%.1f%%", 100 * mean(passes[[nm]]))
        crit <- criteria_string_compact(nm, bounds, tail_metrics, wavelet_active, wavelet_pars)

        .log_info(sprintf("%-12s %10d %8s  %-30s", nm, n_pass, rate, crit))
      }
    }
    .log_info(strrep("-", 75))
  }

  # Status line
  .log_info(sprintf("%s Pool: %d / %d (%.1f%%) | Need: %d | Status: %s",
                   status_icon, pool_size, n_total, pool_pct, target, status))

  if (pool_size >= target) {
    .log_info(strrep("=", 75))
  }

  invisible(NULL)
}

#' Log final summary
#'
#' @description
#' Displays final filtering results.
#'
#' @param pool_size Final pool size
#' @param n_total Total realizations
#' @param n_sampled Number sampled
#' @param relaxation_level Relaxation level reached
#'
#' @keywords internal
#' @export
log_final_summary <- function(pool_size, n_total, n_sampled, relaxation_level) {
  .log_info("")
  .log_info(strrep("=", 75))
  .log_info("FILTERING COMPLETE")
  .log_info(strrep("=", 75))
  .log_info(sprintf("Final pool: %d / %d realizations (%.1f%%)",
                   pool_size, n_total, 100 * pool_size / n_total))
  .log_info(sprintf("Sampled: %d realizations", n_sampled))
  .log_info(sprintf("Relaxation level: %s", relaxation_level))
  .log_info(strrep("=", 75))
}

#' Compact criteria string for a filter
#'
#' @description
#' Creates a compact human-readable string describing current filter criteria.
#'
#' @param filter_name Character filter name
#' @param bounds Bounds environment or list
#' @param tail_metrics Tail metrics list
#' @param wavelet_active Logical
#' @param wavelet_pars Wavelet parameters
#'
#' @return Character string (compact)
#'
#' @keywords internal
#' @export
criteria_string_compact <- function(filter_name, bounds, tail_metrics,
                                    wavelet_active, wavelet_pars) {

  if (filter_name == "mean") {
    return(sprintf("tol = %.4f", bounds$mean))
  }

  if (filter_name == "sd") {
    return(sprintf("tol = %.4f", bounds$sd))
  }

  if (filter_name == "tail_low") {
    return(sprintf("p=%.2f, log.tol=%.4f", bounds$tail.low.p, bounds$tail.tol.log))
  }

  if (filter_name == "tail_high") {
    return(sprintf("p=%.2f, log.tol=%.4f", bounds$tail.high.p, bounds$tail.tol.log))
  }

  if (filter_name == "wavelet") {
    if (!isTRUE(wavelet_active)) return("inactive")
    return(sprintf("sig.frac >= %.2f", bounds$sig.frac))
  }

  "NA"
}

#' Format numeric with specified digits
#'
#' @param x Numeric value
#' @param digits Number of decimal places
#'
#' @return Character string
#'
#' @keywords internal
#' @export
fmt_num <- function(x, digits = 4L) {
  if (is.null(x) || length(x) == 0L || any(!is.finite(x))) return("NA")
  formatC(x, format = "f", digits = digits)
}

#' Format percentage
#'
#' @param x Numeric proportion (0 to 1)
#' @param digits Number of decimal places
#'
#' @return Character string with percent sign
#'
#' @keywords internal
#' @export
fmt_pct <- function(x, digits = 1L) {
  sprintf(paste0("%.", digits, "f%%"), 100 * x)
}

#' Pad string to the right
#'
#' @param x Character string
#' @param width Target width
#'
#' @return Padded character string
#'
#' @keywords internal
#' @export
pad_right <- function(x, width) {
  x <- as.character(x)
  n <- nchar(x, type = "width")
  ifelse(n >= width, substr(x, 1, width), paste0(x, strrep(" ", width - n)))
}


