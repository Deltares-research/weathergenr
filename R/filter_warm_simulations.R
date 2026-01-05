#' Filter and sample WARM simulations by distributional, tail-mass, and wavelet criteria
#'
#' @description
#' Filters a matrix of simulated annual series (WARM realizations) against an observed
#' annual series using fast distributional checks (mean, standard deviation),
#' robust tail-shape checks based on lower/upper tail mass relative to observed
#' quantile thresholds, and an observed-relevant wavelet (GWS) filter.
#'
#' The function applies a tiered relaxation strategy when too few candidates pass.
#' If still insufficient after bounds$relax.max.iter iterations, it falls back
#' to a deterministic selection of the sample.num realizations with the smallest
#' absolute relative mean difference.
#'
#' If make.plots=TRUE, diagnostic plots are returned (time series overlay,
#' relative-difference scatter facets, and wavelet GWS ribbon/lines).
#'
#' @param series.obs Numeric vector. Observed annual series.
#' @param series.sim Numeric matrix. Simulated annual series, with years in rows and
#'   realizations in columns.
#' @param sample.num Integer scalar. Number of realizations to sample from the final
#'   candidate pool.
#' @param seed Optional integer. Random seed for reproducible window sampling and
#'   candidate sampling.
#' @param padding Logical scalar. If TRUE, pads the observed significant-period
#'   set by one index on each side when evaluating presence of observed-relevant signal.
#' @param bounds Named list. Overrides defaults for filtering and relaxation parameters.
#'   Common entries include:
#'   \itemize{
#'     \item mean, sd: relative-difference tolerances.
#'     \item tail.low.p, tail.high.p: quantile probabilities for tail thresholds.
#'     \item tail.tol.log: tolerance on |log(M_sim+eps)-log(M_obs+eps)|.
#'     \item tail.eps: epsilon used inside the log transform.
#'     \item sig.frac, wavelet.region.tol, wavelet.contrast.tol,
#'           wavelet.require_presence, wavelet.presence.frac: wavelet filter controls.
#'     \item relax.*: relaxation settings, including relax.max.iter.
#'   }
#' @param wavelet.pars Named list passed to wavelet_spectral_analysis() for
#'   observed and simulated series. Expected entries:
#'   signif.level, noise.type, period.lower.limit, detrend.
#' @param make.plots Logical scalar. If TRUE, returns diagnostic ggplot objects.
#' @param verbose Logical scalar. If TRUE, prints iteration diagnostics using
#'   logger::log_info().
#'
#' @details
#' \strong{Length harmonization:} If observed and simulated series have different
#' lengths, the function uses n_use = min(length(series.obs), nrow(series.sim)).
#' A random contiguous window of length n_use is selected from the observed
#' series (if needed), and each realization receives an independently sampled window
#' of length n_use (if needed).
#'
#' \strong{Tail mass metrics:} Let t_L and t_H be observed quantile
#' thresholds at tail.low.p and tail.high.p. Define lower-tail deficit
#' mass and upper-tail excess mass (normalized by n_use * robust_scale(obs)),
#' then compare simulations to observations using log-distance with epsilon
#' tail.eps.
#'
#' \strong{Wavelet filter (observed-relevant GWS):} Significant observed periods are
#' identified where observed GWS / signif > 1. These are grouped into contiguous
#' regions. For each realization, regional integrated power and regional-to-background
#' contrast are compared to observed within tolerances, requiring at least
#' sig.frac of regions to pass. Optionally requires presence of signal.
#'
#' \strong{Improvements in this version:}
#' \itemize{
#'   \item wavelet.presence.frac parameter now exposed (was hardcoded to signif.level)
#'   \item GWS always cached for efficiency (was conditional on make.plots)
#'   \item Simplified relaxation logic (removed complex tie-breaking)
#'   \item Extracted helper functions for better modularity and testability
#' }
#'
#' @return A list with elements:
#' \itemize{
#'   \item subsetted: matrix of realizations in the final candidate pool
#'     (original rows from series.sim).
#'   \item sampled: matrix of sample.num sampled realizations
#'     (original rows from series.sim).
#'   \item filter_summary: data.frame with per-filter pass counts and percentages,
#'     and the final relaxation level.
#'   \item diagnostics: list containing length/window metadata, pool indices,
#'     relaxation log, tail-mass details, wavelet diagnostics, and final parameters.
#'   \item plots: NULL or a list of ggplot objects if make.plots=TRUE.
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
#'
#' @export
filter_warm_simulations <- function(series.obs = NULL,
                                    series.sim = NULL,
                                    sample.num = 5,
                                    seed = NULL,
                                    padding = TRUE,
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
    warning("sample.num exceeds number of realizations; returning at most ncol(series.sim).", call. = FALSE)
    sample.num <- n_realizations
  }

  # ---------------------------------------------------------------------------
  # Bounds defaults - IMPROVEMENT: Added wavelet.presence.frac parameter
  # ---------------------------------------------------------------------------
  b_list <- modifyList(list(
    mean = 0.05,
    sd   = 0.05,

    tail.low.p  = 0.10,
    tail.high.p = 0.90,
    tail.tol.log = log(1.10),
    tail.eps = 1e-5,

    sig.frac = 0.80,
    wavelet.region.tol = 0.50,
    wavelet.contrast.tol = 0.30,
    wavelet.min_bg = 1e-12,
    wavelet.require_presence = TRUE,
    wavelet.presence.frac = NULL,  # NEW: NULL means use wavelet.pars$signif.level

    plot.wavelet_q = c(0.50, 0.95),

    relax.mult = 1.5,
    relax.mean.max = 0.25,
    relax.sd.max = 0.25,

    relax.tail.tol.log.max = log(3.0),
    relax.tail.p.step = 0.05,
    relax.tail.p.low.max  = 0.45,
    relax.tail.p.high.min = 0.55,

    relax.wavelet.sig.frac.step = 0.05,
    relax.wavelet.sig.frac.min = 0.40,
    relax.wavelet.region.tol.step = 0.10,
    relax.wavelet.region.tol.max = 1.00,
    relax.wavelet.contrast.tol.step = 0.10,
    relax.wavelet.contrast.tol.max = 1.00,

    relax.max.iter = 10L
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
    log_step("Computing distributional statistics",
             sprintf("mean, sd, tail mass for %d realizations", n_realizations))
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
    log_step("Computing wavelet spectra",
             sprintf("observed + %d simulated series", n_realizations))
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
  # IMPROVEMENT 3: Use helper function for pass vectors
  # ---------------------------------------------------------------------------
  .compute_pass_vectors <- function() {
    compute_pass_vectors(
      rel.diff.mean = rel.diff.mean,
      rel.diff.sd = rel.diff.sd,
      tail_metrics = tail_metrics,
      P_sim_reg = P_sim_reg,
      P_sim_bg = P_sim_bg,
      presence_rpad = presence_rpad,
      wavelet_diag = wavelet_diag,
      bounds = b,
      wavelet_active = wavelet_active,
      wavelet_pars = wavelet.pars,
      n_realizations = n_realizations
    )
  }

  .pool_from_pass <- function(pass_list) {
    active <- names(pass_list)[!vapply(pass_list, is.null, logical(1))]
    if (length(active) == 0) return(seq_len(n_realizations))
    ok <- rep(TRUE, n_realizations)
    for (nm in active) ok <- ok & pass_list[[nm]]
    which(ok)
  }

  # ---------------------------------------------------------------------------
  # Relaxation priority order (used throughout)
  # ---------------------------------------------------------------------------
  RELAX_ORDER <- c("mean", "sd", "tail_low", "tail_high", "wavelet")

  # ---------------------------------------------------------------------------
  # IMPROVEMENT 4: Simplified relaxation (no complex tie-breaking)
  # ---------------------------------------------------------------------------
  .relax_one <- function(filter_name) {
    relax_bounds_one_filter(
      filter_name = filter_name,
      bounds_env = b,
      wavelet_active_env = environment(),
      recompute_tailmass_fn = .recompute_tailmass
    )
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
  # Initial evaluation
  # ---------------------------------------------------------------------------

  passes <- .compute_pass_vectors()
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
    passes <- .compute_pass_vectors()
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
      log_info(sprintf(
        "Fallback activated: closest_mean to guarantee pool size = %d (from %d total).",
        sample.num, n_realizations
      ))
    }
  }

  # Sample from pool
  if (verbose) {
    log_step("Sampling realizations",
             sprintf("selecting %d from pool of %d", sample.num, length(pool)))
  }

  if (!is.null(seed)) set.seed(seed)
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
  # IMPROVEMENT 5: Use helper function for plotting
  # ---------------------------------------------------------------------------
  plots_out <- NULL
  if (make.plots) {
    if (verbose) {
      log_step("Creating diagnostic plots",
               "time series, statistics, and wavelet spectra")
    }

    if (length(pool) < 1L) {
      warning("make.plots=TRUE but final pool is empty; plots=NULL.", call. = FALSE)
    } else {
      plots_out <- create_filter_diagnostic_plots(
        obs.use = obs.use,
        series_sim_for_stats = series_sim_for_stats,
        pool = pool,
        rel.diff.mean = rel.diff.mean,
        rel.diff.sd = rel.diff.sd,
        tail_metrics = tail_metrics,
        power.period = power.period,
        power.obs = power.obs,
        power.signif = power.signif_unmasked,
        gws_cache_mat = gws_cache_mat,
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
      tailmass = list(
        low.p  = b$tail.low.p,
        high.p = b$tail.high.p,
        tol.log = b$tail.tol.log,
        eps = b$tail.eps,
        scale_obs = tail_metrics$scale_obs,
        low.thr  = if (is.finite(tail_metrics$thr_low)) tail_metrics$thr_low else NULL,
        high.thr = if (is.finite(tail_metrics$thr_high)) tail_metrics$thr_high else NULL,
        M_obs_low  = tail_metrics$M_obs_low,
        M_obs_high = tail_metrics$M_obs_high
      ),
      wavelet = wavelet_diag,
      final_params = as.list(b)
    ),
    plots = plots_out
  )
}
