# ==============================================================================
# Helper Functions for filter_warm_simulations()
# ==============================================================================

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
#' and ALWAYS caches GWS for efficiency.
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
#'   \item gws_cache: ALWAYS cached GWS matrix (n_periods x n_realizations)
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
    return_diagnostics = FALSE
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

  # Extract significance curve
  power.signif <- wv_obs$gws_signif
  power.signif_unmasked = wv_obs$gws_signif_unmasked

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
  gws_cache_mat <- NULL

  if (length(periods_sig_core) == 0) {
    warning("No significant periods found in observed GWS (R_obs <= 1 everywhere). Wavelet filter disabled for this run.", call. = FALSE)
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


    gws_cache_mat <- matrix(NA_real_, nrow = length(power.period), ncol = n_realizations)

    # Compute wavelet for each realization -------------------------------------

    for (j in seq_len(n_realizations)) {

      wv <- wavelet_spectral_analysis(
        series_sim_for_stats[, j],
        signif.level = wavelet.pars$signif.level,
        noise.type = wavelet.pars$noise.type,
        period.lower.limit = wavelet.pars$period.lower.limit,
        detrend = isTRUE(wavelet.pars$detrend),
        return_diagnostics = FALSE
      )

      # Regrid and cache
      gj <- gws_regrid(wv, power.period, use_unmasked = FALSE)
      gws_cache_mat[, j] <- gj

      # Background power
      bg <- sum(gj[bg_idx], na.rm = TRUE)
      bg <- max(bg, min_bg)
      P_sim_bg[j] <- bg

      # Regional powers
      for (r in seq_len(n_regions)) {
        P_sim_reg[j, r] <- sum(gj[regions[[r]]], na.rm = TRUE)
      }

      # Presence check
      rpad <- gj[periods_sig_pad] / signif_grid[periods_sig_pad]
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
    gws_cache = gws_cache_mat  # ALWAYS available now
  )
}

#' Compute pass vectors for all filters
#'
#' @description
#' Evaluates which realizations pass each filter criterion.
#'
#' @param rel.diff.mean Relative differences in mean
#' @param rel.diff.sd Relative differences in SD
#' @param tail_metrics List from compute_tailmass_metrics()
#' @param P_sim_reg Matrix of regional powers (n_realizations x n_regions)
#' @param P_sim_bg Vector of background powers (n_realizations)
#' @param presence_rpad Logical vector of presence indicators
#' @param wavelet_diag List of wavelet diagnostics
#' @param bounds Bounds environment or list
#' @param wavelet_active Logical for wavelet filter status
#' @param wavelet_pars List of wavelet parameters
#' @param n_realizations Number of realizations
#'
#' @return Named list of logical vectors (pass$mean, pass$sd, pass$tail_low, etc.)
#'
#' @keywords internal
#' @export
compute_pass_vectors <- function(rel.diff.mean, rel.diff.sd, tail_metrics,
                                 P_sim_reg, P_sim_bg, presence_rpad,
                                 wavelet_diag, bounds, wavelet_active,
                                 wavelet_pars, n_realizations) {

  pass <- list()

  # Basic filters
  pass$mean <- abs(rel.diff.mean) < bounds$mean
  pass$sd   <- abs(rel.diff.sd)   < bounds$sd

  # Tail-mass filters
  pass$tail_low  <- tail_metrics$logdiff_low  <= bounds$tail.tol.log
  pass$tail_high <- tail_metrics$logdiff_high <= bounds$tail.tol.log

  # Wavelet filter
  if (!isTRUE(wavelet_active)) {
    pass$wavelet <- NULL
  } else {
    # Regional power ratio
    P_ratio <- sweep(P_sim_reg, 2, pmax(wavelet_diag$P_obs_reg, bounds$wavelet.min_bg), FUN = "/")

    # Regional contrast
    C_sim_reg <- sweep(P_sim_reg, 1, P_sim_bg, FUN = "/")
    C_ratio <- sweep(C_sim_reg, 2, pmax(wavelet_diag$C_obs_reg, bounds$wavelet.min_bg), FUN = "/")

    # Check tolerances
    P_ok <- (P_ratio >= (1 - bounds$wavelet.region.tol)) & (P_ratio <= (1 + bounds$wavelet.region.tol))
    C_ok <- (C_ratio >= (1 - bounds$wavelet.contrast.tol)) & (C_ratio <= (1 + bounds$wavelet.contrast.tol))

    frac_regions_ok <- rowMeans(P_ok & C_ok)

    # Presence requirement
    if (isTRUE(bounds$wavelet.require_presence)) {
      # IMPROVEMENT: Use exposed parameter instead of hardcoding to signif.level
      presence_frac <- if (!is.null(bounds$wavelet.presence.frac)) {
        bounds$wavelet.presence.frac
      } else {
        as.numeric(wavelet_pars$signif.level)
      }

      cond2 <- apply(
        P_sim_reg, 1,
        function(v) any(v > (presence_frac * pmax(wavelet_diag$P_obs_reg, bounds$wavelet.min_bg)))
      )
      has_signal <- presence_rpad & cond2
    } else {
      has_signal <- rep(TRUE, n_realizations)
    }

    pass$wavelet <- (frac_regions_ok >= bounds$sig.frac) & has_signal
  }

  pass
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
