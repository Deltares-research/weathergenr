#' Continuous Wavelet Spectral Analysis with COI-Aware Significance
#'
#' @description
#' Performs continuous wavelet transform (CWT) analysis using the Morlet mother
#' wavelet following Torrence & Compo (1998), with extensions for:
#' \itemize{
#'   \item cone-of-influence (COI) aware global spectrum estimation,
#'   \item effective sample size (neff)-based significance testing,
#'   \item defensible AR(1) red-noise background estimation,
#'   \item stable and complete signal reconstruction from significant scales.
#' }
#'
#' The function is designed for **detection and diagnostic analysis of dominant
#' time scales**, not for producing orthogonal components. Reconstructed components
#' are scale-localized but not independent.
#'
#' @param variable Numeric vector. Input time series (regularly spaced, no missing
#'   values). Minimum length is 16 observations.
#' @param signif.level Numeric scalar in (0, 1). Significance level for wavelet
#'   power and global wavelet spectrum tests. Default is 0.90.
#' @param noise.type Character. Background noise model for significance testing.
#'   One of \code{"white"} or \code{"red"} (default). For red noise, an AR(1) model
#'   is estimated using Yule-Walker.
#' @param period.lower.limit Numeric scalar. Minimum Fourier period (in time units)
#'   to consider when identifying significant scales. Default is 2.
#' @param detrend Logical. If \code{TRUE}, removes a linear trend **slope only**
#'   (preserves the series mean level) before wavelet analysis. Default is \code{FALSE}.
#' @param mode Character. Computation mode: \code{"fast"} (default) for speed-optimized
#'   filtering with essential outputs only, or \code{"complete"} for comprehensive analysis
#'   including reconstruction and diagnostics. Fast mode returns 9 outputs (~50 KB),
#'   complete mode returns 18+ outputs (~5 MB). Fast mode is ~3x faster by skipping
#'   signal reconstruction and detailed diagnostics.
#' @param return_diagnostics Logical. DEPRECATED. Use \code{mode = "complete"} instead.
#'   If \code{TRUE}, forces complete mode and returns internal diagnostic quantities.
#' @param return_recon_error Logical. If \code{TRUE}, returns scalar reconstruction
#'   error metrics verifying closure of reconstructed components.
#' @param lag1_ci Logical. If \code{TRUE}, computes a bootstrap confidence interval
#'   for the AR(1) coefficient (diagnostic only; does not affect inference).
#' @param lag1_ci_level Numeric scalar in (0, 1). Confidence level for the AR(1)
#'   bootstrap interval. Default is 0.95.
#' @param lag1_boot_n Integer. Number of bootstrap replicates used to estimate the
#'   AR(1) confidence interval. Default is 500.
#' @param seed Optional numeric scalar. Random seed for reproducible bootstrap
#'   diagnostics. Default is \code{NULL}.
#' @param warn_neff Logical. If \code{TRUE} (default), emits a warning when the
#'   effective sample size is small for most scales.
#' @param neff_warn_min Numeric scalar. Threshold below which neff is considered
#'   small. Default is 5.
#' @param neff_warn_frac Numeric scalar in (0, 1). Fraction of scales with
#'   \code{neff < neff_warn_min} required to trigger a warning. Default is 0.60.
#'
#' @details
#' The implementation follows Torrence & Compo (1998) for the Morlet wavelet,
#' with the following methodological refinements:
#' \itemize{
#'   \item two-sided (TC98-consistent) wavelet normalization,
#'   \item COI-masked global wavelet spectrum,
#'   \item effective degrees of freedom based on wavelet decorrelation time,
#'   \item separation of inference curves from plotting curves.
#' }
#'
#' Reconstructed components correspond to **individual wavelet scales** and are
#' additive but not orthogonal. They should not be interpreted as independent
#' stochastic modes.
#'
#' @return A list with outputs depending on mode:
#'
#' \strong{Fast mode (9 outputs):}
#' \describe{
#'   \item{gws}{Numeric vector. COI-masked global wavelet spectrum.}
#'   \item{gws_unmasked}{Numeric vector. Unmasked global wavelet spectrum (plotting only).}
#'   \item{gws_period}{Numeric vector. Fourier periods corresponding to wavelet scales.}
#'   \item{gws_signif}{Numeric vector. COI- and neff-corrected significance threshold (inference).}
#'   \item{gws_signif_unmasked}{Numeric vector. Unmasked significance threshold (plotting only).}
#'   \item{has_significance}{Logical. Indicates whether any significant scales were found.}
#'   \item{signif_periods}{Integer vector. Indices of significant wavelet scales.}
#'   \item{coi}{Numeric vector. Cone of influence in time units.}
#'   \item{power}{Numeric matrix. Wavelet power spectrum (for basic plotting).}
#' }
#'
#' \strong{Complete mode (additional 10+ outputs):}
#' \describe{
#'   \item{power_coi}{Numeric matrix. COI-masked wavelet power.}
#'   \item{sigm}{Numeric matrix. Pointwise significance ratio (unmasked).}
#'   \item{sigm_coi}{Numeric matrix. COI-masked pointwise significance ratio.}
#'   \item{power_signif_coi}{Logical matrix. COI-masked pointwise significance mask.}
#'   \item{wave}{Complex matrix. Wavelet coefficients (scales x time).}
#'   \item{comps}{Numeric matrix. Reconstructed components at significant scales plus residual.}
#'   \item{comps_names}{Character vector. Column names of comps.}
#'   \item{gws_n_coi}{Numeric vector. Number of time points inside the COI per scale.}
#'   \item{gws_neff}{Numeric vector. Effective sample size per scale (masked).}
#'   \item{gws_neff_unmasked}{Numeric vector. Effective sample size per scale (unmasked).}
#'   \item{diagnostics}{List. Returned only if \code{return_diagnostics = TRUE}.}
#'   \item{reconstruction_error}{List. Returned only if \code{return_recon_error = TRUE}.}
#' }
#'
#' @references
#' Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet analysis.
#' \emph{Bulletin of the American Meteorological Society}, 79(1), 61-78.
#'
#' @seealso
#' \code{\link{morlet_wavelet}}, \code{\link{extract_wavelet_components}}
#'
#' @export
wavelet_spectral_analysis <- function(variable,
                                      signif.level = 0.90,
                                      noise.type = "red",
                                      period.lower.limit = 2,
                                      detrend = FALSE,
                                      mode = c("fast", "complete"),
                                      return_diagnostics = TRUE,
                                      return_recon_error = FALSE,
                                      lag1_ci = FALSE,
                                      lag1_ci_level = 0.95,
                                      lag1_boot_n = 500,
                                      seed = NULL,
                                      warn_neff = FALSE,
                                      neff_warn_min = 5,
                                      neff_warn_frac = 0.60) {

  # --- Input Validation ---
  if (!is.numeric(variable)) stop("variable must be numeric")
  if (anyNA(variable)) stop("variable contains missing values")
  if (length(variable) < 16) stop("variable must have at least 16 observations")

  # --- Mode validation and backward compatibility ---
  mode <- match.arg(mode)

  # Backward compatibility: return_diagnostics forces complete mode
  if (isTRUE(return_diagnostics)) {
    mode <- "complete"
    if (!missing(return_diagnostics)) {
      warning("'return_diagnostics' is deprecated. Use mode='complete' instead.",
              call. = FALSE)
    }
  }

  compute_full <- (mode == "complete")

  if (!(noise.type %in% c("white", "red"))) stop("noise.type must be 'white' or 'red'")

  if (!is.numeric(signif.level) || length(signif.level) != 1L ||
      signif.level <= 0 || signif.level >= 1) {
    stop("signif.level must be between 0 and 1")
  }

  if (!is.numeric(period.lower.limit) || length(period.lower.limit) != 1L ||
      period.lower.limit < 0) {
    stop("period.lower.limit must be a non-negative number")
  }

  if (!is.logical(return_diagnostics) || length(return_diagnostics) != 1L) {
    stop("return_diagnostics must be TRUE/FALSE")
  }

  if (!is.logical(return_recon_error) || length(return_recon_error) != 1L) {
    stop("return_recon_error must be TRUE/FALSE")
  }

  if (!is.logical(lag1_ci) || length(lag1_ci) != 1L) stop("lag1_ci must be TRUE/FALSE")
  if (!is.numeric(lag1_ci_level) || length(lag1_ci_level) != 1L ||
      lag1_ci_level <= 0 || lag1_ci_level >= 1) {
    stop("lag1_ci_level must be between 0 and 1")
  }
  if (!is.numeric(lag1_boot_n) || length(lag1_boot_n) != 1L || lag1_boot_n < 50) {
    stop("lag1_boot_n must be >= 50")
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    stop("seed must be NULL or a single number")
  }

  if (!is.logical(warn_neff) || length(warn_neff) != 1L) stop("warn_neff must be TRUE/FALSE")
  if (!is.numeric(neff_warn_min) || length(neff_warn_min) != 1L || neff_warn_min <= 0) {
    stop("neff_warn_min must be > 0")
  }
  if (!is.numeric(neff_warn_frac) || length(neff_warn_frac) != 1L ||
      neff_warn_frac <= 0 || neff_warn_frac >= 1) {
    stop("neff_warn_frac must be between 0 and 1")
  }

  # --- Wavelet Transform Analysis ---
  variable_org <- as.numeric(variable)

  if (isTRUE(detrend)) {
    tt <- seq_along(variable_org)
    fit <- stats::lm(variable_org ~ tt)
    b <- stats::coef(fit)[2]
    if (!is.finite(b)) b <- 0
    variable_org <- variable_org - b * (tt - mean(tt))
  }

  variance1 <- stats::var(variable_org)
  n1 <- length(variable_org)

  sdx <- stats::sd(variable_org)
  if (!is.finite(sdx) || sdx <= 0) stop("variable has zero or non-finite standard deviation")

  variable_mean <- mean(variable_org)

  # Standardize (for transform only)
  variable_std <- (variable_org - variable_mean) / sdx

  # Zero-pad
  base2 <- floor(log2(n1) + 0.4999)
  variable_pad <- c(variable_std, rep(0, (2^(base2 + 1) - n1)))
  n <- length(variable_pad)

  dt <- 1
  dj <- 0.25
  s0 <- 2 * dt
  J <- floor((1 / dj) * log2(n1 * dt / s0))
  scale <- s0 * 2^((0:J) * dj)

  k <- c(0:(floor(n / 2)), -rev(1:floor((n - 1) / 2))) * ((2 * pi) / (n * dt))
  f <- stats::fft(variable_pad)

  params <- morlet_parameters(k0 = 6)
  fourier_factor <- params["fourier_factor"]
  coi_base <- params["coi"]

  wave <- t(sapply(seq_len(J + 1), function(a1) {
    daughter <- morlet_wavelet(k, scale[a1], k0 = 6)
    stats::fft(f * daughter, inverse = TRUE) / n
  }))

  period <- fourier_factor * scale

  coi_indices <- c(
    1e-5,
    seq_len(floor((n1 - 1) / 2)),
    rev(seq_len(ceiling((n1 - 1) / 2))),
    1e-5
  )
  coi <- coi_base * dt * coi_indices[1:n1]

  wave <- wave[, 1:n1, drop = FALSE]
  power <- abs(wave)^2

  # --- COI mask (used by gws + plotting masks) ---
  coi_mask <- outer(period, coi, FUN = "<=")
  n_coi <- rowSums(coi_mask)

  # --- Significance Testing (pointwise) ---
  empir <- c(2, 0.776, 2.32, 0.60)
  dofmin <- empir[1]
  gamma_fac <- empir[3]

  # --- Lag1 estimation (reuse centered variable) ---
  x_centered <- variable_org - mean(variable_org)

  if (noise.type == "white") {
    lag1 <- 0
  } else {
    lag1 <- tryCatch({
      fit <- stats::ar(x_centered, aic = FALSE, order.max = 1, method = "yw")
      as.numeric(fit$ar[1])
    }, error = function(e) {
      stats::cor(x_centered[-length(x_centered)], x_centered[-1])
    })

    if (!is.finite(lag1)) lag1 <- 0
    lag1 <- max(min(lag1, 0.999), -0.999)
  }

  # --- Optional lag1 uncertainty (diagnostic only) ---
  lag1_ci_out <- NULL
  if (lag1_ci && noise.type == "red") {

    # RNG state management
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        old_seed_lag1 <- .Random.seed
        has_seed_lag1 <- TRUE
      } else {
        has_seed_lag1 <- FALSE
      }
      on.exit({
        if (has_seed_lag1) .Random.seed <<- old_seed_lag1
      }, add = TRUE)
      set.seed(seed)
    }

    sig_eps <- stats::sd(x_centered) * sqrt(max(1e-12, 1 - lag1^2))

    phi_boot <- replicate(lag1_boot_n, {
      xb <- as.numeric(stats::arima.sim(n = n1, model = list(ar = lag1), sd = sig_eps))
      xb <- xb - mean(xb)

      ph <- tryCatch({
        fitb <- stats::ar(xb, aic = FALSE, order.max = 1, method = "yw")
        as.numeric(fitb$ar[1])
      }, error = function(e) {
        stats::cor(xb[-length(xb)], xb[-1])
      })

      if (!is.finite(ph)) ph <- 0
      max(min(ph, 0.999), -0.999)
    })

    alpha <- (1 - lag1_ci_level) / 2
    lag1_ci_out <- list(
      level = lag1_ci_level,
      n_boot = lag1_boot_n,
      lag1_hat = lag1,
      lower = as.numeric(stats::quantile(phi_boot, probs = alpha, na.rm = TRUE)),
      upper = as.numeric(stats::quantile(phi_boot, probs = 1 - alpha, na.rm = TRUE))
    )
  }

  freq <- dt / period
  fft_theor <- (1 - lag1^2) / (1 - 2 * lag1 * cos(freq * 2 * pi) + lag1^2)

  # --- COI-masked power (needed for gws in all modes) ---
  power_coi <- power
  power_coi[!coi_mask] <- NA_real_

  # --- Complete-mode-only: pointwise significance ---
  sigm <- NULL
  sigm_coi <- NULL
  power_signif_coi <- NULL

  if (compute_full) {
    chisquare <- stats::qchisq(signif.level, dofmin) / dofmin
    signif <- fft_theor * chisquare
    sigm <- sweep(power, 1, signif, FUN = "/")

    sigm_coi <- sigm
    sigm_coi[!coi_mask] <- NA_real_

    power_signif_coi <- sigm_coi > 1
    power_signif_coi[is.na(sigm_coi)] <- FALSE
  }

  # --- gws (masked for inference + unmasked for plotting) ---
  mean_power_coi <- rowMeans(power_coi, na.rm = TRUE)
  mean_power_coi[!is.finite(mean_power_coi)] <- NA_real_

  gws_unmasked <- variance1 * rowMeans(power)
  gws_unmasked[!is.finite(gws_unmasked)] <- mean(gws_unmasked[is.finite(gws_unmasked)], na.rm = TRUE)
  gws <- variance1 * mean_power_coi

  # --- COI-consistent neff (decorrelation-adjusted; inference curve) ---
  neff <- (n_coi * dt) / (gamma_fac * scale)
  neff[neff < 1] <- NA_real_
  neff[!is.finite(neff)] <- NA_real_

  # --- Short-record / small-neff warnings ---
  if (warn_neff) {
    ok_neff <- is.finite(neff)
    frac_small <- if (any(ok_neff)) mean(neff[ok_neff] < neff_warn_min) else 1

    if (!any(ok_neff) || frac_small >= neff_warn_frac) {
      warning(
        sprintf(
          "Wavelet significance may be unreliable: %.0f%% of scales have neff < %.2f (n=%d). Consider longer series or coarser scales/dj.",
          100 * frac_small, neff_warn_min, n1
        ),
        call. = FALSE
      )
    }
  }

  dof <- dofmin * neff
  dof[is.finite(dof) & dof < dofmin] <- dofmin

  chisq <- rep(NA_real_, length(dof))
  ok <- is.finite(dof) & dof > 0
  chisq[ok] <- stats::qchisq(signif.level, dof[ok]) / dof[ok]

  gws_signif <- rep(NA_real_, length(dof))
  gws_signif[ok] <- fft_theor[ok] * variance1 * chisq[ok]

  # --- UNMASKED gws significance (plotting only) ---
  neff_unmasked <- (n1 * dt) / (gamma_fac * scale)
  neff_unmasked[neff_unmasked < 1] <- 1
  neff_unmasked[!is.finite(neff_unmasked)] <- 1

  dof_unmasked <- dofmin * neff_unmasked
  dof_unmasked[is.finite(dof_unmasked) & dof_unmasked < dofmin] <- dofmin

  chisq_unmasked <- rep(NA_real_, length(dof_unmasked))
  ok_u <- is.finite(dof_unmasked) & dof_unmasked > 0
  chisq_unmasked[ok_u] <- stats::qchisq(signif.level, dof_unmasked[ok_u]) / dof_unmasked[ok_u]
  chisq_unmasked[!ok_u] <- 1

  gws_signif_unmasked <- rep(NA_real_, length(dof_unmasked))
  gws_signif_unmasked[ok_u] <- fft_theor[ok_u] * variance1 * chisq_unmasked[ok_u]
  gws_signif_unmasked[!ok_u] <- fft_theor[!ok_u] * variance1 * chisq_unmasked[!ok_u]

  # --- Identify Significant Periods (use inference curves only) ---
  sig_periods <- which(
    is.finite(gws) & is.finite(gws_signif) &
      (gws > gws_signif) & (period > period.lower.limit)
  )

  out_recon_err <- NULL

  if (length(sig_periods) == 0) {

    signif_periods <- integer(0)

    comps <- matrix(variable_org, ncol = 1)
    colnames(comps) <- "noise"

    if (return_recon_error) {
      out_recon_err <- list(recon_max_abs = 0, recon_rmse = 0)
    }

  } else {

    sig_periods_grp <- split(sig_periods, cumsum(c(1, diff(sig_periods) != 1)))
    signif_periods <- as.integer(unlist(
      lapply(sig_periods_grp, function(x) x[which.max(gws[x])]),
      use.names = FALSE
    ))

    # RECONSTRUCTION: Only in complete mode (skip in fast mode for speed)
    if (compute_full) {
      comps_tb <- extract_wavelet_components(
        wave = wave,
        signif_periods = signif_periods,
        scale = scale,
        dj = dj,
        dt = dt,
        variable_sd = sdx,
        variable_mean = variable_mean,
        Cdelta = 0.776,
        w0_0 = pi^(-1/4),
        include_residual = TRUE
      )

      comps_mat <- as.matrix(comps_tb)

      period_labels <- round(period[signif_periods], 2)
      comp_names <- paste0("period_", period_labels)

      if (ncol(comps_mat) != (length(signif_periods) + 1L)) {
        stop("Internal error: unexpected number of reconstructed component columns.")
      }

      colnames(comps_mat) <- c(comp_names, "noise")
      comps <- comps_mat

      if (return_recon_error) {
        recon_vec <- variable_org - rowSums(comps)
        out_recon_err <- list(
          recon_max_abs = max(abs(recon_vec)),
          recon_rmse = sqrt(mean(recon_vec^2))
        )
      }
    } else {
      comps <- NULL
    }
  }

  signif_periods <- as.integer(signif_periods)
  has_significance <- length(signif_periods) > 0

  # Build return list based on mode
  if (mode == "fast") {
    out <- list(
      gws = gws,
      gws_unmasked = gws_unmasked,
      gws_period = period,
      gws_signif = gws_signif,
      gws_signif_unmasked = gws_signif_unmasked,
      has_significance = has_significance,
      signif_periods = signif_periods,
      coi = coi,
      power = power
    )
  } else {
    out <- list(
      gws = gws,
      gws_unmasked = gws_unmasked,
      gws_period = period,
      gws_signif = gws_signif,
      gws_signif_unmasked = gws_signif_unmasked,
      has_significance = has_significance,
      signif_periods = signif_periods,
      coi = coi,
      power = power,
      power_coi = power_coi,
      sigm = sigm,
      sigm_coi = sigm_coi,
      power_signif_coi = power_signif_coi,
      wave = wave,
      comps = comps,
      comps_names = if (!is.null(comps)) colnames(comps) else NULL,
      gws_n_coi = n_coi,
      gws_neff = neff,
      gws_neff_unmasked = neff_unmasked
    )

    if (return_recon_error && !is.null(comps)) {
      out$reconstruction_error <- out_recon_err
    }

    if (return_diagnostics) {
      out$diagnostics <- list(
        lag1 = lag1,
        lag1_ci = lag1_ci_out,
        variance = variance1,
        n_original = n1,
        n_padded = n,
        dof = dof,
        dof_unmasked = dof_unmasked,
        scale = scale,
        fourier_factor = fourier_factor,
        n_coi = n_coi,
        gws_neff = neff,
        gws_neff_unmasked = neff_unmasked
      )
    }
  }

  out
}


# ==============================================================================
# Wavelet Utility Functions
# ==============================================================================
# This file contains utility functions for wavelet analysis and filtering.
# New additions for filter_warm_simulations support:
# - fill_nearest()
# - extract_signif_curve()
# - gws_regrid()
# ==============================================================================

#' Fill NA values with nearest non-NA value
#'
#' @description
#' Forward and backward fills NA values in a vector using the nearest
#' non-NA value. Useful for handling edge effects in interpolated data.
#'
#' @param v Numeric vector possibly containing NAs
#'
#' @return Numeric vector with NAs filled
#'
#' @keywords internal
#' @export
fill_nearest <- function(v) {
  if (all(is.na(v))) return(v)

  # Forward fill from first non-NA
  if (is.na(v[1])) {
    first <- which(!is.na(v))[1]
    if (!is.na(first)) v[1:(first - 1)] <- v[first]
  }

  # Forward fill
  for (i in 2:length(v)) {
    if (is.na(v[i])) v[i] <- v[i - 1]
  }

  # Backward fill
  for (i in (length(v) - 1):1) {
    if (is.na(v[i])) v[i] <- v[i + 1]
  }

  v
}

#' Extract significance curve from wavelet analysis result
#'
#' @description
#' Searches for GWS significance curve in wavelet analysis output.
#' Tries multiple possible field names.
#'
#' @param wv List output from wavelet_spectral_analysis()
#'
#' @return Numeric vector of significance values, or NULL if not found
#'
#' @keywords internal
#' @export
extract_signif_curve <- function(wv) {
  # Try multiple possible names
  cand <- c("gws_signif", "gws_sig", "signif_gws", "signif_gws_curve",
            "gws_significance", "signif")

  for (nm in cand) {
    if (!is.null(wv[[nm]]) && is.numeric(wv[[nm]])) {
      return(as.numeric(wv[[nm]]))
    }
  }

  NULL
}

#' Regrid GWS to target period vector
#'
#' @description
#' Interpolates GWS from wavelet analysis onto a target period grid.
#' Uses linear interpolation and fills edge NAs.
#'
#' @param wv List output from wavelet_spectral_analysis()
#' @param target_period Numeric vector of target periods
#' @param use_unmasked Logical. If TRUE and available, use gws_unmasked
#'
#' @return Numeric vector of GWS values on target_period grid
#'
#' @keywords internal
#' @export
gws_regrid <- function(wv, target_period, use_unmasked = FALSE) {
  # Get GWS
  g <- if (isTRUE(use_unmasked) && !is.null(wv$gws_unmasked)) {
    wv$gws_unmasked
  } else {
    wv$gws
  }

  # Get period
  p <- wv$gws_period

  # Validate
  if (is.null(g) || !is.numeric(g)) {
    stop("wavelet_spectral_analysis() missing numeric $gws.", call. = FALSE)
  }
  if (is.null(p) || !is.numeric(p)) {
    stop("wavelet_spectral_analysis() missing numeric $gws_period.", call. = FALSE)
  }

  g <- as.numeric(g)
  p <- as.numeric(p)

  # Remove non-finite values
  ok <- is.finite(p) & is.finite(g)
  p <- p[ok]
  g <- g[ok]

  # Handle edge case: too few points
  if (length(g) < 2L) {
    out <- rep(if (length(g) == 1L) g[1] else NA_real_, length(target_period))
    return(fill_nearest(out))
  }

  # Interpolate
  out <- stats::approx(x = p, y = g, xout = target_period, rule = 2)$y
  fill_nearest(out)
}

# ==============================================================================
# Original wavelet_helpers.R functions below
# ==============================================================================

#' Morlet Wavelet in Fourier Domain
#'
#' @description
#' Computes the Morlet wavelet function in Fourier space for a given scale.
#' This is the "daughter" wavelet used in continuous wavelet transform analysis.
#' The function applies normalization following Torrence & Compo (1998).
#'
#' @param k Numeric vector. Wave number vector (angular frequency).
#' @param s Numeric scalar. Scale parameter (inverse of frequency).
#' @param k0 Numeric scalar. Omega0 parameter for the Morlet wavelet (default = 6).
#'   This value balances time and frequency localization.
#'
#' @return Numeric vector (complex). The Morlet wavelet in Fourier space,
#'   with length equal to length(k).
#'
#' @details
#' The Morlet wavelet is defined in Fourier space as:
#' psi(s k) = pi^(-1/4) * sqrt(s k_2 n) * exp(-(s k - k_0)^2 / 2)
#' where the exponential term is applied only to positive frequencies.
#'
#' The normalization factor ensures energy conservation and proper inverse
#' transform reconstruction.
#'
#' @references
#' Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet Analysis.
#' Bulletin of the American Meteorological Society, 79(1), 61-78.
#'
#' @keywords internal
#' @export
morlet_wavelet <- function(k, s, k0 = 6) {

  # Basic validation
  if (!is.numeric(k) || !is.numeric(s) || !is.numeric(k0)) stop("All inputs must be numeric")
  if (length(s) != 1 || s <= 0) stop("'s' must be a positive scalar")
  if (length(k0) != 1 || k0 <= 0) stop("'k0' must be a positive scalar")
  if (length(k) < 2) stop("'k' must have length >= 2")

  if (abs(k[1]) > 1e-10) stop("'k' must start at zero frequency (standard FFT ordering)")
  if (k[2] <= 0) stop("'k[2]' must be positive (frequency spacing)")
  if (any(!is.finite(k))) stop("'k' contains non-finite values (NA, NaN, or Inf)")

  nn <- length(k)
  z <- as.numeric(k > 0)
  expnt <- -((s * k - k0)^2 / 2) * z

  # TC98 Eq. 4 normalization (TWO-SIDED / standard)
  norm <- sqrt(s * k[2]) * (pi^(-0.25)) * sqrt(nn)

  daughter <- norm * exp(expnt) * z
  daughter
}

#' Compute Morlet Wavelet Parameters
#'
#' @description
#' Calculates key parameters for the Morlet wavelet transform, including
#' the Fourier factor (scale-to-period conversion), cone of influence (COI),
#' and minimum degrees of freedom.
#'
#' @param k0 Numeric scalar. Omega0 parameter for the Morlet wavelet (default = 6).
#'   Standard choice is 6, which gives ~6 oscillations per wavelet.
#'
#' @return Named numeric vector with three elements:
#'   fourier_factor: Conversion factor from wavelet scale to Fourier period
#'   coi: Cone of influence e-folding time
#'   dofmin: Minimum degrees of freedom (2 for Morlet)
#'
#' @details
#' The Fourier factor allows conversion between wavelet scale s and
#' equivalent Fourier period T:
#' T = (4*pi / (k_0 + sqrt(2 + k_0^2))) * s
#'
#' The cone of influence (COI) defines the e-folding time for edge effects:
#' COI = fourier_factor / sqrt(2)
#'
#' Values outside the COI are subject to edge artifacts and should be
#' interpreted with caution.
#'
#' @references
#' Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet Analysis.
#' Bulletin of the American Meteorological Society, 79(1), 61-78.
#' See Table 1 and Section 3f.
#'
#' @examples
#' # Standard Morlet wavelet parameters
#' params <- morlet_parameters(k0 = 6)
#' print(params)
#' # fourier_factor: 1.03
#' # coi: 0.73
#' # dofmin: 2
#'
#' @keywords internal
#' @export
morlet_parameters <- function(k0 = 6) {

  # Input validation
  if (!is.numeric(k0) || length(k0) != 1 || k0 <= 0) {
    stop("'k0' must be a positive numeric scalar")
  }

  # Fourier factor: converts scale to period (TC98 Table 1)
  fourier_factor <- (4 * pi) / (k0 + sqrt(2 + k0^2))

  # Cone of influence e-folding time (TC98 Table 1)
  coi <- fourier_factor / sqrt(2)

  # Degrees of freedom for Morlet wavelet (TC98 Table 1)
  dofmin <- 2

  return(c(
    fourier_factor = fourier_factor,
    coi = coi,
    dofmin = dofmin
  ))
}

#' Extract and Reconstruct Wavelet Components from Significant Periods
#'
#' @description
#' Reconstructs time series components from the wavelet transform for
#' user-specified significant periods. Optionally includes a residual component
#' representing all non-significant scales, ensuring complete signal decomposition.
#'
#' @param wave Complex matrix. Wavelet transform coefficients with dimensions
#'   (scales x time). Output from continuous wavelet transform.
#' @param signif_periods Integer vector. Indices of significant period scales
#'   to reconstruct. Each index corresponds to a row in wave.
#' @param scale Numeric vector. Scale values corresponding to rows of wave.
#' @param dj Numeric scalar. Scale resolution parameter used in wavelet transform
#'   (default = 0.25). Smaller values give finer resolution.
#' @param dt Numeric scalar. Time step of the original series (default = 1).
#' @param variable_sd Numeric scalar. Standard deviation of the original time series
#'   before standardization. Used to restore original units.
#' @param variable_mean Numeric scalar. Mean of the original time series
#'   before standardization (default = 0). Added back to noise component to ensure
#'   complete reconstruction: rowSums(output) = original_signal.
#' @param Cdelta Numeric scalar. Reconstruction constant for Morlet wavelet
#'   (default = 0.776). Depends on wavelet type; see Torrence & Compo (1998) Table 2.
#' @param w0_0 Numeric scalar. Normalization constant psi_0(0) = pi^(-1/4)
#'   for Morlet wavelet (default = pi^(-1/4)).
#' @param include_residual Logical. If TRUE (default), includes a "Noise" component
#'   representing the sum of all non-significant scales, ensuring complete
#'   signal decomposition where rowSums(output) equals the reconstructed signal.
#'
#' @return A matrix with columns:
#'   Component_1, Component_2, ...: One column per significant period
#'   Noise: (if include_residual = TRUE) Sum of all non-significant scales
#'   Number of rows equals the length of the original time series.
#'   Property: rowSums(output) approximately equals the original signal.
#'
#' @details
#' The reconstruction formula (Torrence & Compo 1998, Eq. 11) is:
#'
#' x_n = (delta_j * sqrt(delta_t)) / (C_delta * psi_0(0)) *
#'       sum over j of Re(W_n(s_j)) / sqrt(s_j)
#'
#' This function performs a complete wavelet decomposition by:
#' 1. Extracting each significant period as a separate component
#' 2. Summing all remaining (non-significant) scales into a "Noise" component
#' 3. Ensuring: original_signal = sum(significant_components) + Noise
#'
#' This differs from partial reconstruction, which would only sum significant
#' scales and lose information from non-significant scales.
#'
#' @references
#' Torrence, C. and Compo, G.P. (1998). A Practical Guide to Wavelet Analysis.
#' Bulletin of the American Meteorological Society, 79(1), 61-78.
#' See Section 3i and Table 2.
#'
#' @examples
#' \dontrun{
#' # After performing wavelet analysis
#' components <- extract_wavelet_components(
#'   wave = wave_transform,
#'   signif_periods = c(5, 10, 20),  # indices of significant scales
#'   scale = scale_vector,
#'   dj = 0.25,
#'   dt = 1,
#'   variable_sd = sd(original_data),
#'   variable_mean = mean(original_data),
#'   Cdelta = 0.776,
#'   include_residual = TRUE
#' )
#'
#' # Verify complete reconstruction
#' reconstructed <- rowSums(components)
#' plot(original_data)
#' lines(reconstructed, col = "red")
#'
#' # Components include:
#' # - Component_1: First significant period
#' # - Component_2: Second significant period
#' # - Component_3: Third significant period
#' # - Noise: All non-significant scales combined + mean
#' }
#'
#' @keywords internal
#' @export
extract_wavelet_components <- function(wave,
                                       signif_periods,
                                       scale,
                                       dj = 0.25,
                                       dt = 1,
                                       variable_sd,
                                       variable_mean = 0,
                                       Cdelta = 0.776,
                                       w0_0 = pi^(-1/4),
                                       include_residual = TRUE) {

  if (!is.matrix(wave) && !is.array(wave)) stop("'wave' must be a matrix or array")
  if (!is.numeric(signif_periods) || length(signif_periods) == 0) stop("'signif_periods' must be a non-empty numeric vector")
  if (!all(signif_periods == as.integer(signif_periods))) stop("'signif_periods' must contain integer indices")
  if (any(signif_periods < 1) || any(signif_periods > nrow(wave))) stop("'signif_periods' indices must be between 1 and nrow(wave)")
  if (anyDuplicated(signif_periods)) {
    warning("'signif_periods' contains duplicate indices; duplicates will be removed")
    signif_periods <- unique(signif_periods)
  }

  if (!is.numeric(scale) || length(scale) != nrow(wave)) stop("'scale' must be numeric with length equal to nrow(wave)")
  if (!is.numeric(dj) || length(dj) != 1 || dj <= 0) stop("'dj' must be a positive scalar")
  if (!is.numeric(dt) || length(dt) != 1 || dt <= 0) stop("'dt' must be a positive scalar")
  if (!is.numeric(variable_sd) || length(variable_sd) != 1 || variable_sd <= 0) stop("'variable_sd' must be a positive scalar")
  if (!is.numeric(variable_mean) || length(variable_mean) != 1) stop("'variable_mean' must be a numeric scalar")
  if (!is.numeric(Cdelta) || length(Cdelta) != 1 || Cdelta <= 0) stop("'Cdelta' must be a positive scalar")
  if (!is.numeric(w0_0) || length(w0_0) != 1 || w0_0 <= 0) stop("'w0_0' must be a positive scalar")
  if (!is.logical(include_residual) || length(include_residual) != 1) stop("'include_residual' must be a single logical value")

  num_periods <- length(signif_periods)
  n_time <- ncol(wave)
  n_scales <- nrow(wave)

  n_components <- num_periods + (if (include_residual) 1L else 0L)
  COMPS <- matrix(0, nrow = n_time, ncol = n_components)

  # TC98 Eq. 11 reconstruction factor (TWO-SIDED / standard)
  recon_fac <- variable_sd * (dj * sqrt(dt) / (Cdelta * w0_0))

  inv_sqrt_scale <- 1 / sqrt(scale)
  Ww_all <- Re(wave) * inv_sqrt_scale

  for (i in seq_len(num_periods)) {
    j <- signif_periods[i]
    COMPS[, i] <- recon_fac * Ww_all[j, ]
  }

  if (include_residual) {
    nonsig_idx <- setdiff(seq_len(n_scales), signif_periods)

    if (length(nonsig_idx) > 0) {
      COMPS[, num_periods + 1L] <- recon_fac * colSums(Ww_all[nonsig_idx, , drop = FALSE])
    } else {
      COMPS[, num_periods + 1L] <- 0
      warning("All scales are marked as significant; residual component will be zero.")
    }

    COMPS[, num_periods + 1L] <- COMPS[, num_periods + 1L] + variable_mean
  }

  colnames(COMPS) <- c(
    paste0("Component_", seq_len(num_periods)),
    if (include_residual) "Noise" else NULL
  )

  attr(COMPS, "signif_periods") <- signif_periods
  attr(COMPS, "n_scales") <- n_scales
  attr(COMPS, "n_significant") <- num_periods
  attr(COMPS, "reconstruction_complete") <- include_residual

  COMPS
}


#' Wavelet Autoregressive Modeling (WARM)
#'
#' Simulates synthetic time series by modeling each wavelet component (signal or noise)
#' with an ARIMA model, and then summing the simulated components. Optionally enforces
#' variance matching to preserve statistical properties of the original components.
#'
#' @param wavelet.components A list or matrix where each column (or list element) is a
#'   numeric vector corresponding to a wavelet component (low-frequency signal or noise).
#' @param sim.year.num Integer. Desired length (number of years or timesteps) of each
#'   simulated series.
#' @param sim.num Integer. Number of synthetic series to produce. Default: 1000.
#' @param seed Optional. Integer random seed for reproducibility.
#' @param match.variance Logical. If TRUE, rescale simulated components to match the
#'   variance of the original components (default: TRUE). Recommended for preserving
#'   statistical properties.
#' @param variance.tolerance Numeric. Relative tolerance for variance matching
#'   (default: 0.1 = 10\%). Only applies if match.variance = TRUE.
#' @param check.diagnostics Logical. If TRUE, perform basic ARIMA model diagnostics
#'   and issue warnings if models appear inadequate (default: FALSE).
#' @param verbose Logical. If TRUE, emit informative logger::log_info() messages
#'   (default: TRUE). If FALSE, suppresses all logger::log_info() output from this function.
#'
#' @return A matrix of dimension \code{sim.year.num} x \code{sim.num}, where each
#'   column is a synthetic time series realization.
#'
#' @importFrom stats sd simulate
#' @importFrom forecast auto.arima
#' @export
wavelet_arima <- function(wavelet.components = NULL,
                          sim.year.num = NULL,
                          sim.num = 1000,
                          seed = NULL,
                          match.variance = TRUE,
                          variance.tolerance = 0.1,
                          check.diagnostics = FALSE,
                          verbose = TRUE) {

  # ============================================================================
  # Input Validation
  # ============================================================================

  if (is.null(wavelet.components)) {
    stop("Input 'wavelet.components' must not be NULL.")
  }

  if (is.null(sim.year.num)) {
    stop("Input 'sim.year.num' must be specified.")
  }

  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be logical (TRUE/FALSE).")
  }

  if (!.is_int_scalar(sim.year.num) || sim.year.num < 1L) {
    stop("'sim.year.num' must be a positive integer.")
  }
  sim.year.num <- as.integer(sim.year.num)

  if (!.is_int_scalar(sim.num) || sim.num < 1L) {
    stop("'sim.num' must be a positive integer.")
  }
  sim.num <- as.integer(sim.num)

  if (!is.logical(match.variance) || length(match.variance) != 1L) {
    stop("'match.variance' must be logical (TRUE/FALSE).")
  }

  if (!is.numeric(variance.tolerance) || length(variance.tolerance) != 1L ||
      !is.finite(variance.tolerance) || variance.tolerance < 0 || variance.tolerance > 1) {
    stop("'variance.tolerance' must be between 0 and 1.")
  }

  if (!is.logical(check.diagnostics) || length(check.diagnostics) != 1L) {
    stop("'check.diagnostics' must be logical (TRUE/FALSE).")
  }

  # ============================================================================
  # Efficient component list conversion
  # ============================================================================

  if (is.matrix(wavelet.components)) {
    if (!is.numeric(wavelet.components)) {
      stop("Matrix 'wavelet.components' must be numeric.")
    }
    ncomp <- ncol(wavelet.components)
    comp_list <- vector("list", ncomp)
    for (k in seq_len(ncomp)) comp_list[[k]] <- wavelet.components[, k]

  } else if (is.data.frame(wavelet.components)) {
    col_types <- vapply(wavelet.components, is.numeric, logical(1))
    if (!all(col_types)) stop("All columns of 'wavelet.components' must be numeric.")
    ncomp <- ncol(wavelet.components)
    comp_list <- vector("list", ncomp)
    for (k in seq_len(ncomp)) comp_list[[k]] <- wavelet.components[[k]]

  } else if (is.list(wavelet.components)) {
    elem_types <- vapply(wavelet.components, is.numeric, logical(1))
    if (!all(elem_types)) stop("All elements of 'wavelet.components' must be numeric vectors.")
    ncomp <- length(wavelet.components)
    comp_list <- wavelet.components

  } else {
    stop("'wavelet.components' must be a matrix, data.frame, or list of numeric vectors.")
  }

  # --- FIX 2: explicit NA check inside each component vector (before sd/mean)
  na_comp <- vapply(comp_list, function(x) anyNA(x), logical(1))
  if (any(na_comp)) {
    bad <- which(na_comp)
    stop(
      "Missing values detected in wavelet component(s): ",
      paste(bad, collapse = ", "),
      ". Remove/impute NAs before calling wavelet_arima().",
      call. = FALSE
    )
  }

  # ============================================================================
  # RNG State Management
  # ============================================================================

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be NULL or a single finite number.", call. = FALSE)
    }

    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- .Random.seed
      has_seed <- TRUE
    } else {
      has_seed <- FALSE
    }

    on.exit({
      if (has_seed) .Random.seed <<- old_seed
    }, add = TRUE)
  }

  # ============================================================================
  # Pre-allocate output matrix
  # ============================================================================

  output <- matrix(0, nrow = sim.year.num, ncol = sim.num)
  variance_corrections <- logical(ncomp)

  # ============================================================================
  # Component Simulation Loop
  # ============================================================================

  for (k in seq_len(ncomp)) {

    component <- as.numeric(comp_list[[k]])
    n_obs <- length(component)

    comp_sd <- stats::sd(component)
    if (!is.finite(comp_sd)) {
      stop("Component ", k, " has non-finite standard deviation.", call. = FALSE)
    }

    if (comp_sd < 1e-10) {
      warning(
        "Component ", k, " is essentially constant. ",
        "Returning constant values without ARIMA modeling.",
        call. = FALSE
      )
      output <- output + mean(component)
      next
    }

    if (n_obs < 10) {
      warning(
        "Component ", k, " has only ", n_obs,
        " observations. ARIMA modeling may be unreliable.",
        call. = FALSE
      )
    }

    comp_mean <- mean(component)
    centered <- component - comp_mean
    target_sd <- comp_sd

    if (!is.null(seed)) set.seed(as.integer(seed) + k * 1000L)

    MODEL <- tryCatch(
      {
        suppressPackageStartupMessages({
          forecast::auto.arima(
            centered,
            max.p = 2,
            max.q = 2,
            max.P = 0,
            max.Q = 0,
            stationary = TRUE,
            seasonal = FALSE
          )
        })
      },
      error = function(e) {
        warning(
          "ARIMA fitting failed for component ", k, ": ", e$message,
          "\nFalling back to resampling with replacement.",
          call. = FALSE
        )
        NULL
      }
    )

    if (is.null(MODEL)) {
      for (j in seq_len(sim.num)) {
        output[, j] <- output[, j] + sample(component, sim.year.num, replace = TRUE)
      }
      next
    }

    if (check.diagnostics) {
      if (!is.null(MODEL$convergence) && MODEL$convergence != 0) {
        warning("ARIMA model for component ", k, " did not converge properly.", call. = FALSE)
      }

      if (length(MODEL$residuals) > 1) {
        ljung_test <- tryCatch(
          stats::Box.test(MODEL$residuals, type = "Ljung-Box",
                          lag = min(10, length(MODEL$residuals) - 1)),
          error = function(e) NULL
        )

        if (!is.null(ljung_test) && is.finite(ljung_test$p.value) && ljung_test$p.value < 0.05) {
          warning(
            "Component ", k, ": Residuals show significant autocorrelation ",
            "(Ljung-Box p = ", round(ljung_test$p.value, 4), "). ",
            "Model may be inadequate.",
            call. = FALSE
          )
        }
      }
    }

    intercept <- 0
    if ("intercept" %in% names(MODEL$coef)) {
      intercept <- MODEL$coef[["intercept"]]
    } else if ("drift" %in% names(MODEL$coef)) {
      intercept <- MODEL$coef[["drift"]]
    }

    model_sd <- sqrt(MODEL$sigma2)
    needs_variance_check <- isTRUE(match.variance)

    for (j in seq_len(sim.num)) {

      simulated <- as.numeric(stats::simulate(MODEL, sim.year.num, sd = model_sd)) +
        intercept + comp_mean

      if (needs_variance_check) {
        sim_sd <- stats::sd(simulated)

        if (is.finite(sim_sd) && sim_sd > 0) {
          relative_diff <- abs(sim_sd - target_sd) / target_sd
          if (relative_diff > variance.tolerance) {
            sim_mean <- mean(simulated)
            simulated <- comp_mean + (simulated - sim_mean) * (target_sd / sim_sd)
            if (!variance_corrections[k]) variance_corrections[k] <- TRUE
          }
        }
      }

      output[, j] <- output[, j] + simulated
    }
  }

  if (verbose && match.variance && any(variance_corrections) &&
      requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("[WARM] Variance correction applied to component(s)")
  }

  output
}

