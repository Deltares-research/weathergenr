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

  # Detrend slope only (preserve mean level)
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

    if (!is.null(seed)) set.seed(seed)

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
