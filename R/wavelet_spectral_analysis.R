wavelet_spectral_analysis <- function(variable,
                                      signif.level = 0.90,
                                      noise.type = "red",
                                      period.lower.limit = 2,
                                      detrend = FALSE,
                                      return_diagnostics = TRUE,
                                      return_recon_error = FALSE) {

  # --- Input Validation ---
  if (!is.numeric(variable)) stop("variable must be numeric")
  if (anyNA(variable)) stop("variable contains missing values")
  if (length(variable) < 16) stop("variable must have at least 16 observations")
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

  # --- Wavelet Transform Analysis ---
  variable_org <- as.numeric(variable)

  if (detrend) {
    trend <- stats::fitted(stats::lm(variable_org ~ seq_along(variable_org)))
    variable_org <- variable_org - trend
  }

  variance1 <- stats::var(variable_org)
  n1 <- length(variable_org)

  sdx <- stats::sd(variable_org)
  if (!is.finite(sdx) || sdx <= 0) stop("variable has zero or non-finite standard deviation")

  variable_std <- (variable_org - mean(variable_org)) / sdx

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

  period <- as.numeric(fourier_factor * scale)

  coi_indices <- c(
    1e-5,
    seq_len(floor((n1 - 1) / 2)),
    rev(seq_len(ceiling((n1 - 1) / 2))),
    1e-5
  )
  coi <- as.numeric(coi_base * dt * coi_indices[1:n1])

  wave <- wave[, 1:n1, drop = FALSE]
  POWER <- abs(wave)^2

  # --- COI-masked Global Wavelet Spectrum (GWS) + plotting version ---
  coi_mask <- outer(period, coi, FUN = "<=")           # [n_scales x n_time] logical
  n_coi <- rowSums(coi_mask)                           # length = n_scales

  POWER_coi <- POWER
  POWER_coi[!coi_mask] <- NA_real_

  mean_power_coi <- rowMeans(POWER_coi, na.rm = TRUE)
  mean_power_coi[!is.finite(mean_power_coi)] <- NA_real_

  GWS_unmasked <- as.numeric(variance1 * rowMeans(POWER))
  GWS <- as.numeric(variance1 * mean_power_coi)

  # --- Significance Testing ---
  empir <- c(2, 0.776, 2.32, 0.60)
  dofmin <- empir[1]
  gamma_fac <- empir[3]


  # if (noise.type == "white") {
  #   lag1 <- 0
  # } else {
  #   lag1_raw <- stats::cor(variable_org[-n1], variable_org[-1])
  #   bias_correction <- (1 + 2 * lag1_raw) / (n1 - 2)
  #   lag1 <- lag1_raw - bias_correction
  #   lag1 <- max(min(lag1, 0.999), -0.999)
  # }

  if (noise.type == "white") {
    lag1 <- 0
  } else {

    x <- variable_org
    x <- x - mean(x)  # remove mean for AR estimation stability

    # Yule-Walker AR(1) estimate (standard, stable)
    # ar() may fail for very short/degenerate inputs; catch and fall back
    lag1 <- tryCatch({
      fit <- stats::ar(x, aic = FALSE, order.max = 1, method = "yw")
      as.numeric(fit$ar[1])
    }, error = function(e) {
      # Fallback: sample lag-1 autocorrelation (mean-removed)
      stats::cor(x[-length(x)], x[-1])
    })

    # Guard rails
    if (!is.finite(lag1)) lag1 <- 0
    lag1 <- max(min(lag1, 0.999), -0.999)
  }


  freq <- dt / period
  fft_theor <- (1 - lag1^2) / (1 - 2 * lag1 * cos(freq * 2 * pi) + lag1^2)

  chisquare <- stats::qchisq(signif.level, dofmin) / dofmin
  signif <- fft_theor * chisquare
  sigm <- sweep(POWER, 1, signif, FUN = "/")

  # --- Masked but uncorrected significance (COI mask only; plotting clarity) ---
  dof_masked_uncorrected <- dofmin * as.numeric(n_coi)
  dof_masked_uncorrected[is.finite(dof_masked_uncorrected) & dof_masked_uncorrected < dofmin] <- dofmin
  dof_masked_uncorrected[!is.finite(dof_masked_uncorrected)] <- NA_real_

  chisq_masked_uncorrected <- rep(NA_real_, length(dof_masked_uncorrected))
  ok_mu <- is.finite(dof_masked_uncorrected) & dof_masked_uncorrected > 0
  chisq_masked_uncorrected[ok_mu] <-
    stats::qchisq(signif.level, dof_masked_uncorrected[ok_mu]) / dof_masked_uncorrected[ok_mu]

  GWS_signif_masked_uncorrected <- rep(NA_real_, length(dof_masked_uncorrected))
  GWS_signif_masked_uncorrected[ok_mu] <-
    as.numeric(fft_theor[ok_mu] * variance1 * chisq_masked_uncorrected[ok_mu])

  # --- COI-consistent Neff based on wavelet decorrelation time (inference curve) ---
  Neff <- (as.numeric(n_coi) * dt) / (gamma_fac * scale)
  Neff[Neff < 1] <- NA_real_
  Neff[!is.finite(Neff)] <- NA_real_

  dof <- dofmin * Neff
  dof[is.finite(dof) & dof < dofmin] <- dofmin

  chisq <- rep(NA_real_, length(dof))
  ok <- is.finite(dof) & dof > 0
  chisq[ok] <- stats::qchisq(signif.level, dof[ok]) / dof[ok]

  GWS_signif <- rep(NA_real_, length(dof))
  GWS_signif[ok] <- as.numeric(fft_theor[ok] * variance1 * chisq[ok])

  # --- UNMASKED GWS significance (plotting) ---
  # Legacy curve (exact old behavior) kept ONLY for backward-comparison plots
  dof_unmasked_legacy <- n1 - scale
  dof_unmasked_legacy[dof_unmasked_legacy < 1] <- 1
  dof_unmasked_legacy <- dofmin * sqrt(1 + (dof_unmasked_legacy * dt / gamma_fac / scale)^2)
  dof_unmasked_legacy[dof_unmasked_legacy < dofmin] <- dofmin

  chisq_unmasked_legacy <- stats::qchisq(signif.level, dof_unmasked_legacy) / dof_unmasked_legacy
  GWS_signif_unmasked_legacy <- as.numeric(fft_theor * variance1 * chisq_unmasked_legacy)

  # Defensible unmasked curve: Neff based on wavelet decorrelation time over full record
  Neff_unmasked <- (n1 * dt) / (gamma_fac * scale)
  Neff_unmasked[Neff_unmasked < 1] <- NA_real_
  Neff_unmasked[!is.finite(Neff_unmasked)] <- NA_real_

  dof_unmasked <- dofmin * Neff_unmasked
  dof_unmasked[is.finite(dof_unmasked) & dof_unmasked < dofmin] <- dofmin

  chisq_unmasked <- rep(NA_real_, length(dof_unmasked))
  ok_u <- is.finite(dof_unmasked) & dof_unmasked > 0
  chisq_unmasked[ok_u] <- stats::qchisq(signif.level, dof_unmasked[ok_u]) / dof_unmasked[ok_u]

  GWS_signif_unmasked <- rep(NA_real_, length(dof_unmasked))
  GWS_signif_unmasked[ok_u] <- as.numeric(fft_theor[ok_u] * variance1 * chisq_unmasked[ok_u])

  # --- Identify Significant Periods (use inference curves only) ---
  sig_periods <- which(
    is.finite(GWS) & is.finite(GWS_signif) &
      (GWS > GWS_signif) & (period > period.lower.limit))

  has_significance <- length(sig_periods) > 0

  out_recon_err <- NULL

  if (length(sig_periods) == 0) {

    signif_periods <- integer(0)
    COMPS <- matrix(variable_org, ncol = 1)
    colnames(COMPS) <- "NOISE"

    if (return_recon_error) {
      out_recon_err <- list(recon_max_abs = 0, recon_rmse = 0)
    }

  } else {

    sig_periods_grp <- split(sig_periods, cumsum(c(1, diff(sig_periods) != 1)))

    signif_periods <- as.integer(unlist(
      lapply(sig_periods_grp, function(x) x[which.max(GWS[x])]),
      use.names = FALSE))

    comps_mat <- extract_wavelet_components(
      wave = wave,
      signif_periods = signif_periods,
      scale = scale,
      dj = dj,
      dt = dt,
      variable_sd = stats::sd(variable_org),
      variable_mean = 0,
      Cdelta = 0.776,
      w0_0 = pi^(-1/4),
      include_residual = FALSE
    )

    period_labels <- round(period[signif_periods], 2)
    comp_names <- paste0("Period_", period_labels)
    colnames(comps_mat) <- comp_names

    noise <- variable_org - rowSums(comps_mat)
    COMPS <- cbind(comps_mat, NOISE = noise)
    colnames(COMPS) <- c(comp_names, "NOISE")

    if (return_recon_error) {
      recon_vec <- variable_org - rowSums(COMPS)
      out_recon_err <- list(
        recon_max_abs = max(abs(recon_vec)),
        recon_rmse = sqrt(mean(recon_vec^2))
      )
    }
  }

  # Explicit no-significance handling outputs
  signif_periods <- as.integer(signif_periods)
  has_significance <- length(signif_periods) > 0
  significance_status <- if (has_significance) "significant_scales_found" else "no_significant_scales"
  signif_period_values <- if (has_significance) period[signif_periods] else numeric(0)

  out <- list(
    GWS = GWS,
    GWS_unmasked = GWS_unmasked,
    GWS_signif = GWS_signif,
    GWS_signif_masked_uncorrected = GWS_signif_masked_uncorrected,
    GWS_signif_unmasked = GWS_signif_unmasked,
    GWS_signif_unmasked_legacy = GWS_signif_unmasked_legacy,
    GWS_period = period,
    signif_periods = signif_periods,
    wave = wave,
    power = POWER,
    coi = coi,
    sigm = sigm,
    COMPS = COMPS,
    COMPS_names = colnames(COMPS),
    GWS_n_coi = n_coi,
    GWS_Neff = Neff,
    GWS_Neff_unmasked = Neff_unmasked,
    has_significance = has_significance,
    significance_status = significance_status,
    signif_period_values = signif_period_values
  )

  if (return_recon_error) out$reconstruction_error <- out_recon_err

  if (return_diagnostics) {
    out$diagnostics <- list(
      lag1 = lag1,
      variance = variance1,
      n_original = n1,
      n_padded = n,
      dof = dof,                    # inference dof (masked + Neff)
      dof_unmasked = dof_unmasked,  # plotting dof (unmasked + Neff)
      scale = scale,
      fourier_factor = fourier_factor,
      n_coi = n_coi,
      GWS_Neff = Neff,
      GWS_Neff_unmasked = Neff_unmasked
    )
  }

  out
}
