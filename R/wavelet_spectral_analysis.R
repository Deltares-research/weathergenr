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

  # Zero-pad (unchanged; not part of Priority 0)
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
  GWS <- as.numeric(variance1 * rowMeans(POWER))

  # --- Significance Testing ---
  empir <- c(2, 0.776, 2.32, 0.60)
  dofmin <- empir[1]
  gamma_fac <- empir[3]

  if (noise.type == "white") {
    lag1 <- 0
  } else {
    lag1_raw <- stats::cor(variable_org[-n1], variable_org[-1])
    bias_correction <- (1 + 2 * lag1_raw) / (n1 - 2)
    lag1 <- lag1_raw - bias_correction
    lag1 <- max(min(lag1, 0.999), -0.999)
  }

  freq <- dt / period
  fft_theor <- (1 - lag1^2) / (1 - 2 * lag1 * cos(freq * 2 * pi) + lag1^2)

  chisquare <- stats::qchisq(signif.level, dofmin) / dofmin
  signif <- fft_theor * chisquare
  sigm <- sweep(POWER, 1, signif, FUN = "/")

  dof <- n1 - scale
  dof[dof < 1] <- 1
  dof <- dofmin * sqrt(1 + (dof * dt / gamma_fac / scale)^2)
  dof[dof < dofmin] <- dofmin

  chisquare_GWS <- stats::qchisq(signif.level, dof) / dof
  GWS_signif <- as.numeric(fft_theor * variance1 * chisquare_GWS)

  # --- Identify Significant Periods ---
  sig_periods <- which(GWS > GWS_signif & period > period.lower.limit)

  if (length(sig_periods) == 0) {
    signif_periods <- integer(0)

    # Stable contract: COMPS always exists and has nrow == length(variable)
    COMPS <- tibble::tibble(NOISE = variable_org)

    # Optional scalar closure diagnostics (trivially perfect)
    if (return_recon_error) out_recon_err <- list(recon_max_abs = 0, recon_rmse = 0)

  } else {
    sig_periods_grp <- split(sig_periods, cumsum(c(1, diff(sig_periods) != 1)))

    signif_periods <- as.integer(unlist(
      lapply(sig_periods_grp, function(x) x[which.max(GWS[x])]),
      use.names = FALSE
    ))

    # Priority-0 fix: no residual in helper; compute residual exactly once here
    comps_tb <- extract_wavelet_components(
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

    # Rename only component columns (no residual column exists here)
    period_labels <- round(period[signif_periods], 2)
    comp_names <- paste0("Period_", period_labels)
    names(comps_tb) <- comp_names

    # Residual once
    noise <- variable_org - rowSums(as.matrix(comps_tb))
    COMPS <- cbind(comps_tb, NOISE = noise)

    # Optional closure diagnostic
    out_recon_err <- NULL
    if (return_recon_error) {
      recon_vec <- variable_org - rowSums(as.matrix(COMPS[, c(comp_names, "NOISE"), drop = FALSE]))
      out_recon_err <- list(
        recon_max_abs = max(abs(recon_vec)),
        recon_rmse = sqrt(mean(recon_vec^2))
      )
    }


  }

  out <- list(
    GWS = GWS,
    GWS_signif = GWS_signif,
    GWS_period = period,
    signif_periods = signif_periods,
    wave = wave,
    power = POWER,
    coi = coi,
    sigm = sigm,
    COMPS = COMPS
  )

  if (return_recon_error) out$reconstruction_error <- out_recon_err

  if (return_diagnostics) {
    out$diagnostics <- list(
      lag1 = lag1,
      variance = variance1,
      n_original = n1,
      n_padded = n,
      dof = dof,
      scale = scale,
      fourier_factor = fourier_factor
    )
  }

  out
}
