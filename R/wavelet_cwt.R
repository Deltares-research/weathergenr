# ==============================================================================
# Continuous Wavelet Transform (CWT) Analysis
# ==============================================================================
#
# Morlet wavelet-based spectral analysis with COI-aware significance testing.
# Provides visualization and diagnostic capabilities for time-frequency analysis.
#
# For orthogonal component extraction suitable for ARIMA modeling, see
# wavelet_modwt.R which provides the analyze_wavelet_orthogonal() function.
#
# Contents:
# - Morlet wavelet functions (morlet_wavelet, morlet_parameters)
# - Utility functions (fill_nearest, extract_signif_curve, gws_regrid)
# - CWT component extraction (extract_wavelet_components)
# - Main analysis function (analyze_wavelet_spectrum)
#
# References:
# Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet analysis.
# Bulletin of the American Meteorological Society, 79(1), 61-78.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# Morlet Wavelet Functions
# ------------------------------------------------------------------------------

#' Morlet Wavelet in Fourier Domain
#'
#' @description
#' Computes the Morlet wavelet "daughter" in Fourier space for a given scale.
#' Normalization follows Torrence & Compo (1998).
#'
#' @param k Numeric vector. Angular frequencies in FFT ordering.
#' @param scale Numeric scalar. Wavelet scale (inverse frequency), must be positive.
#' @param k0 Numeric scalar. Morlet nondimensional frequency (default: 6).
#'
#' @return Numeric vector (complex). Morlet wavelet in Fourier space.
#'
#' @references
#' Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet analysis.
#' \emph{Bulletin of the American Meteorological Society}, 79(1), 61-78.
#'
#' @seealso \code{\link{morlet_parameters}}, \code{\link{analyze_wavelet_spectrum}}
#'
#' @keywords internal
#' @export
morlet_wavelet <- function(k, scale, k0 = 6) {

 if (!is.numeric(k) || !is.numeric(scale) || !is.numeric(k0)) stop("All inputs must be numeric")
 if (length(scale) != 1 || scale <= 0) stop("'scale' must be a positive scalar")
 if (length(k0) != 1 || k0 <= 0) stop("'k0' must be a positive scalar")
 if (length(k) < 2) stop("'k' must have length >= 2")

 if (abs(k[1]) > 1e-10) stop("'k' must start at zero frequency (standard FFT ordering)")
 if (k[2] <= 0) stop("'k[2]' must be positive (frequency spacing)")
 if (any(!is.finite(k))) stop("'k' contains non-finite values (NA, NaN, or Inf)")

 nn <- length(k)
 z <- as.numeric(k > 0)
 expnt <- -((scale * k - k0)^2 / 2) * z

 norm <- sqrt(scale * k[2]) * (pi^(-0.25)) * sqrt(nn)
 norm * exp(expnt) * z
}


#' Compute Morlet Wavelet Parameters
#'
#' @description
#' Computes key parameters for Morlet wavelets: Fourier factor (scale to period),
#' cone-of-influence (COI) e-folding time, and minimum degrees of freedom.
#'
#' @param k0 Numeric scalar. Morlet nondimensional frequency (default: 6).
#'
#' @return Named numeric vector: \code{fourier_factor}, \code{coi}, \code{dofmin}.
#'
#' @references
#' Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet analysis.
#' \emph{Bulletin of the American Meteorological Society}, 79(1), 61-78.
#'
#' @seealso \code{\link{morlet_wavelet}}, \code{\link{analyze_wavelet_spectrum}}
#'
#' @keywords internal
#' @export
morlet_parameters <- function(k0 = 6) {

 if (!is.numeric(k0) || length(k0) != 1 || k0 <= 0) {
   stop("'k0' must be a positive numeric scalar")
 }

 fourier_factor <- (4 * pi) / (k0 + sqrt(2 + k0^2))
 coi <- fourier_factor / sqrt(2)
 dofmin <- 2

 c(fourier_factor = fourier_factor, coi = coi, dofmin = dofmin)
}


# ------------------------------------------------------------------------------
# Utility Functions
# ------------------------------------------------------------------------------

#' Fill NA Values with the Nearest Non-NA Neighbor
#'
#' @description
#' Fills missing values in a numeric vector by propagating the nearest available
#' non-missing value. This is a deterministic forward-fill followed by a backward-fill.
#' It is mainly used to stabilize curves after interpolation (e.g., regridded spectra).
#'
#' @param x Numeric vector possibly containing \code{NA}.
#'
#' @return Numeric vector with missing values filled. If \code{x} is all \code{NA},
#'   it is returned unchanged.
#'
#' @keywords internal
#' @export
fill_nearest <- function(x) {
 if (all(is.na(x))) return(x)

 if (is.na(x[1])) {
   first <- which(!is.na(x))[1]
   if (!is.na(first)) x[1:(first - 1)] <- x[first]
 }

 for (i in 2:length(x)) {
   if (is.na(x[i])) x[i] <- x[i - 1]
 }

 for (i in (length(x) - 1):1) {
   if (is.na(x[i])) x[i] <- x[i + 1]
 }

 x
}


#' Extract Global-Spectrum Significance Curve
#'
#' @description
#' Extracts a global-spectrum significance curve from a wavelet analysis output list.
#' Intended as a small compatibility helper for downstream code that expects a single
#' numeric vector. If not found, returns \code{NULL}.
#'
#' @param wavelet List output from \code{\link{analyze_wavelet_spectrum}}.
#'
#' @return Numeric vector of significance values, or \code{NULL} if not found.
#'
#' @keywords internal
#' @export
extract_signif_curve <- function(wavelet) {
 cand <- c(
   "gws_signif",
   "gws_sig",
   "signif_gws",
   "signif_gws_curve",
   "gws_significance",
   "signif"
 )

 for (nm in cand) {
   if (!is.null(wavelet[[nm]]) && is.numeric(wavelet[[nm]])) {
     return(as.numeric(wavelet[[nm]]))
   }
 }

 NULL
}


#' Regrid Global Wavelet Spectrum to a Target Period Grid
#'
#' @description
#' Linearly interpolates the global wavelet spectrum from a wavelet analysis output onto
#' a target period grid and fills edge values using nearest-neighbor filling.
#'
#' @param wavelet List output from \code{\link{analyze_wavelet_spectrum}}.
#' @param target_period Numeric vector. Target periods for interpolation.
#' @param use_unmasked Logical. If \code{TRUE} and available, uses \code{wavelet$gws_unmasked};
#'   otherwise uses \code{wavelet$gws}.
#'
#' @return Numeric vector of interpolated GWS values on \code{target_period}.
#'
#' @keywords internal
#' @export
gws_regrid <- function(wavelet, target_period, use_unmasked = FALSE) {

 g <- if (isTRUE(use_unmasked) && !is.null(wavelet$gws_unmasked)) {
   wavelet$gws_unmasked
 } else {
   wavelet$gws
 }

 p <- wavelet$period

 if (is.null(g) || !is.numeric(g)) {
   stop("analyze_wavelet_spectrum() output missing numeric $gws.", call. = FALSE)
 }
 if (is.null(p) || !is.numeric(p)) {
   stop("analyze_wavelet_spectrum() output missing numeric $period.", call. = FALSE)
 }

 g <- as.numeric(g)
 p <- as.numeric(p)

 ok <- is.finite(p) & is.finite(g)
 p <- p[ok]
 g <- g[ok]

 if (length(g) < 2L) {
   out <- rep(if (length(g) == 1L) g[1] else NA_real_, length(target_period))
   return(fill_nearest(out))
 }

 out <- stats::approx(x = p, y = g, xout = target_period, rule = 2)$y
 fill_nearest(out)
}


# ------------------------------------------------------------------------------
# CWT Component Extraction
# ------------------------------------------------------------------------------

#' Reconstruct Wavelet Components from Selected Scales (CWT)
#'
#' @description
#' Reconstructs time-domain components from selected wavelet scales using the
#' Torrence & Compo (1998) inverse transform approximation. Optionally adds a residual
#' component representing non-selected scales to ensure additive closure.
#'
#' Note: CWT-reconstructed components are additive but NOT orthogonal due to scale
#' overlap. For orthogonal decomposition, use \code{\link{extract_modwt_components}}.
#'
#' @param wave Complex matrix (scales x time). Wavelet coefficients.
#' @param signif_periods Integer vector. Scale indices to reconstruct.
#' @param scale Numeric vector. Scale values aligned with \code{nrow(wave)}.
#' @param dj Numeric scalar. Scale resolution used in the transform.
#' @param dt Numeric scalar. Time step of the series.
#' @param series_sd Numeric scalar. Standard deviation of the original series (pre-standardization).
#' @param series_mean Numeric scalar. Mean of the original series; added back to residual.
#' @param Cdelta Numeric scalar. Morlet reconstruction constant (default: 0.776).
#' @param w0_0 Numeric scalar. psi_0(0) constant (default: pi^(-1/4)).
#' @param include_residual Logical. If \code{TRUE}, adds a residual (non-selected) component.
#'
#' @return Numeric matrix with one column per selected component and, if requested,
#'   a final \code{"Noise"} residual column.
#'
#' @references
#' Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet analysis.
#' \emph{Bulletin of the American Meteorological Society}, 79(1), 61-78.
#'
#' @seealso \code{\link{analyze_wavelet_spectrum}}, \code{\link{extract_modwt_components}}
#'
#' @keywords internal
#' @export
extract_wavelet_components <- function(
   wave,
   signif_periods,
   scale,
   dj = 0.25,
   dt = 1,
   series_sd,
   series_mean = 0,
   Cdelta = 0.776,
   w0_0 = pi^(-1/4),
   include_residual = TRUE
) {

 if (!is.matrix(wave) && !is.array(wave)) stop("'wave' must be a matrix or array")
 if (!is.numeric(signif_periods) || length(signif_periods) == 0) stop("'signif_periods' must be non-empty")
 if (!all(signif_periods == as.integer(signif_periods))) stop("'signif_periods' must contain integer indices")
 if (any(signif_periods < 1) || any(signif_periods > nrow(wave))) stop("'signif_periods' out of bounds")
 if (anyDuplicated(signif_periods)) {
   warning("'signif_periods' contains duplicates; duplicates removed", call. = FALSE)
   signif_periods <- unique(signif_periods)
 }

 if (!is.numeric(scale) || length(scale) != nrow(wave)) stop("'scale' must match nrow(wave)")
 if (!is.numeric(dj) || length(dj) != 1 || dj <= 0) stop("'dj' must be positive")
 if (!is.numeric(dt) || length(dt) != 1 || dt <= 0) stop("'dt' must be positive")
 if (!is.numeric(series_sd) || length(series_sd) != 1 || series_sd <= 0) stop("'series_sd' must be positive")
 if (!is.numeric(series_mean) || length(series_mean) != 1) stop("'series_mean' must be scalar")
 if (!is.numeric(Cdelta) || length(Cdelta) != 1 || Cdelta <= 0) stop("'Cdelta' must be positive")
 if (!is.numeric(w0_0) || length(w0_0) != 1 || w0_0 <= 0) stop("'w0_0' must be positive")
 if (!is.logical(include_residual) || length(include_residual) != 1) stop("'include_residual' must be logical scalar")

 num_periods <- length(signif_periods)
 n_time <- ncol(wave)
 n_scales <- nrow(wave)

 n_components <- num_periods + (if (include_residual) 1L else 0L)
 comps <- matrix(0, nrow = n_time, ncol = n_components)

 recon_fac <- series_sd * (dj * sqrt(dt) / (Cdelta * w0_0))

 inv_sqrt_scale <- 1 / sqrt(scale)
 Ww_all <- Re(wave) * inv_sqrt_scale

 for (i in seq_len(num_periods)) {
   j <- signif_periods[i]
   comps[, i] <- recon_fac * Ww_all[j, ]
 }

 if (include_residual) {
   nonsig_idx <- setdiff(seq_len(n_scales), signif_periods)

   if (length(nonsig_idx) > 0) {
     comps[, num_periods + 1L] <- recon_fac * colSums(Ww_all[nonsig_idx, , drop = FALSE])
   } else {
     comps[, num_periods + 1L] <- 0
     warning("All scales are marked significant; residual component is zero.", call. = FALSE)
   }

   comps[, num_periods + 1L] <- comps[, num_periods + 1L] + series_mean
 }

 colnames(comps) <- c(
   paste0("Component_", seq_len(num_periods)),
   if (include_residual) "Noise" else NULL
 )

 attr(comps, "signif_periods") <- signif_periods
 attr(comps, "n_scales") <- n_scales
 attr(comps, "n_significant") <- num_periods
 attr(comps, "reconstruction_complete") <- include_residual

 comps
}


# ------------------------------------------------------------------------------
# Main CWT Analysis Function
# ------------------------------------------------------------------------------

#' Continuous Wavelet Spectral Analysis with COI-Aware Significance
#'
#' @description
#' Performs continuous wavelet transform (CWT) analysis using a Morlet mother wavelet
#' (Torrence & Compo, 1998) and returns wavelet power, cone of influence (COI),
#' global wavelet spectrum (GWS), and significance thresholds. The implementation includes:
#' \itemize{
#'   \item COI-aware global spectrum estimation,
#'   \item effective sample size (\code{neff}) based significance testing,
#'   \item AR(1) red-noise background estimation (Yule-Walker),
#'   \item optional scale-component reconstruction in \code{mode = "complete"}.
#' }
#'
#' This function is intended for detection and diagnostic analysis of dominant
#' time scales. Reconstructed scale components are additive but not orthogonal.
#' For truly orthogonal decomposition suitable for ARIMA modeling, see
#' \code{\link{analyze_wavelet_orthogonal}} which combines CWT diagnostics with
#' MODWT-based orthogonal component extraction.
#'
#' @param series Numeric vector. Input time series (regularly spaced, no missing values).
#'   Minimum length is 16 observations.
#' @param signif Numeric scalar in (0, 1). Significance level for wavelet power and GWS
#'   tests. Default is 0.90.
#' @param noise Character. Background noise model for significance testing. One of
#'   \code{"white"} or \code{"red"} (default). For red noise, an AR(1) model is estimated.
#' @param min_period Numeric scalar. Minimum Fourier period (in time units) used when
#'   identifying significant scales (default: 2).
#' @param detrend Logical. If \code{TRUE}, removes a linear-trend slope only (mean preserved)
#'   prior to analysis (default: \code{FALSE}).
#' @param mode Character. \code{"fast"} (default) returns essential outputs only; \code{"complete"}
#'   additionally returns pointwise significance and (optionally) reconstruction products.
#' @param diagnostics Logical. If \code{TRUE}, returns an additional \code{diagnostics} list.
#'   When \code{mode = "complete"}, diagnostics also include reconstruction error metrics
#'   when reconstructed components are computed.
#' @param lag1_ci Logical. If \code{TRUE}, computes a bootstrap confidence interval for the AR(1)
#'   coefficient (diagnostic only; does not affect inference).
#' @param lag1_ci_level Numeric scalar in (0, 1). Confidence level for the AR(1) bootstrap interval.
#' @param lag1_boot_n Integer. Number of bootstrap replicates for the AR(1) interval (>= 50).
#' @param seed Optional numeric scalar. Random seed used for bootstrap diagnostics.
#' @param warn_neff Logical. If \code{TRUE}, warns when effective sample size is small for
#'   most scales.
#' @param neff_warn_min Numeric scalar. Threshold below which neff is considered small.
#' @param neff_warn_frac Numeric scalar in (0, 1). Fraction of scales with \code{neff < neff_warn_min}
#'   required to trigger a warning.
#'
#' @return A list. Always includes (at minimum):
#' \describe{
#'   \item{gws}{Numeric vector. COI-masked global wavelet spectrum.}
#'   \item{gws_unmasked}{Numeric vector. Unmasked global wavelet spectrum (plotting only).}
#'   \item{period}{Numeric vector. Fourier periods corresponding to wavelet scales.}
#'   \item{gws_signif}{Numeric vector. COI- and neff-corrected GWS significance threshold (inference).}
#'   \item{gws_signif_unmasked}{Numeric vector. Unmasked GWS significance threshold (plotting only).}
#'   \item{has_significance}{Logical. Whether any significant scales were detected.}
#'   \item{signif_periods}{Integer vector. Indices of significant scales retained (one per contiguous band).}
#'   \item{coi}{Numeric vector. Cone of influence in time units.}
#'   \item{power}{Numeric matrix. Wavelet power spectrum (scales x time).}
#' }
#'
#' In \code{mode = "complete"}, additional fields may include:
#' \describe{
#'   \item{power_coi}{Numeric matrix. COI-masked wavelet power.}
#'   \item{sigm}{Numeric matrix. Pointwise significance ratio.}
#'   \item{sigm_coi}{Numeric matrix. COI-masked pointwise significance ratio.}
#'   \item{power_signif_coi}{Logical matrix. COI-masked pointwise significance mask.}
#'   \item{wave}{Complex matrix. Wavelet coefficients (scales x time).}
#'   \item{comps}{Numeric matrix. Reconstructed significant-scale components plus residual (if computed).}
#'   \item{comps_names}{Character vector. Component names.}
#'   \item{gws_n_coi}{Numeric vector. Number of time points inside COI per scale.}
#'   \item{neff}{Numeric vector. Effective sample size per scale (masked).}
#'   \item{neff_unmasked}{Numeric vector. Effective sample size per scale (unmasked).}
#'   \item{diagnostics}{List. Returned only when \code{diagnostics = TRUE}.}
#' }
#'
#' @references
#' Torrence, C., & Compo, G. P. (1998). A practical guide to wavelet analysis.
#' \emph{Bulletin of the American Meteorological Society}, 79(1), 61-78.
#'
#' @seealso \code{\link{morlet_wavelet}}, \code{\link{extract_wavelet_components}},
#'   \code{\link{analyze_wavelet_orthogonal}},
#'   \code{\link{plot_wavelet_power}}, \code{\link{plot_wavelet_global_spectrum}}
#'
#' @export
analyze_wavelet_spectrum <- function(
   series,
   signif = 0.90,
   noise = "red",
   min_period = 2,
   detrend = FALSE,
   mode = c("fast", "complete"),
   diagnostics = FALSE,
   lag1_ci = FALSE,
   lag1_ci_level = 0.95,
   lag1_boot_n = 500,
   seed = NULL,
   warn_neff = FALSE,
   neff_warn_min = 5,
   neff_warn_frac = 0.60
) {

 # --- Input Validation ---
 if (!is.numeric(series)) stop("series must be numeric")
 if (anyNA(series)) stop("series contains missing values")
 if (length(series) < 16) stop("series must have at least 16 observations")

 mode <- match.arg(mode)

 if (!(noise %in% c("white", "red"))) stop("noise must be 'white' or 'red'")

 if (!is.numeric(signif) || length(signif) != 1L || signif <= 0 || signif >= 1) {
   stop("signif must be between 0 and 1")
 }

 if (!is.numeric(min_period) || length(min_period) != 1L || min_period < 0) {
   stop("min_period must be a non-negative number")
 }

 if (!is.logical(diagnostics) || length(diagnostics) != 1L) {
   stop("diagnostics must be TRUE/FALSE")
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

 compute_full <- (mode == "complete")

 # --- Wavelet Transform Analysis ---
 series_org <- as.numeric(series)

 if (isTRUE(detrend)) {
   tt <- seq_along(series_org)
   fit <- stats::lm(series_org ~ tt)
   b <- stats::coef(fit)[2]
   if (!is.finite(b)) b <- 0
   series_org <- series_org - b * (tt - mean(tt))
 }

 variance1 <- stats::var(series_org)
 n1 <- length(series_org)

 sdx <- stats::sd(series_org)
 if (!is.finite(sdx) || sdx <= 0) stop("series has zero or non-finite standard deviation")

 series_mean <- mean(series_org)

 # Standardize (for transform only)
 series_std <- (series_org - series_mean) / sdx

 # Zero-pad
 base2 <- floor(log2(n1) + 0.4999)
 series_pad <- c(series_std, rep(0, (2^(base2 + 1) - n1)))
 n <- length(series_pad)

 dt <- 1
 dj <- 0.25
 s0 <- 2 * dt
 J <- floor((1 / dj) * log2(n1 * dt / s0))
 scale <- s0 * 2^((0:J) * dj)

 k <- c(0:(floor(n / 2)), -rev(1:floor((n - 1) / 2))) * ((2 * pi) / (n * dt))
 f <- stats::fft(series_pad)

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

 # --- COI mask ---
 coi_mask <- outer(period, coi, FUN = "<=")
 n_coi <- rowSums(coi_mask)

 # --- Significance Testing (pointwise) ---
 empir <- c(2, 0.776, 2.32, 0.60)
 dofmin <- empir[1]
 gamma_fac <- empir[3]

 # --- Lag1 estimation ---
 x_centered <- series_org - mean(series_org)

 if (noise == "white") {
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
 if (lag1_ci && noise == "red") {

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
   chisquare <- stats::qchisq(signif, dofmin) / dofmin
   signif_pt <- fft_theor * chisquare
   sigm <- sweep(power, 1, signif_pt, FUN = "/")

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

 # --- COI-consistent neff ---
 neff <- (n_coi * dt) / (gamma_fac * scale)
 neff[neff < 1] <- NA_real_
 neff[!is.finite(neff)] <- NA_real_

 if (warn_neff) {
   ok_neff <- is.finite(neff)
   frac_small <- if (any(ok_neff)) mean(neff[ok_neff] < neff_warn_min) else 1

   if (!any(ok_neff) || frac_small >= neff_warn_frac) {
     warning(
       sprintf(
         "Wavelet significance may be unreliable: %.0f%% of scales have neff < %.2f (n=%s). Consider longer series or coarser scales/dj.",
         100 * frac_small, neff_warn_min, format(n1, big.mark = ",")
       ),
       call. = FALSE
     )
   }
 }

 dof <- dofmin * neff
 dof[is.finite(dof) & dof < dofmin] <- dofmin

 chisq <- rep(NA_real_, length(dof))
 ok <- is.finite(dof) & dof > 0
 chisq[ok] <- stats::qchisq(signif, dof[ok]) / dof[ok]

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
 chisq_unmasked[ok_u] <- stats::qchisq(signif, dof_unmasked[ok_u]) / dof_unmasked[ok_u]
 chisq_unmasked[!ok_u] <- 1

 gws_signif_unmasked <- rep(NA_real_, length(dof_unmasked))
 gws_signif_unmasked[ok_u] <- fft_theor[ok_u] * variance1 * chisq_unmasked[ok_u]
 gws_signif_unmasked[!ok_u] <- fft_theor[!ok_u] * variance1 * chisq_unmasked[!ok_u]

 # --- Identify Significant Periods (use inference curves only) ---
 sig_periods <- which(
   is.finite(gws) & is.finite(gws_signif) &
     (gws > gws_signif) & (period > min_period)
 )

 # --- Reconstruction (only in complete mode); diagnostics include recon error ---
 comps <- NULL
 comps_names <- NULL
 recon_err <- NULL

 if (compute_full) {

   if (length(sig_periods) == 0) {
     signif_periods <- integer(0)
     comps <- matrix(series_org, ncol = 1)
     colnames(comps) <- "noise"

     recon_err <- list(recon_max_abs = 0, recon_rmse = 0)

   } else {

     sig_periods_grp <- split(sig_periods, cumsum(c(1, diff(sig_periods) != 1)))
     signif_periods <- as.integer(unlist(
       lapply(sig_periods_grp, function(x) x[which.max(gws[x])]),
       use.names = FALSE
     ))

     comps_tb <- extract_wavelet_components(
       wave = wave,
       signif_periods = signif_periods,
       scale = scale,
       dj = dj,
       dt = dt,
       series_sd = sdx,
       series_mean = series_mean,
       Cdelta = 0.776,
       w0_0 = pi^(-1/4),
       include_residual = TRUE
     )

     comps_mat <- as.matrix(comps_tb)
     if (ncol(comps_mat) != (length(signif_periods) + 1L)) {
       stop("Internal error: unexpected number of reconstructed component columns.")
     }

     period_labels <- round(period[signif_periods], 2)
     comp_names <- paste0("period_", period_labels)
     colnames(comps_mat) <- c(comp_names, "noise")
     comps <- comps_mat

     recon_vec <- series_org - rowSums(comps)
     recon_err <- list(
       recon_max_abs = max(abs(recon_vec)),
       recon_rmse = sqrt(mean(recon_vec^2))
     )
   }

   comps_names <- if (!is.null(comps)) colnames(comps) else NULL

 } else {
   signif_periods <- integer(0)
 }

 signif_periods <- as.integer(signif_periods)
 has_significance <- length(signif_periods) > 0

 # --- Build return ---
 if (mode == "fast") {
   out <- list(
     gws = gws,
     gws_unmasked = gws_unmasked,
     period = period,
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
     period = period,
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
     comps_names = comps_names,
     gws_n_coi = n_coi,
     neff = neff,
     neff_unmasked = neff_unmasked
   )

   if (diagnostics) {
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
       neff = neff,
       neff_unmasked = neff_unmasked,
       reconstruction_error = recon_err
     )
   }
 }

 out
}
