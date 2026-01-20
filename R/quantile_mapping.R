#' Adjust Daily Precipitation with Gamma Quantile Mapping Under Monthly Scenario Factors
#'
#' @description
#' Applies a month-wise **Gamma quantile mapping (QM)** to daily precipitation to impose
#' prescribed changes in **wet-day intensity** (monthly wet-day mean and variance).
#'
#' The workflow is:
#' \enumerate{
#'   \item Identify "wet-day intensities" as \code{precip > intensity_threshold}.
#'   \item For each calendar month, fit a **baseline Gamma distribution** to wet-day
#'     intensities using \code{fitdistrplus::fitdist()}.
#'   \item For each \strong{simulation year index} \code{y} and month \code{m}, construct a
#'     **target Gamma distribution** by scaling baseline moments using
#'     \code{mean_factor[y, m]} and \code{var_factor[y, m]} (and optionally
#'     \code{scale_var_with_mean}).
#'   \item Map each wet-day intensity through the baseline CDF (\code{pgamma}) and then through
#'     the target inverse CDF (\code{qgamma}).
#' }
#'
#' Values \code{<= intensity_threshold} are treated as "dry" for this function and are
#' returned unchanged. Therefore, this function does \strong{not} change wet-day frequency.
#' If you need changes in wet/dry spell structure, apply a separate occurrence perturbation.
#'
#' @param precip Numeric vector. Daily precipitation (mm/day). Must be \code{>= 0}.
#'   \code{NA} values are allowed and pass through unchanged.
#'
#' @param mean_factor Numeric matrix with dimensions \code{(n_years x 12)}.
#'   Monthly multiplicative factors applied to the \strong{baseline wet-day mean}
#'   for each month. Entry \code{mean_factor[y, m]} is used for simulation year index
#'   \code{y} and calendar month \code{m}.
#'
#' @param var_factor Numeric matrix with dimensions \code{(n_years x 12)}.
#'   Monthly multiplicative factors applied to the \strong{baseline wet-day variance}
#'   for each month. Entry \code{var_factor[y, m]} is used for simulation year index
#'   \code{y} and calendar month \code{m}.
#'
#' @param scale_var_with_mean Logical scalar. If \code{TRUE}, the effective variance factor
#'   used to build target distributions is:
#'   \deqn{var\_factor\_use = var\_factor \times mean\_factor^2}
#'   This tends to keep the coefficient of variation (CV) more stable when mean changes
#'   are applied. Default is \code{FALSE}.
#'
#' @param exaggerate_extremes Logical scalar. If \code{TRUE}, applies a tail-exponent
#'   transform in probability space above \code{extreme_prob_threshold} before mapping
#'   into the target distribution. This amplifies relative changes in the upper tail while
#'   keeping the Gamma-to-Gamma mapping structure. Default is \code{FALSE}.
#'
#' @param extreme_prob_threshold Numeric scalar in \code{(0, 1)}. Defines the start of the
#'   "extreme" tail in baseline probability space. For example, \code{0.95} corresponds to
#'   the top 5\% of wet-day intensities under the baseline month-specific Gamma.
#'   Used only when \code{exaggerate_extremes = TRUE}. Default is \code{0.95}.
#'
#' @param extreme_k Numeric scalar \code{> 0}. Tail exponent controlling tail amplification
#'   when \code{exaggerate_extremes = TRUE}. Values \code{> 1} amplify extremes; values in
#'   \code{(0, 1)} dampen them. Default is \code{1.2}.
#'
#' @param enforce_target_mean Logical scalar. If \code{TRUE}, rescales mapped wet-day values
#'   within each \code{(year index, month)} group so that the wet-day mean matches the
#'   intended target mean \code{baseline_mean(month) * mean_factor[y, m]}.
#'   This is most relevant when \code{exaggerate_extremes = TRUE} because tail amplification
#'   can shift the mean away from the intended target. Default is \code{TRUE}.
#'
#' @param month Integer vector, same length as \code{precip}. Calendar month for each day
#'   (\code{1}--\code{12}).
#'
#' @param year Integer vector, same length as \code{precip}. \strong{Simulation year index}
#'   for each day (\code{1 =} first simulated year, \code{2 =} second, ...).
#'   Must be a contiguous index set \code{1:n_years}. Do not pass calendar years.
#'   \code{max(year)} must equal \code{nrow(mean_factor)} and \code{nrow(var_factor)}.
#'
#' @param intensity_threshold Numeric scalar \code{>= 0}. Defines which values are treated as
#'   wet-day intensities for fitting and mapping: \code{precip > intensity_threshold}.
#'   Values \code{<= intensity_threshold} are returned unchanged. Default is \code{0}.
#'
#' @param fit_method Character scalar. Estimation method passed to
#'   \code{fitdistrplus::fitdist()} for Gamma fitting (e.g., \code{"mme"}, \code{"mle"}).
#'   Default is \code{"mme"}.
#'
#' @param min_events Integer scalar. Minimum number of wet-day intensities required within a
#'   month to fit the baseline Gamma. Months with fewer wet-day values are skipped; wet-day
#'   values in skipped months pass through unchanged. Default is \code{10}.
#'
#' @param validate_output Logical scalar. If \code{TRUE}, replaces any \code{Inf}/\code{NaN}
#'   produced by the mapping with the original values at those positions and clamps any
#'   negative outputs to zero. Default is \code{TRUE}.
#'
#' @param diagnostics Logical scalar. If \code{TRUE}, returns a list containing:
#'   the adjusted series, diagnostics from \code{validate_quantile_mapping()}, and the fitted
#'   objects needed by higher-level workflows (baseline and target Gamma parameters).
#'   Default is \code{FALSE}.
#'
#' @param seed Optional integer. If provided, sets a temporary RNG seed via \code{set.seed()}.
#'   The previous RNG state is restored on exit. Default is \code{NULL}.
#'
#' @param verbose Logical scalar. If \code{TRUE}, prints progress and warnings for skipped
#'   months and failed fits. Default is \code{FALSE}.
#'
#' @return
#' If \code{diagnostics = FALSE}, returns a numeric vector of the same length as \code{precip}.
#' The returned vector has attributes:
#' \itemize{
#'   \item \code{perturbed_months}: months (1--12) successfully fitted and perturbed
#'   \item \code{skipped_months}: months (1--12) skipped due to insufficient data or failed fits
#'   \item \code{n_failed_fits}: number of monthly Gamma fits that failed
#' }
#'
#' If \code{diagnostics = TRUE}, returns a list:
#' \itemize{
#'   \item \code{adjusted}: numeric vector as above (with the same attributes)
#'   \item \code{diagnostics}: output of \code{validate_quantile_mapping()}
#' }
#'
#' @details
#' \strong{Interpretation of mean/variance factors}
#' \itemize{
#'   \item Factors apply to wet-day intensities only (\code{precip > intensity_threshold}).
#'   \item They do not control wet-day frequency; apply a separate occurrence model if required.
#' }
#'
#' @seealso
#' \code{\link[fitdistrplus]{fitdist}}, \code{\link[stats]{pgamma}}, \code{\link[stats]{qgamma}},
#' \code{\link{validate_quantile_mapping}}, \code{\link{diagnose_precip_qm}}
#'
#' @export
adjust_precipitation_qm <- function(
    precip = NULL,
    mean_factor = NULL,
    var_factor = NULL,
    scale_var_with_mean = TRUE,
    exaggerate_extremes = FALSE,
    extreme_prob_threshold = 0.95,
    extreme_k = 1.2,
    enforce_target_mean = TRUE,
    month = NULL,
    year = NULL,
    intensity_threshold = 0,
    fit_method = "mme",
    min_events = 10,
    validate_output = TRUE,
    diagnostics = FALSE,
    seed = NULL,
    verbose = FALSE) {

  # ==========================================================================
  # INPUT VALIDATION
  # ==========================================================================

  if (is.null(precip)) stop("'precip' must not be NULL", call. = FALSE)
  if (is.null(mean_factor)) stop("'mean_factor' must not be NULL", call. = FALSE)
  if (is.null(var_factor)) stop("'var_factor' must not be NULL", call. = FALSE)
  if (is.null(month)) stop("'month' must not be NULL", call. = FALSE)
  if (is.null(year)) stop("'year' must not be NULL", call. = FALSE)

  if (!is.numeric(precip)) stop("'precip' must be numeric", call. = FALSE)
  if (!is.numeric(month)) stop("'month' must be numeric", call. = FALSE)
  if (!is.numeric(year)) stop("'year' must be numeric", call. = FALSE)

  if (!is.numeric(intensity_threshold) || length(intensity_threshold) != 1L ||
      !is.finite(intensity_threshold) || intensity_threshold < 0) {
    stop("'intensity_threshold' must be a single finite numeric >= 0.", call. = FALSE)
  }

  # Allow NA, but disallow NaN/Inf
  bad_inf <- is.infinite(precip) | is.nan(precip)
  if (any(bad_inf, na.rm = TRUE)) {
    stop("'precip' must not contain Inf/NaN (use NA for missing).", call. = FALSE)
  }
  if (any(precip < 0, na.rm = TRUE)) stop("'precip' must be non-negative", call. = FALSE)

  n <- length(precip)
  if (length(month) != n) stop("'month' must have same length as 'precip'", call. = FALSE)
  if (length(year) != n) stop("'year' must have same length as 'precip'", call. = FALSE)

  if (any(!is.finite(month), na.rm = TRUE) || any(!is.finite(year), na.rm = TRUE)) {
    stop("'month' and 'year' must contain only finite values.", call. = FALSE)
  }

  if (any(month < 1 | month > 12, na.rm = TRUE)) stop("'month' must be in 1..12", call. = FALSE)
  if (any(year < 1, na.rm = TRUE)) stop("'year' must contain positive integers", call. = FALSE)

  # Enforce integer year index and contiguity 1..n_years
  year <- as.integer(year)
  uy <- sort(unique(year))
  if (any(is.na(uy))) stop("'year' must not contain NA.", call. = FALSE)
  if (!identical(uy, seq_len(max(uy)))) {
    stop(
      "'year' must be a contiguous simulation-year index 1..n_years. Do not pass calendar years.",
      call. = FALSE
    )
  }

  # Logical scalars
  .is_scalar_lgl <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (!.is_scalar_lgl(scale_var_with_mean)) stop("'scale_var_with_mean' must be TRUE/FALSE", call. = FALSE)
  if (!.is_scalar_lgl(exaggerate_extremes)) stop("'exaggerate_extremes' must be TRUE/FALSE", call. = FALSE)
  if (!.is_scalar_lgl(enforce_target_mean)) stop("'enforce_target_mean' must be TRUE/FALSE", call. = FALSE)
  if (!.is_scalar_lgl(validate_output)) stop("'validate_output' must be TRUE/FALSE", call. = FALSE)
  if (!.is_scalar_lgl(diagnostics)) stop("'diagnostics' must be TRUE/FALSE", call. = FALSE)
  if (!.is_scalar_lgl(verbose)) stop("'verbose' must be TRUE/FALSE", call. = FALSE)

  # min_events
  if (!is.numeric(min_events) || length(min_events) != 1L || !is.finite(min_events) || min_events < 1) {
    stop("'min_events' must be a single positive integer.", call. = FALSE)
  }
  min_events <- as.integer(min_events)

  # fit_method
  if (!is.character(fit_method) || length(fit_method) != 1L || is.na(fit_method) || nchar(fit_method) == 0L) {
    stop("'fit_method' must be a single non-empty character value.", call. = FALSE)
  }

  # mean_factor matrix checks
  if (!is.matrix(mean_factor)) stop("'mean_factor' must be a matrix (n_years x 12)", call. = FALSE)
  if (ncol(mean_factor) != 12L) stop("'mean_factor' must have 12 columns (one per month)", call. = FALSE)
  if (any(!is.finite(mean_factor), na.rm = TRUE)) stop("'mean_factor' must contain only finite values", call. = FALSE)
  if (any(mean_factor <= 0, na.rm = TRUE)) stop("'mean_factor' must contain positive values", call. = FALSE)

  # Determine n_years from year indices
  n_years <- max(year)
  if (nrow(mean_factor) != n_years) {
    stop("'mean_factor' must have ", n_years, " rows (one per simulation year index)", call. = FALSE)
  }

  # Tail amplification checks
  if (!is.numeric(extreme_prob_threshold) || length(extreme_prob_threshold) != 1L ||
      !is.finite(extreme_prob_threshold) || extreme_prob_threshold <= 0 || extreme_prob_threshold >= 1) {
    stop("'extreme_prob_threshold' must be a single numeric in (0, 1)", call. = FALSE)
  }
  if (!is.numeric(extreme_k) || length(extreme_k) != 1L || !is.finite(extreme_k) || extreme_k <= 0) {
    stop("'extreme_k' must be a single positive numeric", call. = FALSE)
  }

  # var_factor matrix checks
  if (!is.matrix(var_factor)) stop("'var_factor' must be a matrix (n_years x 12)", call. = FALSE)
  if (ncol(var_factor) != 12L) stop("'var_factor' must have 12 columns (one per month)", call. = FALSE)
  if (nrow(var_factor) != n_years) stop("'var_factor' must have ", n_years, " rows (one per year)", call. = FALSE)
  if (any(!is.finite(var_factor), na.rm = TRUE)) stop("'var_factor' must contain only finite values", call. = FALSE)
  if (any(var_factor <= 0, na.rm = TRUE)) stop("'var_factor' must contain positive values", call. = FALSE)

  # --------------------------------------------------------------------------
  # Variance factor handling (COMBINE, do not ignore)
  # --------------------------------------------------------------------------
  if (isTRUE(scale_var_with_mean)) {
    var_factor_use <- var_factor * (mean_factor^2)
    if (isTRUE(verbose)) {
      message("scale_var_with_mean=TRUE: using var_factor_use = var_factor * mean_factor^2.")
    }
  } else {
    var_factor_use <- var_factor
  }

  # Optional: keep copy for diagnostics
  if (isTRUE(diagnostics)) precip_org <- precip

  # RNG hygiene: set + restore
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be a single finite numeric/integer if provided", call. = FALSE)
    }
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  # ==========================================================================
  # IDENTIFY WET EVENTS
  # ==========================================================================

  is_wet <- precip > intensity_threshold
  is_wet[is.na(is_wet)] <- FALSE

  if (!any(is_wet)) {
    if (isTRUE(verbose)) message("No wet-day precipitation values found.")
    return(precip)
  }

  precip_wet  <- precip[is_wet]
  month_wet <- month[is_wet]

  # ==========================================================================
  # FIT BASE DISTRIBUTIONS BY MONTH
  # ==========================================================================

  n_wet_by_month <- table(month_wet)
  months_ok <- as.integer(names(n_wet_by_month)[n_wet_by_month >= min_events])
  months_skipped <- setdiff(1:12, months_ok)

  if (length(months_ok) == 0L) {
    if (isTRUE(verbose)) {
      message("No months have at least ", format(min_events, big.mark = ","), " wet-day events.")
    }
    return(precip)
  }

  if (isTRUE(verbose) && length(months_skipped) > 0L) {
    message("Skipping months with insufficient data: ", paste(months_skipped, collapse = ", "), ".")
  }

  base_gamma <- fit_monthly_distributions(
    precip_wet = precip_wet,
    month_wet = month_wet,
    months_ok = months_ok,
    fit_method = fit_method,
    verbose = verbose
  )

  months_fit <- sort(unique(base_gamma$month))
  months_failed_fit <- setdiff(months_ok, months_fit)
  if (length(months_failed_fit) > 0L) {
    months_skipped <- sort(unique(c(months_skipped, months_failed_fit)))
  }
  months_ok <- months_fit

  if (nrow(base_gamma) == 0L || length(months_ok) == 0L) {
    if (isTRUE(verbose)) message("All monthly distribution fits failed; returning original precipitation.")
    return(precip)
  }

  # ==========================================================================
  # COMPUTE TARGET DISTRIBUTION PARAMETERS
  # ==========================================================================

  target_gamma <- compute_target_parameters(
    base_gamma = base_gamma,
    mean_factor = mean_factor,
    var_factor = var_factor_use,
    months_ok = months_ok,
    n_years = n_years
  )

  # ==========================================================================
  # APPLY QUANTILE MAPPING (WET DAYS IN FITTED MONTHS ONLY)
  # ==========================================================================

  idx_perturb <- which(is_wet & (month %in% months_ok))
  if (length(idx_perturb) == 0L) {
    if (isTRUE(verbose)) message("No values to perturb.")
    return(precip)
  }

  precip_out <- precip

  # Ensure NA pass-through
  idx_perturb <- idx_perturb[!is.na(precip_out[idx_perturb])]
  if (length(idx_perturb) == 0L) {
    if (isTRUE(verbose)) message("No non-missing values to perturb.")
    return(precip)
  }

  precip_out[idx_perturb] <- apply_quantile_mapping(
    precip = precip[idx_perturb],
    mon = month[idx_perturb],
    year = year[idx_perturb],
    base_gamma = base_gamma,
    target_gamma = target_gamma,
    exaggerate_extremes = exaggerate_extremes,
    extreme_prob_threshold = extreme_prob_threshold,
    extreme_k = extreme_k,
    enforce_target_mean = enforce_target_mean
  )

  # ==========================================================================
  # VALIDATION AND OUTPUT
  # ==========================================================================

  if (isTRUE(validate_output)) {
    invalid_mask <- !is.finite(precip_out)
    invalid_mask[is.na(invalid_mask)] <- FALSE

    if (any(invalid_mask)) {
      n_invalid <- sum(invalid_mask)
      if (isTRUE(verbose)) {
        warning(
          "Replaced ", format(n_invalid, big.mark = ","),
          " invalid output values (NaN/Inf) with original values.",
          call. = FALSE
        )
      }
      precip_out[invalid_mask] <- precip[invalid_mask]
    }

    neg_mask <- precip_out < 0
    neg_mask[is.na(neg_mask)] <- FALSE
    if (any(neg_mask)) {
      if (isTRUE(verbose)) warning("Negative values produced; clamping to zero.", call. = FALSE)
      precip_out[neg_mask] <- 0
    }
  }

  attr(precip_out, "perturbed_months") <- months_ok
  attr(precip_out, "skipped_months") <- months_skipped
  attr(precip_out, "n_failed_fits") <- attr(base_gamma, "n_failed_fits")

  if (isTRUE(diagnostics)) {

    d <- validate_quantile_mapping(
      precip_org = precip_org,
      precip_adjusted = precip_out,
      month = month,
      year = year,
      mean_factor = mean_factor,
      var_factor = var_factor_use
    )

    if (isTRUE(verbose)) {
      cat("\n")
      print(d)
    }

    return(list(
      adjusted = precip_out,
      diagnostics = d
    ))
  }

  precip_out
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Fit Monthly Gamma Distributions for Wet-Day Precipitation
#'
#' @description
#' Fits a gamma distribution to wet-day precipitation amounts for each month in
#' `months_ok`. Uses `fitdistrplus::fitdist()` and returns month-specific gamma
#' parameters and implied moments. Months that fail to fit are dropped.
#'
#' @param precip_wet Numeric vector of wet-day precipitation values (strictly > 0
#'   recommended). Length must match `month_wet`.
#' @param month_wet Integer vector of months (1-12) corresponding to each element
#'   of `precip_wet`.
#' @param months_ok Integer vector of months to attempt fitting (subset of 1-12).
#' @param fit_method Character scalar. Estimation method forwarded to
#'   `fitdistrplus::fitdist()` (e.g., `"mle"`, `"mme"`).
#' @param verbose Logical. If `TRUE`, emits a warning when a monthly fit fails.
#'
#' @return
#' A `data.frame` with columns:
#' \itemize{
#'   \item `month` Month number.
#'   \item `shape` Gamma shape parameter.
#'   \item `scale` Gamma scale parameter.
#'   \item `mean` Implied mean (`shape * scale`).
#'   \item `var` Implied variance (`shape * scale^2`).
#' }
#' Months with failed fits are removed. The returned data frame has an attribute
#' `n_failed_fits` with the number of failed months.
#'
#' @keywords internal
fit_monthly_distributions <- function(precip_wet,
                                      month_wet,
                                      months_ok,
                                      fit_method,
                                      verbose) {
  n_months <- length(months_ok)
  n_failed <- 0L

  params <- data.frame(
    month = months_ok,
    shape = rep(NA_real_, n_months),
    scale = rep(NA_real_, n_months),
    mean  = rep(NA_real_, n_months),
    var   = rep(NA_real_, n_months),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(months_ok)) {

    m <- months_ok[i]
    precip_m <- precip_wet[month_wet == m]
    precip_m <- precip_m[is.finite(precip_m)]

    # Robustness: MME gamma fit degenerates when variance is ~0
    if (length(precip_m) < 2L || stats::sd(precip_m, na.rm = TRUE) <= sqrt(.Machine$double.eps)) {
      if (isTRUE(verbose)) {
        warning("Month ", m, ": wet-day intensities have ~zero variance; skipping Gamma fit.", call. = FALSE)
      }
      n_failed <- n_failed + 1L
      next
    }

    # Robustness: extremely low discreteness can also destabilize fitting
    if (length(unique(precip_m)) < 3L) {
      if (isTRUE(verbose)) {
        warning("Month ", m, ": too few unique wet-day values; skipping Gamma fit.", call. = FALSE)
      }
      n_failed <- n_failed + 1L
      next
    }


    fit_result <- tryCatch(
      {
        fit <- fitdistrplus::fitdist(precip_m, "gamma", method = fit_method)

        # fitdist() returns shape and rate; convert to scale
        shape <- unname(fit$estimate["shape"])
        rate  <- unname(fit$estimate["rate"])
        scale <- 1 / rate

        mean_val <- shape * scale
        var_val  <- shape * scale^2

        list(
          success = TRUE,
          shape   = shape,
          scale   = scale,
          mean    = mean_val,
          var     = var_val
        )
      },
      error = function(e) {
        if (isTRUE(verbose)) {
          warning(
            "Failed to fit distribution for month ", m, ": ", e$message,
            call. = FALSE
          )
        }
        list(success = FALSE)
      }
    )

    if (isTRUE(fit_result$success)) {
      params$shape[i] <- fit_result$shape
      params$scale[i] <- fit_result$scale
      params$mean[i]  <- fit_result$mean
      params$var[i]   <- fit_result$var
    } else {
      n_failed <- n_failed + 1L
    }
  }

  params <- params[!is.na(params$shape), , drop = FALSE]
  attr(params, "n_failed_fits") <- n_failed

  params
}


#' Compute Target Gamma Parameters from Baseline and Change Factors
#'
#' @description
#' Computes year-specific target gamma distribution parameters for each month by
#' scaling baseline monthly mean and variance using change factors, then
#' converting target moments to gamma `shape` and `scale`.
#'
#' @param base_gamma Data frame produced by `fit_monthly_distributions()`, with
#'   columns `month`, `mean`, and `var` (and typically `shape`, `scale`).
#' @param mean_factor Numeric matrix of multiplicative mean factors with
#'   dimensions `n_years x 12` (year x month).
#' @param var_factor Numeric matrix of multiplicative variance factors with
#'   dimensions `n_years x 12` (year x month).
#' @param months_ok Integer vector of months expected to be represented (kept for
#'   interface consistency; the function uses `base_gamma$month`).
#' @param n_years Integer. Number of years (must match `nrow(mean_factor)` and
#'   `nrow(var_factor)`).
#'
#' @return
#' A list with elements:
#' \itemize{
#'   \item `months` Integer vector of months (from `base_gamma$month`).
#'   \item `shape` Numeric matrix `[n_months x n_years]` of target shapes.
#'   \item `scale` Numeric matrix `[n_months x n_years]` of target scales.
#'   \item `mean` Numeric matrix `[n_months x n_years]` of target means.
#'   \item `var` Numeric matrix `[n_months x n_years]` of target variances.
#' }
#'
#' @keywords internal
compute_target_parameters <- function(base_gamma,
                                      mean_factor,
                                      var_factor,
                                      months_ok,
                                      n_years) {
  n_months <- nrow(base_gamma)

  target_mean  <- matrix(NA_real_, nrow = n_months, ncol = n_years)
  target_var   <- matrix(NA_real_, nrow = n_months, ncol = n_years)
  target_shape <- matrix(NA_real_, nrow = n_months, ncol = n_years)
  target_scale <- matrix(NA_real_, nrow = n_months, ncol = n_years)

  for (i in seq_len(n_months)) {
    m <- base_gamma$month[i]

    target_mean[i, ] <- base_gamma$mean[i] * mean_factor[, m]
    target_var[i, ]  <- base_gamma$var[i]  * var_factor[, m]

    # Gamma moments: mean = shape * scale, var = shape * scale^2
    # => scale = var / mean, shape = mean^2 / var
    target_scale[i, ] <- target_var[i, ] / target_mean[i, ]
    target_shape[i, ] <- (target_mean[i, ]^2) / target_var[i, ]
  }

  list(
    months = base_gamma$month,
    shape  = target_shape,
    scale  = target_scale,
    mean   = target_mean,
    var    = target_var
  )
}


#' Apply Gamma-to-Gamma Quantile Mapping for Precipitation
#'
#' @description
#' Internal helper implementing gamma-based quantile mapping for precipitation:
#' values are mapped to probability space under a baseline monthly gamma, then
#' mapped back to precipitation under a target monthly/year gamma. Optionally
#' applies a tail-probability exaggeration and an enforcement step to match the
#' target mean at the month-year subset level.
#'
#' @param precip Numeric vector of precipitation values (subset, typically wet days).
#' @param mon Integer vector of months (1-12), same length as `precip`.
#' @param year Integer vector of year indices used to select columns in
#'   `target_gamma$shape/scale/mean`. This is an index (1..n_years), not
#'   necessarily the calendar year value.
#' @param base_gamma Data frame with `month`, `shape`, and `scale` columns.
#' @param target_gamma List produced by `compute_target_parameters()`, containing
#'   `shape`, `scale`, and `mean` matrices indexed as `[month_row, year_index]`.
#' @param exaggerate_extremes Logical. If `TRUE`, modifies probabilities above
#'   `extreme_prob_threshold` as `u' = 1 - (1-u)^k`.
#' @param extreme_prob_threshold Numeric in (0,1). Tail threshold `u0`.
#' @param extreme_k Positive numeric. Tail exponent `k`.
#' @param enforce_target_mean Logical. If `TRUE`, rescales mapped values within
#'   each (month, year) subset to match `target_gamma$mean`.
#'
#' @return Numeric vector of mapped precipitation values, same length as `precip`.
#'
#' @keywords internal
apply_quantile_mapping <- function(precip,
                                   mon,
                                   year,
                                   base_gamma,
                                   target_gamma,
                                   exaggerate_extremes = FALSE,
                                   extreme_prob_threshold = 0.95,
                                   extreme_k = 1.2,
                                   enforce_target_mean = TRUE) {
  n <- length(precip)
  result <- precip

  if (n == 0L) {
    return(result)
  }

  keys <- unique(data.frame(mon = mon, year = year, stringsAsFactors = FALSE))
  eps <- 1e-12

  for (i in seq_len(nrow(keys))) {
    m <- keys$mon[i]
    y <- keys$year[i]

    idx <- which(mon == m & year == y)
    if (length(idx) == 0L) {
      next
    }

    i_month <- which(base_gamma$month == m)
    if (length(i_month) == 0L) {
      result[idx] <- precip[idx]
      next
    }

    base_shape <- base_gamma$shape[i_month]
    base_scale <- base_gamma$scale[i_month]

    target_shape <- target_gamma$shape[i_month, y]
    target_scale <- target_gamma$scale[i_month, y]
    target_mean  <- target_gamma$mean[i_month, y]

    if (!is.finite(base_shape) || !is.finite(base_scale) ||
        !is.finite(target_shape) || !is.finite(target_scale) ||
        base_shape <= 0 || base_scale <= 0 ||
        target_shape <= 0 || target_scale <= 0) {
      result[idx] <- precip[idx]
      next
    }

    u <- stats::pgamma(precip[idx], shape = base_shape, scale = base_scale)
    u <- pmin(pmax(u, eps), 1 - eps)

    if (isTRUE(exaggerate_extremes)) {
      u0 <- extreme_prob_threshold
      k  <- extreme_k

      if (is.finite(u0) && u0 > 0 && u0 < 1 && is.finite(k) && k > 0) {
        above <- u > u0
        if (any(above)) {
          ut <- u[above]
          u_new <- 1 - (1 - ut)^k
          u[above] <- pmin(pmax(u_new, eps), 1 - eps)
        }
      }
    }

    mapped <- stats::qgamma(u, shape = target_shape, scale = target_scale)

    if (isTRUE(enforce_target_mean) &&
        is.finite(target_mean) && target_mean > 0) {
      mbar <- mean(mapped, na.rm = TRUE)
      if (is.finite(mbar) && mbar > 0) {
        mapped <- mapped * (target_mean / mbar)
      }
    }

    bad <- !is.finite(mapped)
    if (any(bad)) {
      mapped[bad] <- precip[idx][bad]
    }
    mapped[mapped < 0] <- 0

    result[idx] <- mapped
  }

  result
}
