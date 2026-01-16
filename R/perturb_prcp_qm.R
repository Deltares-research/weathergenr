#' Quantile Mapping for Climate Change Perturbation of Precipitation
#'
#' @description
#' Applies quantile mapping to adjust a daily precipitation time series to reflect
#' changes in monthly mean and variance, consistent with climate change scenarios.
#' The method fits a gamma distribution to observed nonzero precipitation in each
#' calendar month, then modifies the distribution's mean and variance according to
#' supplied monthly/yearly change factors, and maps the original values to their
#' new quantiles in the perturbed distribution.
#'
#' This function is designed for climate stress-testing workflows to impose
#' user-specified changes in mean and variance on daily precipitation data while
#' preserving realistic distributional characteristics.
#'
#' @param prcp Numeric vector. Original daily precipitation values to perturb
#'   (typically for one grid cell). Must be non-negative.
#' @param mean_factor Numeric matrix of mean change factors (multiplicative),
#'   dimension: `n_years` x 12 (year, month). Each entry scales that year-month's mean.
#' @param var_factor Numeric matrix of variance change factors (multiplicative),
#'   dimension: `n_years` x 12 (year, month). Each entry scales that year-month's variance.
#' @param month Integer vector (same length as `prcp`). Calendar month for each
#'   day (1-12).
#' @param year Integer vector (same length as `prcp`). Simulation year index
#'   for each day (1 = first year, etc).
#' @param fit_method Character. Method for fitting the base gamma distribution;
#'   passed to [fitdistrplus::fitdist()]. Default is `"mme"` (method of moments).
#' @param min_events Integer. Minimum number of non-zero precipitation events
#'   required per month for distribution fitting. Default is 10.
#' @param validate_output Logical. If TRUE, checks output for NaN/Inf and replaces
#' @param diagnostics Logical. If TRUE, return diagnostics information in the output.
#' @param seed Integer. Optional random seed used when diagnostics require stochastic elements.
#'   with original values. Default is TRUE.
#' @param verbose Logical. If TRUE, prints diagnostic messages about months that
#'   cannot be perturbed. Default is FALSE.
#'
#' @return
#' Numeric vector, same length as `prcp`. Precipitation time series perturbed
#' according to quantile mapping procedure. Includes attributes:
#' \itemize{
#'   \item `perturbed_months`: Integer vector of months that were successfully perturbed
#'   \item `skipped_months`: Integer vector of months skipped due to insufficient data
#'   \item `n_failed_fits`: Number of distribution fits that failed
#' }
#'
#' @details
#' ## Distribution Fitting
#' The function fits a gamma distribution to non-zero precipitation in each calendar
#' month. Months with fewer than `min_events` non-zero days are skipped. The gamma
#' distribution is parameterized with shape and scale parameters derived from the
#' method of moments or maximum likelihood estimation.
#'
#' ## Quantile Mapping Procedure
#' 1. Fit base gamma distribution to observed non-zero precipitation by month
#' 2. Compute target distribution parameters by scaling base mean and variance
#' 3. Map each original prcp to its quantile in the base distribution
#' 4. Transform to the same quantile in the target distribution
#'
#' ## Limitations
#' - Zero precipitation values remain zero (dry day frequency unchanged)
#' - Requires sufficient non-zero events per month for reliable fitting
#' - Extrapolation beyond observed range may produce unrealistic extreme values
#' - Change factors must be positive; negative or zero values will cause errors
#'
#' @importFrom fitdistrplus fitdist
#' @importFrom stats pgamma qgamma
#'
#' @examples
#' \dontrun{
#' # Example: 2 years of daily data, 5% mean and 10% variance increase
#' set.seed(123)
#' n_days <- 730
#' year <- rep(1:2, each = 365)
#' month_idx <- rep(rep(1:12, times = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)), 2)[1:n_days]
#' daily_precip <- rgamma(n_days, shape = 1, scale = 5)
#'
#' mean_factor <- matrix(1.05, nrow = 2, ncol = 12)
#' var_factor <- matrix(1.10, nrow = 2, ncol = 12)
#'
#' perturbed_precip <- perturb_prcp_qm(
#'   prcp = daily_precip,
#'   mean_factor = mean_factor,
#'   var_factor = var_factor,
#'   month = month_idx,
#'   year = year
#' )
#'
#' # Check which months were perturbed
#' attr(perturbed_precip, "perturbed_months")
#' attr(perturbed_precip, "skipped_months")
#' }
#'
#' @seealso \code{\link[fitdistrplus]{fitdist}}
#' @export
perturb_prcp_qm <- function(
    prcp = NULL,
    mean_factor = NULL,
    var_factor = NULL,
    month = NULL,
    year = NULL,
    fit_method = "mme",
    min_events = 10,
    validate_output = TRUE,
    diagnostics = FALSE,
    seed = NULL,
    verbose = FALSE) {

  # ==========================================================================
  # INPUT VALIDATION
  # ==========================================================================

  # Check for NULL inputs
  if (is.null(prcp)) stop("'prcp' must not be NULL")
  if (is.null(mean_factor)) stop("'mean_factor' must not be NULL")
  if (is.null(var_factor)) stop("'var_factor' must not be NULL")
  if (is.null(month)) stop("'month' must not be NULL")
  if (is.null(year)) stop("'year' must not be NULL")

  # Check types
  if (!is.numeric(prcp)) stop("'prcp' must be numeric")
  if (!is.numeric(mean_factor)) stop("'mean_factor' must be numeric")
  if (!is.numeric(var_factor)) stop("'var_factor' must be numeric")
  if (!is.numeric(month)) stop("'month' must be numeric")
  if (!is.numeric(year)) stop("'year' must be numeric")

  # Check dimensions
  n <- length(prcp)
  if (length(month) != n) {
    stop("'month' must have same length as 'prcp'")
  }
  if (length(year) != n) {
    stop("'year' must have same length as 'prcp'")
  }

  # Check matrix structure
  if (!is.matrix(mean_factor)) {
    stop("'mean_factor' must be a matrix with nrow = n_years, ncol = 12")
  }
  if (!is.matrix(var_factor)) {
    stop("'var_factor' must be a matrix with nrow = n_years, ncol = 12")
  }
  if (ncol(mean_factor) != 12) {
    stop("'mean_factor' must have 12 columns (one per month)")
  }
  if (ncol(var_factor) != 12) {
    stop("'var_factor' must have 12 columns (one per month)")
  }

  # Check prcp ranges
  if (any(prcp < 0, na.rm = TRUE)) {
    stop("'prcp' must be non-negative")
  }
  if (any(month < 1 | month > 12, na.rm = TRUE)) {
    stop("'month' must contain integers between 1 and 12")
  }
  if (any(year < 1, na.rm = TRUE)) {
    stop("'year' must contain positive integers")
  }

  # Validate change factors
  if (any(mean_factor <= 0, na.rm = TRUE)) {
    stop("'mean_factor' must contain positive values")
  }
  if (any(var_factor <= 0, na.rm = TRUE)) {
    stop("'var_factor' must contain positive values")
  }

  # Check matrix dimensions match data
  n_years <- max(year)
  if (nrow(mean_factor) != n_years) {
    stop("'mean_factor' must have ", n_years, " rows (one per year)")
  }
  if (nrow(var_factor) != n_years) {
    stop("'var_factor' must have ", n_years, " rows (one per year)")
  }

  # Store original values if diagnostics requested
  if (diagnostics) {
    prcp_org <- prcp
  }

  # ==========================================================================
  # IDENTIFY NON-ZERO PRECIPITATION EVENTS
  # ==========================================================================

  # Create logical mask for non-zero values
  is_wet <- prcp > 0

  # If no non-zero precipitation, return unchanged
  if (!any(is_wet)) {
    if (verbose) message("No non-zero precipitation values found")
    return(prcp)
  }

  # Extract non-zero values and their months
  prcp_wet <- prcp[is_wet]
  month_wet <- month[is_wet]

  # ==========================================================================
  # FIT BASE DISTRIBUTIONS BY MONTH
  # ==========================================================================

  # Count non-zero events per month
  n_wet_by_month <- table(month_wet)

  # Identify months with sufficient data
  months_ok <- as.integer(names(n_wet_by_month)[n_wet_by_month >= min_events])
  months_skipped <- setdiff(1:12, months_ok)

  if (length(months_ok) == 0) {
    if (verbose) {
      message("No months have sufficient non-zero events (min = ", min_events, ")")
    }
    return(prcp)
  }

  if (verbose && length(months_skipped) > 0) {
    message(
      "Skipping months with insufficient data: ",
      paste(months_skipped, collapse = ", ")
    )
  }

  # Fit distributions for valid months
  base_gamma <- fit_monthly_distributions(
    prcp_wet = prcp_wet,
    month_wet = month_wet,
    months_ok = months_ok,
    fit_method = fit_method,
    verbose = verbose
  )

  # ==========================================================================
  # COMPUTE TARGET DISTRIBUTION PARAMETERS
  # ==========================================================================

  target_gamma <- compute_target_parameters(
    base_gamma = base_gamma,
    mean_factor = mean_factor,
    var_factor = var_factor,
    months_ok = months_ok,
    n_years = n_years
  )

  # ==========================================================================
  # APPLY QUANTILE MAPPING
  # ==========================================================================

  # Identify indices to perturb (non-zero values in valid months)
  idx_perturb <- which(is_wet & (month %in% months_ok))

  if (length(idx_perturb) == 0) {
    if (verbose) message("No values to perturb")
    return(prcp)
  }

  # Apply mapping
  prcp_out <- prcp  # Copy to preserve original

  prcp_out[idx_perturb] <- apply_quantile_mapping(
    prcp = prcp[idx_perturb],
    mon = month[idx_perturb],
    year = year[idx_perturb],
    base_gamma = base_gamma,
    target_gamma = target_gamma
  )

  # ==========================================================================
  # VALIDATION AND OUTPUT
  # ==========================================================================

  # Check for invalid outputs
  if (validate_output) {
    invalid_mask <- !is.finite(prcp_out)
    if (any(invalid_mask)) {
      n_invalid <- sum(invalid_mask)
      if (verbose) {
        warning(
          "Replaced ", n_invalid, " invalid output values (NaN/Inf) with original values"
        )
      }
      prcp_out[invalid_mask] <- prcp[invalid_mask]
    }
  }

  # Add diagnostic attributes
  attr(prcp_out, "perturbed_months") <- months_ok
  attr(prcp_out, "skipped_months") <- months_skipped
  attr(prcp_out, "n_failed_fits") <- attr(base_gamma, "n_failed_fits")

  if (diagnostics) {
    diagnostics <- validate_quantile_mapping(
      prcp_org = prcp_org,
      prcp_adjusted = prcp_out,
      month = month,
      year = year,
      mean_factor = mean_factor,
      var_factor = var_factor
    )

    if (verbose) {
      cat("\n")
      print(diagnostics)
    }

    return(list(
      adjusted = prcp_out,
      diagnostics = diagnostics
    ))
  }

  prcp_out
}



# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Fit gamma distributions for each month
#' @keywords internal
fit_monthly_distributions <- function(prcp_wet, month_wet, months_ok,
                                      fit_method, verbose) {

  n_months <- length(months_ok)
  n_failed <- 0

  # Pre-allocate results
  params <- data.frame(
    month = months_ok,
    shape = numeric(n_months),
    scale = numeric(n_months),
    mean = numeric(n_months),
    var = numeric(n_months),
    stringsAsFactors = FALSE
  )

  # Fit distribution for each month
  for (i in seq_along(months_ok)) {
    m <- months_ok[i]
    prcp_m <- prcp_wet[month_wet == m]

    # Attempt to fit
    fit_result <- tryCatch(
      {
        fit <- fitdistrplus::fitdist(prcp_m, "gamma", method = fit_method)

        # Extract parameters (fitdist returns shape and rate)
        shape <- fit$estimate["shape"]
        rate <- fit$estimate["rate"]
        scale <- 1 / rate

        # Compute moments
        mean_val <- shape * scale
        var_val <- shape * scale^2

        list(
          success = TRUE,
          shape = shape,
          scale = scale,
          mean = mean_val,
          var = var_val
        )
      },
      error = function(e) {
        if (verbose) {
          warning("Failed to fit distribution for month ", m, ": ", e$message)
        }
        list(success = FALSE)
      }
    )

    if (fit_result$success) {
      params$shape[i] <- fit_result$shape
      params$scale[i] <- fit_result$scale
      params$mean[i] <- fit_result$mean
      params$var[i] <- fit_result$var
    } else {
      # Mark as failed - will be excluded from perturbation
      n_failed <- n_failed + 1
      params$shape[i] <- NA
      params$scale[i] <- NA
      params$mean[i] <- NA
      params$var[i] <- NA
    }
  }

  # Remove months with failed fits
  params <- params[!is.na(params$shape), ]

  # Add attribute for diagnostics
  attr(params, "n_failed_fits") <- n_failed

  params
}


#' Compute target distribution parameters from base params and change factors
#' @keywords internal
compute_target_parameters <- function(base_gamma, mean_factor, var_factor,
                                      months_ok, n_years) {

  n_months <- nrow(base_gamma)

  # Pre-allocate arrays (months x years)
  target_mean <- matrix(NA_real_, nrow = n_months, ncol = n_years)
  target_var <- matrix(NA_real_, nrow = n_months, ncol = n_years)
  target_shape <- matrix(NA_real_, nrow = n_months, ncol = n_years)
  target_scale <- matrix(NA_real_, nrow = n_months, ncol = n_years)

  # Vectorized computation across years
  for (i in seq_len(n_months)) {
    m <- base_gamma$month[i]

    # Scale mean and variance by change factors
    target_mean[i, ] <- base_gamma$mean[i] * mean_factor[, m]
    target_var[i, ] <- base_gamma$var[i] * var_factor[, m]

    # Compute gamma parameters from moments
    # For gamma: mean = shape * scale, var = shape * scale^2
    # Therefore: scale = var / mean, shape = mean^2 / var
    target_scale[i, ] <- target_var[i, ] / target_mean[i, ]
    target_shape[i, ] <- target_mean[i, ]^2 / target_var[i, ]
  }

  list(
    months = base_gamma$month,
    shape = target_shape,
    scale = target_scale,
    mean = target_mean,
    var = target_var
  )
}


#' Apply quantile mapping transformation
#' @keywords internal
apply_quantile_mapping <- function(prcp, mon, year, base_gamma, target_gamma) {

  n <- length(prcp)
  result <- numeric(n)

  # Process each unique month-year combination
  # (more efficient than processing each prcp individually)
  month_year_keys <- unique(data.frame(mon = mon, year = year))

  for (i in seq_len(nrow(month_year_keys))) {
    m <- month_year_keys$mon[i]
    y <- month_year_keys$year[i]

    # Find indices for this month-year
    idx <- which(mon == m & year == y)

    # Find parameters for this month
    i_month <- which(base_gamma$month == m)

    if (length(i_month) == 0) next  # Skip if month not in base params

    # Get base distribution parameters
    base_shape <- base_gamma$shape[i_month]
    base_scale <- base_gamma$scale[i_month]

    # Get target distribution parameters
    target_shape <- target_gamma$shape[i_month, y]
    target_scale <- target_gamma$scale[i_month, y]

    # Check for invalid parameters
    if (!is.finite(base_shape) || !is.finite(base_scale) ||
        !is.finite(target_shape) || !is.finite(target_scale) ||
        base_shape <= 0 || base_scale <= 0 ||
        target_shape <= 0 || target_scale <= 0) {
      # Use original values if parameters are invalid
      result[idx] <- prcp[idx]
      next
    }

    # Quantile mapping: original prcp -> quantile -> new prcp
    quantiles <- stats::pgamma(
      prcp[idx],
      shape = base_shape,
      scale = base_scale
    )

    result[idx] <- stats::qgamma(
      quantiles,
      shape = target_shape,
      scale = target_scale
    )
  }



  result
}
