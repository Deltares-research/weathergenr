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
#' @param value Numeric vector. Original daily precipitation values to perturb
#'   (typically for one grid cell). Must be non-negative.
#' @param mean.change Numeric matrix of mean change factors (multiplicative),
#'   dimension: `n_years` x 12 (year, month). Each entry scales that year-month's mean.
#' @param var.change Numeric matrix of variance change factors (multiplicative),
#'   dimension: `n_years` x 12 (year, month). Each entry scales that year-month's variance.
#' @param mon.ts Integer vector (same length as `value`). Calendar month for each
#'   day (1-12).
#' @param year.ts Integer vector (same length as `value`). Simulation year index
#'   for each day (1 = first year, etc).
#' @param fit.method Character. Method for fitting the base gamma distribution;
#'   passed to [fitdistrplus::fitdist()]. Default is `"mme"` (method of moments).
#' @param min.events Integer. Minimum number of non-zero precipitation events
#'   required per month for distribution fitting. Default is 10.
#' @param validate.output Logical. If TRUE, checks output for NaN/Inf and replaces
#'   with original values. Default is TRUE.
#' @param verbose Logical. If TRUE, prints diagnostic messages about months that
#'   cannot be perturbed. Default is FALSE.
#'
#' @return
#' Numeric vector, same length as `value`. Precipitation time series perturbed
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
#' month. Months with fewer than `min.events` non-zero days are skipped. The gamma
#' distribution is parameterized with shape and scale parameters derived from the
#' method of moments or maximum likelihood estimation.
#'
#' ## Quantile Mapping Procedure
#' 1. Fit base gamma distribution to observed non-zero precipitation by month
#' 2. Compute target distribution parameters by scaling base mean and variance
#' 3. Map each original value to its quantile in the base distribution
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
#' year_idx <- rep(1:2, each = 365)
#' month_idx <- rep(rep(1:12, times = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)), 2)[1:n_days]
#' daily_precip <- rgamma(n_days, shape = 1, scale = 5)
#'
#' mean.change <- matrix(1.05, nrow = 2, ncol = 12)
#' var.change <- matrix(1.10, nrow = 2, ncol = 12)
#'
#' perturbed_precip <- quantile_mapping(
#'   value = daily_precip,
#'   mean.change = mean.change,
#'   var.change = var.change,
#'   mon.ts = month_idx,
#'   year.ts = year_idx
#' )
#'
#' # Check which months were perturbed
#' attr(perturbed_precip, "perturbed_months")
#' attr(perturbed_precip, "skipped_months")
#' }
#'
#' @seealso \code{\link[fitdistrplus]{fitdist}}
#' @export
quantile_mapping <- function(
    value = NULL,
    mean.change = NULL,
    var.change = NULL,
    mon.ts = NULL,
    year.ts = NULL,
    fit.method = "mme",
    min.events = 10,
    validate.output = TRUE,
    verbose = FALSE,
    compute.diagnostics = FALSE) {

  # ==========================================================================
  # INPUT VALIDATION
  # ==========================================================================

  # Check for NULL inputs
  if (is.null(value)) stop("'value' must not be NULL")
  if (is.null(mean.change)) stop("'mean.change' must not be NULL")
  if (is.null(var.change)) stop("'var.change' must not be NULL")
  if (is.null(mon.ts)) stop("'mon.ts' must not be NULL")
  if (is.null(year.ts)) stop("'year.ts' must not be NULL")

  # Check types
  if (!is.numeric(value)) stop("'value' must be numeric")
  if (!is.numeric(mean.change)) stop("'mean.change' must be numeric")
  if (!is.numeric(var.change)) stop("'var.change' must be numeric")
  if (!is.numeric(mon.ts)) stop("'mon.ts' must be numeric")
  if (!is.numeric(year.ts)) stop("'year.ts' must be numeric")

  # Check dimensions
  n <- length(value)
  if (length(mon.ts) != n) {
    stop("'mon.ts' must have same length as 'value'")
  }
  if (length(year.ts) != n) {
    stop("'year.ts' must have same length as 'value'")
  }

  # Check matrix structure
  if (!is.matrix(mean.change)) {
    stop("'mean.change' must be a matrix with nrow = n_years, ncol = 12")
  }
  if (!is.matrix(var.change)) {
    stop("'var.change' must be a matrix with nrow = n_years, ncol = 12")
  }
  if (ncol(mean.change) != 12) {
    stop("'mean.change' must have 12 columns (one per month)")
  }
  if (ncol(var.change) != 12) {
    stop("'var.change' must have 12 columns (one per month)")
  }

  # Check value ranges
  if (any(value < 0, na.rm = TRUE)) {
    stop("'value' must be non-negative")
  }
  if (any(mon.ts < 1 | mon.ts > 12, na.rm = TRUE)) {
    stop("'mon.ts' must contain integers between 1 and 12")
  }
  if (any(year.ts < 1, na.rm = TRUE)) {
    stop("'year.ts' must contain positive integers")
  }

  # Validate change factors
  if (any(mean.change <= 0, na.rm = TRUE)) {
    stop("'mean.change' must contain positive values")
  }
  if (any(var.change <= 0, na.rm = TRUE)) {
    stop("'var.change' must contain positive values")
  }

  # Check matrix dimensions match data
  ymax <- max(year.ts)
  if (nrow(mean.change) != ymax) {
    stop("'mean.change' must have ", ymax, " rows (one per year)")
  }
  if (nrow(var.change) != ymax) {
    stop("'var.change' must have ", ymax, " rows (one per year)")
  }

  # Store original values if diagnostics requested
  if (compute.diagnostics) {
    value.original <- value
  }

  # ==========================================================================
  # IDENTIFY NON-ZERO PRECIPITATION EVENTS
  # ==========================================================================

  # Create logical mask for non-zero values
  nonzero.mask <- value > 0

  # If no non-zero precipitation, return unchanged
  if (!any(nonzero.mask)) {
    if (verbose) message("No non-zero precipitation values found")
    return(value)
  }

  # Extract non-zero values and their months
  value.nz <- value[nonzero.mask]
  mon.nz <- mon.ts[nonzero.mask]

  # ==========================================================================
  # FIT BASE DISTRIBUTIONS BY MONTH
  # ==========================================================================

  # Count non-zero events per month
  events.per.month <- table(mon.nz)

  # Identify months with sufficient data
  valid.months <- as.integer(names(events.per.month)[events.per.month >= min.events])
  skipped.months <- setdiff(1:12, valid.months)

  if (length(valid.months) == 0) {
    if (verbose) {
      message("No months have sufficient non-zero events (min = ", min.events, ")")
    }
    return(value)
  }

  if (verbose && length(skipped.months) > 0) {
    message(
      "Skipping months with insufficient data: ",
      paste(skipped.months, collapse = ", ")
    )
  }

  # Fit distributions for valid months
  base.params <- fit_monthly_distributions(
    value.nz = value.nz,
    mon.nz = mon.nz,
    valid.months = valid.months,
    fit.method = fit.method,
    verbose = verbose
  )

  # ==========================================================================
  # COMPUTE TARGET DISTRIBUTION PARAMETERS
  # ==========================================================================

  target.params <- compute_target_parameters(
    base.params = base.params,
    mean.change = mean.change,
    var.change = var.change,
    valid.months = valid.months,
    ymax = ymax
  )

  # ==========================================================================
  # APPLY QUANTILE MAPPING
  # ==========================================================================

  # Identify indices to perturb (non-zero values in valid months)
  perturb.idx <- which(nonzero.mask & (mon.ts %in% valid.months))

  if (length(perturb.idx) == 0) {
    if (verbose) message("No values to perturb")
    return(value)
  }

  # Apply mapping
  value.out <- value  # Copy to preserve original

  value.out[perturb.idx] <- apply_quantile_mapping(
    value = value[perturb.idx],
    mon = mon.ts[perturb.idx],
    year = year.ts[perturb.idx],
    base.params = base.params,
    target.params = target.params
  )

  # ==========================================================================
  # VALIDATION AND OUTPUT
  # ==========================================================================

  # Check for invalid outputs
  if (validate.output) {
    invalid.mask <- !is.finite(value.out)
    if (any(invalid.mask)) {
      n.invalid <- sum(invalid.mask)
      if (verbose) {
        warning(
          "Replaced ", n.invalid, " invalid output values (NaN/Inf) with original values"
        )
      }
      value.out[invalid.mask] <- value[invalid.mask]
    }
  }

  # Add diagnostic attributes
  attr(value.out, "perturbed_months") <- valid.months
  attr(value.out, "skipped_months") <- skipped.months
  attr(value.out, "n_failed_fits") <- attr(base.params, "n_failed_fits")

  if (compute.diagnostics) {
    diagnostics <- validate_quantile_mapping(
      value.original = value.original,
      value.adjusted = value.out,
      mon.ts = mon.ts,
      year.ts = year.ts,
      mean.change = mean.change,
      var.change = var.change
    )

    if (verbose) {
      cat("\n")
      print(diagnostics)
    }

    return(list(
      adjusted = value.out,
      diagnostics = diagnostics
    ))
  }

  value.out
}



# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Fit gamma distributions for each month
#' @keywords internal
fit_monthly_distributions <- function(value.nz, mon.nz, valid.months,
                                      fit.method, verbose) {

  n.months <- length(valid.months)
  n.failed <- 0

  # Pre-allocate results
  params <- data.frame(
    month = valid.months,
    shape = numeric(n.months),
    scale = numeric(n.months),
    mean = numeric(n.months),
    var = numeric(n.months),
    stringsAsFactors = FALSE
  )

  # Fit distribution for each month
  for (i in seq_along(valid.months)) {
    m <- valid.months[i]
    values.m <- value.nz[mon.nz == m]

    # Attempt to fit
    fit.result <- tryCatch(
      {
        fit <- fitdistrplus::fitdist(values.m, "gamma", method = fit.method)

        # Extract parameters (fitdist returns shape and rate)
        shape <- fit$estimate["shape"]
        rate <- fit$estimate["rate"]
        scale <- 1 / rate

        # Compute moments
        mean.val <- shape * scale
        var.val <- shape * scale^2

        list(
          success = TRUE,
          shape = shape,
          scale = scale,
          mean = mean.val,
          var = var.val
        )
      },
      error = function(e) {
        if (verbose) {
          warning("Failed to fit distribution for month ", m, ": ", e$message)
        }
        list(success = FALSE)
      }
    )

    if (fit.result$success) {
      params$shape[i] <- fit.result$shape
      params$scale[i] <- fit.result$scale
      params$mean[i] <- fit.result$mean
      params$var[i] <- fit.result$var
    } else {
      # Mark as failed - will be excluded from perturbation
      n.failed <- n.failed + 1
      params$shape[i] <- NA
      params$scale[i] <- NA
      params$mean[i] <- NA
      params$var[i] <- NA
    }
  }

  # Remove months with failed fits
  params <- params[!is.na(params$shape), ]

  # Add attribute for diagnostics
  attr(params, "n_failed_fits") <- n.failed

  params
}


#' Compute target distribution parameters from base params and change factors
#' @keywords internal
compute_target_parameters <- function(base.params, mean.change, var.change,
                                      valid.months, ymax) {

  n.months <- nrow(base.params)

  # Pre-allocate arrays (months x years)
  target.mean <- matrix(NA_real_, nrow = n.months, ncol = ymax)
  target.var <- matrix(NA_real_, nrow = n.months, ncol = ymax)
  target.shape <- matrix(NA_real_, nrow = n.months, ncol = ymax)
  target.scale <- matrix(NA_real_, nrow = n.months, ncol = ymax)

  # Vectorized computation across years
  for (i in seq_len(n.months)) {
    m <- base.params$month[i]

    # Scale mean and variance by change factors
    target.mean[i, ] <- base.params$mean[i] * mean.change[, m]
    target.var[i, ] <- base.params$var[i] * var.change[, m]

    # Compute gamma parameters from moments
    # For gamma: mean = shape * scale, var = shape * scale^2
    # Therefore: scale = var / mean, shape = mean^2 / var
    target.scale[i, ] <- target.var[i, ] / target.mean[i, ]
    target.shape[i, ] <- target.mean[i, ]^2 / target.var[i, ]
  }

  list(
    months = base.params$month,
    shape = target.shape,
    scale = target.scale,
    mean = target.mean,
    var = target.var
  )
}


#' Apply quantile mapping transformation
#' @keywords internal
apply_quantile_mapping <- function(value, mon, year, base.params, target.params) {

  n <- length(value)
  result <- numeric(n)

  # Process each unique month-year combination
  # (more efficient than processing each value individually)
  month.year.combos <- unique(data.frame(mon = mon, year = year))

  for (i in seq_len(nrow(month.year.combos))) {
    m <- month.year.combos$mon[i]
    y <- month.year.combos$year[i]

    # Find indices for this month-year
    idx <- which(mon == m & year == y)

    # Find parameters for this month
    m.idx <- which(base.params$month == m)

    if (length(m.idx) == 0) next  # Skip if month not in base params

    # Get base distribution parameters
    base.shape <- base.params$shape[m.idx]
    base.scale <- base.params$scale[m.idx]

    # Get target distribution parameters
    target.shape <- target.params$shape[m.idx, y]
    target.scale <- target.params$scale[m.idx, y]

    # Check for invalid parameters
    if (!is.finite(base.shape) || !is.finite(base.scale) ||
        !is.finite(target.shape) || !is.finite(target.scale) ||
        base.shape <= 0 || base.scale <= 0 ||
        target.shape <= 0 || target.scale <= 0) {
      # Use original values if parameters are invalid
      result[idx] <- value[idx]
      next
    }

    # Quantile mapping: original value -> quantile -> new value
    quantiles <- stats::pgamma(
      value[idx],
      shape = base.shape,
      scale = base.scale
    )

    result[idx] <- stats::qgamma(
      quantiles,
      shape = target.shape,
      scale = target.scale
    )
  }



  result
}
