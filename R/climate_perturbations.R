# ==============================================================================
# CLIMATE PERTURBATIONS
# ==============================================================================

#' Apply Monthly Climate Perturbations to Gridded Daily Weather Series
#'
#' @description
#' Applies **monthly climate perturbations** to **daily gridded weather simulations**
#' in a modular and physically interpretable way. The function operates independently
#' on precipitation intensity, precipitation occurrence, temperature, and (optionally)
#' potential evapotranspiration (PET).
#'
#' The perturbation workflow is:
#' \enumerate{
#'   \item Construct a **simulation-year index** from \code{date} (\code{1..n_years}).
#'   \item Apply **precipitation occurrence perturbations** (optional) by modifying
#'     monthly wet-day probabilities within each (year index, month) group.
#'   \item Apply **precipitation intensity perturbations** to wet days using
#'     Gamma-based quantile mapping via \code{\link{adjust_precipitation_qm}}.
#'   \item Apply **temperature perturbations** using monthly additive deltas
#'     (step change or transient ramp).
#'   \item Optionally recompute **PET** from perturbed temperature fields.
#' }
#'
#' This function is intended for **climate stress testing and scenario analysis**,
#' where controlled changes in mean climate, variability, and extremes are required
#' while preserving realistic daily structure.
#'
#' @section Year indexing convention (critical):
#' All precipitation perturbations rely on a **simulation-year index**
#' \code{year_idx = 1..n_years}, derived internally from \code{date}.
#' Calendar years are \strong{not} passed downstream.
#' Any factor matrices supplied as \code{n_years x 12} are indexed using this
#' simulation-year convention.
#'
#' @param data
#' List of data.frames, one per grid cell. Each data.frame must contain:
#' \itemize{
#'   \item \code{prcp} - daily precipitation (mm/day),
#'   \item \code{temp} - daily mean temperature (degC),
#'   \item \code{temp_min} - daily minimum temperature (degC),
#'   \item \code{temp_max} - daily maximum temperature (degC).
#' }
#' If \code{compute_pet = TRUE}, a column \code{pet} is added or overwritten.
#'
#' @param grid
#' Data frame of grid metadata with \code{nrow(grid) == length(data)}.
#' Must include a latitude column: \code{lat} (preferred) or \code{y} (legacy).
#'
#' @param date
#' Date vector with length equal to the number of rows in each grid-cell data.frame.
#' Used to derive month and simulation-year indices.
#'
#' @param prcp_mean_factor
#' Numeric vector of length 12 or numeric matrix \code{(n_years x 12)}.
#' Monthly multiplicative factors applied to **wet-day precipitation mean**
#' during intensity perturbation.
#'
#' @param prcp_var_factor
#' Numeric vector of length 12 or numeric matrix \code{(n_years x 12)}.
#' Monthly multiplicative factors applied to **wet-day precipitation variance**.
#' If \code{scale_var_with_mean = TRUE}, this factor is combined with mean scaling
#' (see below).
#'
#' @param prcp_occurrence_factor
#' Optional numeric vector of length 12 or numeric matrix \code{(n_years x 12)}.
#' Multiplicative factors applied to **monthly wet-day probability**.
#' Values greater than 1 increase wet-day frequency; values less than 1 decrease it.
#'
#' @param prcp_intensity_threshold
#' Numeric scalar \code{>= 0}. Defines wet days as \code{prcp > prcp_intensity_threshold}.
#' This threshold is:
#' \itemize{
#'   \item used to define wet days for occurrence perturbation, and
#'   \item passed to \code{\link{adjust_precipitation_qm}} as its
#'     \code{intensity_threshold}.
#' }
#' Default is \code{0}.
#'
#' @param temp_delta
#' Numeric vector of length 12. Monthly additive temperature deltas (degC) applied
#' to \code{temp}, \code{temp_min}, and \code{temp_max}.
#'
#' @param temp_transient
#' Logical. If \code{TRUE}, temperature deltas ramp linearly from zero to twice
#' \code{temp_delta} across simulation years, preserving the same mean change.
#'
#' @param prcp_transient
#' Logical. If \code{TRUE}, precipitation mean and variance factors ramp linearly
#' across simulation years using the same transient logic.
#'
#' @param prcp_occurrence_transient
#' Logical. If \code{TRUE}, precipitation occurrence factors ramp linearly
#' across simulation years.
#'
#' @param compute_pet
#' Logical. If \code{TRUE}, recompute PET from perturbed temperature using
#' \code{\link{calculate_monthly_pet}}.
#'
#' @param pet_method
#' Character string specifying the PET method passed to
#' \code{\link{calculate_monthly_pet}} (default: \code{"hargreaves"}).
#'
#' @param qm_fit_method
#' Character string specifying the Gamma fitting method passed to
#' \code{\link{adjust_precipitation_qm}} (e.g., \code{"mme"}, \code{"mle"}).
#'
#' @param scale_var_with_mean
#' Logical. If \code{TRUE}, precipitation variance scaling is computed as:
#' \deqn{var\_factor\_use = prcp\_var\_factor \times prcp\_mean\_factor^2}
#' This tends to preserve the coefficient of variation under mean changes.
#'
#' @param exaggerate_extremes
#' Logical. If TRUE, amplify upper-tail wet-day precipitation during the
#' intensity quantile mapping step (forwarded to \code{\link{adjust_precipitation_qm}}).
#'
#' @param extreme_prob_threshold
#' Numeric scalar in (0, 1). Probability threshold defining the start of the tail
#' region for amplification in \code{adjust_precipitation_qm}. For example, 0.95 means
#' the top 5% of baseline wet-day intensities (month-specific).
#'
#' @param extreme_k
#' Numeric scalar > 0. Tail exponent controlling amplification strength in
#' \code{adjust_precipitation_qm}. Values > 1 amplify extremes; values in (0, 1) dampen.
#'
#' @param enforce_target_mean
#' Logical. If TRUE, rescale mapped wet-day values within each (year index, month)
#' group so the wet-day mean matches the intended target mean after tail amplification.
#' Forwarded to \code{\link{adjust_precipitation_qm}}.
#'
#' @param prcp_cap_mm_day
#' Optional numeric scalar. Absolute upper cap (mm/day) applied to precipitation
#' after all perturbations. Intended as a numerical safeguard.
#'
#' @param prcp_floor_mm_day
#' Optional numeric scalar. Minimum precipitation amount (mm/day) applied to wet days
#' after perturbation. Useful to avoid unrealistically small drizzle amounts.
#'
#' @param prcp_cap_quantile
#' Optional numeric scalar in \code{(0, 1)}. Quantile-based cap computed from the
#' perturbed precipitation series and applied as an additional upper bound.
#'
#' @param seed
#' Optional integer. Sets the random seed for reproducible precipitation occurrence
#' perturbations and any stochastic components. The previous RNG state is restored on exit.
#'
#' @param verbose
#' Logical. If \code{TRUE}, prints progress messages and warnings.
#'
#' @return
#' A list of data.frames, one per grid cell, containing the perturbed weather variables.
#' The structure matches the input \code{data} list.
#'
#' @seealso
#' \code{\link{adjust_precipitation_qm}},
#' \code{\link{calculate_monthly_pet}}
#'
#' @export
apply_climate_perturbations <- function(
    data = NULL,
    grid = NULL,
    date = NULL,
    prcp_mean_factor = NULL,
    prcp_var_factor = NULL,
    prcp_occurrence_factor = NULL,
    prcp_intensity_threshold = 0,
    temp_delta = NULL,
    temp_transient = TRUE,
    prcp_transient = TRUE,
    prcp_occurrence_transient = TRUE,
    compute_pet = TRUE,
    pet_method = "hargreaves",
    qm_fit_method = "mme",
    scale_var_with_mean = TRUE,
    exaggerate_extremes = FALSE,
    extreme_prob_threshold = 0.95,
    extreme_k = 1.2,
    enforce_target_mean = TRUE,
    prcp_cap_mm_day = NULL,
    prcp_floor_mm_day = NULL,
    prcp_cap_quantile = NULL,
    seed = NULL,
    verbose = FALSE) {

  # --------------------------------------------------------------------------
  # INPUT VALIDATION (legacy messages preserved)
  # --------------------------------------------------------------------------

  if (is.null(data)) stop("'climate.data' must not be NULL", call. = FALSE)
  if (is.null(grid)) stop("'grid' must not be NULL", call. = FALSE)
  if (is.null(date)) stop("'sim.dates' must not be NULL", call. = FALSE)
  if (is.null(prcp_mean_factor)) stop("'change.factor.precip.mean' must not be NULL", call. = FALSE)
  if (is.null(prcp_var_factor) && !isTRUE(scale_var_with_mean)) {
    stop("'change.factor.precip.variance' must not be NULL", call. = FALSE)
  }
  if (is.null(temp_delta)) stop("'change.factor.temp.mean' must not be NULL", call. = FALSE)

  if (!is.list(data)) stop("'climate.data' must be a list of data frames", call. = FALSE)
  if (!is.data.frame(grid)) stop("'grid' must be a data frame", call. = FALSE)
  if (!inherits(date, "Date")) stop("'sim.dates' must be a Date vector", call. = FALSE)

  n_grid <- length(data)
  if (n_grid != nrow(grid)) {
    stop(
      "Length of 'climate.data' (", n_grid, ") must match ",
      "number of rows in 'grid' (", nrow(grid), ")",
      call. = FALSE
    )
  }

  if (!("lat" %in% names(grid)) && !("y" %in% names(grid))) {
    stop("'grid' must contain a 'y' column (latitude)", call. = FALSE)
  }
  lat_col <- if ("lat" %in% names(grid)) "lat" else "y"

  required_cols <- c("prcp", "temp", "temp_min", "temp_max")
  n_day <- length(date)

  for (i in seq_along(data)) {
    miss <- setdiff(required_cols, names(data[[i]]))
    if (length(miss) > 0) {
      stop("Grid cell ", i, " missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
    }
    if (nrow(data[[i]]) != n_day) {
      stop("Grid cell ", i, " row count does not match 'sim.dates'", call. = FALSE)
    }
  }

  if (!is.numeric(prcp_intensity_threshold) || length(prcp_intensity_threshold) != 1L ||
      !is.finite(prcp_intensity_threshold) || prcp_intensity_threshold < 0) {
    stop("'prcp_intensity_threshold' must be a single finite numeric >= 0", call. = FALSE)
  }

  if (!is.null(prcp_cap_mm_day)) {
    if (!is.numeric(prcp_cap_mm_day) || length(prcp_cap_mm_day) != 1L ||
        !is.finite(prcp_cap_mm_day) || prcp_cap_mm_day <= 0) {
      stop("'prcp_cap_mm_day' must be a single positive finite numeric", call. = FALSE)
    }
  }

  if (!is.null(prcp_floor_mm_day)) {
    if (!is.numeric(prcp_floor_mm_day) || length(prcp_floor_mm_day) != 1L ||
        !is.finite(prcp_floor_mm_day) || prcp_floor_mm_day < 0) {
      stop("'prcp_floor_mm_day' must be a single finite numeric >= 0", call. = FALSE)
    }
  }

  if (!is.null(prcp_cap_quantile)) {
    if (!is.numeric(prcp_cap_quantile) || length(prcp_cap_quantile) != 1L ||
        !is.finite(prcp_cap_quantile) || prcp_cap_quantile <= 0 || prcp_cap_quantile >= 1) {
      stop("'prcp_cap_quantile' must be a single numeric in (0, 1)", call. = FALSE)
    }
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be a single finite numeric/integer", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  # --------------------------------------------------------------------------
  # TIME INDICES (simulation-year index 1..n_years)
  # --------------------------------------------------------------------------

  cal_year <- as.integer(format(date, "%Y"))
  year_idx <- cal_year - min(cal_year) + 1L
  month <- as.integer(format(date, "%m"))
  n_year <- max(year_idx)

  # --------------------------------------------------------------------------
  # HELPERS
  # --------------------------------------------------------------------------

  .as_change_matrix <- function(x, name) {
    if (is.null(dim(x))) {
      if (!is.numeric(x) || length(x) != 12L || any(!is.finite(x))) {
        stop("'", name, "' must be a finite numeric vector of length 12 or a matrix n_years x 12", call. = FALSE)
      }
      matrix(x, nrow = n_year, ncol = 12, byrow = TRUE)
    } else {
      if (!is.matrix(x) || ncol(x) != 12L) {
        stop("'", name, "' must have 12 columns (months)", call. = FALSE)
      }
      if (nrow(x) != n_year) {
        stop("'", name, "' must have nrow == number of years in 'date'", call. = FALSE)
      }
      if (any(!is.finite(x))) {
        stop("'", name, "' contains non-finite values", call. = FALSE)
      }
      x
    }
  }

  # Return n_year x 12 (years x months)
  .ramp_matrix <- function(mat_year_month) {
    out <- matrix(NA_real_, nrow = n_year, ncol = 12)

    for (m in 1:12) {
      end_m <- (mat_year_month[1, m] - 1) * 2 + 1
      out[, m] <- seq(1, end_m, length.out = n_year)
    }

    out
  }

  .apply_prcp_safety_rails <- function(x) {
    if (!is.null(prcp_floor_mm_day)) {
      wet <- x > prcp_intensity_threshold
      wet[is.na(wet)] <- FALSE
      if (any(wet)) x[wet] <- pmax(x[wet], prcp_floor_mm_day)
    }

    cap_mm <- prcp_cap_mm_day

    if (!is.null(prcp_cap_quantile)) {
      cap_q <- suppressWarnings(stats::quantile(x, probs = prcp_cap_quantile, na.rm = TRUE, names = FALSE))
      if (is.finite(cap_q)) {
        cap_mm <- if (is.null(cap_mm)) cap_q else min(cap_mm, cap_q)
      }
    }

    if (!is.null(cap_mm) && is.finite(cap_mm)) {
      x <- pmin(x, cap_mm)
    }

    neg <- x < 0
    neg[is.na(neg)] <- FALSE
    if (any(neg)) x[neg] <- 0

    if (is.numeric(x)) {
      bad <- is.infinite(x) | is.nan(x)
      bad[is.na(bad)] <- FALSE
      if (any(bad)) x[bad] <- 0
    }

    x
  }

  .adjust_occurrence <- function(prcp_vec, mon_vec, yidx_vec, base_gamma, target_gamma, occ_mat) {
    keys <- unique(data.frame(mon = mon_vec, year = yidx_vec, stringsAsFactors = FALSE))

    for (k in seq_len(nrow(keys))) {
      m <- keys$mon[k]
      y <- keys$year[k]

      idx <- which(mon_vec == m & yidx_vec == y)
      if (length(idx) == 0L) next

      wet <- prcp_vec[idx] > prcp_intensity_threshold
      wet[is.na(wet)] <- FALSE

      n_wet <- sum(wet)
      n_tot <- length(idx)

      p_fac <- occ_mat[y, m]
      if (!is.finite(p_fac) || p_fac <= 0) next

      p0 <- if (n_tot > 0L) (n_wet / n_tot) else 0
      p1 <- pmin(1, p_fac * p0)
      n_target <- as.integer(round(p1 * n_tot))
      n_target <- max(0L, min(n_tot, n_target))

      if (n_target == n_wet) next

      i_month <- which(base_gamma$month == m)
      if (length(i_month) != 1L) next

      if (n_target < n_wet) {
        n_remove <- n_wet - n_target
        wet_idx <- which(wet)
        drop_local <- sample(wet_idx, size = n_remove, replace = FALSE)
        prcp_vec[idx[drop_local]] <- 0
      } else {
        n_add <- n_target - n_wet
        dry_idx <- which(!wet)
        if (length(dry_idx) == 0L) next

        add_local <- sample(dry_idx, size = min(n_add, length(dry_idx)), replace = FALSE)

        shape <- target_gamma$shape[i_month, y]
        scale <- target_gamma$scale[i_month, y]
        if (!is.finite(shape) || !is.finite(scale) || shape <= 0 || scale <= 0) next

        new_amt <- stats::rgamma(length(add_local), shape = shape, scale = scale)

        too_small <- new_amt <= prcp_intensity_threshold
        if (any(too_small)) {
          new_amt[too_small] <- prcp_intensity_threshold + 1e-06
        }

        prcp_vec[idx[add_local]] <- new_amt
      }
    }

    prcp_vec
  }

  # --------------------------------------------------------------------------
  # TEMPERATURE DELTAS
  # --------------------------------------------------------------------------

  if (!is.numeric(temp_delta) || length(temp_delta) != 12L || any(!is.finite(temp_delta))) {
    stop("'temp_delta' must be a finite numeric vector of length 12", call. = FALSE)
  }

  # Build n_year x 12
  temp_mat <- matrix(NA_real_, nrow = n_year, ncol = 12)
  for (m in 1:12) {
    temp_mat[, m] <- if (isTRUE(temp_transient)) {
      seq(0, 2 * temp_delta[m], length.out = n_year)
    } else {
      rep(temp_delta[m], n_year)
    }
  }
  temp_day <- temp_mat[cbind(year_idx, month)]

  # --------------------------------------------------------------------------
  # PRECIPITATION FACTORS (intensity)
  # --------------------------------------------------------------------------

  min_factor <- 0.01

  mean_mat_in <- .as_change_matrix(prcp_mean_factor, "prcp_mean_factor")

  if (isTRUE(scale_var_with_mean)) {
    if (!is.null(prcp_var_factor) && isTRUE(verbose)) {
      warning("Ignoring 'prcp_var_factor' because 'scale_var_with_mean = TRUE'.", call. = FALSE)
    }
    var_mat_in <- mean_mat_in^2
  } else {
    var_mat_in <- .as_change_matrix(prcp_var_factor, "prcp_var_factor")
  }

  mean_mat <- if (isTRUE(prcp_transient)) .ramp_matrix(mean_mat_in) else mean_mat_in
  var_mat  <- if (isTRUE(prcp_transient)) .ramp_matrix(var_mat_in)  else var_mat_in

  mean_mat <- pmax(mean_mat, min_factor)
  var_mat  <- pmax(var_mat,  min_factor)

  # --------------------------------------------------------------------------
  # OPTIONAL: OCCURRENCE FACTORS (wet-day frequency)
  # --------------------------------------------------------------------------

  occ_mat <- NULL
  if (!is.null(prcp_occurrence_factor)) {
    occ_in <- .as_change_matrix(prcp_occurrence_factor, "prcp_occurrence_factor")
    occ_mat <- if (isTRUE(prcp_occurrence_transient)) .ramp_matrix(occ_in) else occ_in
    occ_mat <- pmax(occ_mat, min_factor)
  }

  # --------------------------------------------------------------------------
  # APPLY PER GRID
  # --------------------------------------------------------------------------

  out <- vector("list", n_grid)

  for (i in seq_len(n_grid)) {
    cell <- data[[i]]

    qm_out <- adjust_precipitation_qm(
      prcp = cell$prcp,
      mean_factor = mean_mat,
      var_factor = var_mat,
      scale_var_with_mean = FALSE,
      exaggerate_extremes = exaggerate_extremes,
      extreme_prob_threshold = extreme_prob_threshold,
      extreme_k = extreme_k,
      enforce_target_mean = enforce_target_mean,
      month = month,
      year = year_idx,
      intensity_threshold = prcp_intensity_threshold,
      fit_method = qm_fit_method,
      min_events = 10,
      validate_output = TRUE,
      diagnostics = !is.null(occ_mat),
      seed = seed,
      verbose = FALSE
    )

    # adjust_precipitation_qm(diagnostics=TRUE) may still return an atomic vector
    # if it early-exits (e.g., no wet days, no months meet min_events, all fits fail).
    if (is.list(qm_out)) {
      prcp_new <- qm_out$adjusted
      base_gamma <- qm_out$base_gamma
      target_gamma <- qm_out$target_gamma
    } else {
      prcp_new <- qm_out
      base_gamma <- NULL
      target_gamma <- NULL
    }

    # ----------------------------------------------------------------------
    # Optional precipitation occurrence perturbation (requires Gamma params)
    # ----------------------------------------------------------------------
    if (!is.null(occ_mat) && !is.null(base_gamma) && !is.null(target_gamma)) {
      prcp_new <- .adjust_occurrence(
        prcp_vec = prcp_new,
        mon_vec = month,
        yidx_vec = year_idx,
        base_gamma = base_gamma,
        target_gamma = target_gamma,
        occ_mat = occ_mat
      )
    } else if (!is.null(occ_mat) && isTRUE(verbose)) {
      warning(
        "Skipping precipitation occurrence perturbation for grid cell ", i,
        " because Gamma fit outputs are unavailable (likely all-dry or insufficient wet days).",
        call. = FALSE
      )
    }

    prcp_new <- .apply_prcp_safety_rails(prcp_new)
    cell$prcp <- prcp_new

    cell$temp     <- cell$temp     + temp_day
    cell$temp_min <- cell$temp_min + temp_day
    cell$temp_max <- cell$temp_max + temp_day

    if (isTRUE(compute_pet)) {
      cell$pet <- calculate_monthly_pet(
        month = month,
        temp = cell$temp,
        temp_range = cell$temp_max - cell$temp_min,
        lat_deg = grid[[lat_col]][i],
        method = pet_method
      )
      if (is.numeric(cell$pet)) {
        cell$pet[!is.finite(cell$pet)] <- 0
      }
    }

    cell[] <- lapply(cell, function(x) {
      if (is.numeric(x)) {
        x[!is.finite(x)] <- 0
      }
      x
    })

    out[[i]] <- cell
  }

  out
}
