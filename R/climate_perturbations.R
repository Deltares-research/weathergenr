# ==============================================================================
# CLIMATE PERTURBATIONS
# ==============================================================================

#' Apply monthly climate perturbations to gridded daily weather series
#'
#' @description
#' Applies **monthly** climate perturbations to **daily** gridded weather series:
#' \itemize{
#'   \item **Precipitation** is perturbed using quantile mapping with monthly
#'   mean and variance factors via \code{perturb_prcp_qm()}.
#'   \item **Temperature** (\code{temp}, \code{temp_min}, \code{temp_max}) is perturbed
#'   using monthly additive deltas (step or transient).
#'   \item **Potential evapotranspiration (PET)** can be recomputed from perturbed
#'   temperatures using \code{calculate_monthly_pet()}.
#' }
#'
#' Perturbations can be applied as a **step change** (constant in time) or as a
#' **transient change** (linearly ramping over years) while preserving the same
#' *mean* change over the simulation period.
#'
#' @param data List of data.frames, one per grid cell. Each data.frame must contain:
#' \itemize{
#'   \item \code{prcp} (daily precipitation)
#'   \item \code{temp} (daily mean temperature)
#'   \item \code{temp_min} (daily minimum temperature)
#'   \item \code{temp_max} (daily maximum temperature)
#' }
#' If \code{compute_pet = TRUE}, \code{pet} is added or overwritten.
#'
#' @param grid data.frame of grid metadata with \code{nrow(grid) == length(data)}.
#' Must include column \code{lat} (latitude in decimal degrees).
#'
#' @param date Date vector of length \code{nrow(data[[i]])} (identical across cells).
#'
#' @param prcp_mean_factor Numeric vector length 12. Monthly multiplicative factors
#' for precipitation mean.
#'
#' @param prcp_var_factor Numeric vector length 12. Monthly multiplicative factors
#' for precipitation variance.
#'
#' @param temp_delta Numeric vector length 12. Monthly additive temperature deltas (degC)
#' applied to \code{temp}, \code{temp_min}, and \code{temp_max}.
#'
#' @param temp_transient Logical. If \code{TRUE}, temperature deltas ramp linearly
#' from 0 to \code{2 * temp_delta} over years (mean equals \code{temp_delta}).
#'
#' @param prcp_transient Logical. If \code{TRUE}, precipitation factors ramp linearly
#' from 1 to \code{(factor - 1) * 2 + 1} over years.
#'
#' @param compute_pet Logical. If \code{TRUE}, recompute PET using
#' \code{\link{calculate_monthly_pet}}.
#'
#' @param pet_method Character. PET method passed to
#' \code{\link{calculate_monthly_pet}} (default: \code{"hargreaves"}).
#'
#' @param qm_fit_method Character. Distribution-fitting method for
#' \code{perturb_prcp_qm()}.
#'
#' @param verbose Logical. Emit progress logs.
#'
#' @return List of data.frames with perturbed variables.
#'
#' @seealso \code{\link{perturb_prcp_qm}}, \code{\link{diagnose_prcp_qm}},
#'   \code{\link{calculate_monthly_pet}}
#'
#' @export
apply_climate_perturbations <- function(
    data = NULL,
    grid = NULL,
    date = NULL,
    prcp_mean_factor = NULL,
    prcp_var_factor = NULL,
    temp_delta = NULL,
    temp_transient = TRUE,
    prcp_transient = TRUE,
    compute_pet = TRUE,
    pet_method = "hargreaves",
    qm_fit_method = "mme",
    verbose = FALSE) {

  # --------------------------------------------------------------------------
  # INPUT VALIDATION (legacy messages preserved)
  # --------------------------------------------------------------------------

  if (is.null(data)) stop("'climate.data' must not be NULL", call. = FALSE)
  if (is.null(grid)) stop("'climate.grid' must not be NULL", call. = FALSE)
  if (is.null(date)) stop("'sim.dates' must not be NULL", call. = FALSE)
  if (is.null(prcp_mean_factor)) stop("'change.factor.precip.mean' must not be NULL", call. = FALSE)
  if (is.null(prcp_var_factor)) stop("'change.factor.precip.variance' must not be NULL", call. = FALSE)
  if (is.null(temp_delta)) stop("'change.factor.temp.mean' must not be NULL", call. = FALSE)

  if (!is.list(data)) stop("'climate.data' must be a list of data frames", call. = FALSE)
  if (!is.data.frame(grid)) stop("'climate.grid' must be a data frame", call. = FALSE)
  if (!inherits(date, "Date")) stop("'sim.dates' must be a Date vector", call. = FALSE)

  n_grid <- length(data)

  if (n_grid != nrow(grid)) {
    stop(
      "Length of 'climate.data' (", n_grid, ") must match ",
      "number of rows in 'climate.grid' (", nrow(grid), ")",
      call. = FALSE
    )
  }

  if (!"lat" %in% names(grid)) {
    stop("'climate.grid' must contain a 'y' column (latitude)", call. = FALSE)
  }

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

  # --------------------------------------------------------------------------
  # TIME INDICES
  # --------------------------------------------------------------------------

  year <- as.integer(format(date, "%Y"))
  year_idx <- year - min(year) + 1L
  month <- as.integer(format(date, "%m"))
  n_year <- max(year_idx)

  # --------------------------------------------------------------------------
  # HELPERS
  # --------------------------------------------------------------------------

  .as_change_matrix <- function(x) {
    if (is.null(dim(x))) matrix(x, nrow = n_year, ncol = 12) else x
  }

  # --------------------------------------------------------------------------
  # TEMPERATURE DELTAS
  # --------------------------------------------------------------------------

  temp_mat <- if (temp_transient) {
    vapply(1:12, function(m) seq(0, 2 * temp_delta[m], length.out = n_year),
           FUN.VALUE = numeric(n_year))
  } else {
    vapply(1:12, function(m) rep(temp_delta[m], n_year),
           FUN.VALUE = numeric(n_year))
  }
  temp_mat <- .as_change_matrix(temp_mat)
  temp_day <- temp_mat[cbind(year_idx, month)]

  # --------------------------------------------------------------------------
  # PRECIPITATION FACTORS
  # --------------------------------------------------------------------------

  min_factor <- 0.01

  if (prcp_transient) {
    mean_end <- (prcp_mean_factor - 1) * 2 + 1
    var_end  <- (prcp_var_factor  - 1) * 2 + 1

    mean_mat <- vapply(1:12, function(m) seq(1, mean_end[m], length.out = n_year),
                       FUN.VALUE = numeric(n_year))
    var_mat  <- vapply(1:12, function(m) seq(1, var_end[m],  length.out = n_year),
                       FUN.VALUE = numeric(n_year))
  } else {
    mean_mat <- vapply(1:12, function(m) rep(prcp_mean_factor[m], n_year),
                       FUN.VALUE = numeric(n_year))
    var_mat  <- vapply(1:12, function(m) rep(prcp_var_factor[m], n_year),
                       FUN.VALUE = numeric(n_year))
  }

  mean_mat <- pmax(.as_change_matrix(mean_mat), min_factor)
  var_mat  <- pmax(.as_change_matrix(var_mat),  min_factor)

  # --------------------------------------------------------------------------
  # APPLY PER GRID
  # --------------------------------------------------------------------------

  out <- vector("list", n_grid)

  for (i in seq_len(n_grid)) {

    cell <- data[[i]]

    cell$prcp <- perturb_prcp_qm(
      prcp = cell$prcp,
      month = month,
      year = year_idx,
      mean_factor = mean_mat,
      var_factor = var_mat,
      fit_method = qm_fit_method,
      verbose = FALSE
    )

    cell$temp     <- cell$temp     + temp_day
    cell$temp_min <- cell$temp_min + temp_day
    cell$temp_max <- cell$temp_max + temp_day

    if (compute_pet) {
      cell$pet <- calculate_monthly_pet(
        month = month,
        temp = cell$temp,
        temp_range = cell$temp_max - cell$temp_min,
        lat_deg = grid$lat[i],
        method = pet_method
      )
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
