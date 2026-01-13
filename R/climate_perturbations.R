# ==============================================================================
# CLIMATE PERTURBATIONS
# ==============================================================================

#' Apply climate-change perturbations to gridded daily weather data
#'
#' @description
#' Applies monthly climate-change perturbations to gridded daily weather data by:
#' \itemize{
#'   \item adjusting precipitation via quantile mapping with monthly mean/variance factors, and
#'   \item shifting temperature (mean/min/max) via monthly additive deltas.
#' }
#'
#' Perturbations can be applied as a step change (constant over time) or a transient
#' change (linearly ramping over years) while preserving the same *mean* change over
#' the simulation period.
#'
#' @param climate.data List of data.frames, one per grid cell. Each data.frame must
#'   contain columns \code{precip}, \code{temp}, \code{temp_min}, \code{temp_max}.
#'   If \code{calculate.pet = TRUE}, \code{pet} will be (re)computed and added/overwritten.
#' @param climate.grid data.frame of grid metadata. Must have \code{nrow(climate.grid) == length(climate.data)}
#'   and include column \code{y} (latitude in decimal degrees).
#' @param sim.dates Date vector of length \code{nrow(climate.data[[i]])} for all cells.
#' @param change.factor.precip.mean Numeric vector length 12. Monthly multiplicative factors
#'   for precipitation mean (target mean change over the simulation period).
#' @param change.factor.precip.variance Numeric vector length 12. Monthly multiplicative factors
#'   for precipitation variance.
#' @param change.factor.temp.mean Numeric vector length 12. Monthly additive temperature deltas (°C).
#' @param transient.temp.change Logical. If TRUE, temperature deltas ramp linearly from 0 to
#'   \code{2 * change.factor.temp.mean} over years (so the average over years equals the specified delta).
#'   If FALSE, applies step deltas (constant over years).
#' @param transient.precip.change Logical. If TRUE, precipitation factors ramp linearly from 1 to
#'   \code{(factor - 1) * 2 + 1} over years (so the average over years equals the specified factor).
#'   If FALSE, applies step factors (constant over years).
#' @param calculate.pet Logical. If TRUE, recompute PET using \code{\link{pet_hargreaves}} from perturbed temperatures.
#' @param fit.method Character. Distribution fitting method passed to \code{\link{quantile_mapping}}.
#' @param verbose Logical scalar. If TRUE, emit progress logs via \code{.log_info()}.
#'
#' @details
#' Transient precipitation factors are constrained to be at least 0.01 to avoid numerical issues.
#' After perturbations, any \code{Inf} or \code{NaN} values are replaced with 0 (for all columns).
#'
#' @return List of data.frames, same length and row counts as \code{climate.data}, containing perturbed variables.
#'
#' @seealso \code{\link{quantile_mapping}}, \code{\link{pet_hargreaves}}
#'
#' @export
apply_climate_perturbations <- function(
    climate.data = NULL,
    climate.grid = NULL,
    sim.dates = NULL,
    change.factor.precip.mean = NULL,
    change.factor.precip.variance = NULL,
    change.factor.temp.mean = NULL,
    transient.temp.change = TRUE,
    transient.precip.change = TRUE,
    calculate.pet = TRUE,
    fit.method = "mme",
    verbose = FALSE) {

  # ==========================================================================
  # INPUT VALIDATION (keep legacy error messages for tests)
  # ==========================================================================

  if (is.null(climate.data)) stop("'climate.data' must not be NULL", call. = FALSE)
  if (is.null(climate.grid)) stop("'climate.grid' must not be NULL", call. = FALSE)
  if (is.null(sim.dates)) stop("'sim.dates' must not be NULL", call. = FALSE)
  if (is.null(change.factor.precip.mean)) stop("'change.factor.precip.mean' must not be NULL", call. = FALSE)
  if (is.null(change.factor.precip.variance)) stop("'change.factor.precip.variance' must not be NULL", call. = FALSE)
  if (is.null(change.factor.temp.mean)) stop("'change.factor.temp.mean' must not be NULL", call. = FALSE)

  if (!is.list(climate.data)) stop("'climate.data' must be a list of data frames", call. = FALSE)
  if (!is.data.frame(climate.grid)) stop("'climate.grid' must be a data frame", call. = FALSE)
  if (!inherits(sim.dates, "Date")) stop("'sim.dates' must be a Date vector", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L) stop("'verbose' must be logical (TRUE/FALSE)", call. = FALSE)

  ngrids <- length(climate.data)

  if (ngrids != nrow(climate.grid)) {
    stop(
      "Length of 'climate.data' (", ngrids, ") must match ",
      "number of rows in 'climate.grid' (", nrow(climate.grid), ")",
      call. = FALSE
    )
  }

  if (!"y" %in% names(climate.grid)) {
    stop("'climate.grid' must contain a 'y' column (latitude)", call. = FALSE)
  }

  if (length(change.factor.precip.mean) != 12) {
    stop("'change.factor.precip.mean' must have length 12 (one per month)", call. = FALSE)
  }
  if (length(change.factor.precip.variance) != 12) {
    stop("'change.factor.precip.variance' must have length 12 (one per month)", call. = FALSE)
  }
  if (length(change.factor.temp.mean) != 12) {
    stop("'change.factor.temp.mean' must have length 12 (one per month)", call. = FALSE)
  }

  if (any(change.factor.precip.mean <= 0)) stop("'change.factor.precip.mean' must contain positive values", call. = FALSE)
  if (any(change.factor.precip.variance <= 0)) stop("'change.factor.precip.variance' must contain positive values", call. = FALSE)

  required_cols <- c("precip", "temp", "temp_min", "temp_max")
  n_days <- length(sim.dates)

  for (i in seq_along(climate.data)) {
    missing_cols <- setdiff(required_cols, names(climate.data[[i]]))
    if (length(missing_cols) > 0) {
      stop(
        "Grid cell ", i, " is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
    if (nrow(climate.data[[i]]) != n_days) {
      stop(
        "Grid cell ", i, " has ", nrow(climate.data[[i]]), " rows but ",
        "'sim.dates' has length ", n_days,
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # TEMPORAL INDICES
  # ==========================================================================

  year_vec  <- as.integer(format(sim.dates, "%Y"))
  year_ind  <- year_vec - min(year_vec) + 1L
  month_ind <- as.integer(format(sim.dates, "%m"))

  n_years <- max(year_ind)

  .log_info(
    msg = sprintf(
      "Simulation period: %d years (%d-%d), %d days",
      n_years, min(year_vec), max(year_vec), n_days
    ),
    verbose = verbose,
    tag = "PERTURB"
  )

  # ==========================================================================
  # INTERNAL HELPERS
  # ==========================================================================

  .as_change_matrix <- function(mat_or_vec, n_years) {
    # Ensure matrix with nrow = n_years, ncol = 12 (handles n_years == 1)
    if (is.null(dim(mat_or_vec))) {
      mat_or_vec <- matrix(mat_or_vec, nrow = n_years, ncol = 12)
    } else {
      mat_or_vec <- as.matrix(mat_or_vec)
      if (!identical(dim(mat_or_vec), c(n_years, 12L))) {
        # vapply typically returns n_years x 12; this is a safety net
        mat_or_vec <- matrix(as.numeric(mat_or_vec), nrow = n_years, ncol = 12)
      }
    }
    mat_or_vec
  }

  # ==========================================================================
  # COMPUTE TEMPERATURE CHANGE FACTORS
  # ==========================================================================

  if (isTRUE(transient.temp.change)) {
    temp_change_matrix <- vapply(
      1:12,
      function(m) seq(0, change.factor.temp.mean[m] * 2, length.out = n_years),
      FUN.VALUE = numeric(n_years)
    )
  } else {
    temp_change_matrix <- vapply(
      1:12,
      function(m) rep(change.factor.temp.mean[m], n_years),
      FUN.VALUE = numeric(n_years)
    )
  }
  temp_change_matrix <- .as_change_matrix(temp_change_matrix, n_years)

  temp_change_daily <- temp_change_matrix[cbind(year_ind, month_ind)]

  .log_info(
    msg = sprintf(
      "Temperature: %s change (mean monthly delta = %.2f °C)",
      if (isTRUE(transient.temp.change)) "transient" else "step",
      mean(change.factor.temp.mean)
    ),
    verbose = verbose,
    tag = "PERTURB"
  )

  # ==========================================================================
  # COMPUTE PRECIPITATION CHANGE FACTORS
  # ==========================================================================

  min_factor <- 0.01

  if (isTRUE(transient.precip.change)) {

    precip_mean_end <- (change.factor.precip.mean - 1) * 2 + 1
    precip_var_end  <- (change.factor.precip.variance - 1) * 2 + 1

    precip_mean_matrix <- vapply(
      1:12,
      function(m) seq(1, precip_mean_end[m], length.out = n_years),
      FUN.VALUE = numeric(n_years)
    )
    precip_var_matrix <- vapply(
      1:12,
      function(m) seq(1, precip_var_end[m], length.out = n_years),
      FUN.VALUE = numeric(n_years)
    )

    precip_mean_matrix <- .as_change_matrix(precip_mean_matrix, n_years)
    precip_var_matrix  <- .as_change_matrix(precip_var_matrix, n_years)

    precip_mean_matrix[precip_mean_matrix < min_factor] <- min_factor
    precip_var_matrix[precip_var_matrix < min_factor] <- min_factor

  } else {

    precip_mean_matrix <- vapply(
      1:12,
      function(m) rep(change.factor.precip.mean[m], n_years),
      FUN.VALUE = numeric(n_years)
    )
    precip_var_matrix <- vapply(
      1:12,
      function(m) rep(change.factor.precip.variance[m], n_years),
      FUN.VALUE = numeric(n_years)
    )

    precip_mean_matrix <- .as_change_matrix(precip_mean_matrix, n_years)
    precip_var_matrix  <- .as_change_matrix(precip_var_matrix, n_years)
  }

  .log_info(
    msg = sprintf(
      "Precipitation: %s change (mean monthly factor = %.3f)",
      if (isTRUE(transient.precip.change)) "transient" else "step",
      mean(change.factor.precip.mean)
    ),
    verbose = verbose,
    tag = "PERTURB"
  )

  # ==========================================================================
  # APPLY PERTURBATIONS TO EACH GRID CELL
  # ==========================================================================

  .log_info(
    msg = sprintf("Applying perturbations to %d grid cells", ngrids),
    verbose = verbose,
    tag = "PERTURB"
  )

  latitudes <- climate.grid$y

  out <- vector("list", ngrids)

  for (i in seq_len(ngrids)) {

    if (isTRUE(verbose) && (i == 1L || i == ngrids || (i %% 100L) == 0L)) {
      .log_info(
        msg = sprintf("Processing grid %d of %d", i, ngrids),
        verbose = TRUE,
        tag = "PERTURB"
      )
    }

    grid_data <- climate.data[[i]]

    grid_data$precip <- quantile_mapping(
      value = grid_data$precip,
      mon.ts = month_ind,
      year.ts = year_ind,
      mean.change = precip_mean_matrix,
      var.change = precip_var_matrix,
      fit.method = fit.method,
      verbose = FALSE
    )

    grid_data$temp     <- grid_data$temp     + temp_change_daily
    grid_data$temp_min <- grid_data$temp_min + temp_change_daily
    grid_data$temp_max <- grid_data$temp_max + temp_change_daily

    if (isTRUE(calculate.pet)) {
      grid_data$pet <- pet_hargreaves(
        months = month_ind,
        temp = grid_data$temp,
        tdiff = grid_data$temp_max - grid_data$temp_min,
        lat = latitudes[i]
      )
    }

    # Replace Inf/NaN with 0 for numeric columns only
    grid_data[] <- lapply(grid_data, function(col) {
      if (!is.numeric(col)) return(col)
      col[is.infinite(col) | is.nan(col)] <- 0
      col
    })

    out[[i]] <- grid_data
  }

  .log_info(
    msg = "Perturbation complete",
    verbose = verbose,
    tag = "PERTURB"
  )

  out
}
