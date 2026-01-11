#' Apply Climate Change Perturbations to Gridded Weather Data
#'
#' Applies climate change perturbations to gridded daily weather data by modifying
#' precipitation and temperature using monthly change factors. Perturbations can be
#' applied as step changes or gradual transient changes. Both approaches are
#' constructed to yield the same mean change over the simulation period.
#'
#' @param climate.data List of data frames, one per grid cell. Each data frame must
#'   contain columns \code{precip}, \code{temp}, \code{temp_min}, \code{temp_max},
#'   and optionally \code{pet}.
#' @param climate.grid Data frame containing grid cell metadata. Must include a
#'   column \code{y} representing latitude in decimal degrees.
#' @param sim.dates Vector of \code{Date} objects corresponding to the rows of each
#'   element in \code{climate.data}.
#' @param change.factor.precip.mean Numeric vector of length 12. Monthly
#'   multiplicative factors for precipitation mean representing the desired average
#'   change over the simulation period.
#' @param change.factor.precip.variance Numeric vector of length 12. Monthly
#'   multiplicative factors applied to precipitation variance.
#' @param change.factor.temp.mean Numeric vector of length 12. Monthly additive
#'   temperature changes in degrees Celsius representing the desired average change.
#' @param transient.temp.change Logical. If TRUE, temperature changes increase
#'   linearly from zero to twice the specified factor so that the target mean change
#'   is achieved. If FALSE, the full factor is applied uniformly.
#' @param transient.precip.change Logical. If TRUE, precipitation change factors
#'   increase linearly from 1.0 to a calculated endpoint that yields the target mean.
#'   If FALSE, the full factor is applied uniformly.
#' @param calculate.pet Logical. If TRUE, potential evapotranspiration is
#'   recalculated from perturbed temperatures using the Hargreaves method.
#' @param fit.method Character string specifying the distribution fitting method
#'   passed to \code{quantile_mapping}.
#' @param verbose Logical. If TRUE, progress messages are printed.
#'
#' @details
#' Climate change perturbations are applied using the following mechanisms:
#'
#' \itemize{
#'   \item Precipitation is adjusted using quantile mapping with mean and variance scaling.
#'   \item Temperature variables are adjusted using additive monthly shifts.
#'   \item Potential evapotranspiration is recalculated using the Hargreaves equation
#'   when requested.
#' }
#'
#' Step and transient perturbations are formulated to produce the same mean change
#' over the full simulation period. Step changes apply the specified factor
#' uniformly across all years. Transient changes increase linearly from a baseline
#' to an endpoint that yields the same mean effect.
#'
#' For precipitation, transient factors ramp from 1.0 to a computed upper value.
#' For temperature, transient changes ramp from zero to twice the specified mean
#' change. This ensures comparability between step and transient scenarios while
#' representing different temporal pathways.
#'
#' During transient precipitation adjustments, a minimum factor of 0.01 is enforced
#' to prevent numerical instabilities.
#'
#' @return A list of data frames with the same structure as \code{climate.data},
#'   containing perturbed weather variables. Invalid values are replaced with zero.
#'
#' @examples
#' \dontrun{
#' perturbed_step <- apply_climate_perturbations(
#'   climate.data = weather_grids,
#'   climate.grid = grid_meta,
#'   sim.dates = dates,
#'   change.factor.precip.mean = rep(1.2, 12),
#'   change.factor.precip.variance = rep(1.0, 12),
#'   change.factor.temp.mean = rep(2.0, 12),
#'   transient.temp.change = FALSE,
#'   transient.precip.change = FALSE
#' )
#' }
#'
#' @seealso
#' \code{\link{quantile_mapping}}
#'
#' @importFrom logger log_info
#' @import dplyr
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
  # INPUT VALIDATION
  # ==========================================================================

  if (is.null(climate.data)) stop("'climate.data' must not be NULL")
  if (is.null(climate.grid)) stop("'climate.grid' must not be NULL")
  if (is.null(sim.dates)) stop("'sim.dates' must not be NULL")
  if (is.null(change.factor.precip.mean)) stop("'change.factor.precip.mean' must not be NULL")
  if (is.null(change.factor.precip.variance)) stop("'change.factor.precip.variance' must not be NULL")
  if (is.null(change.factor.temp.mean)) stop("'change.factor.temp.mean' must not be NULL")

  if (!is.list(climate.data)) stop("'climate.data' must be a list of data frames")
  if (!is.data.frame(climate.grid)) stop("'climate.grid' must be a data frame")
  if (!inherits(sim.dates, "Date")) stop("'sim.dates' must be a Date vector")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("'verbose' must be logical (TRUE/FALSE)")

  .log <- function(fmt, ...) {
    if (isTRUE(verbose)) logger::log_info(fmt, ...)
    invisible(NULL)
  }

  # Check dimensions
  ngrids <- length(climate.data)

  if (ngrids != nrow(climate.grid)) {
    stop(
      "Length of 'climate.data' (", ngrids, ") must match ",
      "number of rows in 'climate.grid' (", nrow(climate.grid), ")"
    )
  }

  if (!"y" %in% names(climate.grid)) {
    stop("'climate.grid' must contain a 'y' column (latitude)")
  }

  if (length(change.factor.precip.mean) != 12) {
    stop("'change.factor.precip.mean' must have length 12 (one per month)")
  }
  if (length(change.factor.precip.variance) != 12) {
    stop("'change.factor.precip.variance' must have length 12 (one per month)")
  }
  if (length(change.factor.temp.mean) != 12) {
    stop("'change.factor.temp.mean' must have length 12 (one per month)")
  }

  if (any(change.factor.precip.mean <= 0)) stop("'change.factor.precip.mean' must contain positive values")
  if (any(change.factor.precip.variance <= 0)) stop("'change.factor.precip.variance' must contain positive values")

  required_cols <- c("precip", "temp", "temp_min", "temp_max")
  for (i in seq_along(climate.data)) {
    missing_cols <- setdiff(required_cols, names(climate.data[[i]]))
    if (length(missing_cols) > 0) {
      stop(
        "Grid cell ", i, " is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    if (nrow(climate.data[[i]]) != length(sim.dates)) {
      stop(
        "Grid cell ", i, " has ", nrow(climate.data[[i]]), " rows but ",
        "'sim.dates' has length ", length(sim.dates)
      )
    }
  }

  # ==========================================================================
  # TEMPORAL INDICES
  # ==========================================================================

  year_vec  <- as.integer(format(sim.dates, "%Y"))
  year_ind  <- year_vec - min(year_vec) + 1
  month_ind <- as.integer(format(sim.dates, "%m"))

  n_years <- max(year_ind)
  n_days  <- length(sim.dates)

  .log(
    "[PERTURBATION] Simulation period: {n_years} years ({min(year_vec)}-{max(year_vec)}), {n_days} days"
  )

  # ==========================================================================
  # COMPUTE TEMPERATURE CHANGE FACTORS
  # ==========================================================================

  if (transient.temp.change) {
    temp_change_matrix <- sapply(1:12, function(m) {
      seq(0, change.factor.temp.mean[m] * 2, length.out = n_years)
    })
    if (!is.matrix(temp_change_matrix)) {
      temp_change_matrix <- matrix(temp_change_matrix, nrow = 1, ncol = 12)
    }
  } else {
    temp_change_matrix <- sapply(1:12, function(m) {
      rep(change.factor.temp.mean[m], n_years)
    })
    if (!is.matrix(temp_change_matrix)) {
      temp_change_matrix <- matrix(temp_change_matrix, nrow = 1, ncol = 12)
    }
  }

  temp_change_daily <- temp_change_matrix[cbind(year_ind, month_ind)]

  .log(
    "[PERTURBATION] Temperature: {if (isTRUE(transient.temp.change)) 'transient' else 'step'} change (mean = {round(mean(change.factor.temp.mean), 2)} DegC)"
  )

  # ==========================================================================
  # COMPUTE PRECIPITATION CHANGE FACTORS
  # ==========================================================================

  min_factor <- 0.01

  if (transient.precip.change) {

    precip_mean_deltaf <- (change.factor.precip.mean - 1) * 2 + 1
    precip_mean_matrix <- sapply(1:12, function(m) {
      seq(1, precip_mean_deltaf[m], length.out = n_years)
    })
    precip_mean_matrix[precip_mean_matrix < min_factor] <- min_factor
    if (!is.matrix(precip_mean_matrix)) {
      precip_mean_matrix <- matrix(precip_mean_matrix, nrow = 1, ncol = 12)
    }

    precip_var_deltaf <- (change.factor.precip.variance - 1) * 2 + 1
    precip_var_matrix <- sapply(1:12, function(m) {
      seq(1, precip_var_deltaf[m], length.out = n_years)
    })
    precip_var_matrix[precip_var_matrix < min_factor] <- min_factor
    if (!is.matrix(precip_var_matrix)) {
      precip_var_matrix <- matrix(precip_var_matrix, nrow = 1, ncol = 12)
    }

  } else {

    precip_mean_matrix <- sapply(1:12, function(m) rep(change.factor.precip.mean[m], n_years))
    if (!is.matrix(precip_mean_matrix)) {
      precip_mean_matrix <- matrix(precip_mean_matrix, nrow = 1, ncol = 12)
    }

    precip_var_matrix <- sapply(1:12, function(m) rep(change.factor.precip.variance[m], n_years))
    if (!is.matrix(precip_var_matrix)) {
      precip_var_matrix <- matrix(precip_var_matrix, nrow = 1, ncol = 12)
    }
  }

  .log(
    "[PERTURBATION] Precipitation: {if (isTRUE(transient.precip.change)) 'transient' else 'step'} change (mean = {round(mean(change.factor.precip.mean), 3)})"
  )

  # ==========================================================================
  # APPLY PERTURBATIONS TO EACH GRID CELL
  # ==========================================================================

  .log("[PERTURBATION] Applying perturbations to {ngrids} grid cells")

  latitudes <- climate.grid$y

  climate.data <- lapply(seq_len(ngrids), function(i) {

    if (isTRUE(verbose) && (i == 1L || i %% 100L == 0L || i == ngrids)) {
      .log("[PERTURBATION] Processing grid {i} of {ngrids}")
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

    grid_data$temp     <- grid_data$temp + temp_change_daily
    grid_data$temp_min <- grid_data$temp_min + temp_change_daily
    grid_data$temp_max <- grid_data$temp_max + temp_change_daily

    if (calculate.pet) {
      grid_data$pet <- pet_hargreaves(
        months = month_ind,
        temp = grid_data$temp,
        tdiff = grid_data$temp_max - grid_data$temp_min,
        lat = latitudes[i]
      )
    }

    grid_data[] <- lapply(grid_data, function(col) {
      replace(col, is.infinite(col) | is.nan(col), 0)
    })

    grid_data
  })

  .log("[PERTURBATION] Perturbation complete")

  climate.data
}
