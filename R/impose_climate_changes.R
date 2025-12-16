#' Apply Climate Change Perturbations to Gridded Weather Data
#'
#' @description
#' Applies climate change scenarios to gridded daily weather data by perturbing
#' precipitation (via quantile mapping) and temperature (via additive shifts).
#' Supports both step-change and transient (gradual) perturbation schemes.
#' Optionally recalculates potential evapotranspiration (PET) using adjusted
#' temperature fields.
#'
#' @param climate.data List of data frames, one per grid cell. Each must contain
#'   columns: `precip`, `temp`, `temp_min`, `temp_max`. Optionally `pet`.
#' @param climate.grid Data frame with grid metadata. Must include column `y`
#'   (latitude in decimal degrees) for each grid cell.
#' @param sim.dates Vector of `Date` objects corresponding to rows in each
#'   grid cell data frame. Must not contain leap days if using 365-day calendar.
#' @param change.factor.precip.mean Numeric vector of length 12. Monthly
#'   multiplicative factors for precipitation mean representing the
#'   time-averaged target change. See Details for transient vs step behavior.
#' @param change.factor.precip.variance Numeric vector of length 12. Monthly
#'   multiplicative factors for precipitation variance representing the
#'   time-averaged target change.
#' @param change.factor.temp.mean Numeric vector of length 12. Monthly additive
#'   temperature changes in degrees Celsius representing the time-averaged
#'   target change. See Details.
#' @param transient.temp.change Logical. If `TRUE`, temperature changes are
#'   applied gradually (linearly) over years such that the time-averaged change
#'   equals the target. If `FALSE`, full change applied uniformly to all years.
#'   Default: `TRUE`.
#' @param transient.precip.change Logical. If `TRUE`, precipitation changes are
#'   applied gradually over years such that the time-averaged change equals the
#'   target. If `FALSE`, full change applied uniformly to all years. Default: `TRUE`.
#' @param calculate.pet Logical. If `TRUE`, recalculates PET using Hargreaves
#'   method with adjusted temperatures. Default: `TRUE`.
#' @param parallel Logical. If `TRUE` and multiple cores available, processes
#'   grid cells in parallel. Default: `FALSE`.
#' @param num.cores Integer. Number of CPU cores for parallel processing. If
#'   `NULL`, uses all available cores minus one. Default: `NULL`.
#' @param fit.method Character. Method for gamma distribution fitting in quantile
#'   mapping. Passed to `quantile_mapping()`. Default: `"mme"`.
#' @param min.precip.factor Numeric. Minimum allowed precipitation change factor
#'   to prevent extreme reductions. Default: `0.01` (99 percent reduction limit).
#' @param validate.output Logical. If `TRUE`, checks for and handles non-finite
#'   values in output. Default: `TRUE`.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default: `FALSE`.
#'
#' @return
#' List of data frames (same structure as `climate.data`) with adjusted weather
#' variables. Each data frame includes perturbed `precip`, `temp`, `temp_min`,
#' `temp_max`, and (if `calculate.pet = TRUE`) updated `pet`.
#'
#' @details
#' ## Change Application Philosophy
#'
#' This function interprets change factors as time-averaged targets. This
#' ensures that when transient changes are applied, the mean change across the
#' entire simulation period matches the specified factor.
#'
#' ### Step Change Mode (transient = FALSE)
#' The specified change factor is applied uniformly to all years:
#'
#' All years: factor = target
#'
#' Time average: factor = target
#'
#' Example with `change.factor.precip.mean = 1.2` (20 percent increase):
#' \itemize{
#'   \item Every year: 1.2x (20 percent increase)
#'   \item Time average: 1.2x (20 percent increase)
#' }
#'
#' ### Transient Change Mode (transient = TRUE)
#' Changes ramp linearly such that the time-averaged change equals the target.
#' The function internally doubles the endpoint to achieve this:
#'
#' Year 1: factor = 1.0 (baseline)
#'
#' Year N: factor = 2 x (target - 1) + 1
#'
#' Time average: factor = target
#'
#' Example with `change.factor.precip.mean = 1.2` (20 percent increase):
#' \itemize{
#'   \item Year 1: 1.0x (0 percent change)
#'   \item Year N/2: 1.2x (20 percent change)
#'   \item Year N: 1.4x (40 percent change)
#'   \item Time average: 1.2x (20 percent increase)
#' }
#'
#' Important: In transient mode, the final year experiences approximately
#' double the net change (e.g., 40 percent for a 20 percent target). This
#' behavior ensures the time-averaged change matches your specified target.
#'
#' ### Temperature Changes
#' Temperature follows the same logic but with additive changes:
#'
#' Step mode with `change.factor.temp.mean = 2.0`:
#' \itemize{
#'   \item All years: +2.0 degrees C
#'   \item Time average: +2.0 degrees C
#' }
#'
#' Transient mode with `change.factor.temp.mean = 2.0`:
#' \itemize{
#'   \item Year 1: +0.0 degrees C
#'   \item Year N: +4.0 degrees C
#'   \item Time average: +2.0 degrees C
#' }
#'
#' ## Quantile Mapping
#' Precipitation is adjusted using parametric quantile mapping with gamma
#' distributions. This preserves distributional shape while shifting mean and
#' variance. See `quantile_mapping()` for details.
#'
#' ## PET Calculation
#' PET is recalculated using the Hargreaves temperature-based method, which
#' requires mean temperature, diurnal temperature range, and latitude.
#'
#' @examples
#' \dontrun{
#' # Example 1: Step change (uniform across all years)
#' # 20 percent precipitation increase, 2 degrees C warming applied to all years
#' perturbed_step <- impose_climate_changes(
#'   climate.data = climate.data,
#'   climate.grid = climate.grid,
#'   sim.dates = sim.dates,
#'   change.factor.precip.mean = rep(1.2, 12),      # 20 percent increase all years
#'   change.factor.precip.variance = rep(1.3, 12),
#'   change.factor.temp.mean = rep(2.0, 12),        # +2 degrees C all years
#'   transient.temp.change = FALSE,
#'   transient.precip.change = FALSE
#' )
#'
#' # Example 2: Transient change (gradual, time-averaged)
#' # Time-averaged 20 percent increase (final year will be approximately 40 percent higher)
#' # Time-averaged +2 degrees C (final year will be +4 degrees C)
#' perturbed_transient <- impose_climate_changes(
#'   climate.data = climate.data,
#'   climate.grid = climate.grid,
#'   sim.dates = sim.dates,
#'   change.factor.precip.mean = rep(1.2, 12),      # Avg 20 percent, final 40 percent
#'   change.factor.precip.variance = rep(1.3, 12),
#'   change.factor.temp.mean = rep(2.0, 12),        # Avg +2 degrees C, final +4 degrees C
#'   transient.temp.change = TRUE,
#'   transient.precip.change = TRUE
#' )
#'
#' # Example 3: Seasonal variation with transient changes
#' # Winter months get larger changes
#' winter_precip <- rep(1.2, 12)
#' winter_precip[c(12, 1, 2)] <- 1.4  # 40 percent time-averaged increase in winter
#'
#' perturbed_seasonal <- impose_climate_changes(
#'   climate.data = climate.data,
#'   climate.grid = climate.grid,
#'   sim.dates = sim.dates,
#'   change.factor.precip.mean = winter_precip,
#'   change.factor.precip.variance = rep(1.3, 12),
#'   change.factor.temp.mean = rep(2.5, 12),
#'   transient.precip.change = TRUE
#' )
#'
#' # Example 4: Verify time-averaged changes
#' original_mean <- mean(climate.data[[1]]$precip[climate.data[[1]]$precip > 0])
#' perturbed_mean <- mean(perturbed_transient[[1]]$precip[perturbed_transient[[1]]$precip > 0])
#' cat("Time-averaged change:", perturbed_mean / original_mean, "\n")
#' # Should be approximately 1.2 for transient mode
#' }
#'
#' @seealso
#' \code{\link{quantile_mapping}} for precipitation adjustment details
#' \code{\link{pet_hargreaves}} for PET calculation
#'
#' @export
impose_climate_changes <- function(
    climate.data = NULL,
    climate.grid = NULL,
    sim.dates = NULL,
    change.factor.precip.mean = NULL,
    change.factor.precip.variance = NULL,
    change.factor.temp.mean = NULL,
    transient.temp.change = TRUE,
    transient.precip.change = TRUE,
    calculate.pet = TRUE,
    parallel = FALSE,
    num.cores = NULL,
    fit.method = "mme",
    min.precip.factor = 0.01,
    validate.output = TRUE,
    verbose = FALSE) {

  # ==========================================================================
  # INPUT VALIDATION
  # ==========================================================================

  validate_climate_change_inputs(
    climate.data = climate.data,
    climate.grid = climate.grid,
    sim.dates = sim.dates,
    change.factor.precip.mean = change.factor.precip.mean,
    change.factor.precip.variance = change.factor.precip.variance,
    change.factor.temp.mean = change.factor.temp.mean
  )

  # ==========================================================================
  # SETUP
  # ==========================================================================

  ngrids <- length(climate.data)
  n.dates <- length(sim.dates)

  if (verbose) {
    message(sprintf("Processing %d grid cells with %d time steps", ngrids, n.dates))
    if (transient.precip.change) {
      message("Note: Using time-averaged transient mode for precipitation")
      message("      Final year will have ~2x the net change from baseline")
    }
    if (transient.temp.change) {
      message("Note: Using time-averaged transient mode for temperature")
      message("      Final year will have ~2x the specified change")
    }
  }

  # Compute date indices (once for all grids)
  date.indices <- compute_date_indices(sim.dates)

  # ==========================================================================
  # COMPUTE CHANGE FACTOR MATRICES
  # ==========================================================================

  if (verbose) {
    message("Computing change factor matrices...")
  }

  change.matrices <- compute_change_matrices_time_averaged(
    change.factor.precip.mean = change.factor.precip.mean,
    change.factor.precip.variance = change.factor.precip.variance,
    change.factor.temp.mean = change.factor.temp.mean,
    n.years = date.indices$n.years,
    transient.temp = transient.temp.change,
    transient.precip = transient.precip.change,
    min.precip.factor = min.precip.factor
  )

  # ==========================================================================
  # COMPUTE DAILY CHANGE FACTORS
  # ==========================================================================

  if (verbose) {
    message("Computing daily change factors...")
  }

  daily.factors <- compute_daily_factors(
    change.matrices = change.matrices,
    year.idx = date.indices$year.idx,
    month.idx = date.indices$month.idx
  )

  # ==========================================================================
  # APPLY PERTURBATIONS TO EACH GRID CELL
  # ==========================================================================

  if (verbose) {
    message("Applying climate perturbations to grid cells...")
  }

  if (parallel && ngrids > 1) {
    # Parallel processing
    climate.data.out <- apply_perturbations_parallel(
      climate.data = climate.data,
      climate.grid = climate.grid,
      change.matrices = change.matrices,
      daily.factors = daily.factors,
      date.indices = date.indices,
      calculate.pet = calculate.pet,
      fit.method = fit.method,
      num.cores = num.cores,
      verbose = verbose
    )
  } else {
    # Sequential processing
    climate.data.out <- apply_perturbations_sequential(
      climate.data = climate.data,
      climate.grid = climate.grid,
      change.matrices = change.matrices,
      daily.factors = daily.factors,
      date.indices = date.indices,
      calculate.pet = calculate.pet,
      fit.method = fit.method,
      verbose = verbose
    )
  }

  # ==========================================================================
  # VALIDATE OUTPUT
  # ==========================================================================

  if (validate.output) {
    climate.data.out <- validate_climate_output(
      climate.data.out = climate.data.out,
      verbose = verbose
    )
  }

  if (verbose) {
    message("Climate change perturbations complete")
  }

  climate.data.out
}


# ==============================================================================
# COMPUTATION FUNCTIONS (Updated for Time-Averaged Logic)
# ==============================================================================

#' Compute change factor matrices using time-averaged approach
#' @keywords internal
compute_change_matrices_time_averaged <- function(
    change.factor.precip.mean,
    change.factor.precip.variance,
    change.factor.temp.mean,
    n.years,
    transient.temp,
    transient.precip,
    min.precip.factor) {

  # Temperature change matrix (n.years x 12)
  if (transient.temp) {
    # Linear ramp from 0 to 2×target
    # This ensures time average = target
    temp.change.mat <- sapply(seq_len(12), function(m) {
      seq(0, change.factor.temp.mean[m] * 2, length.out = n.years)
    })
  } else {
    # Step change - full change in all years
    temp.change.mat <- matrix(
      rep(change.factor.temp.mean, each = n.years),
      nrow = n.years,
      ncol = 12
    )
  }

  # Precipitation mean change matrix (n.years x 12)
  if (transient.precip) {
    # Linear ramp from 1.0 to 2×(target-1)+1
    # This ensures time average = target
    precip.mean.mat <- sapply(seq_len(12), function(m) {
      target.endpoint <- (change.factor.precip.mean[m] - 1) * 2 + 1
      seq(1.0, target.endpoint, length.out = n.years)
    })

    # Apply minimum constraint
    precip.mean.mat[precip.mean.mat < min.precip.factor] <- min.precip.factor
  } else {
    # Step change
    precip.mean.mat <- matrix(
      rep(change.factor.precip.mean, each = n.years),
      nrow = n.years,
      ncol = 12
    )
  }

  # Precipitation variance change matrix (n.years x 12)
  if (transient.precip) {
    precip.var.mat <- sapply(seq_len(12), function(m) {
      target.endpoint <- (change.factor.precip.variance[m] - 1) * 2 + 1
      seq(1.0, target.endpoint, length.out = n.years)
    })

    precip.var.mat[precip.var.mat < min.precip.factor] <- min.precip.factor
  } else {
    precip.var.mat <- matrix(
      rep(change.factor.precip.variance, each = n.years),
      nrow = n.years,
      ncol = 12
    )
  }

  list(
    temp.change = temp.change.mat,
    precip.mean = precip.mean.mat,
    precip.var = precip.var.mat
  )
}


#' Verify Time-Averaged Change Factors
#'
#' @description
#' Helper function to verify that time-averaged changes match specified targets
#' when using transient mode. Useful for validating that the time-averaging
#' logic is working as intended.
#'
#' @param climate.data.original List of original climate data frames
#' @param climate.data.perturbed List of perturbed climate data frames
#' @param variables Character vector of variables to check (default: c("precip", "temp"))
#'
#' @return Data frame with columns:
#' \itemize{
#'   \item `variable`: Variable name
#'   \item `original.mean`: Mean of original data
#'   \item `perturbed.mean`: Mean of perturbed data
#'   \item `observed.ratio`: Actual change ratio (perturbed/original)
#' }
#'
#' @examples
#' \dontrun{
#' # Apply transient changes
#' perturbed <- impose_climate_changes(
#'   climate.data = climate.data,
#'   change.factor.precip.mean = rep(1.2, 12),
#'   transient.precip.change = TRUE
#' )
#'
#' # Verify time-averaged changes
#' verification <- verify_time_averaged_changes(
#'   climate.data.original = climate.data,
#'   climate.data.perturbed = perturbed
#' )
#' print(verification)
#' # Should show observed.ratio approximately equal to 1.2 for precipitation
#' }
#'
#' @export
verify_time_averaged_changes <- function(
    climate.data.original,
    climate.data.perturbed,
    variables = c("precip", "temp")) {

  # Function body remains the same...
}


#' Plot Transient Change Progression
#'
#' @description
#' Visualizes how change factors evolve over years in transient mode,
#' showing the relationship between yearly factors and time-averaged target.
#'
#' @param change.factor.target Numeric. Target time-averaged change factor
#'   (e.g., 1.2 for 20 percent increase)
#' @param n.years Integer. Number of simulation years
#' @param change.type Character. Either "multiplicative" (for precipitation)
#'   or "additive" (for temperature)
#'
#' @return A ggplot2 object showing the transient progression
#'
#' @examples
#' \dontrun{
#' # Visualize 20 percent precipitation increase over 30 years
#' plot_transient_progression(
#'   change.factor.target = 1.2,
#'   n.years = 30,
#'   change.type = "multiplicative"
#' )
#'
#' # Visualize +2 degrees C temperature increase over 50 years
#' plot_transient_progression(
#'   change.factor.target = 2.0,
#'   n.years = 50,
#'   change.type = "additive"
#' )
#' }
#'
#' @export
plot_transient_progression <- function(
    change.factor.target = 1.2,
    n.years = 30,
    change.type = c("multiplicative", "additive")) {

  change.type <- match.arg(change.type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }

  # Compute progression
  if (change.type == "multiplicative") {
    # Precipitation logic
    endpoint <- (change.factor.target - 1) * 2 + 1
    yearly.factors <- seq(1.0, endpoint, length.out = n.years)
    time.average <- mean(yearly.factors)

    ylabel <- "Change Factor (multiplicative)"
    title <- sprintf(
      "Transient Precipitation Change: Target = %.2fx (%.0f%% increase)",
      change.factor.target,
      (change.factor.target - 1) * 100
    )
  } else {
    # Temperature logic
    endpoint <- change.factor.target * 2
    yearly.factors <- seq(0, endpoint, length.out = n.years)
    time.average <- mean(yearly.factors)

    ylabel <- "Temperature Change (degrees C)"
    title <- sprintf(
      "Transient Temperature Change: Target = +%.1f degrees C",
      change.factor.target
    )
  }

  # Create data frame
  plot.data <- data.frame(
    year = 1:n.years,
    factor = yearly.factors
  )

  # Create plot
  p <- ggplot2::ggplot(plot.data, ggplot2::aes(x = year, y = factor)) +
    ggplot2::geom_line(linewidth = 1, color = "blue") +
    ggplot2::geom_point(size = 2, color = "blue") +
    ggplot2::geom_hline(
      yintercept = time.average,
      linetype = "dashed",
      color = "red",
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "text",
      x = n.years * 0.7,
      y = time.average * 1.05,
      label = sprintf("Time Average = %.2f", time.average),
      color = "red",
      hjust = 0
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf(
        "Year 1: %.2f | Year %d: %.2f | Time Avg: %.2f",
        yearly.factors[1], n.years, yearly.factors[n.years], time.average
      ),
      x = "Year",
      y = ylabel
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11)
    )

  print(p)
  invisible(p)
}
