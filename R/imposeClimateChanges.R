
#' Perturb weather realizations
#'
#' @param climate.data placeholder
#' @param climate.grid placeholder
#' @param sim.dates placeholder
#' @param change.factor.precip.mean placeholder
#' @param change.factor.precip.variance placeholder
#' @param change.factor.temp.mean placeholder
#' @param transient.temp.change placeholder
#' @param transient.precip.change placeholder
#' @param calculate.pet placeholder
#'
#' @return
#' @export
#' @import tibble
#' @import tidyr
imposeClimateChanges <- function(
  climate.data = NULL,
  climate.grid = NULL,
  sim.dates = NULL,
  change.factor.precip.mean = NULL,
  change.factor.precip.variance = NULL,
  change.factor.temp.mean = NULL,
  transient.temp.change = TRUE,
  transient.precip.change = TRUE,
  calculate.pet = TRUE
)

{

  ngrids <- length(climate.data)

  year_vec <- as.numeric(format(sim.dates,"%Y"))
  year_ind <- year_vec - min(year_vec) + 1
  month_ind <- as.numeric(format(sim.dates,"%m"))

  # Temperature change factors #################################################

  if(isTRUE(transient.temp.change)) {

      tempf1 <- sapply(1:12, function(x)
          seq(0, change.factor.temp.mean[x], length.out = max(year_ind)))

  } else {

      tempf1 <- sapply(1:12, function(x)
          rep(change.factor.temp.mean[x], length.out = max(year_ind)))
  }

  tempf2 <- sapply(1:length(sim.dates), function(x) tempf1[year_ind[x], month_ind[x]])

  # Precipitation change factors ###############################################

  if(isTRUE(transient.precip.change)) {

    precip_meanf <- sapply(1:12,
      function(m) seq(1, change.factor.precip.mean[m], length.out = max(year_ind)))
    precip_varf  <- sapply(1:12,
      function(m) seq(1, change.factor.precip.variance[m], length.out = max(year_ind)))

  } else {
    precip_meanf <- sapply(1:12,
      function(m) rep(change.factor.precip.mean[m], max(year_ind)))
    precip_varf  <- sapply(1:12,
      function(m) rep(change.factor.precip.variance[m], max(year_ind)))
  }

  # Apply climate changes (per grid)
  for (x in 1:ngrids) {

    # Perturb daily precipitation using quantile mapping
    climate.data[[x]]$precip <- quantileMapping(
            value = climate.data[[x]]$precip,
            mon.ts = month_ind,
            year.ts = year_ind,
            mean.change = precip_meanf,
            var.change = precip_varf)

    # Perturb temp, temp_min, and temp_max by delta factors
    climate.data[[x]]$temp <- climate.data[[x]]$temp + tempf2
    climate.data[[x]]$temp_min <- climate.data[[x]]$temp_min + tempf2
    climate.data[[x]]$temp_max <- climate.data[[x]]$temp_max + tempf2

    if(isTRUE(calculate.pet)) {
      climate.data[[x]]$pet <- with(climate.data[[x]], hargreavesPet(
          months = month_ind, temp = temp, tdiff = temp_max - temp_min,
          lat = climate.grid$y[x]))
    }

  }

  # Replace possible infinite/NA values with zero
  climate.data <- lapply(1:length(climate.data), function(y)
    do.call(tibble, lapply(climate.data[[y]], function(x) replace(x, is.infinite(x), 0))))

  return(climate.data)
}
