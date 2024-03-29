#' Perturb weather realizations to reflect climate change
#'
#' @param climate.data placeholder
#' @param climate.grid placeholder
#' @param sim.dates placeholder
#' @param change.factor.precip.mean placeholder
#' @param change.factor.precip.variance placeholder
#' @param change.factor.temp.mean placeholder
#' @param transient.temp.change logical value indicating whether temperature changes are applied gradually (transient) or stepwise
#' @param transient.precip.change logical value indicating whether precipitation changes are applied gradually (transient) or stepwise
#' @param calculate.pet shoud pet be calculated? (logical)
#' @param compute.parallel should the function run in parallel mode (logical)
#' @param num.cores number of cores reserved for parallel computing. Default value is maximum cores minus one (omit if compute.parallel is FALSE)
#' @param fit.method placeholder
#'
#' @return
#' @export
#' @import dplyr
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
  calculate.pet = TRUE,
  compute.parallel = FALSE,
  num.cores = NULL,
  fit.method = "mme")

 {

  # Number of grid cells
  ngrids <- length(climate.data)

  # Year and month indices
  year_vec <- as.numeric(format(sim.dates,"%Y"))
  year_ind <- year_vec - min(year_vec) + 1
  month_ind <- as.numeric(format(sim.dates,"%m"))


  ### Define Daily temperature change factors

  # Transient changes (gradual)
  if(isTRUE(transient.temp.change)) {
      tempf1 <- sapply(1:12, function(x)
          seq(0, change.factor.temp.mean[x] * 2, length.out = max(year_ind)))

      # Step changes
      } else {
      tempf1 <- sapply(1:12, function(x)
          rep(change.factor.temp.mean[x], length.out = max(year_ind)))
  }

  tempf2 <- sapply(1:length(sim.dates),
                   function(x) tempf1[year_ind[x], month_ind[x]])

  ### Define Daily precipitation change factors

  # Transient changes (gradual)
  if(isTRUE(transient.precip.change)) {

    # Set constraint on maximum change (i.e., maximum = 99%)
    min_precipf = 0.01

    precip_mean_deltaf <- (change.factor.precip.mean - 1) * 2 + 1
    precip_meanf <- sapply(1:12,
      function(m) seq(1, precip_mean_deltaf[m], length.out = max(year_ind)))
    precip_meanf[precip_meanf < min_precipf] <- min_precipf

    precip_var_deltaf <- (change.factor.precip.variance - 1) * 2 + 1
    precip_varf  <- sapply(1:12,
      function(m) seq(1, precip_var_deltaf[m], length.out = max(year_ind)))
    precip_varf[precip_varf < min_precipf] <- min_precipf

    # Step changes
    } else {
    precip_meanf <- sapply(1:12,
      function(m) rep(change.factor.precip.mean[m], max(year_ind)))
    precip_varf  <- sapply(1:12,
      function(m) rep(change.factor.precip.variance[m], max(year_ind)))
  }

  ##############################################################################

  # Apply climate changes (per grid)
  for (x in 1:ngrids) {

    # Perturb daily precipitation (by quantile mapping)
    climate.data[[x]]$precip <- weathergenr::quantileMapping(
            value = climate.data[[x]]$precip,
            mon.ts = month_ind,
            year.ts = year_ind,
            mean.change = precip_meanf,
            var.change = precip_varf,
            fit.method = "mme")

    # Perturb temp, temp_min, and temp_max (by additive delta factors)
    climate.data[[x]]$temp <- climate.data[[x]]$temp + tempf2
    climate.data[[x]]$temp_min <- climate.data[[x]]$temp_min + tempf2
    climate.data[[x]]$temp_max <- climate.data[[x]]$temp_max + tempf2

    # Calculate adjusted PET based on temperature variables
    if(isTRUE(calculate.pet)) {
      climate.data[[x]]$pet <- with(climate.data[[x]],
                hargreavesPet(months = month_ind,
                              temp = temp,
                              tdiff = temp_max - temp_min,
                              lat = climate.grid$y[x]))
    }

  }

  # Replace possible infinite/NA values with zeros
  climate.data <- lapply(1:length(climate.data), function(y)
    do.call(tibble, lapply(climate.data[[y]],
                           function(x) replace(x, is.infinite(x), 0))))

  return(climate.data)

}

# # Set number of cores for parallel computing
# if(compute.parallel == TRUE) {
#
#   if(is.null(num.cores)) num.cores <- parallel::detectCores()-1
#   cl <- parallel::makeCluster(num.cores)
#   doParallel::registerDoParallel(cl)
#   `%d%` <- foreach::`%dopar%`
#
# } else {
#
#   `%d%` <- foreach::`%do%`
# }

##############################################################################
##############################################################################

# precip <- foreach::foreach(x=seq_len(ngrids)) %d% {
#
#     weathergenr::quantileMapping(
#           value = climate.data[[x]]$precip,
#           mon.ts = month_ind,
#           year.ts = year_ind,
#           mean.change = precip_meanf,
#           var.change = precip_varf,
#           fit.method = fit.method)
# }
#
# if(compute.parallel == TRUE) parallel::stopCluster(cl)
#
# for (x in 1:ngrids) {
#
#   # Perturb temp, temp_min, and temp_max by delta factors
#   climate.data[[x]]$precip <- precip[[x]]
#   climate.data[[x]]$temp   <- climate.data[[x]]$temp + tempf2
#   climate.data[[x]]$temp_min <- climate.data[[x]]$temp_min + tempf2
#   climate.data[[x]]$temp_max <- climate.data[[x]]$temp_max + tempf2
#
#   if(isTRUE(calculate.pet)) {
#     climate.data[[x]]$pet <- with(climate.data[[x]], hargreavesPet(
#         months = month_ind, temp = temp, tdiff = temp_max - temp_min,
#         lat = climate.grid$y[x]))
#   }
#
# }
#
# # Replace possible infinite/NA values with zero
# climate.data <- lapply(1:length(climate.data), function(y)
#       do.call(tibble::tibble, lapply(climate.data[[y]],
#         function(x) replace(x, is.infinite(x) | is.na(x), 0))))
