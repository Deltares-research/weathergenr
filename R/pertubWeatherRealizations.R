
imposeClimateChanges <- function(
  data = NULL,
  climate.grid = NULL,
  sim.dates = NULL,
  change.factor.precip.mean = NULL,
  change.factor.precip.variance = NULL,
  change.factor.temp.mean = NULL,
  change.type.temp = "transient",
  change.type.precip = "transient",
  calculate.pet = TRUE
)

{

  ngrids <- length(data)

  year_vec <- as.numeric(format(sim.dates,"%Y"))
  year_ind <- year_vec - min(year_vec) + 1
  month_ind <- as.numeric(format(sim.dates,"%m"))

  # Temperature change factors
  if (change.type.temp == "transient") {
      tempf1 <- sapply(1:12, function(x)
          seq(0, change.factor.temp.mean[x], length.out = max(year_ind)))
  } else {
      tempf1 <- sapply(1:12, function(x)
          rep(change.factor.temp.mean[x], length.out = max(year_ind)))
  }
  tempf2 <- sapply(1:length(sim.dates), function(x) tempf1[year_ind[x], month_ind[x]])

  # Precipitation change factors
  if(isTRUE(step.change)) {
    mean_a <- sapply(1:12, function(m) seq(1, mean.change[m], length.out = ymax))
    var_a  <- sapply(1:12, function(m) seq(1, var.change[m], length.out = ymax))
  } else {
    mean_a <- sapply(1:12, function(m) rep(mean.change[m], ymax))
    var_a  <- sapply(1:12, function(m) rep(var.change[m], ymax))
  }









  # Precipitation change factors

  # Vector of means and variances (row=years, column=months)
  if(isTRUE(step.change)) {
    mean_a <- sapply(1:12, function(m) seq(1, mean.change[m], length.out = ymax))
    var_a  <- sapply(1:12, function(m) seq(1, var.change[m], length.out = ymax))
  } else {
    mean_a <- sapply(1:12, function(m) rep(mean.change[m], ymax))
    var_a  <- sapply(1:12, function(m) rep(var.change[m], ymax))
  }





  # Apply climate changes (per grid)
  for (x in 1:ngrids) {

    # Perturb daily precipitation using quantile mapping
    data[[x]]$precip <- quantileMapping(
            value = data[[x]]$precip,
            mon.ts = month_ind,
            year.ts = year_ind,
            mean.change = change.factor.precip.mean,
            var.change = change.factor.precip.variance)

    # Perturb temp, temp_min, and temp_max by delta factors
    data[[x]]$temp <- data[[x]]$temp + tempf2
    data[[x]]$temp_min <- data[[x]]$temp_min + tempf2
    data[[x]]$temp_max <- data[[x]]$temp_max + tempf2

    if(isTRUE(calculate.pet)) {
      data[[x]]$pet <- with(data[[x]], hargreavesPet(
          months = month_ind, temp = temp, tdiff = temp_max - temp_min,
          lat = climate.grid$y[x]))
    }

  }

  return(data)
}








# # Write to netcdf
# writeNetcdf(
#     data = rlz_cur,
#     coord.grid = climate.grid,
#     output.path = future_path,
#     origin.date =  sim_dates_d$date[1],
#     calendar.type = "noleap",
#     nc.template.file = nc.template,
#     nc.compression = 4,
#     nc.spatial.ref = "spatial_ref",
#     nc.file.prefix = nc.prefix,
#     nc.file.suffix = NULL
# )
#
#
