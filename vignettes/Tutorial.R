## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message=FALSE)

## ----setup, eval=FALSE--------------------------------------------------------
#  devtools::install_github("tanerumit/gridwegen")

## ----ncfile-------------------------------------------------------------------
library(gridwegen)
ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "gridwegen")
ncdata <- readNetcdf(ncfile)

#Objects stored in the output data
names(ncdata)

## ----ncdata-------------------------------------------------------------------
ncdata$data[[1]]

## ----ncgrids------------------------------------------------------------------
ncdata$grid

## ----ncdates------------------------------------------------------------------
head(ncdata$dates)

## ----stochastic1--------------------------------------------------------------
# Set path to store weather generator results 
output_path <- "C:/testrun/"
variables <- c("precip", "temp", "temp_min", "temp_max")

## ----stochastic2, results='hide', eval = FALSE--------------------------------
#  stochastic_weather <- generateWeatherSeries(
#        weather.data = ncdata$data,
#        weather.grid = ncdata$grid,
#        weather.date = ncdata$date,
#        variable.names = variables,
#        output.path = output_path,
#        month.start = 1,
#        realization.num = 3,
#        warm.variable = "precip",
#        warm.signif.level = 0.90,
#        warm.sample.num = 5000,
#        knn.sample.num = 100,
#        evaluate.model = FALSE,
#        evaluate.grid.num = 20,
#        mc.wet.threshold = 0.2,
#        mc.extreme.quantile = 0.7
#      )

## ----stochastic3, eval = FALSE------------------------------------------------
#  # Resampled dates
#  stochastic_weather$resampled
#  
#  # Date vector
#  head(stochastic_weather$dates)

## ----climchange, eval = FALSE-------------------------------------------------
#  
#  # Temperature changes Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
#  delta_temp_mean <- c(3.0, 3.2, 3.4, 4.0, 4.1, 4.4, 5.0, 3.5, 3.3, 2.9, 2.8, 2.7)
#  
#  # Precipitation changes   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
#  delta_precip_mean     <- c(0.7, 0.7, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
#  delta_precip_variance <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
#  
#  # Select first realization
#  day_order <- match(stochastic_weather$resampled[[1]], ncdata$date)
#  
#  # Obtain stochastic series by re-ordering historical data
#  stochastic_rlz <- lapply(ncdata$data, function(x) x[day_order,])
#  
#  # Apply climate changes to climate data
#  stochastic2 <- imposeClimateChanges(
#      climate.data = stochastic_rlz,
#      climate.grid = ncdata$grid,
#      sim.dates = stochastic_weather$dates,
#      change.factor.precip.mean = delta_precip_mean,
#      change.factor.precip.variance = delta_precip_variance,
#      change.factor.temp.mean = delta_temp_mean,
#      change.type.temp = "transient",
#      change.type.precip = "transient")

## ----climchange2, eval = FALSE------------------------------------------------
#  # Save to netcdf file
#  writeNetcdf(
#      data = stochastic2,
#      coord.grid = ncdata$grid,
#      output.path = output_path,
#      origin.date =  stochastic_weather$dates[1],
#      calendar.type = "noleap",
#      nc.template.file = ncfile,
#      nc.compression = 4,
#      nc.spatial.ref = "spatial_ref",
#      nc.file.prefix = "clim",
#      nc.file.suffix = NULL)

## ----stresstest, eval = FALSE-------------------------------------------------
#  # Bandwith for precipitation mean change
#  precip_mean_cf_min <- rep(0.7,12)
#  precip_mean_cf_max <- rep(1.3,12)
#  precip_variance_cf_min <- rep(1,12)
#  precip_variance_cf_max <- rep(1,12)
#  precip_step_num <- 7
#  
#  # Bandwith for temperature mean change
#  temp_cf_mean_min <- rep(0,12)
#  temp_cf_mean_max <- rep(5,12)
#  temp_cf_step_num <- 6
#  
#  precip_mean_cf_steps <- sapply(1:12, function(m)
#            seq(precip_mean_cf_min[m], precip_mean_cf_max[m],
#                length.out = precip_step_num))
#  
#  precip_variance_cf_steps <- sapply(1:12, function(m)
#            seq(precip_variance_cf_min[m], precip_variance_cf_max[m],
#                length.out = precip_step_num))
#  
#  temp_cf_mean_steps <- sapply(1:12, function(m)
#            seq(temp_cf_mean_min[m], temp_cf_mean_max[m],
#                length.out = temp_cf_step_num))
#  
#  strtest_mat <- tidyr::expand_grid(deltap = 1:precip_step_num,
#            deltat = 1:temp_cf_step_num) %>%
#          mutate(index = 1:n(), .before = 1)
#  
#  delta_precip_mean_mat <- precip_mean_cf_steps[strtest_mat$deltap, ]
#  delta_precip_variance_mat <- precip_variance_cf_steps[strtest_mat$deltap, ]
#  temp_cf_mean_mat <- temp_cf_mean_steps[strtest_mat$deltat, ]
#  
#  smax <- nrow(strtest_mat)
#  
#  # Output folder
#  output_path <- paste0(output_path,"future/")
#  if (!dir.exists(output_path)) dir.create(output_path)
#  
#  # Loop through each stochastic realization
#  for (n in 1:stochastic_num) {
#  
#    # Resample order
#    day_order <- match(stochastic$resampled[[n]], nc_clim$dates)
#  
#    # Obtain stochastic series by re-ordering historical data
#    stochastic_rlz <- lapply(nc_clim$data, function(x) x[day_order,])
#  
#    for (s in 1:smax) {
#  
#      # Apply climate changes to climate data
#      out <- imposeClimateChanges(
#        climate.data = stochastic_rlz,
#        climate.grid = nc_clim$coords,
#        sim.dates = stochastic$dates,
#        change.factor.precip.mean = delta_precip_mean_mat[s,],
#        change.factor.precip.variance = delta_precip_variance_mat[s,],
#        change.factor.temp.mean = temp_cf_mean_mat[s,],
#        change.type.temp = "transient",
#        change.type.precip = "transient")
#  
#      # Save to netcdf file
#      writeNetcdf(
#        data = out,
#        coord.grid = nc_clim$coords,
#        output.path = paste0(output_path,"future/"),
#        origin.date =  stochastic$dates[1],
#        calendar.type = "noleap",
#        nc.template.file = nc_file,
#        nc.compression = 4,
#        nc.spatial.ref = "spatial_ref",
#        nc.file.prefix = "clim",
#        nc.file.suffix = paste0(n,"_", s))
#    }
#  }

