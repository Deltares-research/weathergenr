## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message=FALSE, cache = TRUE)

## ----stochastic2, results='hide', eval = TRUE, cache=TRUE---------------------
stochastic_weather <- generateWeatherSeries(
  weather.data = ncdata$data,
  weather.grid = ncdata$grid,
  weather.date = ncdata$date,
  variable.names = variables,
  variable.labels = variables,
  variable.units = NULL,
  sim.year.num = 20,
  sim.year.start = 2020,
  month.start = 1,
  realization.num = realization_num,
  warm.variable = "precip",
  warm.signif.level = 0.90,
  warm.sample.num = 10000,
  warm.subset.criteria = NULL,
  knn.sample.num = 120,
  mc.wet.quantile= 0.2,
  mc.extreme.quantile = 0.8,
  output.path = output_path,
  seed = 123)

## ----stochastic3, eval = TRUE, cache=TRUE-------------------------------------
# Resampled dates
stochastic_weather$resampled

# Date vector
head(stochastic_weather$dates)
tail(stochastic_weather$dates)

## ----climchange, eval = FALSE-------------------------------------------------
#  # Temperature changes Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
#  delta_temp_mean <- c(3.0, 3.2, 3.4, 4.0, 4.1, 4.4, 5.0, 3.5, 3.3, 2.9, 2.8, 2.7)
#  
#  # Precipitation changes   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
#  delta_precip_mean     <- c(0.7, 0.7, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
#  delta_precip_variance <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
#  
#  # Obtain stochastic series by re-ordering historical data
#  day_order <- match(stochastic_weather$resampled[[1]], ncdata$date)
#  stochastic_rlz <- lapply(ncdata$data, function(x) x[day_order,])
#  
#  # Apply climate changes to climate data
#  stochastic2 <- imposeClimateChanges(
#    climate.data = stochastic_rlz,
#    climate.grid = ncdata$grid,
#    sim.dates = stochastic_weather$dates,
#    change.factor.precip.mean = delta_precip_mean,
#    change.factor.precip.variance = delta_precip_variance,
#    change.factor.temp.mean = delta_temp_mean,
#    transient.temp.change = TRUE,
#    transient.precip.change = TRUE,
#    calculate.pet = TRUE,
#    compute.parallel = TRUE,
#    num.cores = NULL,
#    fit.method = "mme")

## ----climchange2, eval = FALSE------------------------------------------------
#  # Save to netcdf file
#  writeNetcdf(
#    data = stochastic2,
#    coord.grid = ncdata$grid,
#    output.path = output_path,
#    origin.date =  stochastic_weather$dates[1],
#    calendar.type = "noleap",
#    nc.template.file = ncfile,
#    nc.compression = 4,
#    nc.spatial.ref = "spatial_ref",
#    nc.file.prefix = "clim",
#    nc.file.suffix = NULL)

