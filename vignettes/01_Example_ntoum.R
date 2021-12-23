## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message=FALSE)

## ----climdata, eval=FALSE-----------------------------------------------------
#  library(gridwegen)
#  
#  # netcdf file with gridded daily weather data
#  nc_file <- system.file("extdata", "ntoum_era5_data.nc", package = "gridwegen")
#  
#  # names of the dimension variables in the netcdf file
#  nc_dimnames <- list(x = "lon", y = "lat", time = "time")
#  
#  # choose the variables to be used in the weather generator
#  variables <- c("precip", "temp", "temp_min", "temp_max")
#  
#  # define long names and units of the variables specified
#  variable_labels = c("Precipitation", "Mean Temperature", "Minimum Temperature", "Maximum Temperature")
#  variable_units = c("mm/day", "DegC", "DegC", "DegC")
#  
#  # Read-in gridded weather data from netcdf
#  nc_data <- readNetcdf(
#      nc.file = nc_file,
#      nc.dimnames = nc_dimnames,
#      nc.variables = variables,
#      origin.date = as.Date("1981-01-01"),
#      has.leap.days = TRUE)

## ----parameters, eval=FALSE---------------------------------------------------
#  # Climate change perturbations
#  precip_changes <- list()
#  precip_changes$increments <- 2
#  precip_changes$mean <- list()
#  precip_changes$var  <- list()
#  
#  temp_changes <- list()
#  temp_changes$increments <- 3
#  temp_changes$mean <- list()
#  
#  # Monthly precip changes   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
#  precip_changes$mean$min <- c(0.7, 0.7, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
#  precip_changes$mean$max <- c(1.3, 1.3, 1.4, 1.4, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3)
#  precip_changes$var$min <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
#  precip_changes$var$max <- c(1.3, 1.3, 1.5, 1.5, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3)
#  
#  # Monthly temp. changes   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec
#  temp_changes$mean$min <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#  temp_changes$mean$max <- c(3.0, 3.2, 3.4, 4.0, 4.1, 4.4, 5.0, 3.5, 3.3, 2.9, 2.8, 2.7)

## ----run, eval=FALSE, echo=TRUE-----------------------------------------------
#  simulateWeather(
#       climate.data = nc_data$data,
#       climate.grid = nc_data$coords,
#       year.start = 1981,
#       year.num = 40,
#       month.start = 1,
#       variable.names = variables,
#       variable.labels = variable_labels,
#       variable.units = variable_units,
#       sim.year.start = 2020,
#       sim.year.num = 40,
#       realization.num = 3,
#       warm.variable = "precip",
#       warm.signif.level = 0.90,
#       warm.sample.size = 20000,
#       knn.annual.sample.size = 100,
#       wet.state.threshold = 0.1,
#       extreme.state.quantile = 0.8,
#       apply.delta.changes = TRUE,
#       evaluate.model = TRUE,
#       evaluate.grid.num = 20,
#       apply.step.changes = TRUE,
#       delta.precip = precip_changes,
#       delta.temp = temp_changes,
#       output.path = out_path,
#       output.ncfile.template = nc_data,
#       output.ncfile.prefix = "clim_change_rlz",
#       bounds = list(
#          mean = c(0.95,1.05),
#          sd = c(0.90,1.10),
#          min = c(0.90,1.15),
#          max = c(0.85,1.10),
#          power = c(0.50,2.50),
#          nonsignif.threshold = 0.75)
#  )

