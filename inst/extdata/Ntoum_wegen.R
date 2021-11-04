
# Packages needed
library(gridwegen)

# Path to output files
path0 <- "C:/wegentest/"
out_path <- paste0(path0, "TEST2/")

# Read-in gridded weather data from netcdf
nc_path <- system.file('extdata', package = 'gridwegen')
nc_file <- "ntoum.nc"
nc_dimnames <- list(x = "lon", y = "lat", time = "time")
nc_variables <- c("precip", "temp", "temp_min", "temp_max")
origin_date <- as.Date("1981-01-01")
sim_year_start <- 2020
realization_num <- 2

# Climate change perturbations
precip_changes <- list()
precip_changes$increments <- 2
precip_changes$mean <- list()
precip_changes$var  <- list()
precip_changes$mean$min <- c(0.7, 0.7, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
precip_changes$mean$max <- c(1.3, 1.3, 1.4, 1.4, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3)
precip_changes$var$min <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
precip_changes$var$max <- c(1.3, 1.3, 1.5, 1.5, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3)

temp_changes <- list()
temp_changes$increments <- 3
temp_changes$mean <- list()
temp_changes$mean$min <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
temp_changes$mean$max <- c(3.0, 3.2, 3.4, 4.0, 4.1, 4.4, 5.0, 3.5, 3.3, 2.9, 2.8, 2.7)


nc_data <- readNetcdf(
    nc.path = nc_path,
    nc.file = nc_file,
    nc.dimnames = nc_dimnames,
    nc.variables = nc_variables,
    origin.date = origin_date,
    has.leap.days = TRUE)

# Simulate daily realizations of historical weather conditions
out <- simulateWeather(
  output.path = out_path,
  hist.climate = nc_data$data,
  grid.coords = nc_data$coords,
  hist.date.start = origin_date,
  variables = names(nc_data$variables)[1:4],
  sim.year.num = 40,
  sim.year.start = sim_year_start,
  warm.sample.num = 5000,
  realization.num = realization_num,
  return.date.indices = FALSE)


for(n in 1:realization_num) {

  imposeClimateChanges(
      input.data = out[[n]],
      grid.coords = nc_data$coords,
      precip.changes <- precip_changes,
      temp.changes = temp_changes,
      file.suffix = n,
      output.path = out_path,
      sim.year.start = 2020,
      month.start = 1,
      variables = nc_variables,
      nc.dimnames = nc_dimnames,
      save.scenario.matrix = FALSE,
      step.change = TRUE)
}

