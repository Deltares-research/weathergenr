
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

# # Climate change perturbations
# precip_change <- list()
# precip_change$increments <- 3
# precip_change$mean <- list()
# precip_change$var  <- list()
# precip_change$mean$min <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
# precip_change$mean$max <- c(1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3)
# precip_change$var$min <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
# precip_change$var$max <- c(1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3)
#
#
#
#
#
#
#   PARCC$par1 <- list(
#     name = perturb1[[1,2]],
#     change_op <- perturb1[[2,2]],
#     increments = as.numeric(perturb1[[3,2]]),
#     mean = list(min = as.numeric(perturb1[7:18,1:3] %>% pull(2)),
#                 max = as.numeric(perturb1[7:18,1:3] %>% pull(3)),
#                 obs = as.numeric(perturb1[7:18,1:3] %>% pull(1))),
#     var  = list(min = as.numeric(perturb1[22:33,1:3] %>% pull(2)),
#                 max = as.numeric(perturb1[22:33,1:3] %>% pull(3)),
#                 obs = as.numeric(perturb1[22:33,1:3] %>% pull(1)))
#   )
#   PARCC$par2 <- list(
#     name = perturb2[[1,2]],
#     change_op <- perturb2[[2,2]],
#     increments = as.numeric(perturb2[[3,2]]),
#     mean = list(min = as.numeric(perturb2[7:18,1:3] %>% pull(2)),
#                 max = as.numeric(perturb2[7:18,1:3] %>% pull(3)),
#                 obs = as.numeric(perturb2[7:18,1:3] %>% pull(1))),
#     var  = list(min = as.numeric(perturb2[22:33,1:3] %>% pull(2)),
#                 max = as.numeric(perturb2[22:33,1:3] %>% pull(3)),
#                 obs = as.numeric(perturb2[22:33,1:3] %>% pull(1)))
#   )
#
#   PARCC$par1$sind <- 1:PARCC$par1$increments
#   PARCC$par2$sind <- 1:PARCC$par2$increments






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
      coordGrid = nc_data$coords,
      file.suffix = n,
      output.path = paste0(out_path,"future/"),
      sim.year.start = 2020,
      month.start = 1,
      variables = nc_variables,
      nc.dimnames = nc_dimnames,
      change.settings = paste0(nc_path, "/change_factors.xlsx"),
      save.scenario.matrix = FALSE,
      step = TRUE
  )

}

