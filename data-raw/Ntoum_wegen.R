

# Load package
#devtools::install_github("tanerumit/gridwegen@dev")

library(gridwegen)
library(tidyr)
library(dplyr)

# Path to input files and results outputted
path0 <- "C:/Users/taner/OneDrive - Stichting Deltares/_DELTARES/02 Projects/11206634 Gabon/05 Models/wegen/"
nc_path <- paste0(path0, "input/")
nc_file <- "localP_chirpsP_era5T_19810101_chunks.nc"
out_path <- paste0(path0, "TEST/")

# GABON/ ERA5-GRIDDED CLIMATE DATA
nc_data <- readNetcdf(
    in.path = nc_path,
    in.file = nc_file,
    dim.names = list(x = "lon", y = "lat", time = "time"),
    variables = c("precip", "temp", "temp_min", "temp_max"),
    origin.date = as.Date("1981-01-01"),
    leap.year = TRUE)

#Remove extra grids from tidy data table
grid_select <- which(sapply(1:nrow(nc_data$tidy_data), function(x)
  !is.na(mean(nc_data$tidy_data$data[[x]]$temp_min))))[1:20]
climate_tidy <- nc_data$tidy_data[grid_select,] %>% mutate(id = 1:n())

simulateWeather(
  proj.name = "ntoum_test",
  output.dir = out_path,
  climate_tidy = climate_tidy,
  wg.date.begin = as.Date("1981-01-01"),
  wg.vars = c("precip", "temp", "temp_min", "temp_max"),
  wg.var.labs = c("Precipitation", "Avg. Temperature", "Min. Temperature","Max. Temperature"),
  wg.var.units = c("mm/day", "°C", "°C", "°C", "mm/day"),
  warm.var = "precip",
  ymax = 40,
  rmax = 5000,
  nmax = 5,
  nc.dimnames = list(x = "lon", y = "lat", time = "time"),
  validate = TRUE,
  mean.bounds = NULL,
  sdev.bounds = NULL,
  max.bounds  = NULL,
  min.bounds  = NULL,
  power.bounds = NULL,
  nonsig.threshold = NULL
)

# natural variability realizations in the input folder
nvar_filenames <- list.files(paste0(out_path,"historical/"))
nmax <- length(nvar_filenames)

for(n in 2:nmax) {

  imposeClimateChanges(
      proj.name = "ntoum",
      in.path = paste0(out_path,"historical/"),
      in.file = nvar_filenames[n],
      file.suffix = n,
      out.path = paste0(out_path,"future/"),
      sim.date.begin = as.Date("2020-01-01"),
      wg.vars = c("precip", "temp", "temp_min", "temp_max"),
      wg.var.units = c("mm/day", "°C", "°C", "°C", "mm/day"),
      nc.dimnames = list(x = "lon", y = "lat", time = "time"),
      change.settings = paste0(nc_path, "change_factors.xlsx"),
      save.scenario.matrix = FALSE
  )

}








# library(fitdistrplus)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(tibble)
# library(lubridate)
# library(ncdf4)
# library(readxl)
# library(patchwork)
