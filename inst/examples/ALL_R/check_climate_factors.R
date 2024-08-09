
library(readr)
library(weathergenr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Read in extracted ERA5 climate data
fpath <- "C:\\Users\\taner\\Workspace\\test-weathergenr\\data\\Messara2024\\"
ncdata_hist <- paste0(fpath,"extract_historical.nc")
ncdata <- readNetcdf(ncdata_hist)

# Read in resampled date indices
sim_dates <- read_csv(paste0(fpath,"sim_dates.csv"))$x
resampled_dates_ini <- read_csv(paste0(fpath,"resampled_dates.csv"))
resampled_dates <- resampled_dates_ini$rlz_1

# Resample order
day_order <- match(resampled_dates, ncdata$date)

# Obtain stochastic series by re-ordering historical data
stochastic_rlz <- lapply(ncdata$data, function(x) x[day_order,])

n <- 1

# save to NetCDF
weathergenr::writeNetcdf(
  data = stochastic_rlz,
  coord.grid = ncdata$grid,
  output.path = fpath,
  origin.date =  sim_dates[1],
  calendar.type = "noleap",
  nc.template.file = ncdata_hist,
  nc.compression = 4,
  nc.spatial.ref = "spatial_ref",
  nc.file.prefix = "rlz_1",
  nc.file.suffix = NULL
)

# Stochastic weather realization to be perturbed
ffile <- paste0(fpath,"rlz_1_.nc")
rlz_input <- weathergenr::readNetcdf(ffile, leap.days = FALSE)

prcp_mean_changes <- rep(1, 12)
prcp_variance_changes <- rep(1, 12)
temp_mean_changes <- rep(0, 12)

temp_change_transient <- TRUE
precip_change_transient <- TRUE

##### TEST CLIMATE CHANGES
rlz_future <- weathergenr::imposeClimateChanges(
  climate.data = rlz_input$data,
  climate.grid = rlz_input$grid,
  sim.dates = rlz_input$date,
  change.factor.precip.mean = prcp_mean_changes ,
  change.factor.precip.variance = prcp_variance_changes,
  change.factor.temp.mean = temp_mean_changes,
  transient.temp.change = temp_change_transient,
  transient.precip.change = precip_change_transient,
  calculate.pet = TRUE,
  compute.parallel = FALSE,
  num.cores = NULL,
  fit.method = "mme"
)

# # Save to netcdf file
# weathergenr::writeNetcdf(
#   data = rlz_future,
#   coord.grid = rlz_input$grid,
#   output.path = fpath,
#   origin.date =  rlz_input$date[1],
#   calendar.type = "noleap",
#   nc.template.file = ncdata_hist,
#   nc.compression = 4,
#   nc.spatial.ref = "spatial_ref",
#   nc.file.prefix = "rlz_1z",
#   nc.file.suffix = "future"
# )

################################################################################
################################ CHECKS ######################################## 

rlz_hist <- weathergenr::readNetcdf(paste0(fpath,"rlz_1_.nc"), leap.days = FALSE)
#rlz_future <- weathergenr::readNetcdf(paste0(fpath,"rlz_1z_future.nc"), leap.days = FALSE)

grid_num <- 1
temp_hist <- rlz_hist$data[[grid_num]]$temp
temp_future <- rlz_future[[grid_num]]$temp
#temp_future <- rlz_future$data[[grid_num]]$temp
print(c(mean(temp_hist), mean(temp_future)))

precip_hist <- rlz_hist$data[[grid_num]]$precip
precip_future <- rlz_future[[grid_num]]$precip
#precip_future <- rlz_future$data[[grid_num]]$precip
print(c(mean(precip_hist), mean(precip_future)))

dfx <- tibble (hist = precip_hist, future = precip_future) %>%
  mutate(change = future/hist)
View(dfx)
