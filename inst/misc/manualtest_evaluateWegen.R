
library(dplyr)
library(ggplot2)
library(weathergenr)
library(tidyr)

# Path to files
output_path <- "./inst/extdata/"
ncfile <- paste0(output_path, "ntoum_era5_data.nc")
resampled_dates <- read.csv(paste0(output_path, "ntoum_era5_resampled_5traces.csv"), colClasses = "Date")
simulated_dates <- read.csv(paste0(output_path, "ntoum_era5_simdates_5traces.csv"), colClasses = "Date")

# Read-in netcdf file & set other parameters
ncdata <- readNetcdf(ncfile)
realization.num <- ncol(resampled_dates)
evaluate.grid.num <- 20
grids <- ncdata$grid$id
sampleGrids <- sample(grids, size = min(evaluate.grid.num, length(grids)))

# Prepare dail and synthetic weather series
daily_obs <- lapply(ncdata$data[sampleGrids], function(x)
  dplyr::mutate(x, date = ncdata$date, .before = 1))
daily_sim <- vector(mode = "list", length = realization.num)
for (n in 1:realization.num) {

  # Resample order
  day_order <- match(resampled_dates[[n]], ncdata$date)

  # Obtain stochastic series by re-ordering historical data
  daily_sim[[n]] <- lapply(ncdata$data[sampleGrids], function(x)
    x[day_order,] %>% mutate(date = simulated_dates[[1]], .before = 1))
}

# Run evaluateWegen script
results <- evaluateWegen(daily.sim = daily_sim,
              daily.obs = daily_obs,
              output.path = paste0(output_path, "plots"),
              variables = c("precip", "temp", "temp_min", "temp_max"),
              variable.labels = c("precip", "temp", "temp_min", "temp_max"),
              variable.units = NULL,
              realization.num = realization.num,
              wet.quantile = 0.2,
              extreme.quantile = 0.8,
              show.title = TRUE,
              save.plots = TRUE)




################################################################################
# FOR DEBUGGING ONLY

daily.sim = daily_sim
daily.obs = daily_obs
output.path = paste0(output_path, "plots")
variables = c("precip", "temp", "temp_min", "temp_max")
variable.labels = c("precip", "temp", "temp_min", "temp_max")
variable.units = NULL
realization.num = realization.num
wet.quantile = 0.2
extreme.quantile = 0.8
show.title = TRUE
save.plots = TRUE




