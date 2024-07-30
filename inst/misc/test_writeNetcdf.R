

library(readr)
library(weathergenr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Read in extracted ERA5 climate data
ncdata_hist <- paste0("./inst/extdata/ntoum_era5_data.nc")
ncdata <- readNetcdf(ncdata_hist)

# Read in resampled date indices
sim_dates <- read_csv("./inst/extdata/sim_dates.csv")[[1]]
resampled_dates <- read_csv("./inst/extdata/resampled_dates.csv")[[1]]

# Resample order
day_order <- match(resampled_dates, ncdata$date)

# Obtain stochastic series by re-ordering historical data
stochastic_rlz <- lapply(ncdata$data, function(x) x[day_order,])

# save to netcdf
weathergenr::writeNetcdf(
  data = stochastic_rlz,
  coord.grid = ncdata$grid,
  output.path = "./temp/",
  origin.date =  sim_dates[1],
  calendar.type = "noleap",
  nc.template.file = ncdata_hist,
  nc.compression = 4,
  nc.spatial.ref = "spatial_ref",
  nc.file.prefix = "stochastic",
  nc.file.suffix = "_cst_0"
)
