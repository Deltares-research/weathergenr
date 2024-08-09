

library(weathergenr)
library(dplyr)
library(tidyr)
library(ncdf4)

ncdata1 <- readNetcdf(nc.file = "data/ntoum_era5_data.nc")
ncdata2 <- readNetcdf(nc.file = "data/meteo_allvars_regridded.nc")
ncdata3 <- readNetcdf(nc.file = "data/meteo_allvars_regridded_spatialref_nolayer_lonlat.nc")

ncdata3[[1]]
ncdata3[[2]]
ncdata3[[6]]

leap.days = TRUE
omit.empty = TRUE

nc.file = 
spatial.ref = "layer"


nc.file = "data/ntoum_era5_data.nc"
spatial.ref = "spatial_ref"