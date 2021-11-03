## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message=FALSE)

## ----setup--------------------------------------------------------------------
library(gridwegen)

## ----ncdata-------------------------------------------------------------------
nc_data <- readNetcdf(
    nc.path = system.file('extdata', package = 'gridwegen'),
    nc.file = "/ntoum.nc",
    nc.dimnames = list(x = "lon", y = "lat", time = "time"),
    nc.variables = c("precip", "temp", "temp_min", "temp_max"),
    origin.date = as.Date("1981-01-01"),
    has.leap.days = TRUE)

