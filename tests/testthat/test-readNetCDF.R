
library(dplyr)
library(tidyr)
library(ncdf4)

testthat::test_that("Check if netcdf object is correctly created", {

  ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  ncdata <- readNetcdf(nc.file = ncfile, signif.digits = 2)
  expect_equal(length(ncdata), 5)

})
