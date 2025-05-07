
library(dplyr)
library(tidyr)
library(ncdf4)
library(weathergenr)

testthat::test_that("Check if write_netcdf works correctly", {

  signif_digits = 3

  #Read file
  ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  ncdata <- readNetcdf(nc.file = ncfile, leap.days = FALSE,  omit.empty = TRUE,
                       spatial.ref = "spatial_ref", signif.digits = signif_digits)

  # save to new file
  writeNetcdf(
    data = ncdata$data,
    coord.grid = ncdata$grid,
    origin.date =  ncdata$date[1],
    output.path = paste0(tempdir(),"/"),
    calendar.type = "proleptic_gregorian",
    nc.template.file = ncfile,
    nc.compression = 4,
    nc.spatial.ref = "spatial_ref",
    nc.file.prefix = "testfile",
    nc.file.suffix = "7",
    signif.digits = signif_digits
  )

  ncfile2 <- paste0(tempdir(),"/", "testfile_7.nc")
  ncdata2 <- readNetcdf(nc.file = ncfile2, leap.days = FALSE,
                        signif.digits = signif_digits)

  expect_equal(ncdata2$date, ncdata$date)
  expect_equal(ncdata2$data, ncdata$data)
  expect_equal(ncdata2$grid, ncdata$grid)

})
