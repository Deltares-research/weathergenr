library(ncdf4)

testthat::test_that("Check if write_netcdf works correctly", {
  signif_digits <- 3

  # Read file
  ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  ncdata <- read_netcdf(
    nc.file = ncfile, leap.days = FALSE, omit.empty = TRUE,
    spatial.ref = "spatial_ref", signif.digits = signif_digits
  )


  file_prefix <- "testfile"
  file_suffix <- as.integer(Sys.time())
  output_path <- "C:/TEMP/"

  # save to new file
  writeNetcdf(
    data = ncdata$data,
    coord.grid = ncdata$grid,
    origin.date = ncdata$date[1],
    output.path = output_path,
    calendar.type = "proleptic_gregorian",
    nc.template.file = ncfile,
    nc.compression = 4,
    nc.spatial.ref = "spatial_ref",
    nc.file.prefix = file_prefix,
    nc.file.suffix = file_suffix,
    signif.digits = signif_digits
  )

  ncfile2 <- paste0(output_path, file_prefix, "_", file_suffix, ".nc")
  ncdata2 <- read_netcdf(
    nc.file = ncfile2, leap.days = FALSE,
    signif.digits = signif_digits
  )

  expect_equal(ncdata2$date, ncdata$date)
  expect_equal(ncdata2$data, ncdata$data)
  expect_equal(ncdata2$grid, ncdata$grid)
})
