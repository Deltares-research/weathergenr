
test_that("write_netcdf round-trip matches read_netcdf", {

  signif_digits <- 3

  ncfile <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  expect_true(file.exists(ncfile))

  # Read template file (avoid Feb 29 warning)
  ncdata <- weathergenr::read_netcdf(
    nc.file = ncfile,
    leap.days = TRUE,
    omit.empty = TRUE,
    spatial.ref = "spatial_ref",
    signif.digits = signif_digits
  )

  output_path <- tempdir()
  file_prefix <- "testfile"
  file_suffix <- format(Sys.time(), "%Y%m%d%H%M%S")

  # Write
  out_path <- weathergenr::write_netcdf(
    data = ncdata$data,
    coord.grid = ncdata$grid,
    origin.date = as.character(ncdata$date[1]),
    output.path = output_path,
    calendar.type = "proleptic_gregorian",
    nc.template.file = ncfile,
    nc.compression = 4,
    nc.spatial.ref = "spatial_ref",
    nc.file.prefix = file_prefix,
    nc.file.suffix = file_suffix,
    signif.digits = signif_digits
  )

  expect_true(file.exists(out_path))

  # Read written file (avoid Feb 29 warning)
  ncdata2 <- weathergenr::read_netcdf(
    nc.file = out_path,
    leap.days = TRUE,
    signif.digits = signif_digits
  )

  expect_equal(ncdata2$date, ncdata$date)
  expect_equal(ncdata2$data, ncdata$data)
  expect_equal(ncdata2$grid, ncdata$grid)
})




