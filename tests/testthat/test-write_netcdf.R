test_that("write_netcdf round-trip matches read_netcdf", {

  signif_digits <- 3

  template_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  expect_true(file.exists(template_path))

  # Read template file (avoid Feb 29 warning)
  ncdata <- weathergenr::read_netcdf(
    nc_path        = template_path,
    keep_leap_day  = TRUE,
    drop_all_na    = TRUE,
    spatial_ref    = "spatial_ref",
    signif_digits  = signif_digits
  )

  out_dir <- tempdir()
  file_prefix <- "testfile"
  file_suffix <- format(Sys.time(), "%Y%m%d%H%M%S")

  # Write
  out_path <- weathergenr::write_netcdf(
    data          = ncdata$data,
    grid          = ncdata$grid,
    out_dir       = out_dir,
    origin_date   = as.character(ncdata$date[1]),
    calendar      = "proleptic_gregorian",
    template_path = template_path,
    compression   = 4,
    spatial_ref   = "spatial_ref",
    file_prefix   = file_prefix,
    file_suffix   = file_suffix,
    signif_digits = signif_digits,
    verbose       = FALSE
  )

  expect_true(file.exists(out_path))

  # Read written file (avoid Feb 29 warning)
  ncdata2 <- weathergenr::read_netcdf(
    nc_path       = out_path,
    keep_leap_day = TRUE,
    signif_digits = signif_digits
  )

  expect_equal(ncdata2$date, ncdata$date)
  expect_equal(ncdata2$data, ncdata$data)
  expect_equal(ncdata2$grid, ncdata$grid)
})
