library(testthat)

testthat::test_that("Check if netcdf object is correctly created", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")
  ncdata <- read_netcdf(nc.file = nc_path, signif.digits = 2)
  expect_equal(length(ncdata), 5)
})


test_that("read_netcdf returns correct structure for a sample file", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")

  result <- read_netcdf(nc.file = nc_path)
  expect_type(result, "list")
  expect_true(all(c("data", "grid", "date", "dimensions", "attributes") %in% names(result)))
  expect_s3_class(result$grid, "tbl")
  expect_true(is.list(result$data))
  expect_true(inherits(result$date, "Date"))
})

test_that("read_netcdf returns list of correct length when omit.empty = FALSE", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")

  result <- read_netcdf(nc.file = nc_path, omit.empty = FALSE)
  expect_equal(length(result$data), nrow(result$grid))
})

test_that("read_netcdf throws error for missing file", {
  expect_error(read_netcdf("no_such_file.nc"), "does not exist")
})

test_that("read_netcdf returns fewer grids when omit.empty = TRUE", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")

  res_all <- read_netcdf(nc.file = nc_path, omit.empty = FALSE)
  res_omit <- read_netcdf(nc.file = nc_path, omit.empty = TRUE)

  expect_true(length(res_omit$data) <= length(res_all$data))
})

test_that("read_netcdf handles leap.days = FALSE gracefully", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")

  expect_error(
    read_netcdf(nc.file = nc_path, leap.days = FALSE),
    NA # Should not error, even if dates are adjusted
  )
})

test_that("read_netcdf warns or errors for missing spatial reference", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")

  expect_error(read_netcdf(nc.file = nc_path, spatial.ref = "no_such_var"))
})
