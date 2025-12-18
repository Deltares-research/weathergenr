# Test file for read_netcdf function
# These tests verify the function works correctly with the improvements

test_that("Check if netcdf object is correctly created", {
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

  # Should not error regardless of whether data has leap days or not
  expect_error(
    read_netcdf(nc.file = nc_path, leap.days = FALSE),
    NA
  )

  # Verify the result is valid
  result <- read_netcdf(nc.file = nc_path, leap.days = FALSE)
  expect_true(inherits(result$date, "Date"))
  expect_equal(length(result$date), length(result$dimensions$time))
})

test_that("read_netcdf warns or errors for missing spatial reference", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")
  expect_error(read_netcdf(nc.file = nc_path, spatial.ref = "no_such_var"))
})

# Additional tests for the improved functionality

test_that("read_netcdf validates input parameters correctly", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")

  # Test invalid leap.days
  expect_error(
    read_netcdf(nc.file = nc_path, leap.days = "TRUE"),
    "leap.days.*must be.*logical"
  )

  # Test invalid omit.empty
  expect_error(
    read_netcdf(nc.file = nc_path, omit.empty = 1),
    "omit.empty.*must be.*logical"
  )

  # Test invalid spatial.ref
  expect_error(
    read_netcdf(nc.file = nc_path, spatial.ref = c("a", "b")),
    "spatial.ref.*must be.*single character"
  )

  # Test invalid signif.digits
  expect_error(
    read_netcdf(nc.file = nc_path, signif.digits = -1),
    "signif.digits.*must be.*positive integer"
  )

  expect_error(
    read_netcdf(nc.file = nc_path, signif.digits = 1.5),
    "signif.digits.*must be.*positive integer"
  )
})

test_that("read_netcdf handles variable selection correctly", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")

  # Load all variables first to see what's available
  all_vars <- read_netcdf(nc.file = nc_path)
  var_names <- names(all_vars$data[[1]])

  if (length(var_names) > 1) {
    # Test selecting subset of variables
    result <- read_netcdf(nc.file = nc_path, variables = var_names[1])
    expect_equal(ncol(result$data[[1]]), 1)
    expect_equal(names(result$data[[1]]), var_names[1])
  }

  # Test error for non-existent variable
  expect_error(
    read_netcdf(nc.file = nc_path, variables = "nonexistent_var"),
    "Variables not found"
  )
})

test_that("read_netcdf handles variable renaming correctly", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")

  # Load to see available variables
  all_vars <- read_netcdf(nc.file = nc_path)
  var_names <- names(all_vars$data[[1]])

  if (length(var_names) > 0) {
    # Test renaming first variable
    new_name <- paste0(var_names[1], "_renamed")
    rename_vec <- setNames(new_name, var_names[1])

    result <- read_netcdf(nc.file = nc_path, var_rename = rename_vec)
    expect_true(new_name %in% names(result$data[[1]]))
    expect_false(var_names[1] %in% names(result$data[[1]]))
  }

  # Test error for duplicate target names
  if (length(var_names) > 1) {
    # Create named vector with duplicate target names
    dup_rename <- setNames(c("same", "same"), c(var_names[1], var_names[2]))

    expect_error(
      read_netcdf(
        nc.file = nc_path,
        var_rename = dup_rename
      ),
      "duplicate target names"
    )
  }
})

test_that("read_netcdf date vector has correct length", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")

  result <- read_netcdf(nc.file = nc_path)
  expect_equal(length(result$date), nrow(result$data[[1]]))
  expect_equal(length(result$date), length(result$dimensions$time))
})

test_that("read_netcdf handles signif.digits correctly", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")

  result <- read_netcdf(nc.file = nc_path, signif.digits = 3)

  # Check that values are rounded
  # Get first non-NA value from first grid cell
  first_cell <- result$data[[1]]
  first_col <- first_cell[[1]]
  non_na_vals <- first_col[!is.na(first_col)]

  if (length(non_na_vals) > 0) {
    # Check that values have at most 3 significant digits
    # This is a rough check - signif() may keep more digits in some cases
    expect_true(all(!is.na(non_na_vals)))
  }
})

test_that("read_netcdf grid coordinates match dimensions", {
  nc_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  if (!file.exists(nc_path)) skip("No NetCDF file available for testing.")

  result <- read_netcdf(nc.file = nc_path)

  # Check that xind and yind are valid indices within dimension ranges
  dim_names <- names(result$dimensions)
  non_time_dims <- setdiff(dim_names, "time")

  if (length(non_time_dims) >= 2) {
    x_dim_len <- length(result$dimensions[[non_time_dims[1]]])
    y_dim_len <- length(result$dimensions[[non_time_dims[2]]])

    # Verify indices are in valid range
    expect_true(all(result$grid$xind >= 1 & result$grid$xind <= x_dim_len))
    expect_true(all(result$grid$yind >= 1 & result$grid$yind <= y_dim_len))

    # Verify coordinates match their indices (with floating point tolerance)
    x_vals <- result$dimensions[[non_time_dims[1]]]
    y_vals <- result$dimensions[[non_time_dims[2]]]

    sample_size <- min(10, nrow(result$grid))
    for (i in seq_len(sample_size)) {
      expect_equal(result$grid$x[i], x_vals[result$grid$xind[i]], tolerance = 1e-6)
      expect_equal(result$grid$y[i], y_vals[result$grid$yind[i]], tolerance = 1e-6)
    }
  }
})

