# test-read_netcdf.R
# Unit tests for read_netcdf() aligned with current behavior (incl. keep_leap_day warning)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
nc_test_path <- function() {
  system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
}

skip_if_no_netcdf <- function(path) {
  if (!file.exists(path)) skip("No NetCDF file available for testing.")
}

# -----------------------------------------------------------------------------
# Core structure
# -----------------------------------------------------------------------------

test_that("Check if netcdf object is correctly created", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  ncdata <- read_netcdf(nc_path = nc_path, signif_digits = 2)

  # Return is a 5-element list by contract
  expect_type(ncdata, "list")
  expect_true(all(c("data", "grid", "date", "dimensions", "attributes") %in% names(ncdata)))
  expect_equal(length(ncdata), 5)
})

test_that("read_netcdf returns correct structure for a sample file", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  result <- read_netcdf(nc_path = nc_path)

  expect_type(result, "list")
  expect_true(all(c("data", "grid", "date", "dimensions", "attributes") %in% names(result)))
  expect_s3_class(result$grid, "tbl_df")
  expect_true(is.list(result$data))
  expect_true(inherits(result$date, "Date"))
  expect_true(is.list(result$dimensions))
  expect_true(is.list(result$attributes))
})

test_that("read_netcdf returns list of correct length when drop_all_na = FALSE", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  result <- read_netcdf(nc_path = nc_path, drop_all_na = FALSE)

  expect_equal(length(result$data), nrow(result$grid))
})

test_that("read_netcdf throws error for missing file", {
  expect_error(
    read_netcdf("no_such_file.nc"),
    regexp = "File does not exist|does not exist",
    fixed = FALSE
  )
})

test_that("read_netcdf returns fewer-or-equal grids when drop_all_na = TRUE", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  res_all <- read_netcdf(nc_path = nc_path, drop_all_na = FALSE)
  res_omit <- read_netcdf(nc_path = nc_path, drop_all_na = TRUE)

  expect_true(length(res_omit$data) <= length(res_all$data))
  expect_equal(length(res_omit$data), nrow(res_omit$grid))
})

# -----------------------------------------------------------------------------
# keep_leap_day behavior (UPDATED NAMES)
# -----------------------------------------------------------------------------
# Intended behavior:
# - If keep_leap_day=FALSE but Feb 29 exists in derived dates => WARNING, but still returns
#   full date vector matching time dimension (no dropping).

test_that("read_netcdf handles keep_leap_day = FALSE (drops Feb 29 if present)", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  out <- suppressWarnings(
    read_netcdf(nc_path = nc_path, keep_leap_day = FALSE)
  )

  expect_true(inherits(out$date, "Date"))

  # Count leap days in the FULL series (keep_leap_day = TRUE)
  out_full <- suppressWarnings(
    read_netcdf(nc_path = nc_path, keep_leap_day = TRUE)
  )

  expect_equal(length(out_full$date), length(out_full$dimensions$time))

  n_feb29 <- sum(format(out_full$date, "%m-%d") == "02-29")

  # If there are leap days, they should be dropped from out$date
  expect_equal(length(out$date), length(out$dimensions$time) - n_feb29)

  if (n_feb29 > 0) {
    expect_false(any(format(out$date, "%m-%d") == "02-29"))
  }
})


# -----------------------------------------------------------------------------
# Spatial reference
# -----------------------------------------------------------------------------

test_that("read_netcdf errors for missing spatial reference variable", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  expect_error(
    read_netcdf(nc_path = nc_path, spatial_ref = "no_such_var"),
    regexp = "Spatial reference variable|not found",
    fixed = FALSE
  )
})

# -----------------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------------

test_that("read_netcdf validates input parameters correctly", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  expect_error(
    read_netcdf(nc_path = nc_path, keep_leap_day = "TRUE"),
    regexp = "keep_leap_day.*logical|single logical",
    fixed = FALSE
  )

  expect_error(
    read_netcdf(nc_path = nc_path, drop_all_na = 1),
    regexp = "drop_all_na.*logical|single logical",
    fixed = FALSE
  )

  expect_error(
    read_netcdf(nc_path = nc_path, spatial_ref = c("a", "b")),
    regexp = "spatial_ref.*single character",
    fixed = FALSE
  )

  expect_error(
    read_netcdf(nc_path = nc_path, signif_digits = -1),
    regexp = "signif_digits.*positive integer",
    fixed = FALSE
  )

  expect_error(
    read_netcdf(nc_path = nc_path, signif_digits = 1.5),
    regexp = "signif_digits.*positive integer",
    fixed = FALSE
  )
})

# -----------------------------------------------------------------------------
# Variable selection + renaming
# -----------------------------------------------------------------------------

test_that("read_netcdf handles variable selection correctly", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  all_vars <- read_netcdf(nc_path = nc_path)
  var_names <- names(all_vars$data[[1]])

  if (length(var_names) > 1) {
    result <- read_netcdf(nc_path = nc_path, var = var_names[1])
    expect_equal(ncol(result$data[[1]]), 1)
    expect_equal(names(result$data[[1]]), var_names[1])
  }

  expect_error(
    read_netcdf(nc_path = nc_path, var = "nonexistent_var"),
    regexp = "Variables not found in NetCDF|Variables not found",
    fixed = FALSE
  )
})

test_that("read_netcdf handles variable renaming correctly", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  all_vars <- read_netcdf(nc_path = nc_path)
  var_names <- names(all_vars$data[[1]])

  if (length(var_names) > 0) {
    new_name <- paste0(var_names[1], "_renamed")
    rename_vec <- stats::setNames(new_name, var_names[1])

    result <- read_netcdf(nc_path = nc_path, var_name = rename_vec)
    expect_true(new_name %in% names(result$data[[1]]))
    expect_false(var_names[1] %in% names(result$data[[1]]))
  }

  if (length(var_names) > 1) {
    dup_rename <- stats::setNames(c("same", "same"), c(var_names[1], var_names[2]))

    expect_error(
      read_netcdf(nc_path = nc_path, var_name = dup_rename),
      regexp = "duplicate target names|contains duplicate target names",
      fixed = FALSE
    )
  }
})

# -----------------------------------------------------------------------------
# Date + dimension consistency
# -----------------------------------------------------------------------------

test_that("read_netcdf date vector has consistent length", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  result <- read_netcdf(nc_path = nc_path)

  expect_equal(length(result$date), length(result$dimensions$time))
  expect_equal(length(result$date), nrow(result$data[[1]]))
})

# -----------------------------------------------------------------------------
# signif_digits (non-brittle check)
# -----------------------------------------------------------------------------

test_that("read_netcdf handles signif_digits without error and returns numeric columns", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  result <- read_netcdf(nc_path = nc_path, signif_digits = 3)

  cell <- result$data[[1]]
  expect_true(is.data.frame(cell))
  expect_true(all(vapply(cell, is.numeric, logical(1))))
})

# -----------------------------------------------------------------------------
# Grid coordinate consistency
# -----------------------------------------------------------------------------

test_that("read_netcdf grid coordinates match dimensions", {
  nc_path <- nc_test_path()
  skip_if_no_netcdf(nc_path)

  result <- read_netcdf(nc_path = nc_path)

  dim_names <- names(result$dimensions)
  non_time_dims <- setdiff(dim_names, "time")

  if (length(non_time_dims) >= 2) {
    x_dim_len <- length(result$dimensions[[non_time_dims[1]]])
    y_dim_len <- length(result$dimensions[[non_time_dims[2]]])

    expect_true(all(result$grid$xind >= 1 & result$grid$xind <= x_dim_len))
    expect_true(all(result$grid$yind >= 1 & result$grid$yind <= y_dim_len))

    x_vals <- result$dimensions[[non_time_dims[1]]]
    y_vals <- result$dimensions[[non_time_dims[2]]]

    sample_size <- min(10, nrow(result$grid))
    for (i in seq_len(sample_size)) {
      expect_equal(result$grid$x[i], x_vals[result$grid$xind[i]], tolerance = 1e-6)
      expect_equal(result$grid$y[i], y_vals[result$grid$yind[i]], tolerance = 1e-6)
    }
  } else {
    skip("NetCDF does not expose at least two non-time dimensions; grid coordinate test skipped.")
  }
})
