# test-io_zarr.R
# Unit tests for read_zarr() and write_zarr()
# functions under test: read_zarr(), write_zarr()
# relative path: tests/testthat/test-io_zarr.R
# purpose: validate input handling and round-trip behavior for Zarr I/O.

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

skip_if_no_pizzarr <- function() {
  if (!requireNamespace("pizzarr", quietly = TRUE)) {
    skip("Package 'pizzarr' is not available.")
  }
}

# Create a minimal test Zarr store for testing using write_zarr
create_test_zarr <- function(dir, nx = 2, ny = 2, nt = 10, include_leap = FALSE) {
  skip_if_no_pizzarr()

  if (dir.exists(dir)) unlink(dir, recursive = TRUE)

  # Create dates
  if (include_leap) {
    # Include a leap year with Feb 29
    dates <- seq(as.Date("2000-01-01"), by = "day", length.out = nt)
  } else {
    dates <- seq(as.Date("2001-01-01"), by = "day", length.out = nt)
  }

  # Create coordinate values
  x_vals <- seq(10, by = 0.5, length.out = nx)
  y_vals <- seq(50, by = 0.5, length.out = ny)

  # Create grid
  grid <- expand.grid(xind = seq_len(nx), yind = seq_len(ny), KEEP.OUT.ATTRS = FALSE)
  grid$x <- x_vals[grid$xind]
  grid$y <- y_vals[grid$yind]
  grid <- tibble::add_column(tibble::as_tibble(grid), id = seq_len(nrow(grid)), .before = 1)

  # Create data
  ncell <- nx * ny
  set.seed(42)
  data <- lapply(seq_len(ncell), function(i) {
    data.frame(
      precip = runif(nt, 0, 10),
      temp = rnorm(nt, 15, 5)
    )
  })

  # Create attributes
  attributes <- list(
    global = list(Conventions = "CF-1.8", title = "Test Zarr Store"),
    variables = list(
      precip = list(units = "mm", long_name = "precipitation"),
      temp = list(units = "degC", long_name = "temperature")
    )
  )

  # Write using write_zarr
  write_zarr(
    dir = dir,
    data = data,
    grid = grid,
    date = dates,
    attributes = attributes,
    x_dim_name = "lon",
    y_dim_name = "lat"
  )

  # Time values (days since origin) for return
  origin <- as.Date("2000-01-01")
  time_vals <- as.numeric(dates - origin)

  list(
    dir = dir,
    dates = dates,
    x_vals = x_vals,
    y_vals = y_vals,
    time_vals = time_vals,
    nx = nx,
    ny = ny,
    nt = nt
  )
}


# -----------------------------------------------------------------------------
# read_zarr: Input validation tests
# -----------------------------------------------------------------------------

test_that("read_zarr validates dir parameter", {
  expect_error(
    read_zarr(dir = NULL),
    regexp = "dir must be a non-empty character path"
  )

  expect_error(
    read_zarr(dir = ""),
    regexp = "dir must be a non-empty character path"
  )

  expect_error(
    read_zarr(dir = c("a", "b")),
    regexp = "dir must be a non-empty character path"
  )

  expect_error(
    read_zarr(dir = 123),
    regexp = "dir must be a non-empty character path"
  )

  expect_error(
    read_zarr(dir = "/nonexistent/path/to/zarr"),
    regexp = "Directory does not exist"
  )
})

test_that("read_zarr validates keep_leap_day parameter", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  expect_error(
    read_zarr(dir = tmp_dir, keep_leap_day = "TRUE"),
    regexp = "keep_leap_day must be a single logical"
  )

  expect_error(
    read_zarr(dir = tmp_dir, keep_leap_day = c(TRUE, FALSE)),
    regexp = "keep_leap_day must be a single logical"
  )

  expect_error(
    read_zarr(dir = tmp_dir, keep_leap_day = 1),
    regexp = "keep_leap_day must be a single logical"
  )
})

test_that("read_zarr validates drop_all_na parameter", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  expect_error(
    read_zarr(dir = tmp_dir, drop_all_na = "TRUE"),
    regexp = "drop_all_na must be a single logical"
  )

  expect_error(
    read_zarr(dir = tmp_dir, drop_all_na = c(TRUE, FALSE)),
    regexp = "drop_all_na must be a single logical"
  )
})

test_that("read_zarr validates verbose parameter", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  expect_error(
    read_zarr(dir = tmp_dir, verbose = "TRUE"),
    regexp = "verbose must be a single logical"
  )

  expect_error(
    read_zarr(dir = tmp_dir, verbose = c(TRUE, FALSE)),
    regexp = "verbose must be a single logical"
  )
})

test_that("read_zarr validates spatial_ref parameter", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  expect_error(
    read_zarr(dir = tmp_dir, spatial_ref = c("a", "b")),
    regexp = "spatial_ref must be a single character string"
  )

  expect_error(
    read_zarr(dir = tmp_dir, spatial_ref = ""),
    regexp = "spatial_ref must be a single character string"
  )

  expect_error(
    read_zarr(dir = tmp_dir, spatial_ref = 123),
    regexp = "spatial_ref must be a single character string"
  )
})

test_that("read_zarr validates slice_mode parameter", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  expect_error(
    read_zarr(dir = tmp_dir, slice_mode = 2),
    regexp = "slice_mode must be either 0.*or 1"
  )

  expect_error(
    read_zarr(dir = tmp_dir, slice_mode = "1"),
    regexp = "slice_mode must be either 0.*or 1"
  )

  expect_error(
    read_zarr(dir = tmp_dir, slice_mode = c(0, 1)),
    regexp = "slice_mode must be either 0.*or 1"
  )
})

test_that("read_zarr validates signif_digits parameter", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  expect_error(
    read_zarr(dir = tmp_dir, signif_digits = -1),
    regexp = "signif_digits must be a positive integer"
  )

  expect_error(
    read_zarr(dir = tmp_dir, signif_digits = 0),
    regexp = "signif_digits must be a positive integer"
  )

  expect_error(
    read_zarr(dir = tmp_dir, signif_digits = 1.5),
    regexp = "signif_digits must be a positive integer"
  )

  expect_error(
    read_zarr(dir = tmp_dir, signif_digits = "3"),
    regexp = "signif_digits must be a positive integer"
  )
})

test_that("read_zarr validates var_name parameter", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  # Not a named vector
  expect_error(
    read_zarr(dir = tmp_dir, var_name = c("new1", "new2")),
    regexp = "var_name must be a \\*named\\* character vector"
  )

  # Empty names
  expect_error(
    read_zarr(dir = tmp_dir, var_name = stats::setNames("new1", "")),
    regexp = "var_name must be a \\*named\\* character vector"
  )

  # Empty target names
  expect_error(
    read_zarr(dir = tmp_dir, var_name = stats::setNames("", "precip")),
    regexp = "var_name contains empty target names"
  )

  # Duplicate target names
  expect_error(
    read_zarr(dir = tmp_dir, var_name = stats::setNames(c("same", "same"), c("precip", "temp"))),
    regexp = "var_name contains duplicate target names"
  )
})

test_that("read_zarr validates var parameter", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  expect_error(
    read_zarr(dir = tmp_dir, var = 123),
    regexp = "var must be a character vector or NULL"
  )

  expect_error(
    read_zarr(dir = tmp_dir, var = "nonexistent_var"),
    regexp = "Variables not found in Zarr store"
  )
})

test_that("read_zarr errors for missing spatial_ref variable", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  expect_error(
    read_zarr(dir = tmp_dir, spatial_ref = "nonexistent_spatial_ref"),
    regexp = "Spatial reference variable not found"
  )
})


# -----------------------------------------------------------------------------
# read_zarr: Functional tests
# -----------------------------------------------------------------------------

test_that("read_zarr returns correct structure", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  info <- create_test_zarr(tmp_dir)

  result <- read_zarr(dir = tmp_dir)

  expect_type(result, "list")
  expect_true(all(c("data", "grid", "date", "dimensions", "attributes") %in% names(result)))
  expect_equal(length(result), 5)

  expect_s3_class(result$grid, "tbl_df")
  expect_true(is.list(result$data))
  expect_true(inherits(result$date, "Date"))
  expect_true(is.list(result$dimensions))
  expect_true(is.list(result$attributes))
})

test_that("read_zarr returns correct dimensions", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  info <- create_test_zarr(tmp_dir, nx = 3, ny = 4, nt = 15)

  result <- read_zarr(dir = tmp_dir)

  # Grid should have nx * ny cells

expect_equal(nrow(result$grid), info$nx * info$ny)
  expect_equal(length(result$data), info$nx * info$ny)

  # Each data frame should have nt rows
  expect_equal(nrow(result$data[[1]]), info$nt)

  # Date vector should match time dimension
  expect_equal(length(result$date), info$nt)
})

test_that("read_zarr grid coordinates match dimensions", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  info <- create_test_zarr(tmp_dir)

  result <- read_zarr(dir = tmp_dir)

  # Check grid has required columns
  expect_true(all(c("id", "xind", "yind", "x", "y") %in% names(result$grid)))

  # Check indices are valid
  expect_true(all(result$grid$xind >= 1 & result$grid$xind <= info$nx))
  expect_true(all(result$grid$yind >= 1 & result$grid$yind <= info$ny))

  # Check coordinate values match
  for (i in seq_len(min(5, nrow(result$grid)))) {
    expect_equal(result$grid$x[i], info$x_vals[result$grid$xind[i]], tolerance = 1e-6)
    expect_equal(result$grid$y[i], info$y_vals[result$grid$yind[i]], tolerance = 1e-6)
  }
})

test_that("read_zarr handles variable selection", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  # Select single variable
  result <- read_zarr(dir = tmp_dir, var = "precip")
  expect_equal(ncol(result$data[[1]]), 1)
  expect_equal(names(result$data[[1]]), "precip")

  # Select multiple variables
  result2 <- read_zarr(dir = tmp_dir, var = c("precip", "temp"))
  expect_equal(ncol(result2$data[[1]]), 2)
  expect_true(all(c("precip", "temp") %in% names(result2$data[[1]])))
})

test_that("read_zarr handles variable renaming", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  rename_vec <- stats::setNames(c("precipitation", "temperature"), c("precip", "temp"))
  result <- read_zarr(dir = tmp_dir, var_name = rename_vec)

  expect_true("precipitation" %in% names(result$data[[1]]))
  expect_true("temperature" %in% names(result$data[[1]]))
  expect_false("precip" %in% names(result$data[[1]]))
  expect_false("temp" %in% names(result$data[[1]]))
})

test_that("read_zarr errors for var_name with variables not selected", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  # Select only precip but try to rename temp
  expect_error(
    read_zarr(dir = tmp_dir, var = "precip", var_name = stats::setNames("new_temp", "temp")),
    regexp = "var_name contains variables not selected"
  )
})

test_that("read_zarr handles signif_digits", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  result <- read_zarr(dir = tmp_dir, signif_digits = 3)

  cell <- result$data[[1]]
  expect_true(is.data.frame(cell))
  expect_true(all(vapply(cell, is.numeric, logical(1))))
})

test_that("read_zarr handles keep_leap_day = FALSE", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # Create zarr with dates that include Feb 29
  info <- create_test_zarr(tmp_dir, nt = 366, include_leap = TRUE)

  result_full <- read_zarr(dir = tmp_dir, keep_leap_day = TRUE)
  result_no_leap <- read_zarr(dir = tmp_dir, keep_leap_day = FALSE)

  # Count leap days in full result
  n_feb29 <- sum(format(result_full$date, "%m-%d") == "02-29")

  if (n_feb29 > 0) {
    # Dates should be shorter
    expect_equal(length(result_no_leap$date), length(result_full$date) - n_feb29)

    # No Feb 29 in result
    expect_false(any(format(result_no_leap$date, "%m-%d") == "02-29"))

    # Data should also have fewer rows
    expect_equal(nrow(result_no_leap$data[[1]]), nrow(result_full$data[[1]]) - n_feb29)
  }
})

test_that("read_zarr returns attributes", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  result <- read_zarr(dir = tmp_dir)

  expect_true(is.list(result$attributes))
  expect_true("global" %in% names(result$attributes))
  expect_true("variables" %in% names(result$attributes))

  # Check global attributes
  expect_equal(result$attributes$global$Conventions, "CF-1.8")
  expect_equal(result$attributes$global$title, "Test Zarr Store")

  # Check variable attributes
  expect_true("precip" %in% names(result$attributes$variables))
  expect_equal(result$attributes$variables$precip$units, "mm")
})

test_that("read_zarr verbose output works", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  create_test_zarr(tmp_dir)

  expect_message(
    read_zarr(dir = tmp_dir, verbose = TRUE),
    regexp = "Reading variable"
  )
})


# -----------------------------------------------------------------------------
# write_zarr: Input validation tests
# -----------------------------------------------------------------------------

test_that("write_zarr validates dir parameter", {
  expect_error(
    write_zarr(dir = NULL),
    regexp = "dir must be a non-empty character path"
  )

  expect_error(
    write_zarr(dir = ""),
    regexp = "dir must be a non-empty character path"
  )

  expect_error(
    write_zarr(dir = 123),
    regexp = "dir must be a non-empty character path"
  )
})

test_that("write_zarr validates overwrite parameter", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # Create a dir that already exists
  dir.create(tmp_dir, showWarnings = FALSE)

  expect_error(
    write_zarr(dir = tmp_dir, overwrite = FALSE, data = list(), grid = data.frame(), date = Sys.Date()),
    regexp = "Zarr store already exists.*Set overwrite=TRUE"
  )

  expect_error(
    write_zarr(dir = tmp_dir, overwrite = "TRUE"),
    regexp = "overwrite must be a single logical"
  )
})

test_that("write_zarr validates verbose parameter", {
  grid <- data.frame(id = 1, xind = 1, yind = 1, x = 10, y = 50)
  data <- list(data.frame(var1 = 1:5))
  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = 5)

  expect_error(
    write_zarr(dir = tempfile(), data = data, grid = grid, date = dates, verbose = "TRUE"),
    regexp = "verbose must be a single logical"
  )
})

test_that("write_zarr validates compression parameter", {
  grid <- data.frame(id = 1, xind = 1, yind = 1, x = 10, y = 50)
  data <- list(data.frame(var1 = 1:5))
  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = 5)

  expect_error(
    write_zarr(dir = tempfile(), data = data, grid = grid, date = dates, compression = "invalid"),
    regexp = "compression must be one of"
  )
})

test_that("write_zarr validates compression_level parameter", {
  grid <- data.frame(id = 1, xind = 1, yind = 1, x = 10, y = 50)
  data <- list(data.frame(var1 = 1:5))
  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = 5)

  expect_error(
    write_zarr(dir = tempfile(), data = data, grid = grid, date = dates, compression_level = -1),
    regexp = "compression_level must be between 0 and 9"
  )

  expect_error(
    write_zarr(dir = tempfile(), data = data, grid = grid, date = dates, compression_level = 10),
    regexp = "compression_level must be between 0 and 9"
  )
})

test_that("write_zarr validates slice_mode parameter", {
  grid <- data.frame(id = 1, xind = 1, yind = 1, x = 10, y = 50)
  data <- list(data.frame(var1 = 1:5))
  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = 5)

  expect_error(
    write_zarr(dir = tempfile(), data = data, grid = grid, date = dates, slice_mode = 2),
    regexp = "slice_mode must be either 0.*or 1"
  )
})

test_that("write_zarr validates required components", {
  skip_if_no_pizzarr()

  expect_error(
    write_zarr(dir = tempfile(), data = NULL, grid = data.frame(), date = Sys.Date()),
    regexp = "data must be provided"
  )

  expect_error(
    write_zarr(dir = tempfile(), data = list(), grid = NULL, date = Sys.Date()),
    regexp = "grid must be provided"
  )

  expect_error(
    write_zarr(dir = tempfile(), data = list(), grid = data.frame(), date = NULL),
    regexp = "date must be provided"
  )
})

test_that("write_zarr validates data structure", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # Grid with 2 cells but data with 3 elements
  grid <- data.frame(id = 1:2, xind = 1:2, yind = c(1, 1), x = c(10, 11), y = c(50, 50))
  data <- list(data.frame(var1 = 1:5), data.frame(var1 = 1:5), data.frame(var1 = 1:5))
  dates <- seq(as.Date("2000-01-01"), by = "day", length.out = 5)

  expect_error(
    write_zarr(dir = tmp_dir, data = data, grid = grid, date = dates),
    regexp = "data must be a list with one data.frame per grid cell"
  )
})

test_that("write_zarr validates data columns", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  grid <- data.frame(id = 1, xind = 1, yind = 1, x = 10, y = 50)
  data <- list(data.frame())  # Empty data frame with no columns
  dates <- seq(as.Date("2000-01-01"), by = "day", length.out = 5)

  expect_error(
    write_zarr(dir = tmp_dir, data = data, grid = grid, date = dates),
    regexp = "data\\[\\[1\\]\\] must have named columns"
  )
})


# -----------------------------------------------------------------------------
# write_zarr: Functional tests
# -----------------------------------------------------------------------------

test_that("write_zarr creates valid Zarr store", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_write_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # Create test data
  nx <- 2
  ny <- 2
  nt <- 10
  ncell <- nx * ny

  grid <- expand.grid(xind = seq_len(nx), yind = seq_len(ny), KEEP.OUT.ATTRS = FALSE)
  grid$x <- 10 + (grid$xind - 1) * 0.5
  grid$y <- 50 + (grid$yind - 1) * 0.5
  grid <- tibble::add_column(tibble::as_tibble(grid), id = seq_len(nrow(grid)), .before = 1)

  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = nt)

  set.seed(123)
  data <- lapply(seq_len(ncell), function(i) {
    data.frame(
      precip = runif(nt, 0, 10),
      temp = rnorm(nt, 15, 5)
    )
  })

  result <- write_zarr(
    dir = tmp_dir,
    data = data,
    grid = grid,
    date = dates,
    verbose = FALSE
  )

  expect_true(dir.exists(tmp_dir))
  expect_equal(result, tmp_dir)

  # Check zarr structure exists
  expect_true(file.exists(file.path(tmp_dir, ".zgroup")))
})

test_that("write_zarr handles overwrite = TRUE", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_overwrite_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # Create minimal test data
  grid <- data.frame(id = 1, xind = 1, yind = 1, x = 10, y = 50)
  data <- list(data.frame(var1 = 1:5))
  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = 5)

  # Write first time
  write_zarr(dir = tmp_dir, data = data, grid = grid, date = dates)
  expect_true(dir.exists(tmp_dir))

  # Write again with overwrite
  expect_no_error(
    write_zarr(dir = tmp_dir, data = data, grid = grid, date = dates, overwrite = TRUE)
  )
})

test_that("write_zarr handles different compression options", {
  skip_if_no_pizzarr()

  grid <- data.frame(id = 1, xind = 1, yind = 1, x = 10, y = 50)
  data <- list(data.frame(var1 = 1:5))
  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = 5)

  for (comp in c("none", "zlib", "blosc")) {
    if (comp == "blosc" && !requireNamespace("Rarr", quietly = TRUE)) {
      next
    }
    tmp_dir <- tempfile(paste0("zarr_comp_", comp, "_"))
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    expect_no_error(
      write_zarr(dir = tmp_dir, data = data, grid = grid, date = dates, compression = comp)
    )
    expect_true(dir.exists(tmp_dir))
  }
})

test_that("write_zarr handles input_list parameter", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_input_list_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # Create input list mimicking read_zarr output
  grid <- data.frame(id = 1, xind = 1, yind = 1, x = 10, y = 50)
  data <- list(data.frame(var1 = 1:5))
  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = 5)

  input_list <- list(
    data = data,
    grid = grid,
    date = dates,
    dimensions = list(lon = 10, lat = 50, time = 0:4),
    attributes = list(
      global = list(title = "Test"),
      variables = list(var1 = list(units = "mm"))
    )
  )

  expect_no_error(
    write_zarr(dir = tmp_dir, input_list = input_list)
  )
  expect_true(dir.exists(tmp_dir))
})

test_that("write_zarr handles POSIXct dates", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_posix_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  grid <- data.frame(id = 1, xind = 1, yind = 1, x = 10, y = 50)
  data <- list(data.frame(var1 = 1:5))
  dates <- as.POSIXct(seq(as.Date("2001-01-01"), by = "day", length.out = 5), tz = "UTC")

  expect_no_error(
    write_zarr(dir = tmp_dir, data = data, grid = grid, date = dates)
  )
})

test_that("write_zarr verbose output works", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_verbose_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  grid <- data.frame(id = 1, xind = 1, yind = 1, x = 10, y = 50)
  data <- list(data.frame(var1 = 1:5))
  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = 5)

  expect_message(
    write_zarr(dir = tmp_dir, data = data, grid = grid, date = dates, verbose = TRUE),
    regexp = "Writing"
  )
})


# -----------------------------------------------------------------------------
# Round-trip tests
# -----------------------------------------------------------------------------

test_that("write_zarr then read_zarr round-trip preserves data", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_roundtrip_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # Create test data
  nx <- 2
  ny <- 2
  nt <- 10
  ncell <- nx * ny

  grid <- expand.grid(xind = seq_len(nx), yind = seq_len(ny), KEEP.OUT.ATTRS = FALSE)
  grid$x <- 10 + (grid$xind - 1) * 0.5
  grid$y <- 50 + (grid$yind - 1) * 0.5
  grid <- tibble::add_column(tibble::as_tibble(grid), id = seq_len(nrow(grid)), .before = 1)

  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = nt)

  set.seed(456)
  data <- lapply(seq_len(ncell), function(i) {
    data.frame(
      precip = round(runif(nt, 0, 10), 2),
      temp = round(rnorm(nt, 15, 5), 2)
    )
  })

  # Write
  write_zarr(
    dir = tmp_dir,
    data = data,
    grid = grid,
    date = dates,
    fill_value = -9999
  )

  # Read back
  result <- read_zarr(dir = tmp_dir)

  # Check structure
  expect_equal(length(result$data), ncell)
  expect_equal(nrow(result$grid), ncell)
  expect_equal(length(result$date), nt)

  # Check dates match
  expect_equal(result$date, dates)

  # Check grid coordinates
  expect_equal(result$grid$x, grid$x)
  expect_equal(result$grid$y, grid$y)

  # Check data values (with tolerance for float32 precision)
  for (k in seq_len(ncell)) {
    expect_equal(result$data[[k]]$precip, data[[k]]$precip, tolerance = 1e-4)
    expect_equal(result$data[[k]]$temp, data[[k]]$temp, tolerance = 1e-4)
  }
})

test_that("read_zarr output can be passed directly to write_zarr via input_list", {
  skip_if_no_pizzarr()
  tmp_dir1 <- tempfile("zarr_rt1_")
  tmp_dir2 <- tempfile("zarr_rt2_")
  on.exit({
    unlink(tmp_dir1, recursive = TRUE)
    unlink(tmp_dir2, recursive = TRUE)
  }, add = TRUE)

  # Create initial zarr
  info <- create_test_zarr(tmp_dir1)

  # Read it
  data1 <- read_zarr(dir = tmp_dir1)

  # Write it to a new location using input_list
  write_zarr(dir = tmp_dir2, input_list = data1)

  # Read the new one
  data2 <- read_zarr(dir = tmp_dir2)

  # Compare
  expect_equal(length(data2$data), length(data1$data))
  expect_equal(data2$date, data1$date)
  expect_equal(nrow(data2$grid), nrow(data1$grid))
})


# -----------------------------------------------------------------------------
# Edge cases
# -----------------------------------------------------------------------------

test_that("read_zarr handles single grid cell", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_single_cell_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  info <- create_test_zarr(tmp_dir, nx = 1, ny = 1, nt = 5)

  result <- read_zarr(dir = tmp_dir)

  expect_equal(nrow(result$grid), 1)
  expect_equal(length(result$data), 1)
  expect_equal(nrow(result$data[[1]]), 5)
})

test_that("read_zarr handles large grid", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_large_grid_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  info <- create_test_zarr(tmp_dir, nx = 10, ny = 10, nt = 30)

  result <- read_zarr(dir = tmp_dir)

  expect_equal(nrow(result$grid), 100)
  expect_equal(length(result$data), 100)
})

test_that("write_zarr handles NA values correctly", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_na_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  grid <- data.frame(id = 1, xind = 1, yind = 1, x = 10, y = 50)
  data <- list(data.frame(var1 = c(1, NA, 3, NA, 5)))
  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = 5)

  write_zarr(dir = tmp_dir, data = data, grid = grid, date = dates, fill_value = -9999)

  result <- read_zarr(dir = tmp_dir)

  # NA values should be preserved after round-trip (converted to fill_value and back)
  expect_equal(sum(is.na(result$data[[1]]$var1)), 2)
  expect_equal(result$data[[1]]$var1[c(1, 3, 5)], c(1, 3, 5), tolerance = 1e-4)
})

test_that("write_zarr handles custom chunk_size", {
  skip_if_no_pizzarr()
  tmp_dir <- tempfile("zarr_chunk_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  grid <- data.frame(id = 1:4, xind = c(1, 2, 1, 2), yind = c(1, 1, 2, 2),
                     x = c(10, 10.5, 10, 10.5), y = c(50, 50, 50.5, 50.5))
  data <- lapply(1:4, function(i) data.frame(var1 = 1:10))
  dates <- seq(as.Date("2001-01-01"), by = "day", length.out = 10)

  expect_no_error(
    write_zarr(
      dir = tmp_dir,
      data = data,
      grid = grid,
      date = dates,
      chunk_size = list(lon = 1, lat = 1, time = 5)
    )
  )
})
