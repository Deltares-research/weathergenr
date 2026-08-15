


# test-read_netcdf.R
# Unit tests for read_netcdf() aligned with current behavior (incl. keep_leap_day warning)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
nc_test_path <- function() {
  system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
}

skip_if_no_netcdf <- function(path) {
  if (!file.exists(path)) testthat::skip("No NetCDF file available for testing.")
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
    testthat::skip("NetCDF does not expose at least two non-time dimensions; grid coordinate test skipped.")
  }
})


################################################################################
################################################################################


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


################################################################################
# CF calendar handling
################################################################################

test_that(".cf_calendar_kind classifies the CF calendar names", {
  kind <- weathergenr:::.cf_calendar_kind

  # Absent or empty means standard, per CF.
  expect_equal(kind(NULL), "standard")
  expect_equal(kind(NA_character_), "standard")
  expect_equal(kind(""), "standard")

  expect_equal(kind("standard"), "standard")
  expect_equal(kind("proleptic_gregorian"), "standard")
  expect_equal(kind("  Gregorian  "), "standard")
  expect_equal(kind("julian"), "standard")

  expect_equal(kind("noleap"), "noleap")
  expect_equal(kind("365_day"), "noleap")

  expect_equal(kind("360_day"), "unrepresentable")
  expect_equal(kind("all_leap"), "unrepresentable")

  expect_equal(kind("martian"), "unknown")
})

test_that("noleap index conversion round-trips and skips Feb 29", {
  to_date <- weathergenr:::.noleap_offset_to_date

  # 2020 is a Gregorian leap year; on a noleap axis Feb 29 must not appear.
  d <- to_date(as.Date("2020-01-01"), 0:364)
  expect_equal(length(d), 365L)
  expect_equal(d[1], as.Date("2020-01-01"))
  expect_equal(d[365], as.Date("2020-12-31"))
  expect_false(any(format(d, "%m-%d") == "02-29"))

  # Day 59 is 1 March on a noleap calendar, where Gregorian 2020 gives 29 Feb.
  expect_equal(to_date(as.Date("2020-01-01"), 59), as.Date("2020-03-01"))

  # Offsets may run backwards from the origin.
  expect_equal(to_date(as.Date("2020-01-01"), -1), as.Date("2019-12-31"))
  expect_equal(to_date(as.Date("2020-01-01"), -365), as.Date("2019-01-01"))

  # Fractional offsets truncate, matching as.Date() on a POSIXct.
  expect_equal(to_date(as.Date("2020-01-01"), 1.75), as.Date("2020-01-02"))

  # An origin that does not exist on this calendar is rejected.
  expect_error(to_date(as.Date("2020-02-29"), 0), "29 February")
})

test_that(".noleap_date_to_offset inverts .noleap_offset_to_date", {
  origin <- as.Date("2020-01-01")
  offsets <- c(0, 1, 59, 364, 365, 730, 1095)

  d <- weathergenr:::.noleap_offset_to_date(origin, offsets)
  expect_equal(weathergenr:::.noleap_date_to_offset(origin, d), offsets)

  expect_error(
    weathergenr:::.noleap_date_to_offset(origin, as.Date("2020-02-29")),
    "29 February"
  )
})

test_that(".cf_time_to_date honours the calendar attribute", {
  f <- weathergenr:::.cf_time_to_date
  units <- "days since 2020-01-01 00:00:00"

  # The two calendars diverge as soon as a Gregorian leap day is crossed.
  expect_equal(f(units, "standard", 59), as.Date("2020-02-29"))
  expect_equal(f(units, "noleap", 59), as.Date("2020-03-01"))

  # An absent calendar keeps the standard decode.
  expect_equal(f(units, NULL, 59), as.Date("2020-02-29"))

  # Sub-daily units reduce to the same day.
  expect_equal(f("hours since 2020-01-01 00:00:00", "noleap", 24 * 59), as.Date("2020-03-01"))

  # Calendars R's Date class cannot represent are refused, not approximated.
  expect_error(f(units, "360_day", 0), "cannot be represented", fixed = TRUE)
  expect_error(f(units, "all_leap", 0), "cannot be represented", fixed = TRUE)

  # An unrecognised calendar warns and falls back rather than failing the read.
  expect_warning(res <- f(units, "martian", 59), "Unrecognised NetCDF calendar")
  expect_equal(res, as.Date("2020-02-29"))
})

test_that("write_netcdf/read_netcdf round-trip a noleap time axis", {
  template_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  skip_if(!file.exists(template_path), "template NetCDF not available")

  ncdata <- weathergenr::read_netcdf(template_path, keep_leap_day = TRUE, verbose = FALSE)

  # Two 365-day years spanning a Gregorian leap year: this is exactly where the
  # old proleptic-Gregorian decode lost a day.
  expected <- generate_noleap_dates(as.Date("2020-01-01"), 730)
  nt <- length(expected)
  expect_equal(expected[nt], as.Date("2021-12-31"))

  payload <- lapply(1:2, function(i) list(precip = as.numeric(ncdata$data[[i]]$precip[seq_len(nt)])))

  out_dir <- tempfile("nc_noleap_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  out_path <- weathergenr::write_netcdf(
    data          = payload,
    grid          = ncdata$grid[1:2, ],
    out_dir       = out_dir,
    origin_date   = expected[1],
    calendar      = "noleap",
    dates         = expected,
    template_path = template_path,
    verbose       = FALSE
  )

  back <- weathergenr::read_netcdf(out_path, keep_leap_day = TRUE, verbose = FALSE)
  expect_equal(back$date, expected)
  expect_equal(back$date[nt], as.Date("2021-12-31"))
})

test_that("write_netcdf validates calendar and dates against origin_date", {
  template_path <- system.file("extdata", "ntoum_era5_data.nc", package = "weathergenr")
  skip_if(!file.exists(template_path), "template NetCDF not available")

  ncdata <- weathergenr::read_netcdf(template_path, keep_leap_day = TRUE, verbose = FALSE)
  dates <- generate_noleap_dates(as.Date("2020-01-01"), 365)
  nt <- length(dates)
  payload <- lapply(1:2, function(i) list(precip = as.numeric(ncdata$data[[i]]$precip[seq_len(nt)])))
  out_dir <- tempfile("nc_validate_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  wr <- function(...) {
    weathergenr::write_netcdf(
      data = payload, grid = ncdata$grid[1:2, ], out_dir = out_dir,
      template_path = template_path, verbose = FALSE, ...
    )
  }

  expect_error(wr(origin_date = dates[1], calendar = "360_day"), "'calendar' must be one of")
  expect_error(wr(origin_date = dates[1], calendar = "martian"), "'calendar' must be one of")

  # The trap this guards: a conventional epoch instead of the series start.
  expect_error(
    wr(origin_date = as.Date("1970-01-01"), calendar = "noleap", dates = dates),
    "must be the first entry of 'dates'"
  )

  expect_error(
    wr(origin_date = dates[1], calendar = "noleap", dates = dates[-1]),
    "one entry per time step"
  )
})
