# ==============================================================================
# CF calendar handling
# ==============================================================================
#
# A CF time axis is "<unit> since <origin>" plus a `calendar` attribute, and the
# attribute changes what an offset means: on `noleap` the year is always 365
# days, so day 60 after 2020-01-01 is 2 March, where the Gregorian answer is
# 1 March. Decoding without reading the attribute silently loses one day per
# leap year -- which is what this package did to files it had written itself.
#
# The generator runs on a 365-day calendar throughout (see calendar.R), so
# `noleap` is the case that matters here; `all_leap` and `360_day` are rejected
# rather than approximated, because their dates (29 February in a common year,
# 30 February) have no representation in R's Date class at all.

#' Classify a CF calendar attribute
#'
#' @param calendar Character scalar or NULL. Raw `calendar` attribute value.
#' @return One of `"standard"`, `"noleap"`, `"unrepresentable"`, `"unknown"`.
#'   An absent or empty attribute is `"standard"`, per CF.
#' @keywords internal
#' @noRd
.cf_calendar_kind <- function(calendar) {
  if (is.null(calendar) || length(calendar) != 1L) return("standard")
  if (is.na(calendar) || !nzchar(as.character(calendar))) return("standard")

  cal <- tolower(trimws(as.character(calendar)))

  # `julian` diverges from Gregorian only before 1582, which no weather-generator
  # axis reaches; treating it as standard beats refusing to read the file.
  if (cal %in% c("standard", "gregorian", "proleptic_gregorian", "julian")) return("standard")
  if (cal %in% c("noleap", "no_leap", "365_day", "365day")) return("noleap")
  if (cal %in% c("all_leap", "366_day", "366day", "360_day", "360day")) return("unrepresentable")

  "unknown"
}

#' Days elapsed before the first of each month in a 365-day year
#' @keywords internal
#' @noRd
.noleap_month_offsets <- function() {
  c(0L, 31L, 59L, 90L, 120L, 151L, 181L, 212L, 243L, 273L, 304L, 334L)
}

#' Absolute day index on a 365-day calendar
#'
#' Counts days from a fixed origin with no leap adjustment whatsoever, so
#' differences between two indices are exact 365-day-calendar day counts.
#'
#' @param year,month,day Integer vectors.
#' @return Numeric vector of day indices.
#' @keywords internal
#' @noRd
.noleap_ymd_to_index <- function(year, month, day) {
  as.numeric(year) * 365 + .noleap_month_offsets()[month] + (as.numeric(day) - 1)
}

#' Invert [.noleap_ymd_to_index()]
#'
#' @param idx Numeric vector of 365-day-calendar day indices.
#' @return List with integer `year`, `month`, `day` vectors.
#' @keywords internal
#' @noRd
.noleap_index_to_ymd <- function(idx) {
  cum <- .noleap_month_offsets()
  year <- floor(idx / 365)
  doy0 <- idx - year * 365            # 0..364
  month <- findInterval(doy0, cum)    # 1..12
  list(
    year  = as.integer(year),
    month = as.integer(month),
    day   = as.integer(doy0 - cum[month] + 1)
  )
}

#' Convert 365-day-calendar offsets to Date
#'
#' @param origin_date Date scalar. The axis origin. Must not be 29 February,
#'   which does not exist on this calendar.
#' @param offsets_days Numeric vector of offsets in days; fractional values are
#'   floored, matching how `as.Date()` truncates a POSIXct.
#' @return Date vector.
#' @keywords internal
#' @noRd
.noleap_offset_to_date <- function(origin_date, offsets_days) {
  lt <- as.POSIXlt(origin_date)
  if (lt$mon == 1L && lt$mday == 29L) {
    stop("NetCDF time origin is 29 February, which does not exist on a noleap calendar.",
         call. = FALSE)
  }

  base <- .noleap_ymd_to_index(lt$year + 1900L, lt$mon + 1L, lt$mday)
  idx <- floor(base + as.numeric(offsets_days))

  if (any(!is.finite(idx))) {
    stop("NetCDF time axis contains non-finite offsets.", call. = FALSE)
  }

  ymd <- .noleap_index_to_ymd(idx)
  if (any(ymd$year < 0L | ymd$year > 9999L)) {
    stop("NetCDF time axis decodes to a year outside 0-9999 on the noleap calendar.",
         call. = FALSE)
  }

  as.Date(sprintf("%04d-%02d-%02d", ymd$year, ymd$month, ymd$day))
}

#' Convert Dates to 365-day-calendar offsets
#'
#' The inverse of [.noleap_offset_to_date()], used when encoding a time axis.
#'
#' @param origin_date Date scalar.
#' @param dates Date vector. Must contain no 29 February.
#' @return Numeric vector of day offsets from `origin_date`.
#' @keywords internal
#' @noRd
.noleap_date_to_offset <- function(origin_date, dates) {
  d <- as.POSIXlt(dates)
  if (any(d$mon == 1L & d$mday == 29L)) {
    stop("Dates contain 29 February, which does not exist on a noleap calendar.",
         call. = FALSE)
  }

  o <- as.POSIXlt(origin_date)
  .noleap_ymd_to_index(d$year + 1900L, d$mon + 1L, d$mday) -
    .noleap_ymd_to_index(o$year + 1900L, o$mon + 1L, o$mday)
}

#' Decode a CF time axis to Date
#'
#' @param units Character scalar. The `units` attribute, e.g.
#'   `"days since 2020-01-01 00:00:00"`.
#' @param calendar Character scalar or NULL. The `calendar` attribute.
#' @param time_vals Numeric vector of raw time coordinate values.
#' @return Date vector, or NULL when the units string is missing or unsupported.
#' @keywords internal
#' @noRd
.cf_time_to_date <- function(units, calendar, time_vals) {
  if (is.null(units) || length(units) != 1L || is.na(units) || !nzchar(units)) return(NULL)
  if (is.null(time_vals) || !is.numeric(time_vals) || length(time_vals) < 1L) return(NULL)

  m <- regexec("^\\s*(seconds|second|minutes|minute|hours|hour|days|day)\\s+since\\s+(.+?)\\s*$",
               units, ignore.case = TRUE)
  r <- regmatches(units, m)[[1]]
  if (length(r) != 3L) return(NULL)

  mult <- switch(tolower(r[2]),
                 "second" = 1, "seconds" = 1,
                 "minute" = 60, "minutes" = 60,
                 "hour"   = 3600, "hours" = 3600,
                 "day"    = 86400, "days" = 86400,
                 NULL)
  if (is.null(mult)) return(NULL)

  kind <- .cf_calendar_kind(calendar)

  if (identical(kind, "unrepresentable")) {
    stop("NetCDF time uses the '", calendar, "' calendar, whose dates cannot be ",
         "represented by R's Date class (29 February in a common year, or 30 February). ",
         "Convert the file to a standard or noleap calendar before reading.",
         call. = FALSE)
  }

  if (identical(kind, "unknown")) {
    warning("Unrecognised NetCDF calendar '", calendar,
            "'; decoding the time axis as a standard (Gregorian) calendar.",
            call. = FALSE)
    kind <- "standard"
  }

  origin_posix <- suppressWarnings(as.POSIXct(r[3], tz = "UTC"))
  if (is.na(origin_posix)) {
    origin_date <- suppressWarnings(as.Date(r[3]))
    if (is.na(origin_date)) return(NULL)
    origin_posix <- as.POSIXct(origin_date, tz = "UTC")
  }

  if (identical(kind, "standard")) {
    return(as.Date(origin_posix + as.numeric(time_vals) * mult))
  }

  # noleap: count 365-day years rather than adding real days to the origin.
  origin_date <- as.Date(origin_posix)
  origin_frac <- as.numeric(
    difftime(origin_posix, as.POSIXct(origin_date, tz = "UTC"), units = "days")
  )
  .noleap_offset_to_date(origin_date, origin_frac + as.numeric(time_vals) * mult / 86400)
}


#' Read NetCDF variables into tidy data frames
#'
#' Reads one or more variables from a NetCDF file and returns a named list of
#' data.frames (one per variable). Each data.frame contains the coordinate
#' columns (e.g., lon/lat/time) plus a single value column (optionally renamed).
#'
#' Robust to common NetCDF conventions and avoids the
#' "replacement has N rows, data has 0" failure by always building the
#' coordinate table at the correct length before assigning values.
#'
#' @param nc_path Character. Path to NetCDF file.
#' @param var Character vector or NULL. Variables to read. If NULL, reads all
#'   non-coordinate variables (and excludes `spatial_ref` if present).
#' @param var_name NULL, or named character vector.
#'   - If named: names are original variable names and values are new names.
#'   Renaming applies to output list names and the value-column name.
#' @param keep_leap_day Logical. If FALSE and time is convertible to Date, removes Feb 29 rows.
#' @param drop_all_na Logical. If TRUE, drop variables whose values are all NA
#'   after applying _FillValue/missing_value handling.
#' @param spatial_ref Character. Name of a spatial reference variable to ignore
#'   when auto-selecting variables (common in CF/CRS-encoded files).
#' @param signif_digits Integer or NULL. If provided, rounds numeric values via
#'   signif(x, signif_digits).
#' @param verbose Logical. If TRUE, prints basic progress messages.
#'
#' @return Named list of data.frames.
#'
#' @importFrom ncdf4 nc_open nc_close ncvar_get ncatt_get
#' @export
read_netcdf <- function(
    nc_path,
    var = NULL,
    var_name = NULL,
    keep_leap_day = TRUE,
    drop_all_na = TRUE,
    spatial_ref = "spatial_ref",
    signif_digits = NULL,
    verbose = FALSE
) {

  # ---------------------------------------------------------------------------
  # Input validation (match test expectations / regexes)
  # ---------------------------------------------------------------------------
  if (!is.character(nc_path) || length(nc_path) != 1L || !nzchar(nc_path)) {
    stop("nc_path must be a non-empty character path.", call. = FALSE)
  }
  if (!file.exists(nc_path)) {
    stop("File does not exist: ", nc_path, call. = FALSE)
  }

  if (!is.logical(keep_leap_day) || length(keep_leap_day) != 1L) {
    stop("keep_leap_day must be a single logical (TRUE/FALSE).", call. = FALSE)
  }
  if (!is.logical(drop_all_na) || length(drop_all_na) != 1L) {
    stop("drop_all_na must be a single logical (TRUE/FALSE).", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("verbose must be a single logical (TRUE/FALSE).", call. = FALSE)
  }
  if (!is.character(spatial_ref) || length(spatial_ref) != 1L || !nzchar(spatial_ref)) {
    stop("spatial_ref must be a single character string.", call. = FALSE)
  }

  if (!is.null(signif_digits)) {
    if (!is.numeric(signif_digits) || length(signif_digits) != 1L ||
        !is.finite(signif_digits) || signif_digits < 1 || (signif_digits %% 1) != 0) {
      stop("signif_digits must be a positive integer, or NULL.", call. = FALSE)
    }
    signif_digits <- as.integer(signif_digits)
  }

  # var_name: allow NULL or named character mapping old -> new
  rename_map <- NULL
  if (!is.null(var_name)) {
    if (!is.character(var_name) || length(var_name) < 1L) {
      stop("var_name must be a character vector, or NULL.", call. = FALSE)
    }
    nm <- names(var_name)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("var_name must be a *named* character vector: names are original variables.", call. = FALSE)
    }
    if (any(!nzchar(unname(var_name)))) stop("var_name contains empty target names.", call. = FALSE)
    if (any(duplicated(unname(var_name)))) stop("var_name contains duplicate target names.", call. = FALSE)
    rename_map <- var_name
  }

  nc <- ncdf4::nc_open(nc_path)
  on.exit(try(ncdf4::nc_close(nc), silent = TRUE), add = TRUE)

  var_all <- names(nc$var)
  if (length(var_all) == 0L) stop("No variables found in NetCDF.", call. = FALSE)

  # Spatial reference requirement (your tests expect an error if missing)
  if (!is.null(spatial_ref) && nzchar(spatial_ref) && !(spatial_ref %in% var_all)) {
    stop("Spatial reference variable not found: ", spatial_ref, call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  .msg <- function(...) if (isTRUE(verbose)) message(...)

  .mask_missing <- function(nc_obj, var_id, x) {
    fv1 <- ncdf4::ncatt_get(nc_obj, var_id, "_FillValue")
    fv2 <- ncdf4::ncatt_get(nc_obj, var_id, "missing_value")
    fill_vals <- c()

    if (!is.null(fv1) && !isTRUE(fv1$hasatt == FALSE)) fill_vals <- c(fill_vals, fv1$value)
    if (!is.null(fv2) && !isTRUE(fv2$hasatt == FALSE)) fill_vals <- c(fill_vals, fv2$value)

    fill_vals <- unique(as.numeric(fill_vals))
    if (length(fill_vals) > 0L && is.numeric(x)) x[x %in% fill_vals] <- NA_real_
    x
  }

  # Reads the two CF attributes off the time variable and hands them to the
  # file-level decoder, which is shared with write_netcdf()'s encoder so the two
  # directions cannot drift apart.
  .parse_time_to_date <- function(nc_obj, time_var, time_vals) {
    if (is.null(time_var) || is.null(time_vals) || length(time_vals) < 1L) return(NULL)
    if (!is.numeric(time_vals)) return(NULL)

    att <- ncdf4::ncatt_get(nc_obj, time_var, "units")
    if (is.null(att) || isTRUE(att$hasatt == FALSE)) return(NULL)

    cal_att <- ncdf4::ncatt_get(nc_obj, time_var, "calendar")
    calendar <- if (isTRUE(cal_att$hasatt)) as.character(cal_att$value) else NULL

    .cf_time_to_date(
      units     = as.character(att$value),
      calendar  = calendar,
      time_vals = time_vals
    )
  }

  coord_x_alias <- c("lon","longitude","x","rlon")
  coord_y_alias <- c("lat","latitude","y","rlat")
  coord_t_alias <- c("time","Times","t","date")

  # ---------------------------------------------------------------------------
  # Resolve variables to read
  # ---------------------------------------------------------------------------
  coord_like <- tolower(var_all) %in% tolower(c(coord_x_alias, coord_y_alias, coord_t_alias))
  default_var <- var_all[!coord_like]
  default_var <- default_var[default_var != spatial_ref]

  if (is.null(var)) {
    var <- default_var
  } else {
    if (!is.character(var) || length(var) < 1L) {
      stop("var must be a character vector or NULL.", call. = FALSE)
    }
    missing_var <- setdiff(var, var_all)
    if (length(missing_var) > 0L) {
      stop("Variables not found in NetCDF: ", paste(missing_var, collapse = ", "), call. = FALSE)
    }
  }
  if (length(var) == 0L) stop("No data variables selected.", call. = FALSE)

  if (!is.null(rename_map)) {
    bad <- setdiff(names(rename_map), var)
    if (length(bad) > 0L) {
      stop("var_name contains variables not selected: ", paste(bad, collapse = ", "), call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # Determine dimensions from first selected variable
  # ---------------------------------------------------------------------------
  v0 <- nc$var[[var[1]]]
  dim_names <- vapply(v0$dim, function(d) d$name, character(1))
  dim_names_l <- tolower(dim_names)

  time_dim_name <- dim_names[which(dim_names_l %in% tolower(coord_t_alias))[1]]
  if (is.na(time_dim_name) || length(time_dim_name) == 0L) time_dim_name <- NULL

  non_time <- if (!is.null(time_dim_name)) dim_names[dim_names != time_dim_name] else dim_names
  if (length(non_time) < 2L) {
    stop("NetCDF does not expose at least two non-time dimensions.", call. = FALSE)
  }

  x_dim_name <- non_time[which(tolower(non_time) %in% tolower(coord_x_alias))[1]]
  y_dim_name <- non_time[which(tolower(non_time) %in% tolower(coord_y_alias))[1]]
  if (is.na(x_dim_name) || length(x_dim_name) == 0L) x_dim_name <- non_time[1]
  if (is.na(y_dim_name) || length(y_dim_name) == 0L) y_dim_name <- non_time[2]

  get_dim_vals <- function(name) {
    if (name %in% names(nc$dim)) return(nc$dim[[name]]$vals)
    if (name %in% names(nc$var)) return(as.vector(ncdf4::ncvar_get(nc, name)))
    NULL
  }

  x_vals <- get_dim_vals(x_dim_name)
  y_vals <- get_dim_vals(y_dim_name)

  if (is.null(x_vals) || is.null(y_vals) || length(x_vals) < 1L || length(y_vals) < 1L) {
    stop("Could not read spatial coordinate vectors.", call. = FALSE)
  }

  time_vals_raw <- NULL
  if (!is.null(time_dim_name)) time_vals_raw <- get_dim_vals(time_dim_name)
  if (is.null(time_vals_raw)) {
    for (tv in coord_t_alias) {
      if (tv %in% names(nc$var)) {
        time_dim_name <- tv
        time_vals_raw <- as.vector(ncdf4::ncvar_get(nc, tv))
        break
      }
    }
  }
  if (is.null(time_vals_raw) || length(time_vals_raw) < 1L) {
    stop("Could not identify/read time dimension.", call. = FALSE)
  }

  # Parsed exactly once. The date vector and the row mask that drops Feb 29 from
  # every variable must come from the same decode: two parses could disagree
  # about the calendar while still producing matching row counts, which nothing
  # downstream would catch.
  date_full <- .parse_time_to_date(nc, time_dim_name, as.vector(time_vals_raw))
  if (is.null(date_full) || !inherits(date_full, "Date")) {
    stop("Could not convert NetCDF time to Date (missing/unsupported CF units).", call. = FALSE)
  }

  # On a noleap axis this mask is all-TRUE by construction -- Feb 29 cannot be
  # decoded -- so keep_leap_day is simply inert there rather than special-cased.
  leap_keep_idx <- NULL
  date <- date_full
  if (isFALSE(keep_leap_day)) {
    date_lt <- as.POSIXlt(date_full)
    leap_keep_idx <- !(date_lt$mon == 1L & date_lt$mday == 29L)
    date <- date_full[leap_keep_idx]
  }

  # ---------------------------------------------------------------------------
  # Build grid tibble
  # ---------------------------------------------------------------------------
  nx <- length(x_vals)
  ny <- length(y_vals)
  nt <- length(date)  # nolint: object_usage_linter. used in a .log() glue string
  ncell <- nx * ny

  grid_base <- expand.grid(
    xind = seq_len(nx),
    yind = seq_len(ny),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid_base$x <- x_vals[grid_base$xind]
  grid_base$y <- y_vals[grid_base$yind]

  grid <- tibble::as_tibble(grid_base)
  grid <- tibble::add_column(grid, id = seq_len(nrow(grid)), .before = 1)

  # ---------------------------------------------------------------------------
  # Read variables -> matrices [time, cell]
  # ---------------------------------------------------------------------------
  .to_time_by_cell <- function(arr, dim_names, xname, yname, tname, nx, ny, nt) {
    dn <- tolower(dim_names)
    xi <- match(tolower(xname), dn)
    yi <- match(tolower(yname), dn)
    ti <- match(tolower(tname), dn)

    if (is.na(xi)) xi <- which(dn %in% tolower(coord_x_alias))[1]
    if (is.na(yi)) yi <- which(dn %in% tolower(coord_y_alias))[1]
    if (is.na(ti)) ti <- which(dn %in% tolower(coord_t_alias))[1]

    if (is.na(xi) || is.na(yi) || is.na(ti)) {
      stop("Could not map variable dimensions to x/y/time.", call. = FALSE)
    }

    perm <- c(xi, yi, ti)
    extra <- setdiff(seq_along(dim_names), perm)
    if (length(extra) > 0L) {
      extra_sizes <- dim(arr)[extra]
      if (any(extra_sizes != 1L)) {
        stop("Unsupported extra dimensions (not length 1): ", paste(dim_names[extra], collapse = ", "), call. = FALSE)
      }
      perm <- c(perm, extra)
    }

    arr2 <- aperm(arr, perm)
    vec <- as.vector(arr2)

    m <- matrix(vec, nrow = nx * ny, ncol = nt, byrow = FALSE)
    t(m)
  }

  var_mats <- list()
  kept_vars <- character(0)

  # leap_keep_idx was built above from the same decode that produced `date`, so
  # the row mask and the date vector cannot disagree about the calendar.

  for (v in var) {
    .msg("Reading variable: ", v)

    v_obj <- nc$var[[v]]
    v_dim_names <- vapply(v_obj$dim, function(d) d$name, character(1))
    x <- ncdf4::ncvar_get(nc, v)
    x <- .mask_missing(nc, v, x)

    if (is.null(dim(x))) {
      mat <- matrix(as.numeric(x), nrow = length(time_vals_raw), ncol = ncell)
    } else {
      mat <- .to_time_by_cell(x, v_dim_names, x_dim_name, y_dim_name, time_dim_name, nx, ny, length(time_vals_raw))
    }

    # Drop Feb 29 in matrix if we dropped it in date
    if (isFALSE(keep_leap_day)) {
      mat <- mat[leap_keep_idx, , drop = FALSE]
    }

    if (!is.null(signif_digits)) mat <- signif(mat, signif_digits)
    if (isTRUE(drop_all_na) && all(is.na(mat))) next

    out_name <- v
    if (!is.null(rename_map) && v %in% names(rename_map)) out_name <- unname(rename_map[[v]])

    var_mats[[out_name]] <- mat
    kept_vars <- c(kept_vars, out_name)
  }

  data <- vector("list", ncell)
  if (length(kept_vars) == 0L) {
    for (k in seq_len(ncell)) data[[k]] <- data.frame()
  } else {
    var_mats_use <- var_mats[kept_vars]
    n_vars <- length(kept_vars)
    n_time <- nrow(var_mats_use[[1]])
    for (k in seq_len(ncell)) {
      cell_mat <- matrix(NA_real_, nrow = n_time, ncol = n_vars)
      for (v in seq_len(n_vars)) {
        cell_mat[, v] <- as.numeric(var_mats_use[[v]][, k])
      }
      dfk <- as.data.frame(cell_mat, stringsAsFactors = FALSE)
      names(dfk) <- kept_vars
      data[[k]] <- dfk
    }
  }

  dimensions <- list()
  dimensions[[x_dim_name]] <- as.vector(x_vals)
  dimensions[[y_dim_name]] <- as.vector(y_vals)
  dimensions[["time"]] <- as.vector(time_vals_raw)

  global_atts <- list()
  if (length(nc$att) > 0L) {
    for (an in names(nc$att)) global_atts[[an]] <- nc$att[[an]]
  }
  var_atts <- list()
  for (vn in var_all) var_atts[[vn]] <- list()

  attributes <- list(global = global_atts, variables = var_atts)

  list(
    data = data,
    grid = grid,
    date = date,
    dimensions = dimensions,
    attributes = attributes
  )
}

#' Write Gridded Data to a NetCDF File (Template-Based)
#'
#' @description
#' Writes gridded climate/weather time series to a NetCDF file using a template NetCDF
#' as the schema source (dimensions, spatial reference variable, and global attributes).
#'
#' @param data List of length \code{nrow(grid)}. Each element is a named list of
#'   numeric vectors (time series) for each variable to be written. All list elements must
#'   have identical variable names and equal length.
#' @param grid data.frame with at least \code{xind} and \code{yind} integer indices
#'   mapping each list element to an output grid cell.
#' @param out_dir Character. Output directory.
#' @param origin_date Character or Date. \strong{The date of the first time step},
#'   not an arbitrary reference epoch. The time coordinate is written as
#'   \code{0:(nt - 1)}, so the axis is only correct when \code{origin_date} is the
#'   series start -- passing a conventional epoch such as \code{"1970-01-01"} for a
#'   series that begins in 2020 silently relocates it by 50 years. Pass
#'   \code{dates} to have this checked rather than assumed.
#' @param calendar Character. CF calendar for the time axis. One of
#'   \code{"noleap"} (\code{"365_day"}), \code{"standard"}, \code{"gregorian"},
#'   \code{"proleptic_gregorian"}, or \code{"julian"}. Must match the calendar the
#'   data are actually on: \code{\link{generate_weather}} output is 365-day, so
#'   \code{"noleap"} is correct for it, and mislabelling it shifts every date for
#'   any CF-aware reader.
#' @param dates Optional Date vector, one entry per time step, giving the axis the
#'   data are on. When supplied it is validated against \code{nt},
#'   \code{origin_date}, and \code{calendar}, and the time coordinate is computed
#'   from it instead of being assumed contiguous. When \code{NULL} (default) the
#'   axis is assumed to be \code{nt} consecutive days starting at
#'   \code{origin_date}.
#' @param template_path Character. Path to a template NetCDF file.
#' @param compression Integer 0-9. NetCDF4 deflation level.
#' @param spatial_ref Character. Spatial reference variable name in template.
#' @param file_prefix Character. Prefix for output filename.
#' @param file_suffix Character. Optional suffix appended to filename.
#' @param signif_digits Integer. If not \code{NULL}, round values to this many significant digits.
#' @param verbose Logical. If TRUE, emit progress logs via \code{.log()}.
#'
#' @return Invisibly returns the written file path.
#'
#' @export
write_netcdf <- function(
    data = NULL,
    grid = NULL,
    out_dir = NULL,
    origin_date = NULL,
    calendar = "noleap",
    dates = NULL,
    template_path = NULL,
    compression = 4,
    spatial_ref = "spatial_ref",
    file_prefix = "clim_change_rlz",
    file_suffix = "",
    signif_digits = NULL,
    verbose = TRUE
) {

  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  cal_kind <- .cf_calendar_kind(calendar)
  if (!identical(cal_kind, "standard") && !identical(cal_kind, "noleap")) {
    stop(
      "'calendar' must be one of 'noleap', '365_day', 'standard', 'gregorian', ",
      "'proleptic_gregorian' or 'julian'; got '", calendar, "'.",
      call. = FALSE
    )
  }

  tryCatch(
    as.Date(origin_date),
    error = function(e) {
      stop(
        "'origin_date' must be a valid date string (YYYY-MM-DD). Error: ",
        e$message,
        call. = FALSE
      )
    }
  )

  if (is.null(template_path)) {
    stop("'template_path' must be provided.", call. = FALSE)
  }
  if (!file.exists(template_path)) {
    stop("'template_path' does not exist: ", template_path, call. = FALSE)
  }

  if (!is.numeric(compression) || length(compression) != 1L ||
      !is.finite(compression) || compression < 0 || compression > 9) {
    stop("'compression' must be an integer between 0 and 9.", call. = FALSE)
  }
  compression <- as.integer(compression)

  if (is.null(data) || !is.list(data) || length(data) == 0L) {
    stop("'data' must be a non-empty list.", call. = FALSE)
  }
  if (!all(vapply(data, is.list, logical(1)))) {
    stop("'data' must be a list of lists (per-grid payloads).", call. = FALSE)
  }

  if (is.null(grid) || !is.data.frame(grid)) {
    stop("'grid' must be a data.frame.", call. = FALSE)
  }
  if (!all(c("xind", "yind") %in% names(grid))) {
    stop("'grid' must contain columns 'xind' and 'yind'.", call. = FALSE)
  }
  if (nrow(grid) != length(data)) {
    stop("nrow('grid') must match length('data').", call. = FALSE)
  }

  if (is.null(out_dir) || !is.character(out_dir) || length(out_dir) != 1L || !nzchar(out_dir)) {
    stop("'out_dir' must be a non-empty character scalar.", call. = FALSE)
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (!is.null(signif_digits)) {
    if (!is.numeric(signif_digits) || length(signif_digits) != 1L ||
        !is.finite(signif_digits) || signif_digits < 1 || (signif_digits %% 1) != 0) {
      stop("'signif_digits' must be a positive integer, or NULL.", call. = FALSE)
    }
    signif_digits <- as.integer(signif_digits)
  }

  # Determine variables and nt
  var <- names(data[[1]])
  if (length(var) == 0L) stop("No variables found in 'data[[1]]'.", call. = FALSE)

  nt <- length(data[[1]][[1]])
  if (nt < 1L) stop("Time series length 'nt' must be >= 1.", call. = FALSE)

  # Validate all grids consistent
  for (i in seq_along(data)) {
    if (!identical(names(data[[i]]), var)) {
      stop("All grid payloads in 'data' must have identical variable names.", call. = FALSE)
    }
    lens <- vapply(data[[i]], length, integer(1L))
    if (any(lens != nt)) {
      stop("All variable series must have identical length within each grid.", call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # Open template and derive schema
  # ---------------------------------------------------------------------------
  nc_in <- ncdf4::nc_open(template_path)
  on.exit(try(ncdf4::nc_close(nc_in), silent = TRUE), add = TRUE)

  if (!(spatial_ref %in% names(nc_in$var))) {
    stop("Spatial reference variable not found in template: ", spatial_ref, call. = FALSE)
  }

  x_att <- ncdf4::ncatt_get(nc_in, spatial_ref, "x_dim")
  y_att <- ncdf4::ncatt_get(nc_in, spatial_ref, "y_dim")
  if (!isTRUE(x_att$hasatt) || !isTRUE(y_att$hasatt)) {
    stop("Template spatial_ref must have attributes 'x_dim' and 'y_dim'.", call. = FALSE)
  }
  x_dim_name <- x_att$value
  y_dim_name <- y_att$value

  x_vals <- nc_in$dim[[x_dim_name]]$vals
  y_vals <- nc_in$dim[[y_dim_name]]$vals

  nx <- length(x_vals)
  ny <- length(y_vals)

  if (nx < 1L || ny < 1L) stop("Template x/y dimensions are empty.", call. = FALSE)

  .log(
    msg = sprintf(
      "Template loaded: nx=%s, ny=%s, nt=%s",
      format(nx, big.mark = ","),
      format(ny, big.mark = ","),
      format(nt, big.mark = ",")
    ),
    verbose = verbose,
    tag = "IO"
  )

  # ---------------------------------------------------------------------------
  # Output filename
  # ---------------------------------------------------------------------------
  suffix <- if (is.character(file_suffix) && nzchar(file_suffix)) paste0("_", file_suffix) else ""
  nc_file_path <- file.path(out_dir, paste0(file_prefix, suffix, ".nc"))

  # ---------------------------------------------------------------------------
  # Define dimensions
  # ---------------------------------------------------------------------------
  origin <- as.Date(origin_date)
  time_units <- paste0("days since ", format(origin, "%Y-%m-%d"), " 00:00:00")

  # Without `dates`, the axis is assumed to be nt consecutive days from the
  # origin -- correct for everything this package produces, but unverifiable
  # here. With `dates`, the offsets are computed on the declared calendar
  # instead, which is what turns a mislabelled origin into an error.
  if (is.null(dates)) {
    time_vals <- 0:(nt - 1)
  } else {
    if (!inherits(dates, "Date")) {
      stop("'dates' must be a Date vector, or NULL.", call. = FALSE)
    }
    if (length(dates) != nt) {
      stop("'dates' must have one entry per time step: got ", length(dates),
           " for nt = ", nt, ".", call. = FALSE)
    }
    if (anyNA(dates)) {
      stop("'dates' must not contain NA.", call. = FALSE)
    }
    # Value comparison, not identical(): origin_date may arrive as a character,
    # a Date, or a POSIXct-derived Date, and identical() would reject an equal
    # date carrying a different attribute set.
    if (!isTRUE(dates[1] == origin)) {
      stop("'origin_date' (", format(origin), ") must be the first entry of 'dates' (",
           format(dates[1]), "); the time coordinate is written relative to it.",
           call. = FALSE)
    }

    time_vals <- if (identical(cal_kind, "noleap")) {
      .noleap_date_to_offset(origin, dates)
    } else {
      as.numeric(dates - origin)
    }

    if (is.unsorted(time_vals, strictly = TRUE)) {
      stop("'dates' must be strictly increasing on the '", calendar, "' calendar.",
           call. = FALSE)
    }
  }

  dim_time <- ncdf4::ncdim_def(
    name = "time",
    units = time_units,
    vals = time_vals,
    unlim = FALSE
  )

  dim_y <- ncdf4::ncdim_def(name = y_dim_name, units = "", vals = y_vals)
  dim_x <- ncdf4::ncdim_def(name = x_dim_name, units = "", vals = x_vals)

  # ---------------------------------------------------------------------------
  # Define variables (time x y x x)
  # ---------------------------------------------------------------------------
  var_defs <- vector("list", length(var))
  names(var_defs) <- var

  for (v in var) {

    v_units <- ""
    v_long  <- v

    if (v %in% names(nc_in$var)) {
      att_u <- ncdf4::ncatt_get(nc_in, v, "units")$value
      att_l <- ncdf4::ncatt_get(nc_in, v, "long_name")$value
      if (!is.null(att_u) && !is.na(att_u)) v_units <- att_u
      if (!is.null(att_l) && !is.na(att_l)) v_long  <- att_l
    }

    var_defs[[v]] <- ncdf4::ncvar_def(
      name = v,
      units = v_units,
      dim = list(dim_time, dim_y, dim_x),
      missval = NA_real_,
      longname = v_long,
      prec = "float",
      compression = compression
    )
  }

  # ---------------------------------------------------------------------------
  # Define spatial reference variable (copied from template)
  # ---------------------------------------------------------------------------
  sr_template <- nc_in$var[[spatial_ref]]

  sr_prec <- tolower(as.character(sr_template$prec))
  sr_prec <- switch(sr_prec,
                    "int"    = "integer",
                    "int32"  = "integer",
                    "int16"  = "short",
                    "float32"= "float",
                    "float64"= "double",
                    sr_prec)

  if (is.null(sr_prec) || is.na(sr_prec) || !nzchar(sr_prec) ||
      !sr_prec %in% c("short","float","double","integer","char","byte")) {
    sr_prec <- "integer"
  }

  spatial_ref_def <- ncdf4::ncvar_def(
    name = spatial_ref,
    units = "",
    dim = list(),
    missval = NULL,
    longname = spatial_ref,
    prec = sr_prec
  )

  # ---------------------------------------------------------------------------
  # Create file
  # ---------------------------------------------------------------------------
  nc_out <- ncdf4::nc_create(
    nc_file_path,
    vars = c(list(spatial_ref_def), unname(var_defs)),
    force_v4 = TRUE
  )
  on.exit(try(ncdf4::nc_close(nc_out), silent = TRUE), add = TRUE)

  # ---------------------------------------------------------------------------
  # Write spatial reference variable value + attributes from template
  # ---------------------------------------------------------------------------
  sr_val <- ncdf4::ncvar_get(nc_in, spatial_ref)
  try(ncdf4::ncvar_put(nc_out, spatial_ref, sr_val), silent = TRUE)

  sr_atts <- ncdf4::ncatt_get(nc_in, spatial_ref)
  if (length(sr_atts) > 0L) {
    for (an in names(sr_atts)) {
      val <- sr_atts[[an]]
      if (is.atomic(val)) {
        ncdf4::ncatt_put(nc_out, spatial_ref, an, val)
      }
    }
  }

  # Copy global attributes from template when possible
  if (!is.null(nc_in$gatts) && length(nc_in$gatts) > 0L) {
    for (an in names(nc_in$gatts)) {
      val <- nc_in$gatts[[an]]
      if (is.atomic(val)) try(ncdf4::ncatt_put(nc_out, 0, an, val), silent = TRUE)
    }
  }

  # Calendar attribute on time. Not wrapped in try(): every reader depends on
  # this attribute to decode the axis, so a file that silently lost it would be
  # off by a day per leap year for anyone who read it back.
  ncdf4::ncatt_put(nc_out, "time", "calendar", calendar)

  # ---------------------------------------------------------------------------
  # Write data cube per variable (time x y x x)
  # ---------------------------------------------------------------------------
  n_grids <- length(data)

  for (v in var) {

    arr <- array(NA_real_, dim = c(nt, ny, nx))

    for (i in seq_len(n_grids)) {
      xi <- grid$xind[i]
      yi <- grid$yind[i]

      if (!is.finite(xi) || !is.finite(yi)) next
      xi <- as.integer(xi); yi <- as.integer(yi)

      if (xi < 1L || xi > nx || yi < 1L || yi > ny) next

      series <- as.numeric(data[[i]][[v]])
      if (!is.null(signif_digits)) {
        series <- signif(series, digits = signif_digits)
      }

      arr[, yi, xi] <- series
    }

    ncdf4::ncvar_put(nc_out, v, arr)
  }

  .log(
    msg = sprintf(
      "NetCDF written: %s | grids=%s", nc_file_path,
      paste(var, collapse = ", "),
      format(nt, big.mark = ","),
      format(ny, big.mark = ","),
      format(nx, big.mark = ","),
      format(n_grids, big.mark = ",")
    ),
    verbose = verbose,
    tag = "IO"
  )

  invisible(nc_file_path)
}
