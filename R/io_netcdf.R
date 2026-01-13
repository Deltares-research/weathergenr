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
#' @param nc.file Character. Path to NetCDF file.
#' @param variables Character vector or NULL. Variables to read. If NULL, reads all
#'   non-coordinate variables (and excludes `spatial.ref` if present).
#' @param var_rename NULL, or named character vector, or character vector.
#'   - If named: names are original variable names and values are new names.
#'   - If unnamed: must be same length/order as `variables` and will be treated as new names.
#'   Renaming applies to output list names and the value-column name.
#' @param leap.days Logical. If FALSE and time is convertible to Date, removes Feb 29 rows.
#' @param omit.empty Logical. If TRUE, drop variables whose values are all NA
#'   after applying _FillValue/missing_value handling.
#' @param spatial.ref Character. Name of a spatial reference variable to ignore
#'   when auto-selecting variables (common in CF/CRS-encoded files).
#' @param signif.digits Integer or NULL. If provided, rounds numeric values via
#'   signif(x, signif.digits).
#' @param verbose Logical. If TRUE, prints basic progress messages.
#'
#' @return Named list of data.frames.
#'
#' @importFrom ncdf4 nc_open nc_close ncvar_get ncatt_get
#' @export
read_netcdf <- function(
    nc.file,
    variables     = NULL,
    var_rename    = NULL,
    leap.days     = TRUE,
    omit.empty    = TRUE,
    spatial.ref   = "spatial_ref",
    signif.digits = NULL,
    verbose       = FALSE
) {

  # ---------------------------------------------------------------------------
  # Input validation (match test expectations / regexes)
  # ---------------------------------------------------------------------------
  if (!is.character(nc.file) || length(nc.file) != 1L || !nzchar(nc.file)) {
    stop("nc.file must be a non-empty character path.", call. = FALSE)
  }
  if (!file.exists(nc.file)) {
    stop("File does not exist: ", nc.file, call. = FALSE)
  }

  if (!is.logical(leap.days) || length(leap.days) != 1L) {
    stop("leap.days must be a single logical (TRUE/FALSE).", call. = FALSE)
  }

  if (!is.logical(omit.empty) || length(omit.empty) != 1L) {
    stop("omit.empty must be a single logical (TRUE/FALSE).", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("verbose must be a single logical (TRUE/FALSE).", call. = FALSE)
  }
  if (!is.character(spatial.ref) || length(spatial.ref) != 1L || !nzchar(spatial.ref)) {
    stop("spatial.ref must be a single character string.", call. = FALSE)
  }
  if (!is.null(signif.digits)) {
    if (!is.numeric(signif.digits) || length(signif.digits) != 1L ||
        !is.finite(signif.digits) || signif.digits < 1 || (signif.digits %% 1) != 0) {
      stop("signif.digits must be a positive integer, or NULL.", call. = FALSE)
    }
    signif.digits <- as.integer(signif.digits)
  }

  # var_rename: allow NULL or named character mapping old -> new
  rename_map <- NULL
  if (!is.null(var_rename)) {
    if (!is.character(var_rename) || length(var_rename) < 1L) {
      stop("var_rename must be a character vector, or NULL.", call. = FALSE)
    }
    nm <- names(var_rename)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("var_rename must be a *named* character vector: names are original variables.", call. = FALSE)
    }
    if (any(!nzchar(unname(var_rename)))) stop("var_rename contains empty target names.", call. = FALSE)
    if (any(duplicated(unname(var_rename)))) stop("var_rename contains duplicate target names.", call. = FALSE)
    rename_map <- var_rename
  }

  nc <- ncdf4::nc_open(nc.file)
  on.exit(try(ncdf4::nc_close(nc), silent = TRUE), add = TRUE)

  v_all <- names(nc$var)
  if (length(v_all) == 0L) stop("No variables found in NetCDF.", call. = FALSE)

  # Spatial reference requirement (your tests expect an error if missing)
  if (!is.null(spatial.ref) && nzchar(spatial.ref) && !(spatial.ref %in% v_all)) {
    stop("Spatial reference variable not found: ", spatial.ref, call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  .msg <- function(...) if (isTRUE(verbose)) message(...)

  .mask_missing <- function(nc, varname, x) {
    fv1 <- ncdf4::ncatt_get(nc, varname, "_FillValue")
    fv2 <- ncdf4::ncatt_get(nc, varname, "missing_value")
    fill_vals <- c()

    if (!is.null(fv1) && !isTRUE(fv1$hasatt == FALSE)) fill_vals <- c(fill_vals, fv1$value)
    if (!is.null(fv2) && !isTRUE(fv2$hasatt == FALSE)) fill_vals <- c(fill_vals, fv2$value)

    fill_vals <- unique(as.numeric(fill_vals))
    if (length(fill_vals) > 0L && is.numeric(x)) x[x %in% fill_vals] <- NA_real_
    x
  }

  .parse_time_to_date <- function(nc, time_var, time_vals) {
    if (is.null(time_var) || is.null(time_vals) || length(time_vals) < 1L) return(NULL)
    if (!is.numeric(time_vals)) return(NULL)

    att <- ncdf4::ncatt_get(nc, time_var, "units")
    if (is.null(att) || isTRUE(att$hasatt == FALSE)) return(NULL)
    units <- as.character(att$value)
    if (!nzchar(units)) return(NULL)

    m <- regexec("^\\s*(seconds|second|minutes|minute|hours|hour|days|day)\\s+since\\s+(.+?)\\s*$",
                 units, ignore.case = TRUE)
    r <- regmatches(units, m)[[1]]
    if (length(r) != 3L) return(NULL)

    u <- tolower(r[2])
    origin_str <- r[3]

    origin_posix <- suppressWarnings(as.POSIXct(origin_str, tz = "UTC"))
    if (is.na(origin_posix)) {
      origin_date <- suppressWarnings(as.Date(origin_str))
      if (is.na(origin_date)) return(NULL)
      origin_posix <- as.POSIXct(origin_date, tz = "UTC")
    }

    mult <- switch(u,
                   "second"  = 1, "seconds" = 1,
                   "minute"  = 60, "minutes" = 60,
                   "hour"    = 3600, "hours" = 3600,
                   "day"     = 86400, "days" = 86400,
                   NULL
    )
    if (is.null(mult)) return(NULL)

    as.Date(origin_posix + as.numeric(time_vals) * mult)
  }

  coord_x_alias <- c("lon","longitude","x","rlon")
  coord_y_alias <- c("lat","latitude","y","rlat")
  coord_t_alias <- c("time","Times","t","date")

  # ---------------------------------------------------------------------------
  # Resolve variables to read
  # ---------------------------------------------------------------------------
  coord_like <- tolower(v_all) %in% tolower(c(coord_x_alias, coord_y_alias, coord_t_alias))
  default_vars <- v_all[!coord_like]
  default_vars <- default_vars[default_vars != spatial.ref]

  if (is.null(variables)) {
    variables <- default_vars
  } else {
    if (!is.character(variables) || length(variables) < 1L) {
      stop("variables must be a character vector or NULL.", call. = FALSE)
    }
    missing_vars <- setdiff(variables, v_all)
    if (length(missing_vars) > 0L) {
      stop("Variables not found in NetCDF: ", paste(missing_vars, collapse = ", "), call. = FALSE)
    }
  }

  if (length(variables) == 0L) stop("No data variables selected.", call. = FALSE)

  if (!is.null(rename_map)) {
    bad <- setdiff(names(rename_map), variables)
    if (length(bad) > 0L) {
      stop("var_rename contains variables not selected: ", paste(bad, collapse = ", "), call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # Determine dimensions from first selected variable
  # ---------------------------------------------------------------------------
  v0 <- nc$var[[variables[1]]]
  dnames <- vapply(v0$dim, function(d) d$name, character(1))
  dnames_l <- tolower(dnames)

  time_dim_name <- dnames[which(dnames_l %in% tolower(coord_t_alias))[1]]
  if (is.na(time_dim_name) || length(time_dim_name) == 0L) time_dim_name <- NULL

  non_time <- if (!is.null(time_dim_name)) dnames[dnames != time_dim_name] else dnames
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

  date <- .parse_time_to_date(nc, time_dim_name, as.vector(time_vals_raw))
  if (is.null(date) || !inherits(date, "Date")) {
    stop("Could not convert NetCDF time to Date (missing/unsupported CF units).", call. = FALSE)
  }

  if (isFALSE(leap.days) && any(format(date, "%m-%d") == "02-29")) {
    warning("leap.days = FALSE specified but data contains Feb 29; returning full series unchanged.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Build grid tibble
  # ---------------------------------------------------------------------------
  nx <- length(x_vals)
  ny <- length(y_vals)
  nt <- length(date)
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
  .to_time_by_cell <- function(arr, dim_names, xname, yname, tname) {
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

  for (v in variables) {
    .msg("Reading variable: ", v)

    v_obj <- nc$var[[v]]
    v_dim_names <- vapply(v_obj$dim, function(d) d$name, character(1))
    x <- ncdf4::ncvar_get(nc, v)
    x <- .mask_missing(nc, v, x)

    if (is.null(dim(x))) {
      mat <- matrix(as.numeric(x), nrow = nt, ncol = ncell)
    } else {
      mat <- .to_time_by_cell(x, v_dim_names, x_dim_name, y_dim_name, time_dim_name)
    }

    if (!is.null(signif.digits)) mat <- signif(mat, signif.digits)

    if (isTRUE(omit.empty) && all(is.na(mat))) next

    out_name <- v
    if (!is.null(rename_map) && v %in% names(rename_map)) out_name <- unname(rename_map[[v]])

    var_mats[[out_name]] <- mat
    kept_vars <- c(kept_vars, out_name)
  }

  if (length(var_mats) == 0L) kept_vars <- character(0)

  data <- vector("list", ncell)
  if (length(kept_vars) == 0L) {
    for (k in seq_len(ncell)) data[[k]] <- data.frame()
  } else {
    for (k in seq_len(ncell)) {
      dfk <- as.data.frame(lapply(var_mats, function(m) as.numeric(m[, k])),
                           stringsAsFactors = FALSE)
      dfk <- dfk[, kept_vars, drop = FALSE]
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
  for (vn in v_all) var_atts[[vn]] <- list()

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
#' @param data List of length \code{nrow(coord.grid)}. Each element is a named list of
#'   numeric vectors (time series) for each variable to be written. All list elements must
#'   have identical variable names and equal length.
#' @param coord.grid data.frame with at least \code{xind} and \code{yind} integer indices
#'   mapping each list element to an output grid cell.
#' @param output.path Character. Output directory.
#' @param origin.date Character or Date. Origin date used for NetCDF time units
#'   (e.g., \code{"1970-01-01"}).
#' @param calendar.type Character. NetCDF calendar type (e.g., \code{"noleap"}).
#' @param nc.template.file Character. Path to a template NetCDF file.
#' @param nc.compression Integer 0-9. NetCDF4 deflation level.
#' @param nc.spatial.ref Character. Spatial reference variable name in template.
#' @param nc.file.prefix Character. Prefix for output filename.
#' @param nc.file.suffix Character. Optional suffix appended to filename.
#' @param signif.digits Integer. If not \code{NULL}, round values to this many significant digits.
#' @param verbose Logical. If TRUE, emit progress logs via \code{.log_info()}.
#'
#' @return Invisibly returns the written file path.
#'
#' @export
write_netcdf <- function(data = NULL,
                         coord.grid = NULL,
                         output.path = NULL,
                         origin.date = NULL,
                         calendar.type = "noleap",
                         nc.template.file = NULL,
                         nc.compression = 4,
                         nc.spatial.ref = "spatial_ref",
                         nc.file.prefix = "clim_change_rlz",
                         nc.file.suffix = "",
                         signif.digits = NULL,
                         verbose = TRUE) {

  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  tryCatch(
    as.Date(origin.date),
    error = function(e) {
      stop(
        "'origin.date' must be a valid date string (YYYY-MM-DD). ",
        "Error: ", e$message,
        call. = FALSE
      )
    }
  )

  if (is.null(nc.template.file)) {
    stop("'nc.template.file' must be provided.", call. = FALSE)
  }
  if (!file.exists(nc.template.file)) {
    stop("'nc.template.file' does not exist: ", nc.template.file, call. = FALSE)
  }
  if (!is.numeric(nc.compression) || nc.compression < 0 || nc.compression > 9) {
    stop("'nc.compression' must be an integer between 0 and 9.", call. = FALSE)
  }

  if (is.null(data) || !is.list(data) || length(data) == 0L) {
    stop("'data' must be a non-empty list.", call. = FALSE)
  }
  if (!all(vapply(data, is.list, logical(1)))) {
    stop("'data' must be a list of lists (per-grid payloads).", call. = FALSE)
  }

  if (is.null(coord.grid) || !is.data.frame(coord.grid)) {
    stop("'coord.grid' must be a data.frame.", call. = FALSE)
  }
  if (!all(c("xind", "yind") %in% names(coord.grid))) {
    stop("'coord.grid' must contain columns 'xind' and 'yind'.", call. = FALSE)
  }
  if (nrow(coord.grid) != length(data)) {
    stop("nrow('coord.grid') must match length('data').", call. = FALSE)
  }

  if (is.null(output.path) || !is.character(output.path) || length(output.path) != 1L) {
    stop("'output.path' must be a character scalar.", call. = FALSE)
  }
  if (!dir.exists(output.path)) {
    dir.create(output.path, recursive = TRUE, showWarnings = FALSE)
  }

  # Determine variables and nt
  variables <- names(data[[1]])
  if (length(variables) == 0L) stop("No variables found in 'data[[1]]'.", call. = FALSE)

  nt <- length(data[[1]][[1]])
  if (nt < 1L) stop("Time series length 'nt' must be >= 1.", call. = FALSE)

  # Validate all grids consistent
  for (i in seq_along(data)) {
    if (!identical(names(data[[i]]), variables)) {
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
  nc_in <- ncdf4::nc_open(nc.template.file)
  on.exit(try(ncdf4::nc_close(nc_in), silent = TRUE), add = TRUE)

  if (!(nc.spatial.ref %in% names(nc_in$var))) {
    stop("Spatial reference variable not found in template: ", nc.spatial.ref, call. = FALSE)
  }

  x_dim_name <- ncdf4::ncatt_get(nc_in, nc.spatial.ref, "x_dim")$value
  y_dim_name <- ncdf4::ncatt_get(nc_in, nc.spatial.ref, "y_dim")$value

  if (is.null(x_dim_name) || is.null(y_dim_name) || anyNA(c(x_dim_name, y_dim_name))) {
    stop("Template spatial_ref must have attributes 'x_dim' and 'y_dim'.", call. = FALSE)
  }

  x_vals <- nc_in$dim[[x_dim_name]]$vals
  y_vals <- nc_in$dim[[y_dim_name]]$vals

  nx <- length(x_vals)
  ny <- length(y_vals)

  if (nx < 1L || ny < 1L) stop("Template x/y dimensions are empty.", call. = FALSE)

  .log_info(
    msg = sprintf("Template loaded: nx=%d, ny=%d, nt=%d", nx, ny, nt),
    verbose = verbose,
    tag = "IO"
  )

  # ---------------------------------------------------------------------------
  # Output filename
  # ---------------------------------------------------------------------------
  suffix <- if (is.character(nc.file.suffix) && nzchar(nc.file.suffix)) paste0("_", nc.file.suffix) else ""
  nc_file_path <- file.path(output.path, paste0(nc.file.prefix, suffix, ".nc"))

  # ---------------------------------------------------------------------------
  # Define dimensions
  # ---------------------------------------------------------------------------
  time_units <- paste0("days since ", format(as.Date(origin.date), "%Y-%m-%d"), " 00:00:00")

  dim_time <- ncdf4::ncdim_def(
    name = "time",
    units = time_units,
    vals = 0:(nt - 1),
    unlim = FALSE
  )

  dim_y <- ncdf4::ncdim_def(name = y_dim_name, units = "", vals = y_vals)
  dim_x <- ncdf4::ncdim_def(name = x_dim_name, units = "", vals = x_vals)

  # ---------------------------------------------------------------------------
  # Define variables (time x y x x)
  # ---------------------------------------------------------------------------
  # Use template variable metadata when present; otherwise fallback defaults.
  var_defs <- vector("list", length(variables))
  names(var_defs) <- variables

  for (v in variables) {

    # Attempt to inherit units/long_name if exists in template
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
      compression = nc.compression
    )
  }

  # ---------------------------------------------------------------------------
  # Define spatial reference variable (copied from template)
  # ---------------------------------------------------------------------------
  sr_template <- nc_in$var[[nc.spatial.ref]]

  sr_prec <- sr_template$prec
  sr_prec <- tolower(as.character(sr_prec))

  # Map common aliases to ncdf4 accepted precision strings
  sr_prec <- switch(sr_prec,
                    "int"    = "integer",
                    "int32"  = "integer",
                    "int16"  = "short",
                    "float32"= "float",
                    "float64"= "double",
                    sr_prec
  )

  if (is.null(sr_prec) || is.na(sr_prec) || !nzchar(sr_prec) ||
      !sr_prec %in% c("short","float","double","integer","char","byte")) {
    sr_prec <- "integer"
  }

  spatial_ref_def <- ncdf4::ncvar_def(
    name = nc.spatial.ref,
    units = "",
    dim = list(),        # scalar
    missval = NULL,      # IMPORTANT: do not create _FillValue for CRS var
    longname = nc.spatial.ref,
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
  sr_val <- ncdf4::ncvar_get(nc_in, nc.spatial.ref)
  try(ncdf4::ncvar_put(nc_out, nc.spatial.ref, sr_val), silent = TRUE)

  sr_atts <- nc_in$var[[nc.spatial.ref]]$att
  if (!is.null(sr_atts) && length(sr_atts) > 0L) {
    for (an in names(sr_atts)) {
      val <- sr_atts[[an]]
      if (is.atomic(val)) {
        try(ncdf4::ncatt_put(nc_out, nc.spatial.ref, an, val), silent = TRUE)
      }
    }
  }


  # Copy global attributes from template when possible
  if (!is.null(nc_in$gatts) && length(nc_in$gatts) > 0L) {
    for (an in names(nc_in$gatts)) {
      val <- nc_in$gatts[[an]]
      # ncatt_put expects atomic vectors
      if (is.atomic(val)) {
        try(ncdf4::ncatt_put(nc_out, 0, an, val), silent = TRUE)
      }
    }
  }

  # Calendar attribute on time
  try(ncdf4::ncatt_put(nc_out, "time", "calendar", calendar.type), silent = TRUE)

  # ---------------------------------------------------------------------------
  # Write data cube per variable (time x y x x)
  # ---------------------------------------------------------------------------
  n_grids <- length(data)

  for (v in variables) {

    arr <- array(NA_real_, dim = c(nt, ny, nx))

    for (i in seq_len(n_grids)) {
      xi <- coord.grid$xind[i]
      yi <- coord.grid$yind[i]

      if (!is.finite(xi) || !is.finite(yi)) next
      xi <- as.integer(xi); yi <- as.integer(yi)

      if (xi < 1L || xi > nx || yi < 1L || yi > ny) next

      series <- as.numeric(data[[i]][[v]])
      if (!is.null(signif.digits)) {
        series <- signif(series, digits = signif.digits)
      }

      arr[, yi, xi] <- series
    }

    ncdf4::ncvar_put(nc_out, v, arr)
  }

  .log_info(
    msg = sprintf(
      "NetCDF written: %s | vars=%s | dims=%d time x %d y x %d x | grids=%d",
      nc_file_path, paste(variables, collapse = ", "), nt, ny, nx, n_grids
    ),
    verbose = verbose,
    tag = "IO"
  )

  invisible(nc_file_path)
}
