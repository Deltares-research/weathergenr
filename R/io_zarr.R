#' Read Zarr variables into tidy data frames
#'
#' Reads one or more variables from a Zarr store and returns a structured list
#' containing data, grid, date, dimensions, and attributes - matching the output
#' structure of read_netcdf().
#'
#' @param dir Character. Path to Zarr store (directory).
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
#' @param slice_mode Integer. Indexing mode for pizzarr: 1 for R-style (1-based),
#'   0 for Python-style (0-based). Default is 1.
#' @param verbose Logical. If TRUE, prints basic progress messages.
#'
#' @return List with components:
#'   \item{data}{List of data.frames (one per grid cell)}
#'   \item{grid}{Tibble with id, xind, yind, x, y columns}
#'   \item{date}{Vector of Date values}
#'   \item{dimensions}{List of dimension coordinate vectors}
#'   \item{attributes}{List with global and variable attributes}
#'
#' @importFrom pizzarr zarr_open
#' @export
read_zarr <- function(
    dir,
    var = NULL,
    var_name = NULL,
    keep_leap_day = TRUE,
    drop_all_na = TRUE,
    spatial_ref = "spatial_ref",
    signif_digits = NULL,
    slice_mode = 1,
    verbose = FALSE
) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.character(dir) || length(dir) != 1L || !nzchar(dir)) {
    stop("dir must be a non-empty character path.", call. = FALSE)
  }
  if (!dir.exists(dir)) {
    stop("Directory does not exist: ", dir, call. = FALSE)
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

  if (!is.numeric(slice_mode) || length(slice_mode) != 1L || !slice_mode %in% c(0, 1)) {
    stop("slice_mode must be either 0 (Python-style, 0-based) or 1 (R-style, 1-based).", call. = FALSE)
  }
  slice_mode <- as.integer(slice_mode)

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

  # ---------------------------------------------------------------------------
  # Load pizzarr and open zarr store with specified slice mode
  # ---------------------------------------------------------------------------
  if (!requireNamespace("pizzarr", quietly = TRUE)) {
    stop("Package 'pizzarr' is required. Install with: install.packages('pizzarr')", call. = FALSE)
  }

  # Set the slice mode for pizzarr
  old_slice_mode <- getOption("pizzarr.slice_mode")
  on.exit(options(pizzarr.slice_mode = old_slice_mode), add = TRUE)
  options(pizzarr.slice_mode = slice_mode)

  if (verbose) {
    message("Using ", ifelse(slice_mode == 1, "R-style (1-based)", "Python-style (0-based)"), " indexing")
  }

  zgroup <- pizzarr::zarr_open(dir)
  child_dirs <- list.dirs(dir, full.names = FALSE, recursive = FALSE)
  var_all <- child_dirs[file.exists(file.path(dir, child_dirs, ".zarray"))]

  if (length(var_all) == 0L) stop("No variables found in Zarr store.", call. = FALSE)

  # Spatial reference requirement:
  # allow missing default "spatial_ref", but error for explicit non-default requests.
  if (!is.null(spatial_ref) &&
      nzchar(spatial_ref) &&
      !(spatial_ref %in% var_all) &&
      !identical(spatial_ref, "spatial_ref")) {
    stop("Spatial reference variable not found: ", spatial_ref, call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  .msg <- function(...) if (isTRUE(verbose)) message(...)

  .get_zarr_attrs <- function(zobj) {
    attrs <- zobj$get_attrs()
    if (inherits(attrs, "Attributes")) {
      attrs <- attrs$to_list()
    }
    if (is.null(attrs)) return(list())
    if (length(attrs) == 0L) return(list())
    attrs
  }

  .get_item_all <- function(zobj) {
    out <- tryCatch(
      zobj$get_item("..."),
      error = function(e) zobj$get_item()
    )
    if (inherits(out, "NestedArray")) {
      out <- out$as.array()
    }
    out
  }

  .mask_missing <- function(zarray, x) {
    attrs <- .get_zarr_attrs(zarray)
    fill_vals <- c()

    if ("_FillValue" %in% names(attrs)) fill_vals <- c(fill_vals, attrs$`_FillValue`)
    if ("missing_value" %in% names(attrs)) fill_vals <- c(fill_vals, attrs$missing_value)

    fill_vals <- unique(as.numeric(fill_vals))
    if (length(fill_vals) > 0L && is.numeric(x)) x[x %in% fill_vals] <- NA_real_
    x
  }

  .parse_time_to_date <- function(zarray, time_vals) {
    if (is.null(time_vals) || length(time_vals) < 1L) return(NULL)
    if (!is.numeric(time_vals)) return(NULL)

    attrs <- .get_zarr_attrs(zarray)
    if (!"units" %in% names(attrs)) return(NULL)

    units <- as.character(attrs$units)
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
  .find_coord_name <- function(all_names, aliases) {
    idx <- which(tolower(all_names) %in% tolower(aliases))
    if (length(idx) == 0L) return(NULL)
    all_names[idx[1]]
  }

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
      stop("Variables not found in Zarr store: ", paste(missing_var, collapse = ", "), call. = FALSE)
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
  v0_array <- zgroup$get_item(var[1])
  dim_info <- v0_array$get_shape()
  dim_names <- names(dim_info)
  if (is.null(dim_names) || length(dim_names) != length(dim_info) || any(!nzchar(dim_names))) {
    dim_names <- paste0("dim", seq_along(dim_info))
  }
  dim_names_l <- tolower(dim_names)

  coord_x_name <- .find_coord_name(var_all, coord_x_alias)
  coord_y_name <- .find_coord_name(var_all, coord_y_alias)
  coord_t_name <- .find_coord_name(var_all, coord_t_alias)

  time_dim_name <- dim_names[which(dim_names_l %in% tolower(coord_t_alias))[1]]
  if (is.na(time_dim_name) || length(time_dim_name) == 0L) time_dim_name <- coord_t_name
  if (is.null(time_dim_name) && length(dim_names) >= 3L) time_dim_name <- dim_names[length(dim_names)]

  non_time <- if (!is.null(time_dim_name)) dim_names[dim_names != time_dim_name] else dim_names
  if (length(non_time) < 2L) {
    stop("Zarr store does not expose at least two non-time dimensions.", call. = FALSE)
  }

  x_dim_name <- non_time[which(tolower(non_time) %in% tolower(coord_x_alias))[1]]
  y_dim_name <- non_time[which(tolower(non_time) %in% tolower(coord_y_alias))[1]]
  if (is.na(x_dim_name) || length(x_dim_name) == 0L) x_dim_name <- coord_x_name
  if (is.na(y_dim_name) || length(y_dim_name) == 0L) y_dim_name <- coord_y_name
  if (is.null(x_dim_name)) x_dim_name <- non_time[1]
  if (is.null(y_dim_name) || identical(y_dim_name, x_dim_name)) {
    rem <- setdiff(non_time, x_dim_name)
    y_dim_name <- if (length(rem) > 0L) rem[1] else non_time[2]
  }

  get_dim_vals <- function(name) {
    if (zgroup$contains_item(name)) {
      zarray <- zgroup$get_item(name)
      return(as.vector(.get_item_all(zarray)))
    }
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
      if (zgroup$contains_item(tv)) {
        time_dim_name <- tv
        time_vals_raw <- get_dim_vals(tv)
        break
      }
    }
  }
  if (is.null(time_vals_raw) || length(time_vals_raw) < 1L) {
    stop("Could not identify/read time dimension.", call. = FALSE)
  }

  # Parse time to Date
  time_zarray <- zgroup$get_item(time_dim_name)
  date <- .parse_time_to_date(time_zarray, as.vector(time_vals_raw))
  if (is.null(date) || !inherits(date, "Date")) {
    stop("Could not convert Zarr time to Date (missing/unsupported CF units).", call. = FALSE)
  }

  if (isFALSE(keep_leap_day)) {
    date <- date[format(date, "%m-%d") != "02-29"]
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
  .to_time_by_cell <- function(arr, dim_names, xname, yname, tname, nx, ny, nt) {
    if (is.null(dim_names) || length(dim_names) != length(dim(arr)) || any(!nzchar(dim_names))) {
      dim_names <- paste0("dim", seq_along(dim(arr)))
    }
    dn <- tolower(dim_names)
    xi <- if (is.null(xname)) NA_integer_ else match(tolower(xname), dn)
    yi <- if (is.null(yname)) NA_integer_ else match(tolower(yname), dn)
    ti <- if (is.null(tname)) NA_integer_ else match(tolower(tname), dn)

    if (is.na(xi)) xi <- which(dn %in% tolower(coord_x_alias))[1]
    if (is.na(yi)) yi <- which(dn %in% tolower(coord_y_alias))[1]
    if (is.na(ti)) ti <- which(dn %in% tolower(coord_t_alias))[1]
    if (is.na(ti) && length(dim(arr)) >= 3L) ti <- length(dim(arr))
    if (is.na(xi) || is.na(yi)) {
      spare <- setdiff(seq_along(dim(arr)), ti)
      if (is.na(xi) && length(spare) >= 1L) xi <- spare[1]
      if (is.na(yi) && length(spare) >= 2L) yi <- spare[2]
    }

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

  for (v in var) {
    .msg("Reading variable: ", v)

    v_array <- zgroup$get_item(v)
    v_dim_info <- v_array$get_shape()
    v_dim_names <- names(v_dim_info)

    # Read the full array
    x <- as.array(.get_item_all(v_array))
    x <- .mask_missing(v_array, x)

    if (is.null(dim(x))) {
      mat <- matrix(as.numeric(x), nrow = length(time_vals_raw), ncol = ncell)
    } else {
      mat <- .to_time_by_cell(x, v_dim_names, x_dim_name, y_dim_name, time_dim_name, nx, ny, length(time_vals_raw))
    }

    # Drop Feb 29 in matrix if we dropped it in date
    if (isFALSE(keep_leap_day)) {
      time_zarray <- zgroup$get_item(time_dim_name)
      all_dates <- .parse_time_to_date(time_zarray, as.vector(time_vals_raw))
      keep_idx <- format(all_dates, "%m-%d") != "02-29"
      mat <- mat[keep_idx, , drop = FALSE]
    }

    if (!is.null(signif_digits)) mat <- signif(mat, signif_digits)
    if (isTRUE(drop_all_na) && all(is.na(mat))) next

    out_name <- v
    if (!is.null(rename_map) && v %in% names(rename_map)) out_name <- unname(rename_map[[v]])

    var_mats[[out_name]] <- mat
    kept_vars <- c(kept_vars, out_name)
  }

  # ---------------------------------------------------------------------------
  # Build per-cell data frames
  # ---------------------------------------------------------------------------
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

  # ---------------------------------------------------------------------------
  # Collect dimensions and attributes
  # ---------------------------------------------------------------------------
  dimensions <- list()
  dimensions[[x_dim_name]] <- as.vector(x_vals)
  dimensions[[y_dim_name]] <- as.vector(y_vals)
  dimensions[["time"]] <- as.vector(time_vals_raw)

  # Get global attributes (from .zattrs at group level)
  global_atts <- .get_zarr_attrs(zgroup)

  # Get variable attributes
  var_atts <- list()
  for (vn in var_all) {
    if (zgroup$contains_item(vn)) {
      var_atts[[vn]] <- .get_zarr_attrs(zgroup$get_item(vn))
    } else {
      var_atts[[vn]] <- list()
    }
  }

  attributes <- list(global = global_atts, variables = var_atts)

  # ---------------------------------------------------------------------------
  # Return structure matching read_netcdf()
  # ---------------------------------------------------------------------------
  list(
    data = data,
    grid = grid,
    date = date,
    dimensions = dimensions,
    attributes = attributes
  )
}




#' Write tidy data frames to Zarr store
#'
#' Writes gridded time-series data to a Zarr store. Can accept either a
#' structured list (e.g., from read_netcdf/read_zarr) or individual components.
#'
#' @param dir Character. Path to output Zarr store (will be created).
#' @param data List of data.frames (one per grid cell), each with time-series
#'   columns for variables. OR a single data.frame with coordinates.
#' @param grid Tibble/data.frame with columns: id, xind, yind, x, y.
#' @param date Vector of Date or POSIXct values for the time dimension.
#' @param dimensions Named list of dimension coordinate vectors (x, y, time).
#' @param attributes List with $global and $variables sub-lists for metadata.
#' @param input_list Alternative: a single list with components data, grid, date,
#'   dimensions, attributes (as returned by read_netcdf/read_zarr).
#' @param x_dim_name Character. Name for x/longitude dimension (default: "lon").
#' @param y_dim_name Character. Name for y/latitude dimension (default: "lat").
#' @param time_dim_name Character. Name for time dimension (default: "time").
#' @param time_units Character. CF-compliant time units (e.g., "days since 1850-01-01").
#'   If NULL and date is Date, uses "days since YYYY-01-01" from first date.
#' @param fill_value Numeric. Value to use for missing data (default: -9999).
#' @param overwrite Logical. If TRUE, overwrite existing zarr store (default: FALSE).
#' @param compression Character. Compression algorithm: "none", "zlib", "blosc" (default: "zlib").
#' @param compression_level Integer. Compression level 0-9 (default: 5).
#' @param chunk_size Named list. Chunk dimensions (e.g., list(time=365, lat=50, lon=50)).
#'   If NULL, uses automatic chunking.
#' @param slice_mode Integer. Indexing mode for pizzarr: 1 for R-style (1-based),
#'   0 for Python-style (0-based). Default is 1.
#' @param verbose Logical. If TRUE, print progress messages.
#'
#' @return Invisibly returns the path to the created zarr store.
#'
#' @importFrom pizzarr zarr_create_group
#' @export
write_zarr <- function(
    dir,
    data = NULL,
    grid = NULL,
    date = NULL,
    dimensions = NULL,
    attributes = NULL,
    input_list = NULL,
    x_dim_name = "lon",
    y_dim_name = "lat",
    time_dim_name = "time",
    time_units = NULL,
    fill_value = -9999,
    overwrite = FALSE,
    compression = "zlib",
    compression_level = 5,
    chunk_size = NULL,
    slice_mode = 1,
    verbose = FALSE
) {

  # ---------------------------------------------------------------------------
  # Input validation and extraction
  # ---------------------------------------------------------------------------
  if (!is.character(dir) || length(dir) != 1L || !nzchar(dir)) {
    stop("dir must be a non-empty character path.", call. = FALSE)
  }

  if (!is.logical(overwrite) || length(overwrite) != 1L) {
    stop("overwrite must be a single logical (TRUE/FALSE).", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("verbose must be a single logical (TRUE/FALSE).", call. = FALSE)
  }
  if (!is.character(compression) || length(compression) != 1L || !compression %in% c("none", "zlib", "blosc")) {
    stop("compression must be one of: 'none', 'zlib', 'blosc'", call. = FALSE)
  }
  if (!is.numeric(compression_level) || length(compression_level) != 1L ||
      !is.finite(compression_level) || (compression_level %% 1) != 0 ||
      compression_level < 0 || compression_level > 9) {
    stop("compression_level must be between 0 and 9.", call. = FALSE)
  }
  compression_level <- as.integer(compression_level)

  if (dir.exists(dir)) {
    if (!overwrite) {
      stop("Zarr store already exists: ", dir, ". Set overwrite=TRUE to replace.", call. = FALSE)
    }
    unlink(dir, recursive = TRUE)
  }

  # If input_list provided, extract components
  if (!is.null(input_list)) {
    if (!is.list(input_list)) {
      stop("input_list must be a list with components: data, grid, date, dimensions, attributes", call. = FALSE)
    }
    if (is.null(data)) data <- input_list$data
    if (is.null(grid)) grid <- input_list$grid
    if (is.null(date)) date <- input_list$date
    if (is.null(dimensions)) dimensions <- input_list$dimensions
    if (is.null(attributes)) attributes <- input_list$attributes
  }

  # Validate required components
  if (is.null(data)) stop("data must be provided (or in input_list).", call. = FALSE)
  if (is.null(grid)) stop("grid must be provided (or in input_list).", call. = FALSE)
  if (is.null(date)) stop("date must be provided (or in input_list).", call. = FALSE)

  if (!is.numeric(slice_mode) || length(slice_mode) != 1L || !slice_mode %in% c(0, 1)) {
    stop("slice_mode must be either 0 (Python-style, 0-based) or 1 (R-style, 1-based).", call. = FALSE)
  }
  slice_mode <- as.integer(slice_mode)

  # ---------------------------------------------------------------------------
  # Load pizzarr and set slice mode
  # ---------------------------------------------------------------------------
  if (!requireNamespace("pizzarr", quietly = TRUE)) {
    stop("Package 'pizzarr' is required. Install with: install.packages('pizzarr')", call. = FALSE)
  }

  # Set the slice mode for pizzarr
  old_slice_mode <- getOption("pizzarr.slice_mode")
  on.exit(options(pizzarr.slice_mode = old_slice_mode), add = TRUE)
  options(pizzarr.slice_mode = slice_mode)

  .msg <- function(...) if (isTRUE(verbose)) message(...)

  if (verbose) {
    .msg("Using ", ifelse(slice_mode == 1, "R-style (1-based)", "Python-style (0-based)"), " indexing")
  }

  # ---------------------------------------------------------------------------
  # Extract dimensions and prepare coordinate arrays
  # ---------------------------------------------------------------------------
  if (!is.null(dimensions)) {
    x_vals <- dimensions[[x_dim_name]]
    y_vals <- dimensions[[y_dim_name]]
    time_vals_raw <- dimensions[[time_dim_name]]
  } else {
    # Extract from grid and date
    x_vals <- sort(unique(grid$x))
    y_vals <- sort(unique(grid$y))
    time_vals_raw <- seq_along(date)
  }

  nx <- length(x_vals)
  ny <- length(y_vals)
  nt <- length(date)

  .msg("Dimensions: ", x_dim_name, "=", nx, ", ", y_dim_name, "=", ny, ", ", time_dim_name, "=", nt)

  # ---------------------------------------------------------------------------
  # Encode time to numeric
  # ---------------------------------------------------------------------------
  if (is.null(time_units)) {
    if (inherits(date, "Date")) {
      origin_date <- as.Date(paste0(format(date[1], "%Y"), "-01-01"))
      time_units <- paste0("days since ", as.character(origin_date))
      time_vals_raw <- as.numeric(date - origin_date)
    } else if (inherits(date, "POSIXct") || inherits(date, "POSIXlt")) {
      origin_date <- as.POSIXct(paste0(format(date[1], "%Y"), "-01-01 00:00:00"), tz = "UTC")
      time_units <- paste0("seconds since ", format(origin_date, "%Y-%m-%d %H:%M:%S"))
      time_vals_raw <- as.numeric(difftime(date, origin_date, units = "secs"))
    } else {
      stop("date must be Date or POSIXct, or time_units must be specified.", call. = FALSE)
    }
  }

  .msg("Time units: ", time_units)

  # ---------------------------------------------------------------------------
  # Extract variable names from data
  # ---------------------------------------------------------------------------
  if (is.data.frame(data)) {
    # Single data frame - convert to list format
    stop("Single data.frame input not yet implemented. Use list of data.frames per grid cell.", call. = FALSE)
  }

  if (!is.list(data) || length(data) != nrow(grid)) {
    stop("data must be a list with one data.frame per grid cell (length = nrow(grid)).", call. = FALSE)
  }

  var_names <- names(data[[1]])
  if (is.null(var_names) || length(var_names) == 0L) {
    stop("data[[1]] must have named columns for variables.", call. = FALSE)
  }

  .msg("Variables to write: ", paste(var_names, collapse = ", "))

  # ---------------------------------------------------------------------------
  # Convert list of data frames to variable matrices [time, cell]
  # ---------------------------------------------------------------------------
  var_mats <- list()
  ncell <- length(data)

  for (vname in var_names) {
    mat <- matrix(NA_real_, nrow = nt, ncol = ncell)
    for (k in seq_len(ncell)) {
      if (vname %in% names(data[[k]])) {
        mat[, k] <- as.numeric(data[[k]][[vname]])
      }
    }
    var_mats[[vname]] <- mat
  }

  # ---------------------------------------------------------------------------
  # Reshape matrices to [x, y, time] arrays
  # ---------------------------------------------------------------------------
  .mat_to_array <- function(mat, nx, ny, nt) {
    # mat is [time, cell] where cell = x * y (row-major from expand.grid)
    arr <- array(NA_real_, dim = c(nx, ny, nt))

    for (t in seq_len(nt)) {
      slice <- matrix(mat[t, ], nrow = nx, ncol = ny, byrow = FALSE)
      arr[, , t] <- slice
    }
    arr
  }

  var_arrays <- list()
  for (vname in var_names) {
    .msg("Reshaping variable: ", vname)
    arr <- .mat_to_array(var_mats[[vname]], nx, ny, nt)

    # Replace NA with fill_value
    arr[is.na(arr)] <- fill_value

    var_arrays[[vname]] <- arr
  }

  # ---------------------------------------------------------------------------
  # Create Zarr group and write arrays
  # ---------------------------------------------------------------------------
  .msg("Creating Zarr store: ", dir)

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  # Create root group
  zgroup <- pizzarr::zarr_create_group(dir)

  # Set up compression using pizzarr codec classes
  compressor <- NULL
  if (compression == "zlib") {
    compressor <- pizzarr::ZlibCodec$new()
  } else if (compression == "blosc") {
    if (!requireNamespace("Rarr", quietly = TRUE)) {
      stop("Rarr package must be installed to use the Blosc codec. Install with BiocManager::install('Rarr')", call. = FALSE)
    }
    compressor <- pizzarr::BloscCodec$new(cname = "lz4", clevel = compression_level)
  }

  # Determine chunks
  if (is.null(chunk_size)) {
    time_chunk <- min(365, nt)
    x_chunk <- min(50, nx)
    y_chunk <- min(50, ny)
    chunk_size <- list()
    chunk_size[[x_dim_name]] <- x_chunk
    chunk_size[[y_dim_name]] <- y_chunk
    chunk_size[[time_dim_name]] <- time_chunk
  }

  chunks <- c(
    chunk_size[[x_dim_name]],
    chunk_size[[y_dim_name]],
    chunk_size[[time_dim_name]]
  )

  .msg("Chunk size: [", paste(chunks, collapse = ", "), "]")

  .create_dataset <- function(name, data, shape, dtype, chunks) {
    if (is.null(compressor)) {
      zgroup$create_dataset(
        name = name,
        data = data,
        shape = shape,
        dtype = dtype,
        chunks = chunks
      )
    } else {
      zgroup$create_dataset(
        name = name,
        data = data,
        shape = shape,
        dtype = dtype,
        chunks = chunks,
        compressor = compressor
      )
    }
  }

  .set_zarr_attrs <- function(zobj, attrs) {
    if (length(attrs) == 0L) return(invisible(NULL))
    zattrs <- zobj$get_attrs()
    for (aname in names(attrs)) {
      zattrs$set_item(aname, attrs[[aname]])
    }
    invisible(NULL)
  }

  # ---------------------------------------------------------------------------
  # Write coordinate variables
  # ---------------------------------------------------------------------------
  .msg("Writing coordinate: ", x_dim_name)
  x_array <- .create_dataset(
    name = x_dim_name,
    data = as.array(as.numeric(x_vals)),
    shape = nx,
    dtype = "<f8",
    chunks = min(chunk_size[[x_dim_name]], nx)
  )
  .set_zarr_attrs(x_array, list(
    units = "degrees_east",
    long_name = "longitude",
    standard_name = "longitude"
  ))

  .msg("Writing coordinate: ", y_dim_name)
  y_array <- .create_dataset(
    name = y_dim_name,
    data = as.array(as.numeric(y_vals)),
    shape = ny,
    dtype = "<f8",
    chunks = min(chunk_size[[y_dim_name]], ny)
  )
  .set_zarr_attrs(y_array, list(
    units = "degrees_north",
    long_name = "latitude",
    standard_name = "latitude"
  ))

  .msg("Writing coordinate: ", time_dim_name)
  time_array <- .create_dataset(
    name = time_dim_name,
    data = as.array(as.numeric(time_vals_raw)),
    shape = nt,
    dtype = "<f8",
    chunks = min(chunk_size[[time_dim_name]], nt)
  )
  .set_zarr_attrs(time_array, list(
    units = time_units,
    long_name = "time",
    standard_name = "time",
    calendar = "standard"
  ))

  # ---------------------------------------------------------------------------
  # Write data variables
  # ---------------------------------------------------------------------------
  for (vname in var_names) {
    .msg("Writing variable: ", vname)

    var_array <- .create_dataset(
      name = vname,
      data = var_arrays[[vname]],
      shape = c(nx, ny, nt),
      dtype = "<f4",
      chunks = chunks
    )

    # Set variable attributes
    var_attrs <- list(
      `_FillValue` = fill_value,
      missing_value = fill_value
    )

    # Add custom attributes if provided
    if (!is.null(attributes) && !is.null(attributes$variables[[vname]])) {
      custom_attrs <- attributes$variables[[vname]]
      for (aname in names(custom_attrs)) {
        if (!aname %in% names(var_attrs)) {
          var_attrs[[aname]] <- custom_attrs[[aname]]
        }
      }
    }

    .set_zarr_attrs(var_array, var_attrs)
  }

  # ---------------------------------------------------------------------------
  # Write global attributes
  # ---------------------------------------------------------------------------
  global_attrs <- list(
    Conventions = "CF-1.8",
    institution = "Generated by weathergenr",
    history = paste0(Sys.time(), " - Created with write_zarr()"),
    source = "R weathergenr package"
  )

  if (!is.null(attributes) && !is.null(attributes$global)) {
    for (aname in names(attributes$global)) {
      global_attrs[[aname]] <- attributes$global[[aname]]
    }
  }

  .set_zarr_attrs(zgroup, global_attrs)

  .msg("Zarr store written successfully: ", dir)

  invisible(dir)
}
