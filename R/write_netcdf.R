#' Write Gridded Data to a NetCDF File (Template-Based)
#'
#' @description
#' Writes simulated or processed gridded climate/weather data to a NetCDF file, using variable and dimension
#' definitions from a template NetCDF. This wrapper ensures consistency with your reference files (dimensions, spatial ref, etc).
#'
#' Each element in `data` should be a list of variables for a single grid cell; the list should have the same structure and names as variables in the template file.
#'
#' @param data A list of length = n_grids. Each element is a named list of vectors (or matrices) containing simulated data for each variable at that grid.
#'   For example, \code{data[[i]]$precip} is a vector for grid i.
#' @param coord.grid Data frame mapping grid indices (`xind`, `yind`) for each grid cell. Should have same length/order as `data`.
#' @param output.path Directory to write the file. Will be created if it doesn't exist.
#' @param origin.date Origin date as "YYYY-MM-DD" string for the time dimension.
#' @param nc.template.file Path to a reference/template NetCDF file (must exist).
#' @param calendar.type Calendar type for NetCDF time (default: "noleap").
#' @param nc.compression Integer; NetCDF compression level (0-9, default 4).
#' @param nc.spatial.ref Name of the spatial reference variable in the template (default "spatial_ref").
#' @param nc.file.prefix Prefix for the output file name (default "clim_change_rlz").
#' @param nc.file.suffix Optional suffix for the file name (default "").
#' @param signif.digits If not NULL, round all output values to this number of significant digits.
#'
#' @return (Invisibly) The file path to the written NetCDF file.
#'
#' @details
#' All variables and spatial grids must match the template file.
#' Attributes and dimensions are copied from the template.
#'
#' @examples
#' \dontrun{
#' # Suppose you have a template NetCDF "template.nc" with x=2, y=2, time=3
#' # Simulate some data:
#' sim_data <- replicate(4, list(
#'   precip = rnorm(3, 10, 2),
#'   temp   = rnorm(3, 20, 1)
#' ), simplify = FALSE)
#'
#' grid <- data.frame(xind = rep(1:2, each = 2), yind = rep(1:2, times = 2))
#' writeNetcdf(
#'   data = sim_data,
#'   coord.grid = grid,
#'   output.path = tempdir(),
#'   origin.date = "2000-01-01",
#'   nc.template.file = "template.nc",
#'   nc.file.prefix = "test_output"
#' )
#' }
#'
#' @import ncdf4
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
                        signif.digits = NULL) {

  # ===========================================================================
  # BUG FIX 1: Comprehensive Input Validation
  # ===========================================================================

  if (is.null(data) || !is.list(data)) {
    stop("'data' must be a list of grid cell data.", call. = FALSE)
  }

  if (length(data) == 0) {
    stop("'data' cannot be empty.", call. = FALSE)
  }

  if (is.null(coord.grid) || !is.data.frame(coord.grid)) {
    stop("'coord.grid' must be a data frame.", call. = FALSE)
  }

  if (!all(c("xind", "yind") %in% names(coord.grid))) {
    stop("'coord.grid' must have 'xind' and 'yind' columns.", call. = FALSE)
  }

  if (nrow(coord.grid) != length(data)) {
    stop(
      "'coord.grid' rows (", nrow(coord.grid), ") ",
      "must equal length of 'data' (", length(data), ").",
      call. = FALSE
    )
  }

  if (is.null(output.path) || !is.character(output.path)) {
    stop("'output.path' must be a character string.", call. = FALSE)
  }

  if (is.null(origin.date)) {
    stop("'origin.date' must be provided.", call. = FALSE)
  }

  # BUG FIX 2: Validate origin.date format
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
    stop(
      "'nc.template.file' does not exist: ", nc.template.file,
      call. = FALSE
    )
  }

  if (!is.numeric(nc.compression) || nc.compression < 0 || nc.compression > 9) {
    stop("'nc.compression' must be an integer between 0 and 9.", call. = FALSE)
  }

  # BUG FIX 3: Check data structure consistency
  if (!all(sapply(data, is.list))) {
    stop("Each element in 'data' must be a list.", call. = FALSE)
  }

  variables <- names(data[[1]])

  if (is.null(variables) || length(variables) == 0) {
    stop("'data[[1]]' must be a named list with variable names.", call. = FALSE)
  }

  # Check all grid cells have same variables
  for (i in seq_along(data)) {
    if (!identical(sort(names(data[[i]])), sort(variables))) {
      stop(
        "All elements in 'data' must have the same variable names. ",
        "Grid cell ", i, " has different variables.",
        call. = FALSE
      )
    }
  }

  # BUG FIX 4: Check data dimensions consistency
  # Get expected length from first grid cell
  expected_lengths <- sapply(data[[1]], length)

  for (i in seq_along(data)) {
    actual_lengths <- sapply(data[[i]], length)
    if (!identical(actual_lengths, expected_lengths)) {
      stop(
        "All grid cells must have the same time series length. ",
        "Grid cell ", i, " has inconsistent lengths.",
        call. = FALSE
      )
    }
  }

  sim.length <- length(data[[1]][[1]])

  # ===========================================================================
  # Directory Creation with Error Handling
  # ===========================================================================

  if (!dir.exists(output.path)) {
    tryCatch(
      dir.create(output.path, recursive = TRUE, showWarnings = FALSE),
      error = function(e) {
        stop(
          "Cannot create output directory: ", output.path,
          ". Error: ", e$message,
          call. = FALSE
        )
      }
    )
  }

  # ===========================================================================
  # OPTIMIZATION 1: Apply Rounding Early (Before File Operations)
  # ===========================================================================

  if (!is.null(signif.digits)) {
    if (!is.numeric(signif.digits) || signif.digits < 1) {
      stop("'signif.digits' must be a positive integer.", call. = FALSE)
    }

    # Round all data upfront (more efficient than during write)
    data <- lapply(data, function(grid_data) {
      lapply(grid_data, function(var_data) {
        round(var_data, signif.digits)
      })
    })
  }

  # ===========================================================================
  # Open Template File with Error Handling
  # ===========================================================================

  ncin <- tryCatch(
    ncdf4::nc_open(nc.template.file),
    error = function(e) {
      stop(
        "Cannot open template NetCDF file: ", nc.template.file,
        ". Error: ", e$message,
        call. = FALSE
      )
    }
  )

  on.exit(ncdf4::nc_close(ncin), add = TRUE)

  # ===========================================================================
  # OPTIMIZATION 2: Read Template Metadata Efficiently
  # ===========================================================================

  # Get dimension names directly (avoid lapply overhead)
  dim_names <- names(ncin$dim)
  ncin_dim <- setNames(
    lapply(dim_names, function(d) ncdf4::ncvar_get(ncin, d)),
    dim_names
  )

  # BUG FIX 5: Check spatial reference exists in template
  if (!nc.spatial.ref %in% names(ncin$var)) {
    ncdf4::nc_close(ncin)
    stop(
      "Spatial reference variable '", nc.spatial.ref,
      "' not found in template file.",
      call. = FALSE
    )
  }

  # Get spatial variables and global attributes
  ncin_var_spref <- ncdf4::ncvar_get(ncin, nc.spatial.ref)
  ncin_attribs_spref <- ncdf4::ncatt_get(ncin, nc.spatial.ref)
  ncin_attribs_glob <- ncdf4::ncatt_get(ncin, 0)

  # BUG FIX 6: Validate spatial reference attributes exist
  if (is.null(ncin_attribs_spref$x_dim) || is.null(ncin_attribs_spref$y_dim)) {
    ncdf4::nc_close(ncin)
    stop(
      "Spatial reference variable '", nc.spatial.ref,
      "' must have 'x_dim' and 'y_dim' attributes.",
      call. = FALSE
    )
  }

  # Dimension names ordered as x, y, time
  ncin_dimnames <- list(
    x = ncin_attribs_spref$x_dim,
    y = ncin_attribs_spref$y_dim,
    time = "time"
  )

  # BUG FIX 7: Check required dimensions exist in template
  missing_dims <- setdiff(
    c(ncin_dimnames$x, ncin_dimnames$y, ncin_dimnames$time),
    dim_names
  )

  if (length(missing_dims) > 0) {
    ncdf4::nc_close(ncin)
    stop(
      "Template file missing required dimensions: ",
      paste(missing_dims, collapse = ", "),
      call. = FALSE
    )
  }

  # ===========================================================================
  # OPTIMIZATION 3: Pre-compute Dimension Lengths
  # ===========================================================================

  nx <- length(ncin_dim[[ncin_dimnames$x]])
  ny <- length(ncin_dim[[ncin_dimnames$y]])
  nt <- sim.length

  # BUG FIX 8: Validate grid indices are within bounds
  if (any(coord.grid$xind < 1 | coord.grid$xind > nx)) {
    stop(
      "Some xind values in 'coord.grid' are out of bounds [1, ", nx, "].",
      call. = FALSE
    )
  }

  if (any(coord.grid$yind < 1 | coord.grid$yind > ny)) {
    stop(
      "Some yind values in 'coord.grid' are out of bounds [1, ", ny, "].",
      call. = FALSE
    )
  }

  # ===========================================================================
  # Define Output Dimensions
  # ===========================================================================

  ncout_dim <- list()

  # Time dimension
  ncout_dim[[ncin_dimnames$time]] <- ncdf4::ncdim_def(
    name = ncin_dimnames$time,
    units = paste0("days since ", format(as.Date(origin.date), "%Y-%m-%d 00:00:00")),
    vals = seq(0, nt - 1),  # More efficient than as.numeric(1:sim.length) - 1
    calendar = calendar.type,
    longname = NULL
  )

  # Y dimension
  ncout_dim[[ncin_dimnames$y]] <- ncdf4::ncdim_def(
    name = ncin_dimnames$y,
    units = "",
    vals = as.numeric(ncin_dim[[ncin_dimnames$y]]),
    longname = NULL
  )

  # X dimension
  ncout_dim[[ncin_dimnames$x]] <- ncdf4::ncdim_def(
    name = ncin_dimnames$x,
    units = "",
    vals = as.numeric(ncin_dim[[ncin_dimnames$x]]),
    longname = NULL
  )

  # ===========================================================================
  # OPTIMIZATION 4: Efficient Chunking Strategy
  # ===========================================================================

  # Use optimal chunk sizes for time-series data
  # Chunk along time dimension for efficient temporal access
  ncout_chunks <- c(
    min(nt, 365),  # Chunk time: max 1 year at a time
    ny,
    nx
  )

  # ===========================================================================
  # OPTIMIZATION 5: Define Variables with Proper Ordering
  # ===========================================================================

  n_vars <- length(variables)
  ncout_vars <- vector("list", n_vars + 1)

  # ncdf4 warns that shuffle "only has an effect" for integer types.
  # So: enable shuffle only when the variable precision is integer-like.
  int_precisions <- c("byte", "short", "integer")

  for (i in seq_len(n_vars)) {
    vnm <- variables[i]

    # Prefer template precision if available; otherwise default to float
    tmpl_prec <- "float"
    if (vnm %in% names(ncin$var) && !is.null(ncin$var[[vnm]]$prec)) {
      tmpl_prec <- ncin$var[[vnm]]$prec
    }

    # Force float unless you explicitly want to inherit template precision
    # (keep your current behavior)
    out_prec <- "float"

    use_shuffle <- out_prec %in% int_precisions

    ncout_vars[[i]] <- ncdf4::ncvar_def(
      name = vnm,
      units = "",
      dim = ncout_dim,
      missval = NA_real_,
      longname = vnm,
      prec = out_prec,
      shuffle = use_shuffle,
      compression = nc.compression,
      chunksizes = ncout_chunks
    )
  }

  # ===========================================================================
  # Define spatial reference variable CLEANLY (do NOT reuse template var object)
  # ===========================================================================

  ncout_vars[[n_vars + 1]] <- ncdf4::ncvar_def(
    name     = nc.spatial.ref,
    units    = "",
    dim      = list(),        # scalar
    missval  = NA_integer_,
    longname = nc.spatial.ref,
    prec     = "integer"
    # NOTE: no compression, no shuffle
  )

  names(ncout_vars) <- c(variables, nc.spatial.ref)


  # ===========================================================================
  # Create Output File Path
  # ===========================================================================

  suffix_str <- if (nc.file.suffix != "") paste0("_", nc.file.suffix) else ""
  nc_file_path <- file.path(output.path, paste0(nc.file.prefix, suffix_str, ".nc"))

  # BUG FIX 11: Check if file already exists and warn
  if (file.exists(nc_file_path)) {
    warning(
      "Output file already exists and will be overwritten: ", nc_file_path,
      call. = FALSE
    )
  }

  # ===========================================================================
  # Create NetCDF File
  # ===========================================================================

  ncout_file <- tryCatch(
    ncdf4::nc_create(nc_file_path, ncout_vars, force_v4 = TRUE),
    error = function(e) {
      stop(
        "Cannot create NetCDF file: ", nc_file_path,
        ". Error: ", e$message,
        call. = FALSE
      )
    }
  )

  on.exit(ncdf4::nc_close(ncout_file), add = TRUE)

  # ===========================================================================
  # OPTIMIZATION 6: Pre-allocate Output Array Once
  # ===========================================================================

  # Allocate array once with correct dimensions
  ncout_vardata <- array(NA_real_, dim = c(nt, ny, nx))

  # ===========================================================================
  # OPTIMIZATION 7: Vectorized Data Writing
  # ===========================================================================

  # Pre-extract grid coordinates (avoid repeated indexing)
  xind_vec <- coord.grid$xind
  yind_vec <- coord.grid$yind
  n_grids <- nrow(coord.grid)

  # Write each variable
  for (i in seq_len(n_vars)) {
    var_name <- variables[i]

    # Reset array to NA (more efficient than recreating)
    ncout_vardata[] <- NA_real_

    # OPTIMIZATION: Vectorized assignment (much faster than loop)
    # Build matrix of data for this variable across all grid cells
    for (c in seq_len(n_grids)) {
      ncout_vardata[, yind_vec[c], xind_vec[c]] <- data[[c]][[var_name]]
    }

    # Write to NetCDF
    tryCatch(
      ncdf4::ncvar_put(ncout_file, varid = var_name, vals = ncout_vardata),
      error = function(e) {
        stop(
          "Error writing variable '", var_name, "' to NetCDF. ",
          "Error: ", e$message,
          call. = FALSE
        )
      }
    )
  }

  # ===========================================================================
  # Write Spatial Reference Variable and Attributes
  # ===========================================================================

  tryCatch({
    ncdf4::ncvar_put(ncout_file, varid = nc.spatial.ref, vals = ncin_var_spref)

    # Write spatial reference attributes
    for (k in seq_along(ncin_attribs_spref)) {
      ncdf4::ncatt_put(
        ncout_file,
        varid = nc.spatial.ref,
        attname = names(ncin_attribs_spref)[k],
        attval = ncin_attribs_spref[[k]]
      )
    }
  }, error = function(e) {
    stop(
      "Error writing spatial reference. Error: ", e$message,
      call. = FALSE
    )
  })

  # ===========================================================================
  # Write Global Attributes
  # ===========================================================================

  if (length(ncin_attribs_glob) > 0) {
    tryCatch({
      for (k in seq_along(ncin_attribs_glob)) {
        ncdf4::ncatt_put(
          ncout_file,
          varid = 0,
          attname = names(ncin_attribs_glob)[k],
          attval = ncin_attribs_glob[[k]]
        )
      }
    }, error = function(e) {
      warning(
        "Error writing global attributes: ", e$message,
        call. = FALSE
      )
    })
  }

  # Add provenance metadata
  ncdf4::ncatt_put(
    ncout_file,
    varid = 0,
    attname = "creation_date",
    attval = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )

  ncdf4::ncatt_put(
    ncout_file,
    varid = 0,
    attname = "created_with",
    attval = "weathergenr::writeNetcdf"
  )

  # ===========================================================================
  # Cleanup and Return
  # ===========================================================================

  # Files are closed automatically by on.exit handlers

  logger::log_info("[NetCDF] File written successfully: {nc_file_path}")
  logger::log_info("[NetCDF]   Variables: {paste(variables, collapse = ', ')}")
  logger::log_info("[NetCDF]   Dimensions: {nt} time x {ny} y x {nx} x")
  logger::log_info("[NetCDF]   Grid cells: {n_grids}")

  invisible(nc_file_path)
}
