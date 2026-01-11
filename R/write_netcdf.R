#' Write Gridded Data to a NetCDF File (Template-Based)
#'
#' @description
#' Writes gridded climate/weather time series to a NetCDF file using a template NetCDF
#' as the schema source (dimensions, spatial reference variable, and global attributes).
#' The function expects a list of grid-cell payloads (each containing identical variable
#' names and equal-length time series) and places each cell into the output raster cube
#' according to \code{coord.grid} indices. Output is written as NetCDF-4 with optional
#' compression and optional early rounding of values.
#'
#' @details
#' The template file provided via \code{nc.template.file} is opened to:
#' \itemize{
#'   \item validate that the spatial reference variable \code{nc.spatial.ref} exists;
#'   \item read \code{x_dim} and \code{y_dim} attributes from the spatial reference variable
#'         to determine the x/y dimension names in the template;
#'   \item copy the spatial reference variable values and its attributes into the output file;
#'   \item copy template global attributes into the output file (best-effort);
#'   \item define output x/y dimension coordinate values based on template dimension values.
#' }
#'
#' A \code{time} dimension is defined as \code{0:(nt - 1)} with units
#' \code{"days since <origin.date> 00:00:00"} and the specified \code{calendar.type}.
#' Each variable in \code{data[[1]]} is written as a 3D array with dimension order
#' \code{time x y x x}. Grid cells are inserted at \code{[ , yind, xind ]} locations.
#'
#' The output filename is constructed as:
#' \code{file.path(output.path, paste0(nc.file.prefix, optional_suffix, ".nc"))}
#' where \code{optional_suffix} is \code{paste0("_", nc.file.suffix)} when
#' \code{nc.file.suffix != ""}.
#'
#' @param data List of length \code{nrow(coord.grid)}. Each element is a named list of
#'   numeric vectors (time series) for each variable to be written. All list elements must
#'   have identical variable names and identical time series lengths.
#' @param coord.grid Data frame mapping grid cells to template indices. Must contain
#'   integer columns \code{xind} and \code{yind} (1-based indices) and have
#'   \code{nrow(coord.grid) == length(data)}. Indices must be within the template x/y bounds.
#' @param output.path Character string. Directory where the NetCDF file will be written.
#'   Created if it does not exist.
#' @param origin.date Date or character coercible to \code{Date}. Used to define the time units
#'   (days since origin date at 00:00:00).
#' @param calendar.type Character string. NetCDF calendar attribute for the time dimension
#'   (e.g., \code{"noleap"}). Passed to \code{ncdf4::ncdim_def()}.
#' @param nc.template.file Character string. Path to an existing template NetCDF file used to
#'   infer dimensions, coordinate values, spatial reference, and global attributes.
#' @param nc.compression Integer in \code{[0, 9]}. NetCDF-4 deflate compression level applied
#'   to output variables.
#' @param nc.spatial.ref Character string. Name of the spatial reference variable in the
#'   template to copy to output (default \code{"spatial_ref"}). The variable must contain
#'   attributes \code{x_dim} and \code{y_dim}.
#' @param nc.file.prefix Character string. Output file prefix (base name without extension).
#' @param nc.file.suffix Character string. Optional suffix appended to the file name after an
#'   underscore. If \code{""}, no suffix is added.
#' @param signif.digits Optional integer. If provided, all values in \code{data} are rounded
#'   using \code{round(x, signif.digits)} prior to writing.
#' @param verbose Logical scalar. If \code{TRUE}, progress messages are emitted using
#'   \code{logger::log_info()} with a \code{[NetCDF]} tag.
#'
#' @return Invisibly returns the full path to the written NetCDF file (character scalar).
#'
#' @section File structure:
#' \itemize{
#'   \item Dimensions: \code{time}, template \code{y}, template \code{x}
#'   \item Variables: one variable per name in \code{names(data[[1]])}, plus \code{nc.spatial.ref}
#'   \item Global attributes: copied from template (best-effort) plus \code{creation_date} and
#'         \code{created_with}
#' }
#'
#' @section Side effects:
#' Creates directories (if needed), writes a NetCDF file, and may overwrite an existing output
#' file with the same name (a warning is issued).
#'
#' @examples
#' \dontrun{
#' # Minimal example (schematic)
#' coord.grid <- data.frame(xind = c(1, 2), yind = c(1, 1))
#' data <- list(
#'   list(prcp = rnorm(365), tavg = rnorm(365)),
#'   list(prcp = rnorm(365), tavg = rnorm(365))
#' )
#'
#' out_file <- write_netcdf(
#'   data = data,
#'   coord.grid = coord.grid,
#'   output.path = "out",
#'   origin.date = "1990-01-01",
#'   calendar.type = "noleap",
#'   nc.template.file = "template.nc",
#'   nc.file.prefix = "clim_change_rlz",
#'   nc.file.suffix = "run01",
#'   nc.compression = 4,
#'   signif.digits = 3,
#'   verbose = TRUE
#' )
#' }
#'
#' @import ncdf4
#' @importFrom logger log_info
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

  # ===========================================================================
  # Input validation
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
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be logical (TRUE/FALSE).", call. = FALSE)
  }

  .log <- function(fmt, ...) {
    if (isTRUE(verbose)) logger::log_info(fmt, ...)
    invisible(NULL)
  }

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

  if (!all(vapply(data, is.list, logical(1)))) {
    stop("Each element in 'data' must be a list.", call. = FALSE)
  }

  variables <- names(data[[1]])
  if (is.null(variables) || length(variables) == 0) {
    stop("'data[[1]]' must be a named list with variable names.", call. = FALSE)
  }

  for (i in seq_along(data)) {
    if (!identical(sort(names(data[[i]])), sort(variables))) {
      stop(
        "All elements in 'data' must have the same variable names. ",
        "Grid cell ", i, " has different variables.",
        call. = FALSE
      )
    }
  }

  expected_lengths <- vapply(data[[1]], length, integer(1))
  for (i in seq_along(data)) {
    actual_lengths <- vapply(data[[i]], length, integer(1))
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
  # Directory creation
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
  # Round early (optional)
  # ===========================================================================

  if (!is.null(signif.digits)) {
    if (!is.numeric(signif.digits) || signif.digits < 1) {
      stop("'signif.digits' must be a positive integer.", call. = FALSE)
    }
    data <- lapply(data, function(grid_data) {
      lapply(grid_data, function(var_data) round(var_data, signif.digits))
    })
  }

  # ===========================================================================
  # Open template NetCDF
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

  dim_names <- names(ncin$dim)
  ncin_dim <- setNames(
    lapply(dim_names, function(d) ncdf4::ncvar_get(ncin, d)),
    dim_names
  )

  if (!nc.spatial.ref %in% names(ncin$var)) {
    ncdf4::nc_close(ncin)
    stop(
      "Spatial reference variable '", nc.spatial.ref,
      "' not found in template file.",
      call. = FALSE
    )
  }

  ncin_var_spref <- ncdf4::ncvar_get(ncin, nc.spatial.ref)
  ncin_attribs_spref <- ncdf4::ncatt_get(ncin, nc.spatial.ref)
  ncin_attribs_glob <- ncdf4::ncatt_get(ncin, 0)

  if (is.null(ncin_attribs_spref$x_dim) || is.null(ncin_attribs_spref$y_dim)) {
    ncdf4::nc_close(ncin)
    stop(
      "Spatial reference variable '", nc.spatial.ref,
      "' must have 'x_dim' and 'y_dim' attributes.",
      call. = FALSE
    )
  }

  ncin_dimnames <- list(
    x = ncin_attribs_spref$x_dim,
    y = ncin_attribs_spref$y_dim,
    time = "time"
  )

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

  nx <- length(ncin_dim[[ncin_dimnames$x]])
  ny <- length(ncin_dim[[ncin_dimnames$y]])
  nt <- sim.length

  if (any(coord.grid$xind < 1 | coord.grid$xind > nx)) {
    stop("Some xind values in 'coord.grid' are out of bounds [1, ", nx, "].", call. = FALSE)
  }
  if (any(coord.grid$yind < 1 | coord.grid$yind > ny)) {
    stop("Some yind values in 'coord.grid' are out of bounds [1, ", ny, "].", call. = FALSE)
  }

  # ===========================================================================
  # Define output dims/vars
  # ===========================================================================

  ncout_dim <- list()

  ncout_dim[[ncin_dimnames$time]] <- ncdf4::ncdim_def(
    name = ncin_dimnames$time,
    units = paste0("days since ", format(as.Date(origin.date), "%Y-%m-%d 00:00:00")),
    vals = seq(0, nt - 1),
    calendar = calendar.type,
    longname = NULL
  )

  ncout_dim[[ncin_dimnames$y]] <- ncdf4::ncdim_def(
    name = ncin_dimnames$y,
    units = "",
    vals = as.numeric(ncin_dim[[ncin_dimnames$y]]),
    longname = NULL
  )

  ncout_dim[[ncin_dimnames$x]] <- ncdf4::ncdim_def(
    name = ncin_dimnames$x,
    units = "",
    vals = as.numeric(ncin_dim[[ncin_dimnames$x]]),
    longname = NULL
  )

  ncout_chunks <- c(min(nt, 365), ny, nx)

  n_vars <- length(variables)
  ncout_vars <- vector("list", n_vars + 1)

  int_precisions <- c("byte", "short", "integer")

  for (i in seq_len(n_vars)) {
    vnm <- variables[i]

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

  ncout_vars[[n_vars + 1]] <- ncdf4::ncvar_def(
    name     = nc.spatial.ref,
    units    = "",
    dim      = list(),
    missval  = NA_integer_,
    longname = nc.spatial.ref,
    prec     = "integer"
  )

  names(ncout_vars) <- c(variables, nc.spatial.ref)

  # ===========================================================================
  # Output path + create file
  # ===========================================================================

  suffix_str <- if (nc.file.suffix != "") paste0("_", nc.file.suffix) else ""
  nc_file_path <- file.path(output.path, paste0(nc.file.prefix, suffix_str, ".nc"))

  if (file.exists(nc_file_path)) {
    warning("Output file already exists and will be overwritten: ", nc_file_path, call. = FALSE)
  }

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
  # Write variables
  # ===========================================================================

  ncout_vardata <- array(NA_real_, dim = c(nt, ny, nx))

  xind_vec <- coord.grid$xind
  yind_vec <- coord.grid$yind
  n_grids <- nrow(coord.grid)

  for (i in seq_len(n_vars)) {
    var_name <- variables[i]
    ncout_vardata[] <- NA_real_

    for (c in seq_len(n_grids)) {
      ncout_vardata[, yind_vec[c], xind_vec[c]] <- data[[c]][[var_name]]
    }

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
  # Spatial reference + attributes
  # ===========================================================================

  tryCatch({
    ncdf4::ncvar_put(ncout_file, varid = nc.spatial.ref, vals = ncin_var_spref)

    for (k in seq_along(ncin_attribs_spref)) {
      ncdf4::ncatt_put(
        ncout_file,
        varid = nc.spatial.ref,
        attname = names(ncin_attribs_spref)[k],
        attval = ncin_attribs_spref[[k]]
      )
    }
  }, error = function(e) {
    stop("Error writing spatial reference. Error: ", e$message, call. = FALSE)
  })

  # ===========================================================================
  # Global attributes + provenance
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
      warning("Error writing global attributes: ", e$message, call. = FALSE)
    })
  }

  ncdf4::ncatt_put(
    ncout_file, varid = 0,
    attname = "creation_date",
    attval = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )

  ncdf4::ncatt_put(
    ncout_file, varid = 0,
    attname = "created_with",
    attval = "weathergenr::writeNetcdf"
  )

  # ===========================================================================
  # Final logs (gated)
  # ===========================================================================

  .log("[NetCDF] File written successfully: {nc_file_path}")
  .log("[NetCDF] Variables: {paste(variables, collapse = ', ')}")
  .log("[NetCDF] Dimensions: {nt} time x {ny} y x {nx} x")
  .log("[NetCDF] Grid cells: {n_grids}")

  invisible(nc_file_path)
}
