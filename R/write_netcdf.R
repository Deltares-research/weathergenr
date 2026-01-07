#' Write Gridded Data to a NetCDF File (Template-Based)
#'
#' @description
#' Writes simulated or processed gridded climate/weather data to a NetCDF file, using variable and dimension
#' definitions from a template NetCDF. This wrapper ensures consistency with your reference files (dimensions, spatial ref, etc).
#'
#' @param verbose Logical. If TRUE, progress messages are logged via
#'   \code{logger::log_info("[NetCDF] ...")}.
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
