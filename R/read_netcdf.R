#' Read a gridded NetCDF file into tidy weather format (fast + memory efficient)
#'
#' @description
#' Efficiently reads 3D (x, y, time) gridded variables from a NetCDF file and
#' returns:
#' \itemize{
#'   \item per-grid-cell tidy time series (list of tibbles),
#'   \item grid coordinates and indices,
#'   \item a Date vector,
#'   \item dimension metadata,
#'   \item variable and global attributes.
#' }
#'
#' Variable behaviour:
#' \itemize{
#'   \item \code{variables = NULL} -> load all data variables (except \code{spatial.ref}).
#'   \item \code{variables = c("var1","var2")} -> load only these variables.
#'   \item \code{var_rename = c(old1 = "new1", old2 = "new2")} -> rename a subset
#'         of the loaded variables in the output (columns, attributes, dimnames).
#'         All names in \code{var_rename} must be among the selected variables.
#' }
#'
#' @param nc.file Character. Path to the NetCDF file.
#' @param variables Character vector of variable names to load. If \code{NULL},
#'   all variables except \code{spatial.ref} are loaded. Default is \code{NULL}.
#' @param var_rename Named character vector giving new names for a subset of the
#'   loaded variables, e.g. \code{c(precip = "prcp", temp = "T")}. Old names must
#'   be a subset of \code{variables} (or of all variables if \code{variables=NULL}).
#'   Default is \code{NULL}.
#' @param leap.days Logical. Whether the data include leap days. If \code{FALSE},
#'   Feb 29 is removed from the date vector. Default is \code{TRUE}.
#' @param omit.empty Logical. If \code{TRUE}, drop grid cells where all values
#'   are NA for all variables and all time steps. Default is \code{TRUE}.
#' @param spatial.ref Character. Name of the spatial reference variable in the
#'   NetCDF file (used to get x/y dimension names). Default is \code{"spatial_ref"}.
#' @param signif.digits Integer. If not \code{NULL}, values are rounded
#'   to this number of significant digits. Default is \code{NULL}.
#'
#' @return A list with five components:
#' \describe{
#'   \item{data}{List of tibbles, one per grid cell. Each tibble has rows
#'         representing time steps and columns representing variables.}
#'   \item{grid}{Tibble with columns \code{id, xind, yind, x, y} describing
#'         the spatial grid cells.}
#'   \item{date}{Vector of \code{Date} objects with length equal to the number
#'         of time steps.}
#'   \item{dimensions}{Named list of dimension vectors (e.g., x, y, time) as
#'         extracted from the NetCDF file.}
#'   \item{attributes}{Named list containing variable attributes (one element
#'         per variable using final renamed names), spatial reference attributes,
#'         and global attributes.}
#' }
#'
#' @examples
#' \dontrun{
#' # Load all variables from a NetCDF file
#' weather_data <- read_netcdf("climate_data.nc")
#'
#' # Load only specific variables and rename them
#' weather_data <- read_netcdf(
#'   nc.file = "climate_data.nc",
#'   variables = c("precipitation", "temperature"),
#'   var_rename = c(precipitation = "prcp", temperature = "temp")
#' )
#'
#' # Load without leap days
#' weather_data <- read_netcdf("climate_data.nc", leap.days = FALSE)
#'
#' # Round values and keep only non-empty cells
#' weather_data <- read_netcdf(
#'   nc.file = "climate_data.nc",
#'   signif.digits = 4,
#'   omit.empty = TRUE
#' )
#' }
#'
#' @import ncdf4
#' @import tibble
#' @export
read_netcdf <- function(
    nc.file,
    variables     = NULL,
    var_rename    = NULL,
    leap.days     = TRUE,
    omit.empty    = TRUE,
    spatial.ref   = "spatial_ref",
    signif.digits = NULL) {

  ## ---------------------------------------------------------------------------
  ## Input validation
  ## ---------------------------------------------------------------------------
  if (!file.exists(nc.file)) {
    stop("File does not exist: ", nc.file)
  }

  if (!is.logical(leap.days) || length(leap.days) != 1L) {
    stop("`leap.days` must be a single logical value (TRUE or FALSE).")
  }

  if (!is.logical(omit.empty) || length(omit.empty) != 1L) {
    stop("`omit.empty` must be a single logical value (TRUE or FALSE).")
  }

  if (!is.character(spatial.ref) || length(spatial.ref) != 1L) {
    stop("`spatial.ref` must be a single character string.")
  }

  if (!is.null(signif.digits)) {
    if (!is.numeric(signif.digits) || length(signif.digits) != 1L ||
        signif.digits < 1 || signif.digits != round(signif.digits)) {
      stop("`signif.digits` must be a single positive integer.")
    }
  }

  ## ---------------------------------------------------------------------------
  ## Open file
  ## ---------------------------------------------------------------------------
  nc_in <- ncdf4::nc_open(nc.file)
  on.exit(ncdf4::nc_close(nc_in))

  ## Check file has variables
  if (length(nc_in$var) == 0L) {
    stop("NetCDF file contains no variables.")
  }

  ## Helper: extract dimension names of a variable
  ncGetDimNames <- function(f, v) {
    dims <- f$var[[v]]$dim
    vapply(dims, function(x) x$name, character(1L))
  }

  ## ---------------------------------------------------------------------------
  ## Dimensions
  ## ---------------------------------------------------------------------------
  dim_names <- names(nc_in$dim)

  if (length(dim_names) == 0L) {
    stop("NetCDF file contains no dimensions.")
  }

  nc_dim <- lapply(nc_in$dim, function(d) d$vals)
  names(nc_dim) <- dim_names

  ## ---------------------------------------------------------------------------
  ## Spatial reference
  ## ---------------------------------------------------------------------------
  if (!spatial.ref %in% names(nc_in$var)) {
    stop("Spatial reference variable '", spatial.ref, "' not found in NetCDF.")
  }

  spref_attr <- ncdf4::ncatt_get(nc_in, spatial.ref)

  if (is.null(spref_attr$x_dim) || is.null(spref_attr$y_dim)) {
    stop("Spatial reference variable '", spatial.ref,
         "' must have 'x_dim' and 'y_dim' attributes.")
  }

  dim_map <- list(
    x    = spref_attr$x_dim,
    y    = spref_attr$y_dim,
    time = "time"
  )

  ## ---------------------------------------------------------------------------
  ## Time -> Date
  ## ---------------------------------------------------------------------------
  time_units <- ncdf4::ncatt_get(nc_in, "time")$units
  origin_date <- tryCatch(
    as.Date(sub("^\\D+", "", time_units)),
    error = function(e) stop("Time units do not contain a valid origin date.")
  )

  if (!"time" %in% names(nc_dim)) {
    stop("Time dimension 'time' not found in NetCDF.")
  }

  time_vals <- nc_dim[["time"]]

  if (leap.days) {
    datev <- origin_date + time_vals
  } else {
    # User indicates data should have no leap days (365-day calendar)
    all_dates <- origin_date + time_vals

    # Check if any Feb 29s exist in the computed dates
    has_feb29 <- any(format(all_dates, "%m-%d") == "02-29")

    if (has_feb29) {
      # Data has leap days but user specified leap.days = FALSE
      # This is a mismatch - warn but use all dates to match time dimension
      warning("leap.days = FALSE specified but data contains Feb 29 dates. ",
              "Using all dates from file to match time dimension. ",
              "Consider setting leap.days = TRUE if data includes leap years.")
    }

    datev <- all_dates
  }

  ## ---------------------------------------------------------------------------
  ## Variable selection
  ## ---------------------------------------------------------------------------
  all_data_vars <- setdiff(names(nc_in$var), spatial.ref)

  if (length(all_data_vars) == 0L) {
    stop("No data variables found in NetCDF (excluding '", spatial.ref, "').")
  }

  if (is.null(variables)) {
    # Load all data variables
    nc_var_names_original <- all_data_vars
  } else {
    if (!is.character(variables)) {
      stop("`variables` must be NULL or a character vector of variable names.")
    }
    missing_vars <- setdiff(variables, all_data_vars)
    if (length(missing_vars) > 0L) {
      stop("Variables not found in NetCDF: ",
           paste(missing_vars, collapse = ", "))
    }
    nc_var_names_original <- variables
  }

  ## ---------------------------------------------------------------------------
  ## Variable renaming (robust subset)
  ## ---------------------------------------------------------------------------
  if (!is.null(var_rename)) {
    if (!is.character(var_rename) || is.null(names(var_rename))) {
      stop("`var_rename` must be a named character vector, e.g. c(old = 'new').")
    }

    rename_from <- names(var_rename)

    # All names in var_rename must be among the selected variables
    invalid_rename <- setdiff(rename_from, nc_var_names_original)
    if (length(invalid_rename) > 0L) {
      stop(
        "`var_rename` contains names not in selected variables: ",
        paste(invalid_rename, collapse = ", ")
      )
    }

    # Check for duplicate target names
    rename_to <- unname(var_rename)
    if (any(duplicated(rename_to))) {
      stop("`var_rename` contains duplicate target names: ",
           paste(rename_to[duplicated(rename_to)], collapse = ", "))
    }

    nc_var_names_final <- nc_var_names_original
    idx <- match(rename_from, nc_var_names_original)
    nc_var_names_final[idx] <- rename_to

    # Warn if renaming creates conflicts with other variable names
    conflicts <- intersect(rename_to, setdiff(nc_var_names_final, rename_to))
    if (length(conflicts) > 0L) {
      warning("Variable renaming creates name conflicts with existing variables: ",
              paste(conflicts, collapse = ", "))
    }

  } else {
    nc_var_names_final <- nc_var_names_original
  }

  nv <- length(nc_var_names_original)

  ## ---------------------------------------------------------------------------
  ## Read variables into 4D array: (x, y, time, var)
  ## ---------------------------------------------------------------------------
  first_var <- nc_var_names_original[1L]
  raw_first <- ncdf4::ncvar_get(nc_in, first_var)
  dim_names_first <- ncGetDimNames(nc_in, first_var)

  perm_first <- match(c(dim_map$x, dim_map$y, dim_map$time), dim_names_first)
  if (any(is.na(perm_first))) {
    stop(
      "Variable '", first_var,
      "' does not have expected dimensions (must include x, y, time)."
    )
  }

  arr_first <- aperm(raw_first, perm_first)
  if (!is.null(signif.digits)) {
    arr_first <- signif(arr_first, signif.digits)
  }

  nx <- dim(arr_first)[1L]
  ny <- dim(arr_first)[2L]
  nt <- dim(arr_first)[3L]

  var_array <- array(
    NA_real_,
    dim = c(nx, ny, nt, nv),
    dimnames = list(NULL, NULL, NULL, nc_var_names_final)
  )
  var_array[,,,1L] <- arr_first

  if (nv > 1L) {
    for (k in 2L:nv) {
      vname <- nc_var_names_original[k]
      raw_k <- ncdf4::ncvar_get(nc_in, vname)
      dim_names_k <- ncGetDimNames(nc_in, vname)
      perm_k <- match(c(dim_map$x, dim_map$y, dim_map$time), dim_names_k)
      if (any(is.na(perm_k))) {
        stop(
          "Variable '", vname,
          "' does not have expected dimensions (must include x, y, time)."
        )
      }
      arr_k <- aperm(raw_k, perm_k)
      if (!is.null(signif.digits)) {
        arr_k <- signif(arr_k, signif.digits)
      }
      dim(arr_k) <- c(nx, ny, nt)
      var_array[,,,k] <- arr_k
    }
  }

  ## ---------------------------------------------------------------------------
  ## Attributes (read by original names, stored under final names)
  ## ---------------------------------------------------------------------------
  var_attribs <- lapply(nc_var_names_original, function(v) {
    ncdf4::ncatt_get(nc_in, v)
  })
  names(var_attribs) <- nc_var_names_final

  global_attr <- ncdf4::ncatt_get(nc_in, 0)

  ## ---------------------------------------------------------------------------
  ## Grid indices from first var/time slice
  ## ---------------------------------------------------------------------------
  slice0 <- var_array[,,1L,1L, drop = FALSE]
  dim(slice0) <- c(nx, ny)

  idx <- which(!is.na(slice0), arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    warning("No non-NA grid cells found in first time step of first variable. ",
            "Checking all time steps...")
    # Check across all time steps for the first variable
    all_slices <- var_array[,,,1L, drop = FALSE]
    dim(all_slices) <- c(nx, ny, nt)
    has_data <- apply(all_slices, c(1L, 2L), function(x) any(!is.na(x)))
    idx <- which(has_data, arr.ind = TRUE)
    if (nrow(idx) == 0L) {
      stop("No non-NA grid cells found in any time step of first variable.")
    }
  }

  x_vals <- nc_dim[[dim_map$x]]
  y_vals <- nc_dim[[dim_map$y]]

  grid_df <- tibble::tibble(
    id   = seq_len(nrow(idx)),
    xind = idx[, 1L],
    yind = idx[, 2L],
    x    = x_vals[idx[, 1L]],
    y    = y_vals[idx[, 2L]]
  )

  ## ---------------------------------------------------------------------------
  ## Build list of tibbles (time series per grid cell) - optimized
  ## ---------------------------------------------------------------------------
  n_cells <- nrow(grid_df)
  data_list <- vector("list", n_cells)
  indices <- as.matrix(grid_df[, c("xind", "yind")])

  for (i in seq_len(n_cells)) {
    # Direct extraction - faster than drop=FALSE + dim()
    cell_vals <- matrix(
      var_array[indices[i, 1L], indices[i, 2L], , ],
      nrow = nt,
      ncol = nv,
      dimnames = list(NULL, nc_var_names_final)
    )

    # Check if empty during creation (if omit.empty = TRUE)
    if (omit.empty && all(is.na(cell_vals))) {
      data_list[[i]] <- NULL  # Mark for removal
    } else {
      data_list[[i]] <- tibble::as_tibble(cell_vals)
    }
  }

  ## ---------------------------------------------------------------------------
  ## Remove empty cells if requested
  ## ---------------------------------------------------------------------------
  if (omit.empty) {
    keep <- !vapply(data_list, is.null, logical(1L))
    if (!any(keep)) {
      warning("All grid cells are empty (all NA values). Returning empty result.")
    }
    data_list <- data_list[keep]
    grid_df <- grid_df[keep, , drop = FALSE]
    grid_df$id <- seq_len(nrow(grid_df))
  }

  ## ---------------------------------------------------------------------------
  ## Return
  ## ---------------------------------------------------------------------------
  list(
    data = data_list,
    grid = grid_df,
    date = datev,
    dimensions = nc_dim,
    attributes = c(
      var_attribs,
      stats::setNames(list(spref_attr), spatial.ref),
      global = list(global_attr)
    )
  )
}
