
#' Climate Change perturbations function
#'
#' @param output.path To be completed...
#' @param sim.year.start To be completed...
#' @param variables To be completed...
#' @param variable.units To be completed...
#' @param nc.dimnames To be completed...
#' @param change.settings To be completed...
#' @param input.path Placeholder.
#' @param input.data Placeholder.
#' @param file.suffix Placeholder.
#' @param save.scenario.matrix Placeholder.
#'
#' @return
#' @export
#' @import dplyr
#' @import ncdf4
#' @importFrom utils write.csv
#' @importFrom readxl read_excel
#' @importFrom tidyr expand_grid
imposeClimateChanges <- function(
  input.data = NULL,
  grid.coords = NULL,
  precip.changes = NULL,
  temp.changes = NULL,
  file.suffix = "",
  file.prefix = "clim_change_rlz",
  output.path = NULL,
  sim.year.start = NULL,
  month.start = 1,
  variables = NULL,
  variable.units = NULL,
  nc.dimnames = NULL,
  save.scenario.matrix = TRUE,
  step.change = TRUE)

{


  if(is.null(variable.units)) {variable.units <- rep("", length(variables))}

  #:::::::::::::::::::::: CLIMATE CHANGE SETTINGS ::::::::::::::::::::::::::::::

  precip.changes$mean$steps <- sapply(1:12, function(m)
        seq(precip.changes$mean$min[m], precip.changes$mean$max[m],
            length.out = precip.changes$increments))

  precip.changes$var$steps <- sapply(1:12, function(m)
        seq(precip.changes$var$min[m], precip.changes$var$max[m],
            length.out = precip.changes$increments))

  temp.changes$mean$steps <- sapply(1:12, function(m)
        seq(temp.changes$mean$min[m], temp.changes$mean$max[m],
            length.out = temp.changes$increments))

  scn_mat_index <- tidyr::expand_grid(precip_ind = 1:precip.changes$increments,
        temp_ind = 1:temp.changes$increments) %>%
      mutate(ind = 1:n(), .before = 1)

  smax <- nrow(scn_mat_index)

  ##############################################################################

  if(isTRUE(save.scenario.matrix)) {
    write.csv(x = scn_mat_index, file = paste0(output.path, "scenario_matrix.csv"))
  }

  message(cat("\u2713", "|", "scenario matrix created: ", smax, "scenarios in total"))

  #::::::::::::::::::::::: TEMPLATE FOR WRITING TO NETCDF ::::::::::::::::::::::

  # Check starting month !!!!!!!!!!!!!!!! is resampleDates start from wyear start or jan1?

  # Number of grids
  grids  <- grid.coords$id
  ngrids <- length(grids)

  # Date indices
  sim_year_end <- sim.year.start + nrow(input.data[[1]])/365
  date_sim <- seq(as.Date(paste(sim.year.start,"-1-01",sep="")),
    as.Date(paste(sim_year_end,"-12-31",sep="")), by="day")

  sim_dates_d <- tibble(date = date_sim) %>%
    mutate(year = as.numeric(format(date_sim,"%Y")),
           month = as.numeric(format(date_sim,"%m")),
           day = as.numeric(format(date_sim,"%d"))) %>%
    filter(year >= sim.year.start+1 & year <= sim_year_end) %>%
    filter(month!=2 | day!=29)

  year_series <- sim_dates_d$year
  month_series <- sim_dates_d$month
  year_index <- year_series - min(year_series) + 1
  year_num <- length(unique(year_series))


  #Create output directory if doesn't exist
  if (!dir.exists(output.path)) {dir.create(output.path)}

  # Loop through each scenario
  for (s in 1:smax) {

    # New object to store the results
    rlz <- input.data

    # Current perturbation scenario for each variable
    perturb_precip_mean <- precip.changes$mean$steps[scn_mat_index$precip_ind[s],]
    perturb_precip_var <- precip.changes$var$steps[scn_mat_index$precip_ind[s],]
    perturb_temp_mean <- temp.changes$mean$steps[scn_mat_index$temp_ind[s],]

    temp_delta_factors <- sapply(1:12, function(x)
      seq(0, perturb_temp_mean[x], length.out = year_num))

    temp_deltas <- sapply(1:length(year_series), function(x)
      temp_delta_factors[year_index[x], month_series[x]])

    # Loop through each grid cell
    for (x in 1:ngrids) {

      # Perturb daily precipitation using quantile mapping
      rlz[[x]]$precip <- quantileMapping(value = rlz[[x]]$precip,
              mon.ts = month_series, year.ts = year_index,
              mean.change = perturb_precip_mean,
              var.change = perturb_precip_var,
              step.change = step.change)

      # Perturb temp, temp_min, and temp_max by delta factors
      rlz[[x]]$temp <- rlz[[x]]$temp + temp_deltas
      rlz[[x]]$temp_min <- rlz[[x]]$temp_min + temp_deltas
      rlz[[x]]$temp_max <- rlz[[x]]$temp_max + temp_deltas

      # Calculate PET from temp, temp_min, temp_max
      rlz[[x]]$pet <- with(rlz[[x]], hargreavesPet(
        months = month_series, temp = temp, tdiff = temp_max - temp_min,
        lat = grid.coords$y[x]))
    }

    # Write to netcdf
    writeNetcdf(
        nc.temp = nc_data,
        data = rlz,
        coord.grid = grid.coords,
        output.path = paste0(output.path,"future/"),
        nc.dimnames = list(x = "lon", y = "lat", time = "time"),
        origin.date =  sim_dates_d$date[1],
        calendar.type = "no leap",
        variables = c(variables, "pet"),
        variable.units = c(variable.units, "mm/day"),
        file.prefix = file.prefix,
        file.suffix = paste0(file.suffix,"_",s)
    )
  }

  message(cat("\u2713", "|", "Results outputted to", s, "netcdf files at: ", out_path))

}

