
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
  coordGrid = NULL,
  file.suffix = "",
  file.prefix = "clim_change_rlz",
  output.path = NULL,
  sim.year.start = NULL,
  month.start = 1,
  variables = NULL,
  variable.units = NULL,
  nc.dimnames = NULL,
  change.settings = NULL,
  save.scenario.matrix = TRUE,
  step = TRUE)

{


  if(is.null(variable.units)) {variable.units <- rep("", length(variables))}

  #:::::::::::::::::::::: CLIMATE CHANGE SETTINGS ::::::::::::::::::::::::::::::

  perturb1 <- readxl::read_excel(change.settings, sheet = "par1", range = "B3:D36")
  perturb2 <- readxl::read_excel(change.settings, sheet = "par2", range = "B3:D36")

  PARCC <- list(par1 = list(), par2 = list())

  PARCC$par1 <- list(
    name = perturb1[[1,2]],
    change_op <- perturb1[[2,2]],
    increments = as.numeric(perturb1[[3,2]]),
    mean = list(min = as.numeric(perturb1[7:18,1:3] %>% pull(2)),
                max = as.numeric(perturb1[7:18,1:3] %>% pull(3)),
                obs = as.numeric(perturb1[7:18,1:3] %>% pull(1))),
    var  = list(min = as.numeric(perturb1[22:33,1:3] %>% pull(2)),
                max = as.numeric(perturb1[22:33,1:3] %>% pull(3)),
                obs = as.numeric(perturb1[22:33,1:3] %>% pull(1)))
  )
  PARCC$par2 <- list(
    name = perturb2[[1,2]],
    change_op <- perturb2[[2,2]],
    increments = as.numeric(perturb2[[3,2]]),
    mean = list(min = as.numeric(perturb2[7:18,1:3] %>% pull(2)),
                max = as.numeric(perturb2[7:18,1:3] %>% pull(3)),
                obs = as.numeric(perturb2[7:18,1:3] %>% pull(1))),
    var  = list(min = as.numeric(perturb2[22:33,1:3] %>% pull(2)),
                max = as.numeric(perturb2[22:33,1:3] %>% pull(3)),
                obs = as.numeric(perturb2[22:33,1:3] %>% pull(1)))
  )

  PARCC$par1$sind <- 1:PARCC$par1$increments
  PARCC$par2$sind <- 1:PARCC$par2$increments

  for (i in 1:length(PARCC)) {

    PARCC[[i]]$mean$steps <- sapply(1:12, function(m)
      seq(PARCC[[i]]$mean$min[m], PARCC[[i]]$mean$max[m], length.out = PARCC[[i]]$increments)) %>%
      rbind(rep(NA,12))

    PARCC[[i]]$var$steps <-   sapply(1:12, function(m)
      seq(PARCC[[i]]$var$min[m], PARCC[[i]]$var$max[m], length.out = PARCC[[i]]$increments)) %>%
      rbind(rep(NA,12))
  }

  # Scenario matrix
  scn_mat <- tidyr::expand_grid(par1 = 1:PARCC[[1]]$increments, par2 = 1:PARCC[[2]]$increments) %>%
    mutate(id = 1:n(), .before = 1) %>%
    mutate(precip_change_mean = PARCC[[1]]$mean$steps[par1,1], temp_change_mean = PARCC[[2]]$mean$steps[par2,1]) %>%
    mutate(precip_change_variance = PARCC[[1]]$var$steps[par1,1], temp_change_variance = PARCC[[2]]$var$steps[par2,1])

  if(isTRUE(save.scenario.matrix)) {
    write.csv(x = scn_mat, file = paste0(output.path, "scenario_matrix.csv"))
  }

  smax <- nrow(scn_mat)

  message(cat("\u2713", "|", "scenario matrix created: ", smax, "scenarios in total"))

  #::::::::::::::::::::::: TEMPLATE FOR WRITING TO NETCDF ::::::::::::::::::::::

  # Number of grids
  grids  <- coordGrid$id
  ngrids <- length(grids)

  # Date indices
  sim_year_end <- sim.year.start + nrow(input.data[[1]])/365
  date_sim <- seq(as.Date(paste(sim.year.start,"-1-01",sep="")),
    as.Date(paste(sim_year_end,"-12-31",sep="")), by="day")

  sim_dates_d <- tibble(year = as.numeric(format(date_sim,"%Y")),
          wyear = getWaterYear(date_sim, month.start),
          month = as.numeric(format(date_sim,"%m")),
          day = as.numeric(format(date_sim,"%d"))) %>%
      filter(month!=2 | day!=29) %>%
      filter(wyear >= sim.year.start+1 & wyear <= sim_year_end) %>%
      mutate(date = as.Date(paste(year, month, day, sep = "-")), .before=1)

  year_series <- as.numeric(format(sim_dates_d$date,"%Y"))
  month_series <- as.numeric(format(sim_dates_d$date,"%m"))
  year_index <- year_series - min(year_series) + 1
  year_num <- length(unique(year_series))


  #Create output directory if doesn't exist
  if (!dir.exists(output.path)) {dir.create(output.path)}

  # Loop through each scenario
  for (s in 1:smax) {

    # New object to store the results
    rlz <- input.data

    # Current perturbation scenario for each variable
    perturb_precip <- list(mean = PARCC[[1]]$mean$steps[scn_mat$par1[s],],
           var = PARCC[[1]]$var$steps[scn_mat$par1[s],])

    perturb_temp <- list(mean = PARCC[[2]]$mean$steps[scn_mat$par2[s],])

    temp_delta_factors <- sapply(1:12, function(x)
      seq(0, perturb_temp$mean[x], length.out = year_num))

    temp_deltas <- sapply(1:length(year_series), function(x)
      temp_delta_factors[year_index[x], month_series[x]])

    # Loop through each grid cell
    for (x in 1:ngrids) {

      # Perturb daily precipitation using quantile mapping
      rlz[[x]]$precip <- quantileMapping(value = rlz[[x]]$precip,
        mon.ts = month_series, year.ts = year_index, par = perturb_precip, step.change = step)

      # Perturb temp, temp_min, and temp_max by delta factors
      rlz[[x]]$temp <- rlz[[x]]$temp + temp_deltas
      rlz[[x]]$temp_min <- rlz[[x]]$temp_min + temp_deltas
      rlz[[x]]$temp_max <- rlz[[x]]$temp_max + temp_deltas

      # Calculate PET from temp, temp_min, temp_max
      rlz[[x]]$pet <- with(rlz[[x]], hargreavesPet(
        months = month_series, temp = temp, tdiff = temp_max - temp_min,
        lat = coordGrid$y[x]))
    }

    # Write to netcdf
    writeNetcdf(
        nc.temp = nc_data,
        data = rlz,
        coord.grid = coordGrid,
        output.path = paste0(output.path,"future/"),
        nc.dimnames = list(x = "lon", y = "lat", time = "time"),
        origin.date =  as.Date(paste0(sim.year.start,"-01-01")),
        calendar.type = "no leap",
        variables = c(variables, "pet"),
        variable.units = c(variable.units, "mm/day"),
        file.prefix = file.prefix,
        file.suffix = s
    )
  }

  message(cat("\u2713", "|", "Results outputted to", s, "netcdf files at: ", out_path))

}

