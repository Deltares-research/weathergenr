#' Batch create climate response surfaces from CST Toolbox applications
#'
#' @param str.data placeholder
#' @param gcm.data placeholder
#' @param project.dir placeholder
#' @param variable.y.breaks placeholder
#' @param variable.x.breaks placeholder
#' @param gcm.future.period placeholder
#' @param metrics placeholder
#' @param metric.direction placeholder
#' @param metric.labels placeholder
#'
#' @return
#' @export
#' @import dplyr
#' @import ggplot2
#'
#' @examples
batchClimateSurface <- function(
  project.dir = NULL,
  str.data = NULL,
  gcm.data = NULL,
  gcm.legend = TRUE,
  gcm.future.period = "near",
  metrics = c("mean", "min", "max", "q95", "Q7day_max", "Q7day_min",
    "drymonth_mean", "wetmonth_mean"),
  metric.labels = c(bquote(bold('Mean flow '(m^3/s))),
     bquote(bold('Minimum daily flow '(m^3/s))),
     bquote(bold('Maximum daily flow '(m^3/s))),
     bquote(bold('95th quantile flow '(m^3/s))),
     bquote(bold('7-Day maximum flow '(m^3/s))),
     bquote(bold('7-Day minimum flow '(m^3/s))),
     bquote(bold('Mean flow during driest month '(m^3/s))),
     bquote(bold('Mean flow during wettest month '(m^3/s)))),
  metric.direction = c(1,1,0,1,0,1,1,0),
  ...)

{

  if(is.null(project.dir)) project.dir <- paste0(getwd(), "/climSurfaces")


  if(!is.null(gcm.data)) {
    gcm_df <- gcm.data %>% filter(horizon %in% gcm.future.period)

  } else {
    gcm_df <- NULL
  }

  # Post-process input data
  locations <- str.data %>% select(-statistic, -tavg, -prcp) %>% colnames()

  # Path of output plots
  dir.create(project.dir, showWarnings = F)

  # Loop through each location
  for (l in 1:length(locations)) {

    # Create a directory for current location
    dir.create(paste0(project.dir, locations[l]),  showWarnings = FALSE)

    ### SET CURRENT METRIC
    for (m in 1:length(metrics)) {

      p <- climateSurface(
          str.data = str.data %>% filter(statistic == metrics[m]),
          gcm.data = gcm_df,
          variable.x = "prcp",
          variable.y = "tavg",
          variable.z = locations[l],
          plot.title = metric.labels[[m]],
          failure.direction = metric.direction[m],
          ...)

      # Save to file
      ggsave(paste0(project.dir, locations[l],"/", locations[l], "_",
                    metrics[m], ".png"),
             width = 6.63, height = 6.63, dpi = 700, scale = 1)

   } #metrics close
  } #locations close

  message(paste0("All plots saved to: \n", project.dir))

}

