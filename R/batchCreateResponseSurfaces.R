#' Batch create climate response surfaces from CST Toolbox applications
#'
#' @param str.data placeholder
#' @param gcm.data placeholder
#' @param project.dir placeholder
#' @param tavg.steps placeholder
#' @param prcp.steps placeholder
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
batchCreateResponseSurfaces <- function(
  str.data = NULL,
  gcm.data = NULL,
  project.dir = NULL,
  tavg.steps = NULL,
  prcp.steps = NULL,
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
  metric.direction = c(1,1,0,1,0,1,1,0)
)

{

  if(is.null(project.dir)) project.dir <- paste0(getwd(), "/climSurfaces")

  # Plot based on GCM information
  if(!is.null(gcm.future.period)) {
    gcm_transparency <- 0.75
    gcm_legend <- TRUE

  } else {
    gcm_legend <- FALSE
    gcm_transparency <- 0
  }

  # Post-process input data
  locations <- str.data %>% select(-statistic, -tavg, -prcp) %>% colnames()

  # Get incremental steps for each uncertain weather variable
  if(is.null(tavg.steps)) tavg.steps <- sort(unique(str.data$tavg))
  if(is.null(prcp.steps)) prcp.steps <- sort(unique(str.data$prcp))

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
          gcm.data = gcm.data %>% filter(horizon %in% gcm.future.period),
          variable.x = "prcp",
          variable.y = "tavg",
          variable.z = locations[l],
          threshold.z = NULL,
          plot.title = metric.labels[[m]],
          variable.y.label = expression(Delta~"Temperature"),
          variable.x.label = expression(Delta~"Precipitation"),
          failure.direction = metric.direction[m],
          gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
          gcm.bivariate.dist = FALSE,
          gcm.transparency = gcm_transparency,
          gcm.legend = gcm_legend,
          variable.z.min = NULL,
          variable.z.max = NULL,
          variable.z.bin = 15,
          variable.z.min.legend = NULL,
          variable.z.max.legend = NULL,
          variable.z.bin.legend = 11,
          variable.x.breaks = prcp.steps,
          variable.y.breaks = tavg.steps,
          text.scale = 0.6
       )

      # Save to file
      ggsave(paste0(project.dir, locations[l],"/", locations[l], "_",
                    metrics[m], ".png"),
             width = 6.7, height = 6, dpi = 700, scale = 1)

   } #metrics close
  } #locations close

  message(paste0("All plots saved to: \n", project.dir))

}

