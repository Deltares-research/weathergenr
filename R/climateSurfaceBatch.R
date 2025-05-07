#' Batch create climate response surfaces from CST Toolbox applications
#'
#' @param str.data placeholder
#' @param gcm.data placeholder
#' @param save.dir placeholder
#' @param gcm.future.period placeholder
#' @param metrics placeholder
#' @param metric.direction placeholder
#' @param metric.labels placeholder
#' @param gcm.legend
#' @param locations
#' @param location.labels
#' @param relative.results
#' @param ...
#'
#' @return
#' @export
#' @import dplyr
#' @import ggplot2
#'
#' @examples
climateSurfaceBatch <- function(
  save.dir = NULL,
  str.data = NULL,
  gcm.data = NULL,
  gcm.legend = TRUE,
  gcm.future.period = "near",
  metrics = NULL,
  locations = NULL,
  location.labels = NULL,
  metric.labels = NULL,
  metric.direction = NULL,
  relative.results = TRUE,
  ...)

{

  if(is.null(save.dir)) save.dir <- paste0(getwd(), "/climSurfaces")

  if(!is.null(gcm.data)) {
    gcm_df <- gcm.data %>% filter(horizon %in% gcm.future.period)

  } else {
    gcm_df <- NULL
  }

  if(isTRUE(relative.results)) {

    strdata <- str.data %>%
      pivot_longer(cols = starts_with("Q_"), names_to = "location", values_to = "value") %>%
      group_by(statistic, location) %>%
      mutate(value = value/.data[["value"]][which(.data$tavg == 0 & .data$prcp == 0)] * 100 - 100) %>%
      pivot_wider(names_from = "location", values_from = "value")

  } else {
    strdata <- str.data
  }

  # Path of output plots
  dir.create(save.dir, showWarnings = F)

  # Loop through each location
  for (l in 1:length(locations)) {

    # Create a directory for current location
    dir.create(paste0(save.dir, locations[l]),  showWarnings = FALSE)

    ### SET CURRENT METRIC
    for (m in 1:length(metrics)) {

      # Set plot title
      if(!is.null(location.labels)) {
        plot_title <- paste0(location.labels[l], " | ", metric.labels[[m]])
      } else {
        plot_title <-  metric.labels[[m]]
      }

      p <- climateSurface(
          str.data = strdata %>% filter(statistic == metrics[m]),
          gcm.data = gcm_df,
          variable.x = "prcp",
          variable.y = "tavg",
          variable.z = locations[l],
          plot.title = plot_title,
          failure.direction = metric.direction[m],
          ...)

#     if(isTRUE(show.hist.value)) {
#
#        str_hist <- str.data %>% filter(prcp == 0 & tavg == 0 & statistic == metrics[m]) %>%
#          select(tavg, prcp, value = locations[l])
#
#        p <- p + geom_point(data = str_hist, size = 3, shape = 8, color = "chartreuse2", stroke = 2) +
#          geom_label(aes(label = value), data = str_hist, nudge_x = 0.1, nudge_y = 0.2, fill = "white", color = "black")
#      }

      # Save to file
      ggsave(paste0(save.dir, locations[l],"/", locations[l], "_",
                    metrics[m], ".png"),
             width = 6.63, height = 6.63, dpi = 700, scale = 1)

   } #metrics close
  } #locations close

  message(paste0("All plots saved to: \n", paste0(getwd(), save.dir)))

}

