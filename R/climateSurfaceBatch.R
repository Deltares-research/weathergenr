#' Batch Generation of Climate Response Surface Plots
#'
#' @description
#' Automatically generates and saves climate response surface plots for multiple locations and metrics,
#' based on stress-test results from the CST Toolbox and projections from General Circulation Models (GCMs).
#' This function loops over user-specified metrics and locations, creating and saving a set of response surface figures to disk,
#' leveraging the `climateSurface()` plotting function.
#'
#' Typically used in automated analysis pipelines for climate impact and risk studies.
#'
#' @param save.dir Character. Directory path where plots will be saved. If `NULL`, defaults to a subdirectory \code{climSurfaces} in the working directory.
#' @param str.data Data frame. Stress-test data, typically containing surface response variables, metric columns, and covariates (`prcp`, `tavg`, etc.).
#' @param gcm.data Data frame. Downscaled GCM projections for plotting scenario points/ellipses (must include columns matching those in `climateSurface()`).
#' @param gcm.legend Logical. If `TRUE`, include a legend for GCM scenarios in the plot.
#' @param gcm.future.period Character or vector. Which GCM future horizons to include (e.g., "near", "far").
#' @param metrics Character vector. Names of performance metrics/statistics to plot (must match `statistic` column in `str.data`).
#' @param locations Character vector. Names of locations or columns in `str.data` for which to create plots.
#' @param location.labels Character vector. Human-readable labels for locations (for plot titles). Optional.
#' @param metric.labels Character vector. Human-readable labels for metrics (for plot titles). Optional.
#' @param metric.direction Integer vector. Indicates failure direction for each metric (e.g., 1 = high values are worse, -1 = low values are worse).
#' @param relative.results Logical. If `TRUE` (default), express metric values as percent change from baseline; if `FALSE`, use raw values.
#' @param ... Additional arguments passed to [climateSurface()].
#'
#' @return
#' Invisibly returns `NULL`. Side effect: saves PNG files for all locations and metrics to the specified directory.
#'
#' @details
#' For each location and metric, the function:
#'   - Optionally expresses results as percent change from the baseline (centered at zero change).
#'   - Filters and relabels the data for plotting.
#'   - Calls [climateSurface()] to generate the plot.
#'   - Saves the plot as a high-resolution PNG in a location-specific subfolder.
#'
#' Plot filenames are constructed as: \code{<save.dir>/<location>/<location>_<metric>.png}.
#'
#' @import dplyr
#' @import ggplot2
#'
#' @seealso \code{\link{climateSurface}}
#'
#' @examples
#' \dontrun{
#' # Suppose you have stress-test and GCM data frames `str_data` and `gcm_data`
#' metrics <- c("Reliability", "Resilience")
#' locations <- c("Q_A", "Q_B")
#' metric.labels <- c("System Reliability", "System Resilience")
#' location.labels <- c("Catchment A", "Catchment B")
#' metric.direction <- c(1, 1)
#'
#' climateSurfaceBatch(
#'   save.dir = tempdir(),
#'   str.data = str_data,
#'   gcm.data = gcm_data,
#'   metrics = metrics,
#'   locations = locations,
#'   location.labels = location.labels,
#'   metric.labels = metric.labels,
#'   metric.direction = metric.direction
#' )
#' }
#' @export

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

