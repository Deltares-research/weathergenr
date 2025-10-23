#' Quantify Plausible Ranges of Impact metrics Based on GCM Projections
#'
#' @description
#' Calculates the plausible range of climate impact metrics (e.g., system performance, failure probability)
#' for a set of locations and metrics, based on the overlap between user-provided climate response surfaces (`str.data`)
#' and the distribution of future climate change as projected by General Circulation Models (GCMs, `gcm.data`).
#' For each metric and location, the function extracts the range of metric values falling within GCM ellipses
#' (e.g., 50% and 95% probability contours) in delta-precipitation/delta-temperature space.
#'
#' The main application is to support "plausibility" or "stress-testing" analyses in climate risk assessments.
#'
#' @param str.data Data frame. Climate response surface with columns for x (precip), y (temp), and z (metric of interest)
#'   plus `statistic` and location identifier. Typically produced by a climate stress-test workflow.
#' @param gcm.data Data frame. GCM projections, must include columns for `prcp`, `tavg`, `scenario`, and optionally `location`.
#' @param gcm.scenario.list Character vector. GCM scenario names to include (default: c("ssp126", "ssp245", "ssp370", "ssp585")).
#' @param clevel.list Numeric vector. List of probability levels for GCM ellipses (e.g., c(0.5, 0.95)).
#' @param metric.list Character vector. List of metric names (e.g., c("reliability", "failure_prob")) to summarize.
#' @param metric.labs Character vector. Labels for metrics (used in output; defaults to metric.list).
#' @param location.list Character vector. locations or system units to summarize (must match column in str.data/gcm.data).
#'
#' @return
#' A data frame (tibble) with the following columns:
#'   - `location`: location name.
#'   - `metric`: metric name.
#'   - `Baseline`: metric value at baseline (no climate change).
#'   - Additional columns for each probability level in `clevel.list` (e.g., "CL:50%", "CL:95%"), containing
#'     the plausible range ("min to max") of metric values within the GCM ellipse for each location and metric.
#'
#' @details
#' For each location/metric, interpolates the response surface, finds points within the specified GCM probability ellipse
#' (in climate space), and records the range of metric values within that region. Useful for reporting GCM-informed
#' plausible ranges of change for each system/metric.
#'
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr expand_grid unite pivot_wider
#' @importFrom akima interp
#' @importFrom sp point.in.polygon
#'
#' @examples
#' \dontrun{
#' # Example data setup (must match your workflow's response surfaces & GCM format)
#' str.data <- data.frame(
#'   prcp = rep(seq(-20, 20, 5), each = 9),
#'   tavg = rep(seq(-2, 6, 1), times = 9),
#'   z = runif(81, 70, 100),
#'   statistic = "reliability",
#'   location = "SiteA"
#' )
#' gcm.data <- data.frame(
#'   prcp = rnorm(100, 0, 10),
#'   tavg = rnorm(100, 2, 2),
#'   scenario = sample(c("ssp126", "ssp245", "ssp370", "ssp585"), 100, replace = TRUE)
#' )
#' results <- GCMplausiblity(
#'   str.data = str.data,
#'   gcm.data = gcm.data,
#'   metric.list = c("reliability"),
#'   location.list = c("SiteA")
#' )
#' print(results)
#' }
GCMplausiblity <- function(
    str.data = NULL,
    gcm.data = NULL,
    gcm.scenario.list = c("ssp126", "ssp245", "ssp370", "ssp585"),
    clevel.list = c(0.50, 0.95),
    metric.list = NULL,
    metric.labs = NULL,
    location.list = NULL) {


  # Custom function for interpolation between grid cells
  gridInterpolate <- function(x, y, z = NULL, resolution = 100, ...) {

    # Interpolation for three-dimensional array
    if (is.null(z)) {
      z <- rep(0, length(x))
    }

    z <- data.frame(z)

    df1 <- lapply(seq_len(ncol(z)), function(i) {
      akima::interp(x, y, z[, i],
        xo = seq(min(x), max(x), length = resolution),
        yo = seq(min(y), max(y), length = resolution)
      )
    }, ...)

    df2 <- do.call("cbind", lapply(df1, function(x) c(x$z)))
    df3 <- cbind(expand.grid(x = df1[[1]]$x, y = df1[[1]]$y), df2)
  }

  # Surpress warnings
  options(warn = -1)

  # if (is.null(metric.labs)) metric.labs <- metric.list

  clevel_labs <- paste0("CL:", clevel.list * 100, "%")
  plausDF <- expand_grid(location = location.list, clevel = clevel.list, metric = metric.list) %>%
    mutate(baseline = 0, mean = NA, low = 0, high = 0)

  # Filter selected GCM scenarios
  gcm.data <- gcm.data %>%
    filter(scenario %in% gcm.scenario.list) %>%
    rename(x = prcp, y = tavg)


  for (x in 1:nrow(plausDF)) {


    strDF_ini <- str.data %>%
      filter(statistic == plausDF$metric[x]) %>%
      select(x = prcp, y = tavg, z)

    bindex <- which(strDF_ini$x == 0 & strDF_ini$y == 0)
    if (length(bindex) == 0) stop("No baseline point (x=0, y=0) found in str.data for metric ", plausDF$metric[x])

    strDF <- strDF_ini #%>% mutate(z = z / strDF_ini[[bindex, "z"]] * 100 - 100)

    strDF_interp <- gridInterpolate(strDF$x, strDF$y, strDF$z) %>%
      as_tibble() %>%
      rename(x = 1, y = 2, z = 3)

    # Extract points from the response surface
    p <- ggplot(strDF_interp, aes(x, y)) +
      geom_point(aes())
    points <- ggplot_build(p)$data[[1]]

    # Extract ellipse points from GCMs
    p1 <- p + stat_ellipse(data = gcm.data, aes(x, y), level = plausDF$clevel[x], type = "norm")
    ell <- ggplot_build(p1)$data[[2]]

    # Find intersecting points on the ellipse
    con <- which(as.logical(sp::point.in.polygon(points$x, points$y, ell$x, ell$y)))

    if (length(con) == 0) {
      plausDF$low[x] <- NA
      plausDF$high[x] <- NA
    } else {
      plausDF$low[x] <- min(strDF_interp[con, ]$z)
      plausDF$high[x] <- max(strDF_interp[con, ]$z)
    }

    plausDF$baseline[x] <- strDF_ini$z[which(strDF_ini$x == 0 & strDF_ini$y == 0)]

    #Find mean
    x1 <- mean(gcm.data$x)
    y1 <- mean(gcm.data$y)
    distances <- sqrt((strDF_interp$x - x1)^2 + (strDF_interp$y - y1)^2)
    closest_index <- which.min(distances)
    plausDF$mean[x] <- strDF_interp$z[closest_index]

  }

  plausDF %>%
    select(location, metric, clevel, baseline, mean, low, high) %>%
    pivot_wider(
      id_cols = c(location, metric, baseline, mean),
      names_from = clevel,
      values_from = c("low", "high"),
      names_glue = "{.value}_{clevel}")

}
