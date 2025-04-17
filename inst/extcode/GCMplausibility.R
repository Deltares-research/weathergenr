

#' Estimate GCM plausibility ranges
#'
#' @param str.data placeholder
#' @param gcm.data placeholder
#' @param clevel.list placeholder
#' @param metric.list placeholder
#' @param metric.labs placeholder
#' @param location.list placeholder
#'
#' @returns
#' @import sp
#' @import dplyr
#' @import tidyr
#' @import akima
#' @import ggplot2
#' @export
#'
#' @examples
GCMplausiblity <- function(
    str.data = NULL,
    gcm.data = NULL,
    gcm.scenario.list = c("ssp126", "ssp245", "ssp370", "ssp585"),
    clevel.list = c(0.50, 0.95),
    metric.list = NULL,
    metric.labs = NULL,
    location.list = NULL)

{

  gridInterpolate <- function(x, y, z = NULL, resolution = 100, ...) {

    # Interpolation for three-dimensional array
    if (is.null(z)) {z <- rep(0, length(x))}

    z <- data.frame(z)

    df1 <- lapply(seq_len(ncol(z)), function(i) akima::interp(x, y, z[, i],
                                                              xo = seq(min(x), max(x), length = resolution),
                                                              yo = seq(min(y), max(y), length = resolution)), ...)

    df2 <- do.call("cbind", lapply(df1, function(x) c(x$z)))
    df3 <- cbind(expand.grid(x = df1[[1]]$x, y = df1[[1]]$y), df2)

  }


  # Surpress warnings
  options(warn=-1)

  if (is.null(metric.labs)) metric.labs <- metric.list

  clevel_labs <- paste0("CL:",clevel.list*100,"%")
  plausDF <- expand_grid(Location = location.list, clevel = clevel.list, Metric = metric.list) %>%
    mutate(baseline = 0, min = 0, max = 0)

  #Post-process gcm data
  gcmDF <- gcm.data %>% filter(scenario %in% gcm.scenario.list) %>% na.omit()

  # Find most-likely GCM
  ind1 <- with(gcmDF, which.min(abs(x - median(x))))
  ind2 <- with(gcmDF, which.min(abs(y - median(y))))
  median_x <- mean(gcmDF$x[c(ind1, ind2)])
  median_y <- mean(gcmDF$y[c(ind1, ind2)])

  for (x in 1:nrow(plausDF)) {

    strDF_ini <- str.data %>% filter(statistic == plausDF$Metric[x]) %>%
      select(x = prcp, y = tavg, z = plausDF$Location[x])

    bindex <- which(strDF_ini$x == 0 & strDF_ini$y == 0)
    strDF <- strDF_ini %>% mutate(z = z/strDF_ini[[bindex,"z"]] * 100 - 100)

    strDF_interp <- gridInterpolate(strDF$x, strDF$y, strDF$z) %>% as_tibble() %>% rename(x=1, y=2, z=3)

    # Extract points from the response surface
    p <- ggplot(strDF_interp, aes(x, y)) + geom_point(aes())
    points <- ggplot_build(p)$data[[1]]

    # Extract ellipse points from GCMs
    p1 <- p + stat_ellipse(data = gcmDF, aes(x, y), level=plausDF$clevel[x], type = "norm")
    ell <- ggplot_build(p1)$data[[2]]

    # Find intersecting points on the ellipse
    con <- which(as.logical(sp::point.in.polygon(points$x, points$y, ell$x, ell$y)))
    plausDF$Baseline[x] <- strDF_ini$z[which(strDF_ini$x == 0 & strDF_ini$y == 0)]
    plausDF$min[x] <- min(strDF_interp[con,]$z)
    plausDF$max[x] <- max(strDF_interp[con,]$z)

    #Find closest point to median GCM
    which.min(abs(strDF_interp$x - median_x))
    which.min(abs(strDF_interp$y - median_y))

  }

  df2 <- plausDF %>%
    mutate(min = paste0(round(min), "%"),
           max = paste0(round(max), "%")) %>%
    mutate(clevel = factor(clevel, levels = clevel.list, labels = clevel_labs)) %>%
    unite(range, c(min, max), sep = " to ") %>%
    pivot_wider(names_from = clevel, values_from = range,
                id_cols = c(Location, Metric, Baseline))

  return(df2)

}
