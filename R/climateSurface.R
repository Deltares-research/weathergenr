#' Climate Response Surface Generator
#'
#' @param str.data placeholder
#' @param gcm.data placeholder
#' @param variable.x placeholder
#' @param variable.y placeholder
#' @param variable.z placeholder
#' @param threshold.z placeholder
#' @param variable.x.label placeholder
#' @param variable.y.label placeholder
#' @param plot.title placeholder
#' @param failure.direction placeholder
#' @param gcm.bivariate.dist placeholder
#' @param gcm.transparency placeholder
#' @param gcm.legend placeholder
#' @param variable.z.min placeholder
#' @param variable.z.max placeholder
#' @param variable.z.bin placeholder
#' @param variable.z.min.legend placeholder
#' @param variable.z.max.legend placeholder
#' @param variable.z.bin.legend placeholder
#' @param variable.x.breaks placeholder
#' @param variable.y.breaks placeholder
#' @param text.scale placeholder
#' @param gcm.scenario.list placeholder
#'
#' @return
#' @export
#' @import ggplot2
#' @import dplyr
#'
#' @examples
climateSurface <- function(
    str.data = NULL,
    gcm.data = NULL,
    variable.x = NULL,
    variable.y = NULL,
    variable.z = NULL,
    threshold.z = NULL,
    plot.title = "change in variable",
    variable.x.label = "change in temperature",
    variable.y.label = "change in precipitation",
    failure.direction = 1,
    gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
    gcm.bivariate.dist = FALSE,
    gcm.transparency = 0.85,
    gcm.legend = TRUE,
    variable.z.min = NULL,
    variable.z.max = NULL,
    variable.z.bin = 15,
    variable.z.min.legend = NULL,
    variable.z.max.legend = NULL,
    variable.z.bin.legend = 11,
    variable.x.breaks = NULL,
    variable.y.breaks = NULL,
    text.scale = 0.8)

  {

    ## GGplot themes
    gg_theme_surface <- function(size = 18 * text.scale) {
      theme_light() %+replace%
      theme(
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = size+2, hjust = 0, margin = margin(0,0,-10,0)),
        axis.ticks = element_line(colour = "gray60"),
        axis.title = element_text(size = size),
        axis.text = element_text(size = size-2),
        legend.position="top",
        legend.direction="horizontal",
        legend.title.position = "top",
        legend.text = element_text(size = size-2),
        legend.box.margin=margin(-2,-2,-2,-2)
      )
    }

    gg_theme_blank <- theme(axis.ticks = element_blank(),
                            panel.background = element_blank(),
                            axis.text.x      = element_blank(),
                            axis.text.y      = element_blank(),
                            axis.title.x     = element_blank(),
                            axis.title.y     = element_blank())

    empty <- ggplot() + geom_point(aes(1,1), colour="white") + gg_theme_blank

    # Specify x, y breaks
    if(is.null(variable.x.breaks)) variable.x.breaks  <- unique(str.data[[variable.x]])
    if(is.null(variable.y.breaks)) variable.y.breaks  <- unique(str.data[[variable.y]])

    # Specify z range and breaks
    if(is.null(variable.z.min)) variable.z.min <- min(str.data[[variable.z]])
    if(is.null(variable.z.max)) variable.z.max <- max(str.data[[variable.z]])
    z_breaks <- seq(variable.z.min, variable.z.max, length.out = variable.z.bin+1)

    # Specify z legend
    if(is.null(variable.z.min.legend)) variable.z.min.legend  <- variable.z.min
    if(is.null(variable.z.max.legend)) variable.z.max.legend  <- variable.z.max
    if(is.null(variable.z.bin.legend)) variable.z.bin.legend  <- variable.z.bin

    z_legend_breaks <- pretty(z_breaks, variable.z.bin.legend)
    z_legend_limits <- range(z_legend_breaks)

    # Set threshold to mean if not specified
    if(is.null(threshold.z)) {
      threshold.z  <- str.data %>% filter(tavg == 0, prcp == 0) %>%
        pull(variable.z) %>% mean()
    }

    bin_num <- length(z_legend_breaks) - 1
    mid_bin <- findInterval(threshold.z, z_legend_breaks)
    colpal <- vector("character", length(z_legend_breaks)-1)
    bin_num_lw <- mid_bin-1
    bin_num_up <- length(colpal) - mid_bin

    # Set color palette based on direction of increasing performance
    if (failure.direction == 1) {
      colpal[1:bin_num_lw] <- colorRampPalette(c("#FF3300", "#FEE5D9"))(bin_num_lw)
      colpal[(mid_bin+1):length(colpal)]  <- colorRampPalette(c("#EFF3FF", "#0033FF"))(bin_num_up)
    } else {
      colpal[1:bin_num_lw] <- colorRampPalette(c("#0033FF", "#EFF3FF"))(bin_num_lw)
      colpal[(mid_bin+1):length(colpal)]  <- colorRampPalette(c("#FEE5D9","#FF3300"))(bin_num_up)
    }

    colpal[[mid_bin]] <- "white"

    # Core climate response surface
    p <- ggplot(str.data, aes(x = .data[[variable.x]], y = .data[[variable.y]])) +
      # Define theme
      gg_theme_surface() +
      # Place z dimension
      geom_contour_filled(aes(z = .data[[variable.z]],
            fill = after_stat(level_mid)), breaks = z_legend_breaks) +
      # Place threshold line
      geom_contour(aes(z = .data[[variable.z]]),
                   breaks = threshold.z, color = "black", linewidth = 1) +
      # Set x,y, and fill scales
      scale_x_continuous(expand = c(0, 0), breaks = variable.x.breaks, labels = ~ paste0(.x, "%")) +
      scale_y_continuous(expand = c(0, 0), breaks = variable.y.breaks, labels = ~ paste0(.x, "°C")) +
      scale_fill_gradientn(colors = colpal,
                           breaks = z_legend_breaks,
                           limits = z_legend_limits,
                           guide = guide_colorbar(barwidth=25, show.limits=TRUE, ticks.colour = "black",
                                              barheight = 1.30*text.scale, order = 1,
                                              draw.ulim = TRUE, draw.llim = TRUE)) +
      # Set labs
      labs(x = variable.x.label,y = variable.y.label,
        color = "Climate\nProjections", fill = "", title = plot.title)


      ######## GCM Dots ########################################################

      if(!is.null(gcm.data)) {

        gcm_scenario_color <- c("ssp126" = "#003466", "ssp245"	="#f69320",
                                "ssp370"	="#df0000", "ssp585"	="#980002")

        gcm_period_shape <- c("near" = 1, "far" = 4)

        p <- p + geom_point(mapping = aes(x = .data[[variable.x]], y = .data[[variable.y]],
            color = scenario, shape = horizon), data = gcm.data,
            stroke = 1.5, size = 2*text.scale, alpha = gcm.transparency) +
            scale_color_manual(values = gcm_scenario_color) +
            scale_shape_manual(values = gcm_period_shape) +
            labs(shape = "Time\nHorizon")


        # Draw GCM legend
        if(isTRUE(gcm.legend)) {

            # Set legend for GCM color
            p <- p +  guides(color = guide_legend(order = 2, position = "right", direction = "vertical"),
                             shape = guide_legend(order = 3, position = "right", direction = "vertical"))


          } else {

            p <- p +  guides(color = guide_legend(order = 2, position = "right",
                                                  direction = "vertical",
                                                  override.aes = list(alpha = 0),
                                                  theme = theme(legend.title = element_text(color = "transparent"),
                                                                legend.text = element_text(color = "transparent"))),
                             shape = guide_legend(order = 3, position = "right",
                                                  direction = "vertical",
                                                  override.aes = list(alpha = 0),
                                                  theme = theme(legend.title = element_text(color = "transparent"),
                                                                legend.text = element_text(color = "transparent"))))
          }

        if(isTRUE(gcm.bivariate.dist)) {

          p = p + stat_ellipse(aes(x = .data[[variable.x]], y = .data[[variable.y]]), gcm.data, type = "norm", level=0.50) +
            stat_ellipse(aes(x = .data[[variable.x]], y = .data[[variable.y]]), gcm.data, type = "norm", level=0.75) +
            stat_ellipse(aes(x = .data[[variable.x]], y = .data[[variable.y]]), gcm.data, type = "norm", level=0.90) +
            stat_ellipse(aes(x = .data[[variable.x]], y = .data[[variable.y]]), gcm.data, type = "norm", level=0.95)
        }

      } #gcm-data close

    return(p)

  }

