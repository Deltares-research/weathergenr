
#' Climate Response Surface Generator
#'
#' @param str.data
#' @param gcm.data
#' @param variable.x
#' @param variable.y
#' @param variable.z
#' @param threshold.z
#' @param variable.x.label
#' @param variable.y.label
#' @param plot.title
#' @param failure.direction
#' @param gcm.bivariate.dist
#' @param gcm.transparency
#' @param gcm.legend
#' @param variable.z.min
#' @param variable.z.max
#' @param variable.z.bin
#' @param variable.z.min.legend
#' @param variable.z.max.legend
#' @param variable.z.bin.legend
#' @param variable.x.breaks
#' @param variable.y.breaks
#' @param text.scale
#' @param gcm.scenario.list
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
    failure.direction = "low",
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
        plot.title = element_text(size = size+2, hjust = 0),
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
    z_breaks <- seq(variable.z.min, variable.z.max, length.out = variable.z.bin)

    # Specify z legend
    if(is.null(variable.z.min.legend)) variable.z.min.legend  <- variable.z.min
    if(is.null(variable.z.max.legend)) variable.z.max.legend  <- variable.z.max
    if(is.null(variable.z.bin.legend)) variable.z.bin.legend  <- variable.z.bin

    # Set threshold to mean if not specified
    if(is.null(threshold.z)) {
      threshold.z  <- str.data %>% filter(tavg == 0, prcp == 0) %>%
        pull(variable.z) %>% mean()
    }

    bin_num <- length(z_breaks) - 1

    mid_bin <- findInterval(threshold.z, z_breaks)
    colpal <- vector("character", length(z_breaks)-1)
    colpal[[mid_bin]] <- "white"
    bin_num_lw <- mid_bin-1
    bin_num_up <- length(colpal) - mid_bin

    if (failure.direction == "low") {
      colpal[1:bin_num_lw] <- colorRampPalette(c("#99000D", "#FEE5D9"))(bin_num_lw)
      colpal[(mid_bin+1):length(colpal)]  <- colorRampPalette(c("#EFF3FF", "#004B88"))(bin_num_up)
    } else {
      colpal[1:bin_num_lw] <- colorRampPalette(c("#004B88", "#EFF3FF"))(bin_num_lw)
      colpal[(mid_bin+1):length(colpal)]  <- colorRampPalette(c("#FEE5D9","#99000D"))(bin_num_up)
    }

    # Core climate response surface
    p <- ggplot(str.data, aes(x = .data[[variable.x]], y = .data[[variable.y]])) +
      # Define theme
      gg_theme_surface() +
      # Place z dimension
      geom_contour_filled(aes(z = .data[[variable.z]],
            fill = after_stat(level_mid)), breaks = z_breaks) +
      # Place threshold line
      geom_contour(aes(z = .data[[variable.z]]),
                   breaks = threshold.z, color = "black", size = 1) +
      # Set x,y, and fill scales
      scale_x_continuous(expand = c(0, 0), breaks = variable.x.breaks, labels = ~ paste0(.x, "%")) +
      scale_y_continuous(expand = c(0, 0), breaks = variable.y.breaks, labels = ~ paste0(.x, "°C")) +
      scale_fill_gradientn(colors = colpal,
                           breaks = pretty(z_breaks, variable.z.bin.legend),
                           limits = range(pretty(z_breaks, variable.z.bin.legend)),
                           guide = guide_colorbar(barwidth=25, show.limits=TRUE, ticks.colour = "white",
                                              barheight = 1*text.scale, order = 1,
                                              draw.ulim = TRUE, draw.llim = TRUE)) +
      # Set labs
      labs(x = variable.x.label,y = variable.y.label,
        color = "Climate\nProjections", fill = "", title = plot.title)

      ######## GCM Dots
      gcm_scenario_color <- c("ssp126" = "#003466", "ssp245"	="#f69320",
          "ssp370"	="#df0000", "ssp585"	="#980002")

      if(!is.null(gcm.data)) {

        p <- p + geom_point(mapping = aes(x = .data[[variable.x]], y = .data[[variable.y]],
            color = scenario), data = gcm.data, shape = 1,
            stroke = 1.5, size = 2*text.scale, alpha = gcm.transparency) +
            scale_color_manual(values = gcm_scenario_color)

          if(isTRUE(gcm.legend)) {

            # Set legend for GCM color
            p <- p +  guides(color = guide_legend(order = 2, position = "right",
                                                  direction = "vertical"))

          } else {

            # Set legend for GCM color
            p <- p +  guides(color = guide_legend(order = 2, position = "right",
                                                  direction = "vertical",
                                                  override.aes = list(alpha = 0),
                                                  theme = theme(legend.title = element_text(color = "transparent"),
                                                                legend.text = element_text(color = "transparent"))
            ))

          }


      }

    return(p)

  }


# #Plot marginal distributions
# if(gcm.marginal.dist == TRUE & !is.null(gcm.data)) {
#
#     p_top  <- ggplot(gcm.data, aes(x = .data[[variable.x]]))+
#       scale_fill_manual(values = gcm_scenario_color) +
#       scale_color_manual(values = gcm_scenario_color) +
#       geom_density(aes(fill = scenario, color = scenario), alpha = 0.4, position="identity") +
#       scale_x_continuous(expand = c(0, 0), limits = range(variable.x.breaks)) +
#       gg_theme_blank + guides(fill="none", color = "none")
#
#     p_left  <- ggplot(gcm.data, aes(y = .data[[variable.y]]))+
#       scale_fill_manual(values = gcm_scenario_color) +
#       scale_color_manual(values = gcm_scenario_color) +
#       geom_density(aes(fill = scenario, color = scenario), alpha = 0.4, position="identity") +
#       scale_y_continuous(expand = c(0, 0), limits = range(variable.y.breaks)) +
#       gg_theme_blank + guides(fill="none", color = "none") +
#       scale_x_reverse()
#
#     p <- empty + p_top + p_left + p + plot_layout(ncol = 2, widths = c(1, 7), heights = c(1,7))
#
#   }

