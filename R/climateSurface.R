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
#' @param gcm.bvnorm.levels placeholder
#' @param gcm.transparency placeholder
#' @param gcm.legend placeholder
#' @param variable.z.min placeholder
#' @param variable.z.max placeholder
#' @param contour.num placeholder
#' @param variable.x.breaks placeholder
#' @param variable.y.breaks placeholder
#' @param text.scale placeholder
#' @param gcm.scenario.list placeholder
#' @param variable.z.breaks placeholder
#' @param multi.panel placeholder
#' @param panel.variable placeholder
#' @param panel.variable.levels placeholder
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
    variable.x.label = expression(Delta~"Precipitation"),
    variable.y.label = expression(Delta~"Temperature"),
    failure.direction = 1,
    gcm.scenario.list = c("rcp26", "rcp45", "rcp60", "rcp85"),
    gcm.bvnorm.levels = FALSE,
    gcm.transparency = 0.75,
    gcm.legend = TRUE,
    variable.z.breaks = NULL,
    variable.x.breaks = NULL,
    variable.y.breaks = NULL,
    contour.num = 15,
    text.scale = 0.6,
    multi.panel = FALSE,
    panel.variable = NULL,
    panel.variable.levels = NULL)

  {

  ############## GENERAL SETTINGS/GGPLOT TEMPLATES #############################

  # Surpress warnings
  options(warn=-1)

  # General ggplot template
  gg_theme_surface <- function(size = 18 * text.scale) {
    theme_bw() %+replace%
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = size+3, hjust = 0, margin = margin(0,0,-15,0)),
      axis.ticks = element_line(colour = "gray60"),
      axis.title = element_text(size = size +1),
      axis.text = element_text(size = size-2),
      strip.background =element_rect(fill="white"),
      strip.text = element_text(size = size),
      panel.spacing = unit(1, "lines"),
      legend.position="top",
      legend.direction="horizontal",
      legend.title.position = "top",
      legend.text = element_text(size = size-2),
      legend.justification.top = "left",
      legend.justification.right = "left",
      legend.box.margin=margin(-2,-2,-2,-2),
      aspect.ratio = 1)

    }

  gg_theme_blank <- theme(axis.ticks = element_blank(),
                          panel.background = element_blank(),
                          axis.text.x      = element_blank(),
                          axis.text.y      = element_blank(),
                          axis.title.x     = element_blank(),
                          axis.title.y     = element_blank())

  empty <- ggplot() + geom_point(aes(1,1), colour="white") + gg_theme_blank

  # CALCULATE PLOTTING PARAMETERS ##############################################

  # Vulnerability threshold
  if(is.null(threshold.z)) {
    threshold.z  <- str.data %>% filter(tavg == 0, prcp == 0) %>%
      pull(variable.z) %>% mean()
  }

  # Multipanel option
  if(isTRUE(multi.panel)) {
    str.data[[panel.variable]] <- factor(str.data[[panel.variable]], levels = panel.variable.levels)
  }

  # x, y breaks
  if(is.null(variable.x.breaks)) variable.x.breaks  <- unique(str.data[[variable.x]])
  if(is.null(variable.y.breaks)) variable.y.breaks  <- unique(str.data[[variable.y]])

  # Z range and breaks
  if(is.null(variable.z.breaks)) {
    variable.z.breaks <- pretty(c(min(str.data[[variable.z]]), max(str.data[[variable.z]])), contour.num)
  }
  z_breaks <- pretty(range(variable.z.breaks), contour.num)

  # Color palette
  color1a <- "#df0000"; color2a <- "#0033FF"
  color1b <- "#FEE5D9"; color2b <- "#FFFFFF"

  bin_num <- length(z_breaks) - 1
  mid_bin <- findInterval(threshold.z - 1e-5, z_breaks)
  bin_num_lw <- mid_bin
  bin_num_up <- bin_num - mid_bin
  colpal <- vector("character", bin_num)

  if (failure.direction == 1) {
    colpal[1:bin_num_lw] <- colorRampPalette(c(color1a, color1b))(bin_num_lw)
    colpal[(mid_bin+1):length(colpal)]  <- colorRampPalette(c(color2b, color2a))(bin_num_up)
  } else {
    colpal[1:bin_num_lw] <- colorRampPalette(c(color2a, color2b))(bin_num_lw)
    colpal[(mid_bin+1):length(colpal)]  <- colorRampPalette(c(color1b,color1a))(bin_num_up)
  }

  # Create Basic response surface
  p <- ggplot(str.data, aes(x = .data[[variable.x]], y = .data[[variable.y]])) +
    # Define theme
    gg_theme_surface() +
    {if(isTRUE(multi.panel)) facet_wrap(~get(panel.variable), ncol = 3)} +
    # Place z dimension
    geom_contour_filled(aes(z = .data[[variable.z]],
          fill = after_stat(level_mid)), breaks = z_breaks) +
    # Place threshold line
    geom_contour(aes(z = .data[[variable.z]]),
                 breaks = threshold.z, color = "black", linewidth = 0.7) +
    # Set x,y, scales
    scale_x_continuous(expand = c(0, 0), breaks = variable.x.breaks, labels = ~ paste0(.x, "%")) +
    scale_y_continuous(expand = c(0, 0), breaks = variable.y.breaks, labels = ~ paste0(.x, "\u00B0", "C")) +

    # Set fill scale
    scale_fill_gradientn(colors = colpal, breaks = variable.z.breaks, limits = range(variable.z.breaks),
                         guide = guide_coloursteps(barwidth = 25,
                                                show.limits=TRUE,
                                                barheight = 1.50*text.scale, order = 1,
                                                frame.linewidth = 0.1,
                                                frame.colour = "gray60",
                                                ticks.colour = "gray60",
                                                ticks.linewidth = 0.5,
                                                draw.ulim = TRUE,
                                                draw.llim = TRUE)) +

    # Additional fine-tuning
    coord_cartesian(xlim=range(variable.x.breaks),
                    ylim = range(variable.y.breaks),
                    expand = FALSE) +
    # Set labs
    labs(x = variable.x.label,y = variable.y.label,
      color = "GCM\nScenarios", fill = "", title = plot.title)

  if(isTRUE(multi.panel)) p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 0.7))


  ######## GCM INFORMATION #####################################################

    if(!is.null(gcm.data)) {

      gcm.data <- gcm.data %>%
        filter(get(variable.x) > min(variable.x.breaks) &
               get(variable.x) < max(variable.x.breaks) &
               get(variable.y) > min(variable.y.breaks) &
               get(variable.y) < max(variable.y.breaks))

      gcm_scenario_color <- c("ssp126" = "#003466", "ssp245"	="#f69320",
                              "ssp370"	="#df0000", "ssp585"	="#980002")

      gcm_period_shape <- c("near" = 1, "far" = 4)

      p <- p + geom_point(mapping = aes(x = .data[[variable.x]], y = .data[[variable.y]],
          color = scenario, shape = horizon), data = gcm.data,
          stroke = 1.5, size = 3*text.scale, alpha = gcm.transparency) +
          scale_color_manual(values = gcm_scenario_color) +
          scale_shape_manual(values = gcm_period_shape) +
          labs(shape = "Time\nHorizon")

      # Draw GCM legend
      if(isTRUE(gcm.legend)) {

          # Set legend for GCM color
          p <- p +  guides(color = guide_legend(order = 2, position = "right", direction = "vertical"),
                           shape = guide_legend(order = 3, position = "right", direction = "vertical"))
        } else {
          p <- p + guides(shape = "none", color = "none")
       }

      if(!is.null(gcm.bvnorm.levels)) {

        for (s in 1:length(gcm.bvnorm.levels)) {
          p = p + stat_ellipse(aes(x = .data[[variable.x]], y = .data[[variable.y]]),
                               gcm.data, level=gcm.bvnorm.levels[s],
                         type = "norm", color = "gray30", linetype = "dashed", size = 0.5)
        }
      }

    } #gcm-data close

  return(p)

}

