

################################################################################

# Comparison figures

# 1) Monthly statistics (mean, s.deviation, skew) x num.variables
# 2) Intersite correlations x num.variables
# 3) Cross-correlations x num.variables

#' Function to assess the weather generator skill
#'
#' @param daily.sim A matrix of daily time-series of weather variables
#' @param daily.obs A vector of daily time-series of observed weather series
#' @param output.path  A character string to define the path of resulting netcdf file
#' @param variables A vector object specifying the names of weather variables
#' @param variable.labels A vector object specifying the names of weather variable labels
#' @param variable.units A vector object specifying the names of weather variable units
#' @param realization.num number of weather realizations
#' @param wet.quantile precip. quantile for calculating monthly wet state threshold
#' @param extreme.quantile precip. quantile for calculating monthly extremely wet state threshold
#' @param show.title PLACEHOLDER
#'
#' @return
#' @export
#' @import stats utils ggplot2 dplyr
evaluateWegen <- function(
  daily.sim = NULL,
  daily.obs = NULL,
  output.path = NULL,
  variables = NULL,
  variable.labels = NULL,
  variable.units = NULL,
  realization.num = NULL,
  wet.quantile = 0.2,
  extreme.quantile = 0.8,
  show.title = TRUE)

{

  #Workaround for rlang warning
  year <- mon <- day <- precip <- sd <- variable <- value <- id_variable <- 0
  wet_count <- dry_count <- wet <- dry <- id_variable1 <- id_variable2 <- 0
  rlz <- id1 <- id2 <- variable1 <- variable2 <- type <- Observed <- Stochastic <- 0
  wet_th <- 0

  # General options
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)

  # General parameters
  stat_level <- c("mean", "sd", "skewness")
  stat_label <- c("mean", "standard dev.", "skewness")
  name_obs <- "Observed"
  name_st <- "Simulated"
  nsgrids <- length(daily.sim[[1]])
  num_stats <- length(stat_level)
  num_vars <- length(variables)
  var_combs <- expand_grid(var1=variables, var2 = variables)
  num_var_combs <- choose(num_vars, 2)

  # Historical statistics ######################################################

  hist_year_num <- nrow(daily.obs[[1]])/365

  hist_daily_tidy <- lapply(daily.obs, "[", c("date", variables)) %>%
    bind_rows(.id = "id") %>%
    mutate(year = as.numeric(format(date,"%Y")),
       mon = as.numeric(format(date,"%m")),
       day = as.numeric(format(date,"%d")),
       id = as.numeric(id)) %>%
    select(id, year, mon, day, {{variables}})

  # Calculate dry, wet, and extremely wet day thresholds
  mc_thresholds <- hist_daily_tidy %>% group_by(mon) %>%
    group_by(year, mon, day) %>%
    summarize(precip = mean(precip)) %>%
    group_by(mon) %>%
    summarize(wet_th = stats::quantile(precip, wet.quantile, names = F),
              extreme_th = stats::quantile(precip, extreme.quantile, names = F))

  # Monthly stats per variable
  hist_stats <- hist_daily_tidy %>%
    group_by(id, mon) %>%
    dplyr::summarize(across({{variables}},
      list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = value,-id, -year, -mon) %>%
    separate(variable, c("variable","stat"), sep = ":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label)) %>%
    mutate(type = name_obs)

  # Wet and dry days per month
  hist_wetdry_days <- hist_daily_tidy %>%
    left_join(mc_thresholds, by = "mon") %>%
    group_by(id, mon) %>%
    summarize(wet_count = length(which(precip >= wet_th))/hist_year_num,
              dry_count = length(which(precip < wet_th))/hist_year_num) %>%
    gather(key = stat, value = value, wet_count:dry_count) %>%
    mutate(variable = "precip", type = name_obs) %>%
    select(id, mon, variable, stat, Observed = value)

  # Wet and dry spell lengths per month
  hist_wetdry_spells <- hist_daily_tidy %>%
    left_join(mc_thresholds, by = "mon") %>%
    group_by(id, mon) %>%
    summarize(dry = averageSpellLength(precip, threshold = wet_th, below = TRUE),
              wet = averageSpellLength(precip, threshold = wet_th, below = FALSE)) %>%
    group_by(id, mon) %>%
    summarize(dry = mean(dry),
              wet = mean(wet)) %>%
    gather(key = stat, value = value, wet:dry) %>%
    mutate(variable = "precip", type = name_obs) %>%
    select(id, mon, variable, stat, Observed = value)

  # Intersite/crossite correlations
  hist_allcor_ini <- hist_daily_tidy %>%
    gather(key = variable, value = value, -id, -year, -mon, -day) %>%
    unite(id_variable, c("id","variable"),sep = ":") %>%
    spread(id_variable, value) %>%
    dplyr::select(-year, -mon, -day) %>%
    stats::cor(.$value) %>% as_tibble()

  hist_allcor <- hist_allcor_ini %>%
    mutate(id_variable1 = colnames(hist_allcor_ini), .before = 1) %>%
    gather(key = id_variable2, value = value, -id_variable1) %>%
    separate(id_variable1, c("id1","variable1"), sep =":") %>%
    separate(id_variable2, c("id2","variable2"), sep =":") %>%
    rename(Observed = value)

  hist_stats_aavg <- hist_daily_tidy %>%
    group_by(year, mon) %>%
    summarize(across({{variables}},
      list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = value, -mon, -year) %>%
    separate(variable, c("variable","stat"), sep=":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))  %>%
    mutate(type = name_obs)

  hist_stats_aavg_mon <- hist_stats_aavg %>%
    group_by(mon, variable, stat, type) %>%
    summarize(value = mean(value))


  # Simulated statistics ######################################################

  sim_stats <- vector("list", realization.num)
  sim_stats_aavg <- vector("list", realization.num)
  sim_allcor <- vector("list", realization.num)
  sim_wetdry_days <-vector("list", realization.num)
  sim_wetdry_spells <- vector("list", realization.num)

  sim_year_num <- nrow(daily.sim[[1]][[1]])/365

  #### Calculate statistics per realization
  for (n in 1:realization.num) {

    sim_daily_tidy <- lapply(daily.sim[[n]], "[", c("date", variables)) %>%
      bind_rows(.id = "id") %>%
      mutate(year = as.numeric(format(date,"%Y")),
         mon = as.numeric(format(date,"%m")),
         day = as.numeric(format(date,"%d")),
         id = as.numeric(id)) %>%
      select(id, year, mon, day, {{variables}})

    # Grid-based stats
    sim_stats[[n]] <- sim_daily_tidy %>%
      group_by(id, mon) %>%
      dplyr::summarize(across({{variables}},
        list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
      gather(key = variable, value = value,-id, -year, -mon) %>%
      separate(variable, c("variable","stat"), sep = ":") %>%
      mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

    # Area-averaged stats
    sim_stats_aavg[[n]] <- sim_daily_tidy %>%
      group_by(year, mon) %>%
      summarize(across({{variables}},
        list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
      gather(key = variable, value = value, -mon, -year) %>%
      separate(variable, c("variable","stat"), sep=":") %>%
      mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

    # Intersite/Intervariable correlations
    sim_allcor_ini <- sim_daily_tidy %>%
      gather(key = variable, value = value, -id, -year, -mon, -day) %>%
      unite(id_variable, c("id","variable"),sep = ":") %>%
      spread(id_variable, value) %>%
      dplyr::select(-year, -mon, -day) %>%
      stats::cor(.$value) %>% as_tibble()

    sim_allcor[[n]] <- sim_allcor_ini %>%
      mutate(id_variable1 = colnames(sim_allcor_ini), .before = 1) %>%
      gather(key = id_variable2, value = value, -id_variable1) %>%
      separate(id_variable1, c("id1","variable1"), sep =":") %>%
      separate(id_variable2, c("id2","variable2"), sep =":") #%>%
      #filter(variable1 == variable2) %>%
      #filter(id1 != id2) %>%
      #unite(id, c("id1","id2"), sep=":") %>%
      #filter(id %in% id_combs)

    sim_wetdry_days[[n]] <- sim_daily_tidy %>%
        left_join(mc_thresholds, by = "mon") %>%
        group_by(id, mon) %>%
        summarize(wet_count = length(which(precip >= wet_th))/sim_year_num,
                  dry_count = length(which(precip < wet_th))/sim_year_num) %>%
      gather(key = stat, value = value, wet_count:dry_count) %>%
      mutate(variable = "precip", type = name_st) %>%
      select(id, mon, variable, stat, Simulated = value)


    sim_wetdry_spells[[n]] <- sim_daily_tidy %>%
      left_join(mc_thresholds, by = "mon") %>%
      group_by(id, mon) %>%
      summarize(dry = averageSpellLength(precip, threshold = wet_th, below = TRUE),
                wet = averageSpellLength(precip, threshold = wet_th, below = FALSE)) %>%
      group_by(id, mon) %>%
      summarize(dry = mean(dry),
                wet = mean(wet)) %>%
      gather(key = stat, value = value, wet:dry) %>%
      mutate(variable = "precip", type = name_st) %>%
      select(id, mon, variable, stat, Simulated = value)

    # sim_wetdry_spells_aavg[[n]] <- sim_wetdry_spells[[n]] %>%
    #   group_by(mon, stat, type) %>%
    #   summarize(value = mean(value)) %>%
    #     mutate(type = name_st)
    #
    # sim_wetdry_days_aavg[[n]] <- sim_wetdry_days[[n]] %>%
    #   group_by(mon, stat, type) %>%
    #   summarize(value = mean(value)) %>%
    #   mutate(type = name_st)

  }

  sim_stats <- bind_rows(sim_stats, .id = "rlz")
  sim_stats_aavg <- bind_rows(sim_stats_aavg, .id = "rlz")
  sim_wetdry_days <- bind_rows(sim_wetdry_days, .id = "rlz")
  sim_wetdry_spells <- bind_rows(sim_wetdry_spells, .id = "rlz")
  sim_allcor <- bind_rows(sim_allcor, .id = "rlz") %>% rename(Simulated = value)

  #### MERGE STATISTICS #############################################

  var_combs <- apply(combn(variables, 2),2, paste, collapse=":")
  id_combs  <- apply(combn(1:nsgrids, 2),2, paste, collapse=":")

  # Cross-correlation across sites
  stats_allcor <- sim_allcor %>%
    left_join(hist_allcor, by = c("id1","variable1", "id2", "variable2"))

  stats_crosscor <- stats_allcor %>%
    filter(variable1 == variable2) %>%
    filter(id1 != id2) %>%
    unite(id, c("id1","id2"), sep=":") %>%
    filter(id %in% id_combs)

  stats_intercor <- stats_allcor %>%
    filter(id1 == id2) %>%
    filter(variable1 != variable2) %>%
    unite(id, c("id1","id2"),sep=":") %>%
    unite(variable, c("variable1", "variable2"),sep=":") %>%
    filter(variable %in% var_combs)

  # Wet and Dry days
  stats_wetdry_days <- sim_wetdry_days %>%
    left_join(hist_wetdry_days, by = c("id","mon","variable","stat")) %>%
    mutate(stat = factor(stat, levels = c("dry_count", "wet_count"), labels = c("Dry", "Wet")))

  stats_wetdry_spells <- sim_wetdry_spells %>%
    left_join(hist_wetdry_spells, by = c("id","mon","variable","stat")) %>%
    mutate(stat = factor(stat, levels = c("dry", "wet"), labels = c("Dry", "Wet")))

  stats_wetdry_spells_aavg <- stats_wetdry_spells %>%
    group_by(rlz, mon, variable, stat) %>%
    summarize(Simulated = mean(Simulated), Observed = mean(Observed))

  stats_wetdry_days_aavg <- stats_wetdry_days %>%
    group_by(rlz, mon, variable, stat) %>%
    summarize(Simulated = mean(Simulated), Observed = mean(Observed))



  #:::::::::::::::::::::::::::: PLOTS ::::::::::::::::::::::::::::::::::::::::::

  base_plot_length <- 5
  base_font_size <- 12
  show.title <- TRUE
  alpha_val <- 0.4

  # Common theme for the plots
  theme_wgplots <- theme_bw(base_size = base_font_size) +
    theme(plot.title = element_text(size = base_font_size, face = "bold"),
      plot.subtitle = element_text(size = 9))


  ### Wet and dry spell lengths ++++++++++++++++++++++++++++++++++++++++++++++++
  xy_breaks <- pretty(unlist(stats_wetdry_spells[,c(6,7)]), 5)

  p <- ggplot(stats_wetdry_spells, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, size = 0.8) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 1) +
    facet_wrap(stat ~ ., ncol = 2) +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    xlab(name_obs) + ylab(name_st)

    if(isTRUE(show.title)) {
      p <- p +  labs(title =  "Average length of dry and wet spells per month and grid cell (days)",
        subtitle = "Range and median of all stochastic simulations are shown against the observed values")
    }

  ggsave(file.path(output.path,"average_length_drywet_spells.png"),
        height = base_plot_length, width = base_plot_length*1.75)


  ### Wet dry days +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  xy_breaks <- pretty(unlist(stats_wetdry_days[,c(6,7)]), 5)

  p <- ggplot(stats_wetdry_days, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, size = 0.8) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 1) +
    facet_wrap(stat ~ ., ncol = 2) +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    xlab(name_obs) + ylab(name_st)

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Average number of dry and wet days per month and grid cell",
                   subtitle = "Range and median of all stochastic simulations are shown against the observed values")
  }

  ggsave(file.path(output.path,"average_number_drywet_days.png"),
         height = base_plot_length, width = base_plot_length*1.75)

  # Average dry spell lengths ++++++++++++++++++++++++++++++++++++++++++++++++++
  p <- ggplot(stats_wetdry_spells_aavg, aes(x = as.factor(mon), y = Simulated)) +
    theme_wgplots +
    facet_wrap(stat ~ ., ncol = 2, scales = "free_y") +
    stat_summary(fun.max = max, fun.min = min,
       geom = "linerange", alpha = 0.3, linewidth = 1.5) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    geom_point(aes(y = Observed), color = "blue", size = 2) +
    geom_line(aes(y = Observed, group = 1), color = "blue") +
    labs(x="Month", y = "Days")

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Average length of dry and wet spells per month",
                   subtitle = "Stochastic simulation range and median is shown against the observed values (blue color)")
  }

  ggsave(file.path(output.path,"monthly_average_length_drywet_spells.png"),
         height = base_plot_length*1.5, width = base_plot_length*1.75)


  # Average dry spell number ++++++++++++++++++++++++++++++++++++++++++++++++++
  p <- ggplot(stats_wetdry_days_aavg, aes(x = as.factor(mon), y = Simulated)) +
    theme_wgplots +
    facet_wrap(stat ~ ., ncol = 2, scales = "free_y") +
    stat_summary(fun.max = max, fun.min = min,
                 geom = "linerange", alpha = 0.3, linewidth = 1.5) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    geom_point(aes(y = Observed), color = "blue", size = 2) +
    geom_line(aes(y = Observed, group = 1), color = "blue") +
    labs(x="Month", y = "Days")


  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Average number of dry and wet spells per month",
                   subtitle = "Stochastic simulation range and median is shown against the observed values (blue color)")
  }

  ggsave(file.path(output.path,"monthly_average_number_drywet_days.png"),
         height = base_plot_length*1.5, width = base_plot_length*1.75)


  ### Cross site correlations ++++++++++++++++++++++++++++++++++++++++++++++++++
  xy_breaks <- pretty(unlist(stats_crosscor[,c(5,6)]), 3)

  p <- ggplot(stats_crosscor, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, size = 0.8) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 1) +
    facet_wrap(variable1 ~ ., ncol = 2, scales = "free") +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    xlab(name_obs) + ylab(name_st)


  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Cross-site correlations of all daily variables",
                   subtitle = "Range and median values from all stochastic simulations are shown against the observed values for each grid cell.\nCorrelations are calculated over the entire period.")
  }

  ggsave(file.path(output.path,"daily_crossite_correlations.png"),
         height = base_plot_length*1.5, width = base_plot_length*1.75)


  ### Inter-site correlations ++++++++++++++++++++++++++++++++++++++++++++++++++
  xy_breaks <- pretty(unlist(stats_intercor[,c(4,5)]), 5)

  p <- ggplot(stats_intercor, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, size = 0.8) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 1) +
    facet_wrap(variable ~ ., ncol = 3, scales = "free") +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    xlab(name_obs) + ylab(name_st)

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Inter-site correlations between each pair of variables",
                   subtitle = "Range and median values from all stochastic simulations are shown against the observed values for each grid cell.\nCorrelations are calculated over the entire period.")
  }

  ggsave(file.path(output.path,"daily_intersite_correlations.png"),
         height = base_plot_length*1.5, width = base_plot_length*1.75)


  ### Monthly statistics per variable +++++++++++++++++++++++++++++++++++++++++++
  for (v in 1:length(variables)) {

    ##### MONTHLY CYCLE STATISTICS
    plot_cols <- setNames(c("blue", "black"), c(name_obs, name_st))
    dat <- sim_stats_aavg %>% filter(variable == variables[v] & !is.nan(value))

    p <- ggplot(dat, aes(x = as.factor(mon), y=value)) +
      theme_bw(base_size = base_font_size + 1) +
      theme(plot.title = element_text(size = base_font_size+2, face = "bold"),
            plot.subtitle = element_text(size = base_font_size)) +
      geom_boxplot(color = plot_cols[2], alpha = alpha_val) +
      facet_wrap(~ stat, scales = "free", ncol = 2) +
      stat_summary(fun="mean",  alpha = alpha_val,
                   aes(color=names(plot_cols)[2]),  geom="point") +
      geom_point(data = filter(hist_stats_aavg_mon, variable == variables[v]),
                 aes(color=names(plot_cols)[1]), size = 2.5) +
      scale_color_manual("", values=plot_cols) +
      labs(x = "", y = "") +
      theme(legend.position = c(0.875, 0.40),
            legend.background = element_rect(fill = "white", color = NA),
            legend.text=element_text(size=base_font_size))

    if(isTRUE(show.title)) {
      p <- p +  labs(title =  paste0(variable.labels[v],": monthly variability"),
        subtitle = "Monthly means from all stochastic simulations compared to the observed values.\nVariability range is calculated accross each grid cell.")
    }

    ggsave(file.path(output.path, paste0("monthly_variability_", variables[v],".png")),
           height = base_plot_length*2, width = base_plot_length*1.75)

  }

  ###### DAILY STATISTICS ######################################################

  sim_sts <- sim_stats %>% rename(!!name_st:= value)
  obs_sts <- hist_stats %>% rename(!!name_obs:= value) %>% select(-type)

  daily_stats <- sim_sts %>% left_join(obs_sts, by = c("id","mon","variable","stat")) %>%
    mutate(variable = factor(variable)) %>%
    mutate(stat = factor(stat, levels = stat_label)) %>%
    arrange(variable, stat)

  p <- ggplot(filter(daily_stats, stat == "mean"),
              aes(x = .data[[name_obs]], y = .data[[name_st]])) +
    theme_wgplots +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, size = 0.8) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 1) +
    geom_abline(color = "blue") +
    labs(x = name_obs, y = name_st) +
    facet_wrap(variable ~ ., scales = "free", nrow = 2)

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Mean of daily variables calculated across each month and grid cell",
                   subtitle = "Range from all stochastic simulations are shown against the observed values")
  }

  ggsave(file.path(output.path,"daily_mean_all_variables.png" ),
         height = base_plot_length*2, width = base_plot_length*1.75)


  p <- ggplot(filter(daily_stats, stat == "standard dev."),
              aes(x = .data[[name_obs]], y = .data[[name_st]])) +
    theme_wgplots +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, size = 0.8) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 1) +
    geom_abline(color = "blue") +
    labs(x = name_obs, y = name_st) +
    facet_wrap(variable ~ ., scales = "free", nrow = 2)

   if(isTRUE(show.title)) {
     p <- p +  labs(title =  "Standard deviation of daily variables calculated across each month and grid cell",
                    subtitle = "Range from all stochastic simulations are shown against the observed values")
   }

  ggsave(file.path(output.path,"daily_stdev_all_variables.png" ),
      height = base_plot_length*2, width = base_plot_length*1.75)

}







