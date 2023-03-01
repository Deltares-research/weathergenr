#' Compare observed and simulated climate statistics
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
  # General parameters #########################################################

  #Workaround for rlang warning
  year <- mon <- day <- precip <- sd <- variable <- value <- id_variable <- 0
  Wet <- Dry <- wet <- dry <- id_variable1 <- id_variable2 <- 0
  rlz <- id1 <- id2 <- variable1 <- variable2 <- type <- Observed <- Stochastic <- 0
  wet_th <- 0

  # General options
  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)

  # General parameters
  stat_level <- c("mean", "sd", "skewness")
  stat_label <- c("mean", "standard dev.", "skewness")

  num_stats <- length(stat_level)
  num_vars <- length(variables)
  var_combs <- expand_grid(var1=variables, var2 = variables)
  num_var_combs <- choose(num_vars, 2)


  nsgrids <- length(daily.sim[[1]])
  message(cat(as.character(Sys.time()), "- Comparison accross",nsgrids,"grid cells"))

  # Calculate observed climate statistics ######################################################

  message(cat(as.character(Sys.time()), "- Calculating historical climate statistics"))

  # Number of years in the historical record
  hist_year_num <- nrow(daily.obs[[1]])/365

  # Clean-up historical data
  hist_daily_tidy <- lapply(daily.obs, "[", c("date", variables)) %>%
    bind_rows(.id = "id") %>%
    mutate(year = as.numeric(format(date,"%Y")),
       mon = as.numeric(format(date,"%m")),
       day = as.numeric(format(date,"%d")),
       id = as.numeric(id)) %>%
    select(id, year, mon, day, all_of(variables))

  # Calculate summary statistics
  hist_stats <- hist_daily_tidy %>%
    group_by(id, mon) %>%
    dplyr::summarize(across(all_of(variables),
                            list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = Observed, -id, -year, -mon) %>%
    separate(variable, c("variable","stat"), sep = ":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  # Calculate dry, wet, and extremely wet day thresholds
  mc_thresholds <- hist_daily_tidy %>% group_by(mon) %>%
    group_by(year, mon, day) %>%
    summarize(precip = mean(precip)) %>%
    group_by(mon) %>%
    summarize(wet_th = stats::quantile(precip, wet.quantile, names = F),
              extreme_th = stats::quantile(precip, extreme.quantile, names = F))

  # Calculate wet and dry days per month
  hist_wetdry_days <- hist_daily_tidy %>%
    left_join(mc_thresholds, by = "mon") %>%
    group_by(id, mon) %>%
    summarize(Wet = length(which(precip >= wet_th))/hist_year_num,
              Dry = length(which(precip < wet_th))/hist_year_num) %>%
    gather(key = stat, value = Observed, Wet:Dry) %>%
    mutate(variable = "precip", .after = mon)

  # Wet and dry spell lengths per month
  hist_wetdry_spells <- hist_daily_tidy %>%
    left_join(mc_thresholds, by = "mon") %>%
    group_by(id, mon) %>%
    summarize(dry = averageSpellLength(precip, threshold = wet_th, below = TRUE),
              wet = averageSpellLength(precip, threshold = wet_th, below = FALSE)) %>%
    group_by(id, mon) %>%
    summarize(Dry = mean(dry),
              Wet = mean(wet)) %>%
    gather(key = stat, value = Observed, Wet:Dry) %>%
    mutate(variable = "precip", .after = mon)

  # Intersite/crossite correlations
  hist_allcor_ini <- hist_daily_tidy %>%
    gather(key = variable, value = value, -id, -year, -mon, -day) %>%
    unite(id_variable, c("id","variable"),sep = ":") %>%
    spread(id_variable, value) %>%
    dplyr::select(-year, -mon, -day) %>%
    stats::cor(.$value) %>% as_tibble()

  hist_allcor <- hist_allcor_ini %>%
    mutate(id_variable1 = colnames(hist_allcor_ini), .before = 1) %>%
    gather(key = id_variable2, value = Observed, -id_variable1) %>%
    separate(id_variable1, c("id1","variable1"), sep =":") %>%
    separate(id_variable2, c("id2","variable2"), sep =":")

  # Historical stats averaged over grid/cells
  hist_stats_aavg <- hist_daily_tidy %>% group_by(year, mon) %>%
    summarize(across(all_of(variables),
      list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = Observed, -mon, -year) %>%
    separate(variable, c("variable","stat"), sep=":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  # Historical stats averaged over grid/cells & months
  hist_stats_aavg_mon <- hist_stats_aavg %>% group_by(mon, variable, stat) %>%
    summarize(Observed = mean(Observed))


  # Calculate simulated climate statistics ######################################################

  message(cat(as.character(Sys.time()), "- Calculating simulated climate statistics"))

  # Initialize lists to store the results
  sim_stats <- vector("list", realization.num)
  sim_stats_aavg <- vector("list", realization.num)
  sim_allcor <- vector("list", realization.num)
  sim_wetdry_days <-vector("list", realization.num)
  sim_wetdry_spells <- vector("list", realization.num)

  # Total years of simulation
  sim_year_num <- nrow(daily.sim[[1]][[1]])/365

  #### Calculate statistics per realization
  for (n in 1:realization.num) {

    # Simulated series in tidy format
    sim_daily_tidy <- lapply(daily.sim[[n]], "[", c("date", variables)) %>%
      bind_rows(.id = "id") %>%
      mutate(year = as.numeric(format(date,"%Y")),
         mon = as.numeric(format(date,"%m")),
         day = as.numeric(format(date,"%d")),
         id = as.numeric(id)) %>%
      select(id, year, mon, day, all_of(variables))

    # Calculate summary statistics per grid cell
    sim_stats[[n]] <- sim_daily_tidy %>%
      group_by(id, mon) %>%
      dplyr::summarize(across(all_of(variables),
        list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
      gather(key = variable, value = Simulated,-id, -year, -mon) %>%
      separate(variable, c("variable","stat"), sep = ":") %>%
      mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

    # Area-averaged stats
    sim_stats_aavg[[n]] <- sim_daily_tidy %>%
      group_by(year, mon) %>%
      summarize(across(all_of(variables),
        list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
      gather(key = variable, value = Simulated, -mon, -year) %>%
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
      gather(key = id_variable2, value = Simulated, -id_variable1) %>%
      separate(id_variable1, c("id1","variable1"), sep =":") %>%
      separate(id_variable2, c("id2","variable2"), sep =":") #%>%

    sim_wetdry_days[[n]] <- sim_daily_tidy %>%
        left_join(mc_thresholds, by = "mon") %>%
        group_by(id, mon) %>%
        summarize(Wet = length(which(precip >= wet_th))/sim_year_num,
                  Dry = length(which(precip < wet_th))/sim_year_num) %>%
      gather(key = stat, value = value, Wet:Dry) %>%
      mutate(variable = "precip") %>%
      select(id, mon, variable, stat, Simulated = value)

    sim_wetdry_spells[[n]] <- sim_daily_tidy %>%
      left_join(mc_thresholds, by = "mon") %>%
      group_by(id, mon) %>%
      summarize(Dry = averageSpellLength(precip, threshold = wet_th, below = TRUE),
                Wet = averageSpellLength(precip, threshold = wet_th, below = FALSE)) %>%
      group_by(id, mon) %>%
      summarize(Dry = mean(Dry),Wet = mean(Wet)) %>%
      gather(key = stat, value = Simulated, Wet:Dry) %>%
      mutate(variable = "precip", .after = mon)

  }

  sim_stats <- bind_rows(sim_stats, .id = "rlz")
  sim_stats_aavg <- bind_rows(sim_stats_aavg, .id = "rlz")
  sim_wetdry_days <- bind_rows(sim_wetdry_days, .id = "rlz")
  sim_wetdry_spells <- bind_rows(sim_wetdry_spells, .id = "rlz")
  sim_allcor <- bind_rows(sim_allcor, .id = "rlz")


  # Merge results #################################

  var_combs <- apply(combn(variables, 2),2, paste, collapse=":")
  id_combs  <- apply(combn(1:nsgrids, 2),2, paste, collapse=":")

  # Combined summary statistics
  daily_stats <- sim_stats %>%
    left_join(hist_stats, by = c("id","mon","variable","stat")) %>%
    mutate(variable = factor(variable)) %>%
    mutate(stat = factor(stat, levels = stat_label))

  # Cross-correlation across sites
  stats_allcor <- sim_allcor %>%
    left_join(hist_allcor, by = c("id1","variable1", "id2", "variable2"))

  # Cross-site correlations
  stats_intersite_cor <- stats_allcor %>%
    filter(variable1 == variable2) %>%
    filter(id1 != id2) %>%
    unite(id, c("id1","id2"), sep=":") %>%
    filter(id %in% id_combs)

  # Intersite correlations
  stats_cross_cor <- stats_allcor %>%
    filter(id1 == id2) %>%
    filter(variable1 != variable2) %>%
    unite(id, c("id1","id2"),sep=":") %>%
    unite(variable, c("variable1", "variable2"),sep=":") %>%
    filter(variable %in% var_combs)

  # Wet and Dry days across all grids
  stats_wetdry_days <- sim_wetdry_days %>%
    left_join(hist_wetdry_days, by = c("id","mon","variable","stat")) %>%
    mutate(stat = factor(stat, levels = c("Dry", "Wet")))
  stats_wetdry_spells <- sim_wetdry_spells %>%
    left_join(hist_wetdry_spells, by = c("id","mon","variable","stat")) %>%
    mutate(stat = factor(stat, levels = c("Dry", "Wet")))

  # Wet dry spells and days (area-averaged)
  stats_wetdry_spells_aavg <- stats_wetdry_spells %>%
    group_by(rlz, mon, variable, stat) %>%
    summarize(Simulated = mean(Simulated), Observed = mean(Observed))
  stats_wetdry_days_aavg <- stats_wetdry_days %>%
    group_by(rlz, mon, variable, stat) %>%
    summarize(Simulated = mean(Simulated), Observed = mean(Observed))

  # Plot results ###################################################################

  message(cat(as.character(Sys.time()), "- Preparing comparison plots"))

  # Various plotting parameters
  plot_length <- 5
  font_size <- 12
  title_size <- font_size + 2
  alpha_val <- 0.4

  # Common theme for the plots
  theme_wgplots <- theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 14),
          plot.subtitle = element_text(size = 10))


  # 1) Summary statistics for all climate variables

  p <- ggplot(filter(daily_stats, stat == "mean"),
              aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 0.8) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 1) +
    geom_abline(color = "blue") +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", nrow = 2)

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Mean of daily variables calculated across each month and grid cell",
                   subtitle = "Range from all stochastic simulations are shown against the observed values")
  }

  ggsave(file.path(output.path,"daily_mean_all_variables.png" ),
         height = plot_length*2, width = plot_length*1.75)


  p <- ggplot(filter(daily_stats, stat == "standard dev."),
              aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 0.8) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 1) +
    geom_abline(color = "blue") +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", nrow = 2)

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Standard deviation of daily variables calculated across each month and grid cell",
                   subtitle = "Range from all stochastic simulations are shown against the observed values")
  }

  ggsave(file.path(output.path,"daily_stdev_all_variables.png" ),
         height = plot_length*2, width = plot_length*1.75)



# 2) Wet and dry spell statistics

  xy_breaks <- pretty(unlist(stats_wetdry_spells[,c(6,7)]), 5)

  p <- ggplot(stats_wetdry_spells, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 0.8) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 1) +
    facet_wrap(stat ~ ., ncol = 2) +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    xlab("Observed") + ylab("Simulated")

    if(isTRUE(show.title)) {
      p <- p +  labs(title =  "Average length of dry and wet spells per month and grid cell (days)",
        subtitle = "Range and median of all stochastic simulations are shown against the observed values")
    }

  ggsave(file.path(output.path,"avg_length_drywet_spells.png"),
        height = plot_length, width = plot_length*1.75)


  #3) Average number of wet and dry days

  xy_breaks <- pretty(unlist(stats_wetdry_days[,c(6,7)]), 5)

  p <- ggplot(stats_wetdry_days, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 0.8) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 1) +
    facet_wrap(stat ~ ., ncol = 2) +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    xlab("Observed") + ylab("Simulated")

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Average number of dry and wet days per month across all grid cells",
                   subtitle = "Stochastic simulation range is shown against the observed values")
  }

  ggsave(file.path(output.path,"avg_number_drywet_days.png"),
         height = plot_length, width = plot_length*1.75)


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
                   subtitle = "Stochastic simulation range is shown against the observed values (blue color).")
  }

  ggsave(file.path(output.path,"monthly_drywet_spells.png"),
         height = plot_length*1.5, width = plot_length*1.75)


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
                   subtitle = "Stochastic simulation range is shown against the observed values (blue color).")
  }

  ggsave(file.path(output.path,"monthly_drywet_days.png"),
         height = plot_length*1.5, width = plot_length*1.75)

  #4) Inter-site correlations
  xy_breaks <- pretty(unlist(stats_intersite_cor[,c(5,6)]), 3)

  p <- ggplot(stats_intersite_cor, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 1.5) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    facet_wrap(variable1 ~ ., ncol = 2, scales = "free") +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    xlab("Observed") + ylab("Simulated")

  if(isTRUE(show.title)) {
    p <- p +  labs(title = "Intersite correlations for all variables",
                   subtitle = "Stochastic simulation range is shown against the observed values for each grid cell.\nCorrelations are calculated over the entire period.")
  }

  ggsave(file.path(output.path,"intersite_correlations.png"),
         height = plot_length*1.75, width = plot_length*1.75)


  #5) Cross correlations
  xy_breaks <- pretty(unlist(stats_cross_cor[,c(4,5)]), 5)

  p <- ggplot(stats_cross_cor, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 1.5) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    facet_wrap(variable ~ ., ncol = 3, scales = "free") +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    xlab("Observed") + ylab("Simulated")

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Cross correlations between each pair of variables",
                   subtitle = "Stochastic simulation range is shown against the observed values for each grid cell.\nCorrelations are calculated over the entire period.")
  }

  ggsave(file.path(output.path,"cross_correlations.png"),
         height = plot_length*1.75, width = plot_length*1.75)


  #6) Monthly statistics per variable
  for (v in 1:length(variables)) {

    ##### MONTHLY CYCLE STATISTICS
    plot_cols <- setNames(c("blue", "black"), c("Observed", "Simulated"))
    dat <- sim_stats_aavg %>% filter(variable == variables[v] & !is.nan(value))

    p <- ggplot(dat, aes(x = as.factor(mon), y=Simulated)) +
      theme_wgplots +
      geom_boxplot(color = plot_cols[2], alpha = alpha_val) +
      facet_wrap(~ stat, scales = "free", ncol = 2) +
      stat_summary(fun="mean",  alpha = alpha_val,
                   aes(color=names(plot_cols)[2]),  geom="point") +
      geom_point(data = filter(hist_stats_aavg_mon, variable == variables[v]),
                 aes(color=names(plot_cols)[1], y = Observed), size = 2.5) +
      scale_color_manual("", values=plot_cols) +
      labs(x = "", y = "") +
      theme(legend.position = c(0.875, 0.40),
            legend.background = element_rect(fill = "white", color = NA),
            legend.text=element_text(size=12),
            plot.title = element_text(size = 14),
            plot.subtitle = element_text(size = 12))

    if(isTRUE(show.title)) {
      p <- p +  labs(title =  paste0("Monthly variability of ", variable.labels[v]),
        subtitle = "Monthly means from all stochastic simulations compared to the observed values.\nVariability range is averaged accross all grid cells.")
    }

    ggsave(file.path(output.path, paste0("monthly_variability_", variables[v],".png")),
           height = plot_length*2, width = plot_length*1.75)

  }

}





