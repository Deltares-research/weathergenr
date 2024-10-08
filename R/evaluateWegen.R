#' @title Compare observed and simulated climate statistics
#'
#' @description
#' A short description...
#'
#' @details
#' Additional details...
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
#' @param save.plots PLACEHOLDER
#'
#' @return returns a list of ggplot2 objects and figures saved to the desired output path
#' @export
#' @import ggplot2 dplyr
#' @importFrom stats median sd setNames
#' @importFrom utils combn
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
  show.title = TRUE,
  save.plots = TRUE)

{
  # General parameters #########################################################

  #Workaround for rlang warning
  year <- mon <- day <- precip <- sd <- variable <- value <- id_variable <- 0
  Wet <- Dry <- wet <- dry <- id_variable1 <- id_variable2 <- 0
  rlz <- id1 <- id2 <- variable1 <- variable2 <- type <- Observed <- Stochastic <- 0
  wet_th <- Simulated <- 0

  # Set variable names and labels
  if(is.null(variable.labels)) variable.labels <- variables
  if(is.null(variable.units)) variable.units <- rep("", length(variables))

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

  message(cat(as.character(format(Sys.time(),'%H:%M:%S')), "- Comparison accross", nsgrids, "grid cells"))

  # Calculate observed climate statistics ######################################################

  message(cat(as.character(format(Sys.time(),'%H:%M:%S')), "- Calculating historical trace statistics"))

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

  # Calculate summary statistics (per month)
  hist_stats_season <- hist_daily_tidy %>%
    group_by(id, mon) %>%
    dplyr::summarize(across(all_of(variables),
                            list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = Observed, -id, -year, -mon) %>%
    separate(variable, c("variable","stat"), sep = ":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  # Calculate summary statistics (per month x per year)
  hist_stats_mon_aavg <- hist_daily_tidy %>%
    group_by(year, mon) %>%
    dplyr::summarize(across(all_of(variables),
                            list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = Observed, -year, -mon) %>%
    separate(variable, c("variable","stat"), sep = ":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  hist_stats_annual_aavg <- hist_daily_tidy %>%
    group_by(year, mon, day) %>%
    summarize_at(vars(all_of(variables)), mean) %>%
    group_by(year) %>%
    summarize(across(all_of(variables),
                     list(mean=mean, min=min, max=max),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = Observed, -year) %>%
    separate(variable, c("variable","stat"), sep=":") %>%
    mutate(year = year - min(year) + 1)

  # Historical stats averaged over grid/cells & months
  hist_stats_season_aavg <- hist_stats_mon_aavg %>% group_by(mon, variable, stat) %>%
    summarize(Observed = mean(Observed))

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

  # Calculate simulated climate statistics #####################################

  message(cat(as.character(format(Sys.time(),'%H:%M:%S')), "- Calculating synthetic trace statistics"))

  # Initialize lists to store the results
  sim_stats_season <- vector("list", realization.num)
  sim_stats_mon_aavg <- vector("list", realization.num)
  sim_stats_annual_aavg <- vector("list", realization.num)
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
    sim_stats_season[[n]] <- sim_daily_tidy %>%
      group_by(id, mon) %>%
      dplyr::summarize(across(all_of(variables),
        list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
      gather(key = variable, value = Simulated,-id, -year, -mon) %>%
      separate(variable, c("variable","stat"), sep = ":") %>%
      mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

    # Area-averaged stats
    sim_stats_mon_aavg[[n]] <- sim_daily_tidy %>%
      group_by(year, mon) %>%
      summarize(across(all_of(variables),
        list(mean=mean, sd=stats::sd, skewness=e1071::skewness),.names = "{.col}:{.fn}")) %>%
      gather(key = variable, value = Simulated, -mon, -year) %>%
      separate(variable, c("variable","stat"), sep=":") %>%
      mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

    sim_stats_annual_aavg[[n]] <- sim_daily_tidy %>%
      group_by(year, mon, day) %>%
      summarize_at(vars(all_of(variables)), mean) %>%
      group_by(year) %>%
      summarize(across(all_of(variables),
                       list(mean=mean, min=min, max=max),.names = "{.col}:{.fn}")) %>%
      gather(key = variable, value = Simulated, -year) %>%
      separate(variable, c("variable","stat"), sep=":")

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

  sim_stats_season <- bind_rows(sim_stats_season, .id = "rlz")
  sim_stats_mon_aavg <- bind_rows(sim_stats_mon_aavg, .id = "rlz")
  sim_stats_annual_aavg <- bind_rows(sim_stats_annual_aavg, .id = "rlz") %>%
    mutate(year = year - min(year) + 1)
  sim_allcor <- bind_rows(sim_allcor, .id = "rlz")
  sim_wetdry_days <- bind_rows(sim_wetdry_days, .id = "rlz")
  sim_wetdry_spells <- bind_rows(sim_wetdry_spells, .id = "rlz")

  # Merge results ##############################################################

  var_combs <- apply(combn(variables, 2),2, paste, collapse=":")
  id_combs  <- apply(combn(1:nsgrids, 2),2, paste, collapse=":")

  # Combined summary statistics
  daily_stats_season <- sim_stats_season %>%
    left_join(hist_stats_season, by = c("id","mon","variable","stat")) %>%
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


  ##############################################################################

  # Plot results

  message(cat(as.character(format(Sys.time(),'%H:%M:%S')), "- Preparing comparison plots"))

  # Common plotting parameters
  pl_size <- 8  # 4 for each row/column
  font_size <- 12
  title_size <- font_size + 2
  alpha_val <- 0.4
  pl_sub <- "Value range and median values from all simulations are shown against the observed"

  theme_wgplots <- theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 14),
          plot.subtitle = element_text(size = 10))

  plot_cols <- setNames(c("blue3", "gray40"), c("Observed", "Simulated"))
  plots <- list()

  # 1) Daily mean statistics for all variables
  dummy_gg <- daily_stats_season %>%
    filter(stat == "mean") %>%
    group_by(variable) %>%
    summarize(minval = min(Simulated, Observed)*1,
              maxval = max(Simulated, Observed)*1) %>%
    pivot_longer(cols = minval:maxval, names_to = "type", values_to = "value") %>%
    select(variable, Observed = value, Simulated = value)

  p <- ggplot(filter(daily_stats_season, stat == "mean"),
              aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_point(data = dummy_gg, color = NA) +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 1.5) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    geom_abline(color = "blue") +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", nrow = 2)

  if(isTRUE(show.title))
    p <- p + labs(title =  "Daily means for all grid cell and months", subtitle = pl_sub)

  if(isTRUE(save.plots))
    ggsave(file.path(output.path,"daily_means.png"), height=pl_size, width=pl_size)

  plots$daily_means <- p

  # 2) Daily standard deviations for all variables
  dummy_gg <- daily_stats_season %>%
    filter(stat == "standard dev.") %>%
    group_by(variable) %>%
    summarize(minval = min(Simulated, Observed)*1,
              maxval = max(Simulated, Observed)*1) %>%
    pivot_longer(cols = minval:maxval, names_to = "type", values_to = "value") %>%
    select(variable, Observed = value, Simulated = value)


  p <- ggplot(filter(daily_stats_season, stat == "standard dev."),
              aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_point(data = dummy_gg, color = NA) +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 1.5) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    geom_abline(color = "blue") +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", nrow = 2)

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Daily standard deviations for all grid cell and months", subtitle=pl_sub)
  }

  if(isTRUE(save.plots))
    ggsave(file.path(output.path,"daily_stdev.png"), height=pl_size, width=pl_size)

  plots$daily_sd <- p

# 3) Wet and dry spell statistics
  dummy_gg <- stats_wetdry_spells %>%
    group_by(stat) %>%
    summarize(minval = min(Simulated, Observed)*1,
              maxval = max(Simulated, Observed)*1) %>%
    pivot_longer(cols = minval:maxval, names_to = "type", values_to = "value") %>%
    select(stat, Observed = value, Simulated = value)

  p <- ggplot(stats_wetdry_spells, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    geom_point(data = dummy_gg, color = NA) +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 1.5) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    facet_wrap(stat ~ ., ncol = 2, scales = "free") +
    xlab("Observed") + ylab("Simulated")

    if(isTRUE(show.title)) {
      p <- p +  labs(title =  "Average dry and wet spell length per month, across all grid cells", subtitle = pl_sub)
    }

  if(isTRUE(save.plots))
    ggsave(file.path(output.path,"drywet_spell_length.png"), height = pl_size/2, width = pl_size)

  plots$spell_lengths <- p

  #4) Average number of wet and dry days
  dummy_gg <- stats_wetdry_days %>%
    group_by(stat) %>%
    summarize(minval = min(Simulated, Observed)*1,
              maxval = max(Simulated, Observed)*1) %>%
    pivot_longer(cols = minval:maxval, names_to = "type", values_to = "value") %>%
    select(stat, Observed = value, Simulated = value)

  p <- ggplot(stats_wetdry_days, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    geom_point(data = dummy_gg, color = NA) +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 1.5) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    facet_wrap(stat ~ ., ncol = 2, scales = "free") +
    xlab("Observed") + ylab("Simulated")

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  "Average number of dry and wet days per month accross all grid cells", subtitle = pl_sub)
  }

  if(isTRUE(save.plots))
    ggsave(file.path(output.path,"drywet_days_number.png"), height = pl_size/2, width = pl_size)

  plots$spell_duration <- p

  #5) CROSS-GRID CORRELATIONS
  dummy_gg <- stats_intersite_cor %>%
    group_by(variable1) %>%
    summarize(minval = min(Simulated, Observed)*1,
              maxval = max(Simulated, Observed)*1) %>%
    pivot_longer(cols = minval:maxval, names_to = "type", values_to = "value") %>%
    select(variable1, Observed = value, Simulated = value)

  p <- ggplot(stats_intersite_cor, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    geom_point(data = dummy_gg, color = NA) +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 1.5) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    facet_wrap(variable1 ~ ., ncol = 2, scales = "free") +
    xlab("Observed") + ylab("Simulated")

  if(isTRUE(show.title)) {
    p <- p +  labs(title = "Cross-grid correlations",
                   subtitle = paste0(pl_sub,"\nCorrelations are calculated over daily series"))
  }

  if(isTRUE(save.plots))
    ggsave(file.path(output.path,"crossgrid_correlations.png"), height = pl_size, width = pl_size)

  plots$crossgrid_cor <- p

  #6) INTERGRID CORRELATIONS
  dummy_gg <- stats_cross_cor %>%
    group_by(variable) %>%
    summarize(minval = min(Simulated, Observed)*1,
              maxval = max(Simulated, Observed)*1) %>%
    pivot_longer(cols = minval:maxval, names_to = "type", values_to = "value") %>%
    select(variable, Observed = value, Simulated = value)

  p <- ggplot(stats_cross_cor, aes(x = Observed, y = Simulated)) +
    theme_wgplots +
    geom_abline(color = "blue") +
    geom_point(data = dummy_gg, color = NA) +
    stat_summary(geom = "linerange", fun.max = max, fun.min = min,
                 alpha = alpha_val, linewidth = 1.5) +
    stat_summary(fun = "median", geom = "point", alpha = alpha_val, size = 2) +
    facet_wrap(variable ~ ., ncol = 3, scales = "free") +
    xlab("Observed") + ylab("Simulated")

  if(isTRUE(show.title)) {
    p <- p +  labs(title = "Inter-grid correlations",
                   subtitle = paste0(pl_sub,"\nCorrelations are calculated over daily series"))
  }

  if(isTRUE(save.plots))
  ggsave(file.path(output.path,"intergrid_correlations.png"), height = pl_size, width = pl_size*1.25)

  plots$intergrid_cor <- p

  #7) Monthly statistics per variable
  for (v in 1:length(variables)) {

    dat <- sim_stats_mon_aavg %>%
      filter(variable == variables[v]) %>%
      mutate(type = "Simulated") %>%
      rename(value = Simulated)

    dat2 <- hist_stats_mon_aavg %>%
      mutate(rlz = "0", .before = year) %>%
      filter(variable == variables[v] & !is.nan(value)) %>%
      mutate(type = "Observed") %>%
      rename(value = Observed)

    datx <- bind_rows(dat, dat2)

    p <- ggplot(datx, aes(x = as.factor(mon), y = value, fill = type, color = type)) +
      theme_wgplots +
      geom_boxplot(alpha = 0.2) +
      facet_wrap(~ stat, scales = "free", ncol = 2) +
      scale_fill_manual("", values=plot_cols) +
      scale_color_manual("", values=plot_cols) +
      stat_summary(fun="mean",  size = 3, geom="point", position = position_dodge(0.8), shape = 18) +
      labs(x = "", y = "") +
      theme(legend.position = c(1, 0),
            legend.justification = c(1, 0),
            legend.background = element_rect(fill = "white", color = NA),
            legend.text=element_text(size=12),
            plot.title = element_text(size = 14),
            plot.subtitle = element_text(size = 12)) +
      scale_x_discrete(labels = substr(month.name, 1,1))

    if(isTRUE(show.title)) {
      p <- p +  labs(title =  paste0("Monthly patterns for ", variable.labels[v]),
        subtitle = paste0(pl_sub,"\nResults are averaged accross all grid cells."))
    }

    if(isTRUE(save.plots))
      ggsave(file.path(output.path, paste0("monthly_patterns_", variables[v],".png")),
            height = pl_size, width = pl_size)

    plots[[paste0("annual_pattern_", variables[v])]] <- p
  }

  # 8) Monthly mean cycle
  sim_stats_season_aavg <- sim_stats_season %>%
    group_by(rlz, mon, variable) %>%
    filter(stat == "mean") %>%
    summarize(value = mean(Simulated)) %>%
    mutate(type = "Simulated")

  hist_stats_season_aavg2 <- hist_stats_season %>%
    group_by(mon, variable) %>%
    filter(stat == "mean") %>%
    summarize(value = mean(Observed)) %>%
  mutate(type = "Observed")

  p <- ggplot(sim_stats_season_aavg, aes(x = as.factor(mon), y = value)) +
    theme_wgplots +
    facet_wrap(~ variable, scales = "free") +
    geom_line(aes(group = rlz, color = rlz), alpha = 0.8) +
    geom_line(data = hist_stats_season_aavg2, color = "black", group = 1, size = 1.25) +
    scale_x_discrete(labels = substr(month.name, 1,1)) +
    labs(x = "", y = "")

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  paste0("Annual cycles of variables"),
                   subtitle = paste0(pl_sub,"\nResults are averaged accross each month"))
  }

  if(isTRUE(save.plots))
    ggsave(file.path(output.path, paste0("monthly_cycle.png")), height = pl_size, width = pl_size+1)

  plots$annual_cycle <- p

  # Annual precip means as time-series
  sim_annual_aavg_precip <- sim_stats_annual_aavg %>% filter(stat == "mean") %>% filter(variable == "precip")
  hist_annual_aavg_precip <- hist_stats_annual_aavg %>% filter(stat == "mean") %>% filter(variable == "precip")

  p <- ggplot(sim_annual_aavg_precip, aes(x = year)) +
    theme_wgplots +
    geom_line(aes(y = Simulated, group = rlz), color = "gray30", alpha = 0.4) +
    geom_point(aes(y = Simulated, group = rlz), size = 0.5, color = "gray30", alpha = 0.3) +
    geom_line(aes(y = Observed),
              data = hist_annual_aavg_precip, color = "blue", group = 1) +
    geom_point(aes(y = Observed), size = 0.5, alpha = 0.3,
              data = hist_annual_aavg_precip, color = "blue", group = 1) +
    labs(x = "serial year", y = "mm/day")

  if(isTRUE(show.title)) {
    p <- p +  labs(title =  paste0("Annual mean precipitation"))
  }

  if(isTRUE(save.plots))
    ggsave(file.path(output.path, paste0("annual_precip.png")), height = pl_size/1.9, width = pl_size)

  plots$annual_mean <- p

  return(plots)
}
