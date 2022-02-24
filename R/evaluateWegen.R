

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
#'
#' @return
#' @export
#' @import utils ggplot2 dplyr
evaluateWegen <- function(
  daily.sim = NULL,
  daily.obs = NULL,
  output.path = NULL,
  variables = NULL,
  variable.labels = NULL,
  variable.units = NULL,
  realization.num = NULL,
  wet.quantile = 0.2,
  extreme.quantile = 0.8
)

  {

  #Workaround for rlang warning
  year <- mon <- day <- precip <- sd <- variable <- value <- id_variable <- 0
  wet_count <- dry_count <- wet <- dry <- id_variable1 <- id_variable2 <- 0
  rlz <- id1 <- id2 <- variable1 <- variable2 <- type <- Observed <- Simulated <- 0
  hist_stats_mon <- wet_th <- 0

  options(dplyr.summarise.inform = FALSE)
  options(tidyverse.quiet = TRUE)

  stat_level <- c("mean", "sd", "skewness")
  stat_label <- c("Mean", "Standard Deviation", "Skewness")

  nsgrids <- length(daily.sim[[1]])

  num_stats <- length(stat_level)
  num_vars <- length(variables)
  var_combs <- expand_grid(var1=variables, var2 = variables)
  num_var_combs <- choose(num_vars, 2)





  # Scatterplot of simulated and historical daily means all locations (precip + temp)
  # Scatterplot of daily cross correlations all locations
  # Scatterplot of intersite correlations all locations

  # Monthly patterns (precip + temp)
  # Monthly patterns mean, sd, skew (precip + temp)
  # Wet and dry spell stats
  # Wet and dry day stats

  # Emprical density plot?
  # What else???



  mc_thresholds <- hist_stats_ini %>% group_by(mon) %>%
    group_by(year, mon, day) %>%
    summarize(precip = mean(precip)) %>%
    group_by(mon) %>%
    summarize(wet_th = stats::quantile(precip, wet.quantile, names = F),
              extreme_th = stats::quantile(precip, extreme.quantile, names = F))



  # Historical statistics ######################################################

  hist_year_num <- nrow(daily.obs[[1]])/365

  hist_stats_ini <- lapply(daily.obs, "[", c("date", variables)) %>%
    bind_rows(.id = "id") %>%
    mutate(year = as.numeric(format(date,"%Y")),
       mon = as.numeric(format(date,"%m")),
       day = as.numeric(format(date,"%d")),
       id = as.numeric(id)) %>%
    select(id, year, mon, day, {{variables}})

  # Monthly stats per variable
  hist_stats <- hist_stats_ini %>%
    group_by(id, mon) %>%
    dplyr::summarize(across({{variables}},
      list(mean="mean", sd="sd", skewness=e1071::"skewness"),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = value,-id, -year, -mon) %>%
    separate(variable, c("variable","stat"), sep = ":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label)) %>%
    mutate(type = "Observed")

  # Area-averaged stats per variable
  hist_stats_aavg <- hist_stats_ini %>%
    group_by(year, mon) %>%
    summarize(across({{variables}},
      list(mean="mean", sd="sd", skewness=e1071::"skewness"),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = value, -mon, -year) %>%
    separate(variable, c("variable","stat"), sep=":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))  %>%
    mutate(type = "Observed")

  # Intersite/Intervariable correlations
  hist_stats_icor <- hist_stats_ini %>%
    gather(key = variable, value = value, -id, -year, -mon, -day) %>%
    unite(id_variable, c("id","variable"),sep = ":") %>%
    spread(id_variable, value) %>%
    dplyr::select(-year, -mon, -day) %>%
    stats::cor(.$value) %>% as_tibble()

  hist_intercor <- hist_stats_icor %>%
    mutate(id_variable1 = colnames(hist_stats_icor), .before = 1) %>%
    gather(key = id_variable2, value = value, -id_variable1) %>%
    separate(id_variable1, c("id1","variable1"), sep =":") %>%
    separate(id_variable2, c("id2","variable2"), sep =":") %>%
    mutate(type = "Observed")

  hist_wetdry_days <- hist_stats_ini %>%
    left_join(mc_thresholds, by = "mon") %>%
      group_by(id, mon) %>%
      summarize(wet_count = length(which(precip >= wet_th))/sim_year_num,
                dry_count = length(which(precip < wet_th))/sim_year_num) %>%
      gather(key = stat, value = value, wet_count:dry_count) %>%
      mutate(type = "Observed")

  hist_wetdry_spells <- hist_stats_ini %>%
    left_join(mc_thresholds, by = "mon") %>%
    group_by(id, mon) %>%
    summarize(dry = averageSpellLength(precip, threshold = wet_th, below = TRUE),
              wet = averageSpellLength(precip, threshold = wet_th, below = FALSE)) %>%
    group_by(id, mon) %>%
    summarize(dry = mean(dry),
              wet = mean(wet)) %>%
    gather(key = stat, value = value, wet:dry) %>%
    mutate(type = "Observed")

  # Per variable plots
  hist_stats_aavg_mon <- hist_stats_aavg %>%
    group_by(mon, variable, stat, type) %>%
    summarize(value = mean(value))



  # Simulated statistics ######################################################

  sim_stats <- list()
  sim_stats_aavg <- list()
  sim_intercor <- list()
  sim_wetdry_days <-list()
  sim_wetdry_spells <- list()

  sim_year_num <- nrow(daily.sim[[1]][[1]])/365


  #### Calculate statistics per realization

  for (n in 1:realization.num) {

    sim_stats_ini <- lapply(daily.sim[[n]], "[", c("date", variables)) %>%
      bind_rows(.id = "id") %>%
      mutate(year = as.numeric(format(date,"%Y")),
         mon = as.numeric(format(date,"%m")),
         day = as.numeric(format(date,"%d")),
         id = as.numeric(id)) %>%
      select(id, year, mon, day, {{variables}})

    # Grid-based stats
    sim_stats[[n]] <- sim_stats_ini %>%
      group_by(id, mon) %>%
      dplyr::summarize(across({{variables}},
        list(mean="mean", sd="sd", skewness=e1071::"skewness"),.names = "{.col}:{.fn}")) %>%
      gather(key = variable, value = value,-id, -year, -mon) %>%
      separate(variable, c("variable","stat"), sep = ":") %>%
      mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

    # Area-averaged stats
    sim_stats_aavg[[n]] <- sim_stats_ini %>%
      group_by(year, mon) %>%
      summarize(across({{variables}},
        list(mean="mean", sd="sd", skewness=e1071::"skewness"),.names = "{.col}:{.fn}")) %>%
      gather(key = variable, value = value, -mon, -year) %>%
      separate(variable, c("variable","stat"), sep=":") %>%
      mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

    # Intersite/Intervariable correlations
    sim_stats_icor <- sim_stats_ini %>%
      gather(key = variable, value = value, -id, -year, -mon, -day) %>%
      unite(id_variable, c("id","variable"),sep = ":") %>%
      spread(id_variable, value) %>%
      dplyr::select(-year, -mon, -day) %>%
      stats::cor(.$value) %>% as_tibble()

    sim_intercor[[n]] <- sim_stats_icor %>%
      mutate(id_variable1 = colnames(sim_stats_icor), .before = 1) %>%
      gather(key = id_variable2, value = value, -id_variable1) %>%
      separate(id_variable1, c("id1","variable1"), sep =":") %>%
      separate(id_variable2, c("id2","variable2"), sep =":") %>%
      mutate(type = "Simulated")

    sim_wetdry_days[[n]] <- sim_stats_ini %>%
        left_join(mc_thresholds, by = "mon") %>%
        group_by(id, mon) %>%
        summarize(wet_count = length(which(precip >= wet_th))/sim_year_num,
                  dry_count = length(which(precip < wet_th))/sim_year_num) %>%
      mutate(type = "Simulated") %>%
      gather(key = stat, value = value, wet_count:dry_count)

    sim_wetdry_spells[[n]] <- sim_stats_ini %>%
      left_join(mc_thresholds, by = "mon") %>%
      group_by(id, mon) %>%
      summarize(dry = averageSpellLength(precip, threshold = wet_th, below = TRUE),
                wet = averageSpellLength(precip, threshold = wet_th, below = FALSE)) %>%
      group_by(id, mon) %>%
      summarize(dry = mean(dry),
                wet = mean(wet)) %>%
      gather(key = stat, value = value, wet:dry) %>%
      mutate(type = "Simulated")

  }

  sim_stats <- bind_rows(sim_stats, .id = "rlz")
  sim_stats_aavg <- bind_rows(sim_stats_aavg, .id = "rlz")
  sim_wetdry_days <- bind_rows(sim_wetdry_days, .id = "rlz")
  sim_wetdry_spells <- bind_rows(sim_wetdry_spells, .id = "rlz")
  sim_intercor <- bind_rows(sim_intercor, .id = "rlz")



  #### Calculate median statistics

  sim_stats_median <- sim_stats %>%
    group_by(id, mon, variable, stat) %>%
    summarize(value = stats::median(value)) %>%
    mutate(type = "Simulated")

  sim_wetdry_days_meadian <- sim_wetdry_days %>%
    gather(key = stat, value = value, wet_count:dry_count) %>%
    group_by(id, mon, stat) %>%
    summarize(value = stats::median(value)) %>%
    mutate(type = "Simulated")

  sim_wetdry_spells_median <- sim_wetdry_spells %>%
    group_by(id, mon, stat) %>%
    summarize(value = stats::median(value)) %>%
    mutate(type = "Simulated")

  sim_intercor_median <- sim_intercor %>%
    group_by(id1, variable1, id2, variable2) %>%
    summarize(value = stats::median(value)) %>%
    mutate(type = "Simulated")


  #### MERGE STATISTICS #############################################

  var_combs <- apply(combn(variables, 2),2, paste, collapse=":")
  id_combs  <- apply(combn(1:nsgrids, 2),2, paste, collapse=":")

  # Cross-correlation accross sites
  stats_cross <- bind_rows(hist_intercor, sim_intercor_median) %>%
    filter(variable1 == variable2) %>%
    filter(id1 != id2) %>%
    unite(id, c("id1","id2"),sep=":") %>%
    filter(id %in% id_combs) %>%
    spread(key = type, value = value) %>%
    dplyr::select(id, variable = variable1, Observed, Simulated) %>%
    mutate(variable = factor(variable, variables, variable.labels))

  # Intercorrelation
  stats_inter <- bind_rows(hist_intercor, sim_intercor_median) %>%
    filter(id1 == id2) %>%
    filter(variable1 != variable2) %>%
    unite(id, c("id1","id2"),sep=":") %>%
    unite(variable, c("variable1", "variable2"),sep=":") %>%
    filter(variable %in% var_combs) %>%
    spread(key = type, value = value)

  # Wet and Dry days
  stats_wetdry_days <- hist_wetdry_days %>%
    bind_rows(sim_wetdry_days_meadian) %>%
    spread(type, value) %>%
    mutate(stat = factor(stat, levels = c("dry_count", "wet_count"),
           labels = c("Number of Dry Days", "Number of Wet Days")))

  stats_wetdry_spells <- hist_wetdry_spells %>%
    bind_rows(sim_wetdry_spells_median) %>%
    spread(type, value) %>%
    mutate(stat = factor(stat, levels = c("dry", "wet"),
           labels = c("Dry spells", "Wet spells")))


  #:::::::::::::::::::::::::::: PLOTS ::::::::::::::::::::::::::::::::::::::::::

  base_len <- 5

  ### Wet dry spells
  xy_breaks <- pretty(unlist(stats_wetdry_spells[,c(4,5)]), 5)

  p <- ggplot(stats_wetdry_spells, aes(x = Observed, y = Simulated)) +
    theme_bw(base_size = 12) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(color = "blue") +
    facet_wrap(stat ~ ., ncol = 2) +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    labs(x = "Observed", y = "Simulated")

  ggsave(paste0(output.path,"monthly_statistics_wetdry_spells.png"),
        height = base_len, width = base_len*2)

  ### Wet dry days
  xy_breaks <- pretty(unlist(stats_wetdry_days[,c(4,5)]), 5)

  p <- ggplot(stats_wetdry_days, aes(x = Observed, y = Simulated)) +
    theme_bw(base_size = 12) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(color = "blue") +
    facet_wrap(stat ~ ., ncol = g.wdth1/base_len) +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    labs(x = "Observed", y = "Simulated")

  ggsave(paste0(output.path,"monthly_statistics_wetdry_days.png"),
         height = base_len, width = base_len*2)

  ### Cross site correlations
  xy_breaks <- pretty(unlist(stats_cross[,c(3,4)]), 5)

  p <- ggplot(stats_cross, aes(x = Observed, y = Simulated)) +
    theme_bw(base_size = 11) +
    ggtitle("Cross-site correlations") +
    geom_abline(color = "blue") +
    geom_point(alpha = 0.6, size = 1.5) +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", ncol = g.wdth1/base_len) +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks)

  ggsave(paste0(output.path,"monthly_statistics_crosscorrelation.png"),
         height = base_len*1.5, width = base_len*1.5)

  ### Intersite correlations
  xy_breaks <- pretty(unlist(stats_inter[,c(3,4)]), 5)

  p <- ggplot(stats_inter, aes(x = Observed, y = Simulated)) +
    theme_bw(base_size = 12) +
    ggtitle("Inter-site correlations") +
    geom_abline(color = "blue") +
    geom_point(alpha = 0.6, size = 1.5) +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", ncol = 3) +
    scale_x_continuous(limits = range(xy_breaks), breaks = xy_breaks) +
    scale_y_continuous(limits = range(xy_breaks), breaks = xy_breaks)

  ggsave(paste0(output.path,"monthly_statistics_intercorrelation.png" ),
         width = base_len*2, height = base_len*1.5)


  for (v in 1:length(variables)) {

    ##### MONTHLY CYCLE STATISTICS
    dat <- sim_stats_aavg %>% filter(variable == variables[v] & !is.nan(value))

    p <- ggplot(dat, aes(x = as.factor(mon), y=value)) +
      theme_bw(base_size = 12) +
      ggtitle(variable.labels[v]) +
      geom_boxplot() +
      stat_summary(fun="mean", color="blue", aes(color="Simulated mean",  geom="point")) +
      facet_wrap(~ stat, scales = "free", ncol = ceiling(num_stats/2)) +
      geom_point(data = filter(hist_stats_aavg_mon, variable == variables[v]),
                 aes(color="Observed mean",  geom="point"), size = 2) +
      scale_color_manual("", values=c("Observed mean"="red", "Simulated mean"="blue")) +
      labs(x = "", y = "") +
       theme(legend.position = c(0.875, 0.45),
        legend.background = element_rect(fill = "white", color = NA),
        legend.text=element_text(size=13))

    ggsave(paste0(output.path,"monthly_statistics_", variables[v],".png" ),
           height = base_len*1.5, width = base_len*1.5)

    #### DAILY STATISTICS ######################################################

    stats_df1v <- bind_rows(hist_stats_mon, sim_stats_median) %>%
      filter(variable == variables[v]) %>%
      spread(key = type, value = value) %>%
      filter(stat != "Skewness")

    p <- ggplot(stats_df1v, aes(x = Observed, y = Simulated)) +
      theme_light(base_size = 12) +
      ggtitle(variable.labels[v]) +
      geom_point(alpha = 0.4, size = 1.5) +
      geom_abline(color = "blue") +
      labs(x = "Observed", y = "Simulated") +
      facet_wrap(stat ~ ., scales = "free", ncol = ceiling(num_stats/2))

    ggsave(paste0(output.path,"daily_statistics_", variables[v],".png" ),
           height = base_len, width = base_len*1.5)

  }

}







