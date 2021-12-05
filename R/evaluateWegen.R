

################################################################################

# Comparison figures

# 1) Monthly statistics (mean, s.deviation, skew) x num.variables
# 2) Intersite correlations x num.variables
# 3) Cross-correlations x num.variables

#' Function to assess the weather generator skill
#'
#' @param daily.sim A matrix of daily time-series of weather variables
#' @param daily.obs A vector of daily time-series of observed weather series
#' @param output.path  A character string to define the path of resulting netcdf file.
#' @param variables A vector object specifying the names of weather variables.
#' @param variable.labels A vector object specifying the names of weather variable labels.
#' @param variable.units A vector object specifying the names of weather variable units.
#' @param realization.num Placeholder
#'
#' @return
#' @export
#' @import utils ggplot2 dplyr
#' @importFrom stats cor median
#' @importFrom e1071 skewness
#' @importFrom rlang .data
evaluateWegen <- function(
  daily.sim = NULL,
  daily.obs = NULL,
  output.path = NULL,
  variables = NULL,
  variable.labels = NULL,
  variable.units = NULL,
  realization.num = NULL)

  {

  #Workaround for rlang warning
  year <- mon <- day <- precip <- sd <- variable <- value <- id_variable <- 0
  wet_count <- dry_count <- wet <- dry <- id_variable1 <- id_variable2 <- 0
  rlz <- id1 <- id2 <- variable1 <- variable2 <- type <- Observed <- Simulated <- 0


  stat_level <- c("mean", "sd", "skewness")
  stat_label <- c("Mean", "Standard Deviation", "Skewness")

  nsgrids <- length(daily.sim[[1]])
  variable_labels2 <- paste0(variable.labels, " (", variable.units, ")")

  num_stats <- length(stat_level)
  num_vars <- length(variables)
  var_combs <- expand_grid(var1=variables, var2 = variables)
  num_var_combs <- choose(num_vars, 2)
  sim_year_num <- nrow(daily.sim[[1]][[1]])/365

  sim_stats <- NULL
  sim_stats_aavg <- NULL
  sim_stats_icor <- NULL
  sim_wetdry_days <-NULL
  sim_wet_spells <- NULL
  sim_dry_spells <- NULL

  for (n in 1:realization.num) {

    # Calculate for each simulated trace
    daily_sim_tbl <- daily.sim[[1]] %>% bind_rows(.id = "id")
    daily_sim_tbl$year <- as.numeric(format(daily_sim_tbl$date,"%Y"))
    daily_sim_tbl$mon <- as.numeric(format(daily_sim_tbl$date,"%m"))
    daily_sim_tbl$day <- as.numeric(format(daily_sim_tbl$date,"%d"))
    daily_sim_tbl$id <- as.numeric(daily_sim_tbl$id)

    # Calculate wet spells
    sim_wet_spells <- bind_rows(sim_wet_spells,
      lapply(1:nsgrids, function(x)
         table(calculateSpellLength(daily.sim[[n]][[x]]$precip, below = FALSE)) %>%
         tibble(length = as.numeric(names(.)), wet = .)) %>%
      bind_rows(.id = "id") %>%
      mutate(rlz = n, .before = 1))

    # Calculate dry spells
    sim_dry_spells <- bind_rows(sim_dry_spells,
      lapply(1:nsgrids, function(x)
      table(calculateSpellLength(daily.sim[[n]][[x]]$precip, below = TRUE)) %>%
         tibble(length = as.numeric(names(.)), dry = .)) %>%
      bind_rows(.id = "id") %>%
      mutate(rlz = n, .before = 1))

    # Calculate wet and dry days
    sim_wetdry_days <- bind_rows(sim_wetdry_days,
      daily_sim_tbl %>%
      select(id, year, mon, day, precip) %>%
      group_by(id, mon) %>%
      summarize(wet_count = length(precip[precip!=0])/sim_year_num,
                dry_count = length(precip[precip==0])/sim_year_num) %>%
      mutate(rlz = n, .before = 1))

    # Grid-based statistics (id:sample grids, rlz = realizations)
    sim_stats <- bind_rows(sim_stats,
       daily_sim_tbl %>%
       select(id, year, mon, {{variables}}) %>%
       group_by(id, year, mon) %>%
       summarize(across({{variables}}, list(mean="mean", sd="sd", skewness=e1071::"skewness"),
         .names = "{.col}:{.fn}")) %>%
       gather(key = variable, value = value, -id, -year, -mon) %>%
       separate(variable, c("variable","stat"), sep =":") %>%
       mutate(rlz = n, .before = 1))

    # Area-averaged statistics
    sim_stats_aavg <- bind_rows(sim_stats_aavg,
        daily_sim_tbl %>%
        select(id, year, mon, {{variables}}) %>%
        group_by(year, mon) %>%
        summarize(across({{variables}}, list(mean="mean", sd="sd", skewness=e1071::"skewness"),
          .names = "{.col}:{.fn}")) %>%
        gather(key = variable, value = value, -year, -mon) %>%
        separate(variable, c("variable","stat"), sep =":") %>%
        mutate(rlz = n, .before = 1)
    )

    # Intersite/cross-site correlations
    sim_stats_icor <- bind_rows(sim_stats_icor,
        daily_sim_tbl %>%
        dplyr::select(id, year, mon, day, {{variables}}) %>%
        gather(key = variable, value = value, -id, -year, -mon, -day) %>%
        unite(id_variable, c("id","variable"), sep = ":") %>%
        spread(id_variable, value) %>%
        dplyr::select(-year, -mon, -day) %>%
        cor(.$value) %>% as_tibble() %>%
        mutate(rlz = n, .before = 1)
    )

  }

  sim_stats <- sim_stats %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  sim_stats_aavg <- sim_stats_aavg %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  hist_stats_ini <- do.call("rbind", lapply(1:length(daily.obs),
    function(x) mutate(daily.obs[[x]], id = x))) %>%
    mutate(year = as.numeric(format(date,"%Y")),
           mon = as.numeric(format(date,"%m")),
           day = as.numeric(format(date,"%d")),
           id = as.numeric(id))

  hist_stats <- hist_stats_ini %>%
    group_by(id, year, mon) %>%
    dplyr::summarize(across({{variables}},
      list(mean="mean", sd="sd", skewness=e1071::"skewness"),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = value,-id, -year, -mon) %>%
    separate(variable, c("variable","stat"), sep = ":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  hist_stats_aavg <- hist_stats_ini %>%
    group_by(mon) %>%
    summarize(across({{variables}},
      list(mean="mean", sd="sd", skewness=e1071::"skewness"),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = value,-mon) %>%
    separate(variable, c("variable","stat"), sep=":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  hist_stats_icor <- hist_stats_ini %>%
    dplyr::select(id, year, mon, day, {{variables}}) %>%
    gather(key = variable, value = value, -id, -year, -mon, -day) %>%
    unite(id_variable, c("id","variable"),sep = ":") %>%
    spread(id_variable, value) %>%
    dplyr::select(-year, -mon, -day) %>%
    cor(.$value) %>% as_tibble()

  hist_wetdry_days <- hist_stats_ini %>%
      select(id, year, mon, day, precip) %>%
      group_by(id, mon) %>%
      summarize(wet_count = length(precip[precip!=0])/sim_year_num,
                 dry_count = length(precip[precip==0])/sim_year_num) %>%
    mutate(type = "Observed") %>%
    gather(key = stat, value = value, wet_count:dry_count)

  #sim_wetdry_spells
  hist_dry_spells <- lapply(1:nsgrids, function(x)
      table(calculateSpellLength(daily.obs[[x]]$precip, below = TRUE)) %>%
         tibble(length = as.numeric(names(.)), dry = as.numeric(.))) %>%
      bind_rows(.id = "id") %>%
    select(-.)

  hist_wet_spells <- lapply(1:nsgrids, function(x)
     table(calculateSpellLength(daily.obs[[x]]$precip, below = FALSE)) %>%
         tibble(length = as.numeric(names(.)), wet = as.numeric(.))) %>%
      bind_rows(.id = "id") %>%
    select(-.)

  hist_wetdry_spells <- hist_dry_spells %>%
    full_join(hist_wet_spells, by = c("id", "length")) %>%
    gather(key = stat, value = value, wet:dry) %>%
    mutate(type = "Observed")

  # :::::::::: COMPARE MONTHLY STATISTICS ::::::::::::::::::::::::::::::::::::::

  sim_stats_mon_median <- sim_stats %>%
    group_by(id, mon, variable, stat) %>%
    summarize(value = median(value)) %>%
    mutate(type = "Simulated")

  hist_stats_mon <- hist_stats %>%
    group_by(id, mon, variable, stat) %>%
    summarize(value = median(value)) %>%
    mutate(type = "Observed")

  sim_wetdry_days_meadian <- sim_wetdry_days %>%
    gather(key = stat, value = value, wet_count:dry_count) %>%
    group_by(id, mon, stat) %>%
    summarize(value = median(value)) %>%
    mutate(type = "Simulated")

  sim_wet_spells_median <- sim_wet_spells %>%
    group_by(id, length) %>%
    summarize(wet = median(wet))

  sim_dry_spells_median <- sim_dry_spells %>%
    group_by(id, length) %>%
    summarize(dry = median(dry))

  sim_wetdry_spells_median <- sim_wet_spells_median %>%
    full_join(sim_dry_spells_median, by = c("id", "length")) %>%
    gather(key = stat, value = value, wet:dry) %>%
    mutate(type = "Simulated")

  ### Intersite/cross-site correlations ::::::::::::::::::::::::::::::::::::::::

  hist_intercor <- hist_stats_icor %>%
    mutate(id_variable1 = colnames(hist_stats_icor), .before = 1) %>%
    gather(key = id_variable2, value = value, -id_variable1) %>%
    separate(id_variable1, c("id1","variable1"), sep =":") %>%
    separate(id_variable2, c("id2","variable2"), sep =":") %>%
    mutate(type = "Observed")

  sim_intercor_median <- sim_stats_icor %>%
    mutate(id_variable1 = rep(colnames(sim_stats_icor)[-1], realization.num), .before = 2) %>%
    gather(key = id_variable2, value = value, -rlz, -id_variable1) %>%
    separate(id_variable1, c("id1","variable1"),sep=":") %>%
    separate(id_variable2, c("id2","variable2"),sep=":") %>%
    group_by(id1, variable1, id2, variable2) %>%
    summarize(value = median(value)) %>%
    mutate(type = "Simulated")

  var_combs <- apply(combn(variables, 2),2, paste, collapse=":")
  id_combs  <- apply(combn(1:nsgrids, 2),2, paste, collapse=":")

  stats_cross <- bind_rows(hist_intercor, sim_intercor_median) %>%
    filter(variable1 == variable2) %>%
    filter(id1 != id2) %>%
    unite(id, c("id1","id2"),sep=":") %>%
    filter(id %in% id_combs) %>%
    spread(key = type, value = value) %>%
    dplyr::select(id, variable = variable1, Observed, Simulated) %>%
    mutate(variable = factor(variable, variables, variable.labels))

  stats_inter <- bind_rows(hist_intercor, sim_intercor_median) %>%
    filter(id1 == id2) %>%
    filter(variable1 != variable2) %>%
    unite(id, c("id1","id2"),sep=":") %>%
    unite(variable, c("variable1", "variable2"),sep=":") %>%
    filter(variable %in% var_combs) %>%
    spread(key = type, value = value)

  stats_wetdry_days <- bind_rows(hist_wetdry_days, sim_wetdry_days_meadian) %>%
    spread(type, value) %>%
    mutate(stat = factor(stat, levels = c("dry_count", "wet_count"),
           labels = c("Number of Dry Days", "Number of Wet Days")))

  stats_wetdry_spells <- bind_rows(hist_wetdry_spells, sim_wetdry_spells_median) %>%
    spread(type, value) %>%
    mutate(stat = factor(stat, levels = c("dry", "wet"),
           labels = c("Dry spells", "Wet spells")))


  #:::::::::::::::::::::::::::: PLOTS ::::::::::::::::::::::::::::::::::::::::::

  base_len <- 5

  g.hght1 <- ifelse(num_vars > 2, base_len*2, base_len)
  g.wdth1 <- ifelse(num_vars > 1, base_len*2, base_len)

  g.hght2 <- if(num_var_combs > 2) base_len*2 else if (num_var_combs > 4) base_len*3 else base_len
  g.wdth2 <- if(num_var_combs > 1) base_len*2 else if (num_var_combs > 6) base_len*3 else base_len

  ### Wet dry spells
  p <- ggplot(stats_wetdry_spells, aes(x = Observed/sim_year_num, y = Simulated/sim_year_num)) +
    theme_light(base_size = 12) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(color = "blue", size = 1) +
    facet_wrap(stat ~ ., scales = "free", ncol = g.wdth1/base_len) +
    labs(x = "Observed", y = "Simulated")

  ggsave(paste0(output.path,"monthly_stats_wet_dry_spells.png"),
        height = base_len, width = base_len*2)

  ### Wet dry days
  p <- ggplot(stats_wetdry_days, aes(x = Observed, y = Simulated)) +
    theme_light(base_size = 12) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(color = "blue", size = 1) +
    facet_wrap(stat ~ ., scales = "free", ncol = g.wdth1/base_len) +
    labs(x = "Observed", y = "Simulated")

  ggsave(paste0(output.path,"monthly_stats_wet_dry_days.png"),
         height = base_len, width = base_len*2)

  ### Cross site correlations
  p <- ggplot(stats_cross, aes(x = Observed, y = Simulated)) +
    theme_light(base_size = 12) +
    ggtitle("Cross-site correlations") +
    geom_abline(color = "blue", size = 1) +
    geom_point(alpha = 0.6, size = 1.5) +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", ncol = g.wdth1/base_len) +
    scale_x_continuous(breaks = pretty(c(-0.5,1), n=5), limits = c(-0.5, 1)) +
    scale_y_continuous(breaks = pretty(c(-0.5,1), n=5), limits = c(-0.5, 1))

  ggsave(paste0(output.path,"monthly_stats_all_crosssite.png"),
         height = g.hght1, width = g.wdth1)

  ### Intersite correlations
  p <- ggplot(stats_inter, aes(x = Observed, y = Simulated)) +
    theme_light(base_size = 12) +
    ggtitle("Inter-site correlations") +
    geom_abline(color = "blue", size = 1) +
    geom_point(alpha = 0.6, size = 1.5) +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", ncol = g.wdth2/base_len) +
    scale_x_continuous(breaks = pretty(c(-0.5,1), n=5), limits = c(-0.5, 1)) +
    scale_y_continuous(breaks = pretty(c(-0.5,1), n=5), limits = c(-0.5, 1))

  ggsave(paste0(output.path,"monthly_stats_all_intersite.png" ),
         width = g.wdth2, height = g.hght2)

  for (v in 1:length(variables)) {

    ##### MONTHLY CYCLE STATISTICS
    dat <- sim_stats_aavg %>% filter(variable == variables[v] & !is.nan(value))
    p <- ggplot(dat, aes(x = as.factor(mon), y=value)) +
      theme_light(base_size = 12) +
      ggtitle(variable_labels2[v]) +
      geom_boxplot() +
      stat_summary(fun="mean", color="blue", aes(color="Simulated nmean",  geom="point")) +
      facet_wrap(~ stat, scales = "free", ncol = ceiling(num_stats/2)) +
      geom_point(data = filter(hist_stats_aavg, variable == variables[v]),
                 aes(color="Observed mean",  geom="point"), size = 2) +
      scale_color_manual("", values=c("Observed mean"="red", "Simulated mean"="blue")) +
      labs(x = "", y = "") +
       theme(legend.position = c(0.875, 0.45),
        legend.background = element_rect(fill = "white", color = NA),
        legend.text=element_text(size=13))



    ggsave(paste0(output.path,"monthly_stats_", variables[v],".png" ),
           height = base_len*2, width = base_len*2)

    #### DAILY STATISTICS ######################################################

    stats_df1v <- bind_rows(hist_stats_mon, sim_stats_mon_median) %>%
      filter(variable == variables[v]) %>%
      spread(key = type, value = value) %>%
      filter(stat != "Skewness")

    p <- ggplot(stats_df1v, aes(x = Observed, y = Simulated)) +
      theme_light(base_size = 12) +
      ggtitle(variable.labels[v]) +
      geom_point(alpha = 0.4, size = 1.5) +
      geom_abline(color = "blue", size = 1) +
      labs(x = "Observed", y = "Simulated") +
      facet_wrap(stat ~ ., scales = "free", ncol = ceiling(num_stats/2))

    ggsave(paste0(output.path,"daily_stats_", variables[v],".png" ),
           height = base_len, width = g.wdth2)

  }

}

