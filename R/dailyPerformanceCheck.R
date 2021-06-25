

################################################################################

# Comparison figures

# 1) Monthly statistics (mean, s.deviation, skew) x num.variables
# 2) Intersite correlations x num.variables
# 3) Cross-correlations x num.variables



#' Function to assess the weather generator skill
#'
#' @param daily.sim A matrix of daily time-series of weather variables
#' @param daily.obs A vector of daily time-series of observed weather series
#' @param out.path  A character string to define the path of resulting netcdf file.
#' @param variables A vector object specifying the names of weather variables.
#' @param variable.labels A vector object specifying the names of weather variable labels.
#' @param variable.units A vector object specifying the names of weather variable units.
#'
#' @return
#' @export
#'
dailyPerformanceCheck <- function(
  daily.sim = NULL,
  daily.obs = NULL,
  out.path = NULL,
  variables = NULL,
  variable.labels = NULL,
  variable.units = NULL)

  {

  nsgrids <- length(daily.sim[[1]])

  sim_stats <- NULL
  sim_stats_aavg <- NULL
  sim_stats_icor <- NULL

  stat_level <- c("mean", "sd", "skewness")
  stat_label <- c("Mean", "Standard Deviation", "Skewness")

  num_stats <- length(stat_level)
  num_vars <- length(variables)
  var_combs <- expand_grid(var1=variables, var2=variables)
  crossing(var1=variables, var2=variables)
  num_var_combs <- choose(num_vars,2)

  # Calculate for each simulated trace
  for (n in 1:nmax) {

    # Grid-based statistics (id:sample grids, rlz = realizations)
    sim_stats <- bind_rows(sim_stats,
       daily.sim[[n]] %>%
       bind_rows(.id = "id") %>%
       mutate(year = year(date), mon = month(date), id = as.numeric(id)) %>%
       group_by(id, year, mon) %>%
       summarize(across({{variables}}, list(mean=mean, sd=sd, skewness=skewness), .names = "{.col}:{.fn}")) %>%
       gather(key = variable, value = value, -id, -year, -mon) %>%
       separate(variable, c("variable","stat"), sep =":") %>%
       add_column(rlz = n, .before = 1)
    )

    # Area-averaged statistics
    sim_stats_aavg <- bind_rows(sim_stats_aavg,
        daily.sim[[n]] %>%
        bind_rows(.id = "id") %>%
        mutate(year = year(date), mon = month(date)) %>%
        group_by(year, mon) %>%
        summarize(across({{variables}}, list(mean=mean, sd=sd, skewness=skewness), .names = "{.col}:{.fn}")) %>%
        gather(key = variable, value = value, -year, -mon) %>%
        separate(variable, c("variable","stat"), sep =":") %>%
        add_column(rlz = n, .before = 1)
    )

    # Intersite/cross-site correlations
    sim_stats_icor <- bind_rows(sim_stats_icor,
        daily.sim[[n]] %>%
        bind_rows(.id = "id") %>%
        mutate(year = year(date), mon = month(date), day = day(date), id = as.numeric(id)) %>%
        dplyr::select(id, year, mon, day, {{variables}}) %>%
        gather(key = variable, value = value, -id, -year, -mon, -day) %>%
        unite(id_variable, c("id","variable"), sep = ":") %>%
        spread(id_variable, value) %>%
        dplyr::select(-year, -mon, -day) %>%
        cor(.$value) %>% as_tibble() %>%
        add_column(rlz = n, .before = 1)
    )

  }

  sim_stats <- sim_stats %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  sim_stats_aavg <- sim_stats_aavg %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  hist_stats <- lapply(1:length(daily.obs), function(x) mutate(daily.obs[[x]], id = x)) %>%
    do.call("rbind", .) %>%
    mutate(year = year(date), mon = month(date), id = as.numeric(id)) %>%
    group_by(id, year, mon) %>%
    summarize(across({{variables}}, list(mean=mean, sd=sd, skewness=skewness),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = value,-id, -year, -mon) %>%
    separate(variable, c("variable","stat"), sep = ":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  hist_stats_aavg <- lapply(1:length(daily.obs), function(x) mutate(daily.obs[[x]], id = x)) %>%
    do.call("rbind", .) %>%
    mutate(mon = month(date)) %>%
    group_by(mon) %>%
    summarize(across({{variables}}, list(mean=mean, sd=sd, skewness=skewness),.names = "{.col}:{.fn}")) %>%
    gather(key = variable, value = value,-mon) %>%
    separate(variable, c("variable","stat"), sep=":") %>%
    mutate(stat = factor(stat, levels = stat_level, labels = stat_label))

  hist_stats_icor <- lapply(1:length(daily.obs), function(x) mutate(daily.obs[[x]], id = x)) %>%
    do.call("rbind", .) %>%
    mutate(year = year(date), mon = month(date), day = day(date)) %>%
    dplyr::select(id, year, mon, day, {{variables}}) %>%
    gather(key = variable, value = value, -id, -year, -mon, -day) %>%
    unite(id_variable, c("id","variable"),sep = ":") %>%
    spread(id_variable, value) %>%
    dplyr::select(-year, -mon, -day) %>%
    cor(.$value) %>% as_tibble()

  # :::::::::: COMPARE MONTHLY STATISTICS ::::::::::::::::::::::::::::::::::::::

  sim_stats_mon_median <- sim_stats %>%  group_by(id, mon, variable, stat) %>%
    summarize(value = median(value)) %>% mutate(type = "Simulated")

  hist_stats_mon <- hist_stats %>% group_by(id, mon, variable, stat) %>%
    summarize(value = median(value)) %>% mutate(type = "Observed")

  ### Intersite/cross-site correlations ::::::::::::::::::::::::::::::::::::::::

  hist_intercor <- hist_stats_icor %>%
    mutate(id_variable1 = colnames(hist_stats_icor), .before = 1) %>%
    gather(key = id_variable2, value = value, -id_variable1) %>%
    separate(id_variable1, c("id1","variable1"), sep =":") %>%
    separate(id_variable2, c("id2","variable2"), sep =":") %>%
    mutate(type = "Observed")

  sim_intercor_median <- sim_stats_icor %>%
    mutate(id_variable1 = rep(colnames(sim_stats_icor)[-1], nmax), .before = 2) %>%
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


  #:::::::::::::::::::::::::::: PLOTS ::::::::::::::::::::::::::::::::::::::::::

  ### Cross site correlations
  p <- ggplot(stats_cross, aes(x = Observed, y = Simulated)) +
    theme_light(base_size = 12) +
    ggtitle("Cross-site correlations") +
    geom_abline(color = "blue", size = 1) +
    geom_point(alpha = 0.6, size = 1.5) +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", ncol = ceiling(num_vars/2)) +
    scale_x_continuous(breaks = scales::pretty_breaks(5), limits = c(0, 1)) +
    scale_y_continuous(breaks = scales::pretty_breaks(5), limits = c(0, 1))

  ggsave(paste0(out.path,"monthly_correlation_crosssite.png" ),
         height = 4 * ceiling(num_vars/2), width = 4 *  ifelse(num_vars > 1, 2, 1))

  ### Intersite correlations
  p <- ggplot(stats_inter, aes(x = Observed, y = Simulated)) +
    theme_light(base_size = 12) +
    ggtitle("Inter-site correlations") +
    geom_abline(color = "blue", size = 1) +
    geom_point(alpha = 0.6, size = 1.5) +
    labs(x = "Observed", y = "Simulated") +
    facet_wrap(variable ~ ., scales = "free", ncol = ceiling(num_var_combs/2)) +
    scale_x_continuous(breaks = scales::pretty_breaks(5), limits =c(0,1)) +
    scale_y_continuous(breaks = scales::pretty_breaks(5), limits =c(0,1))

  ggsave(paste0(out.path,"monthly_correlation_intersite.png" ),
         width = 4 * ceiling(num_var_combs/2), height = 4 *  ifelse(num_var_combs > 1, 2, 1))


  for (v in 1:length(variables)) {

    ##### MONTHLY CYCLE STATISTICS
    p <- ggplot(filter(sim_stats_aavg, variable == variables[v]),
                aes(x = as.factor(mon), y=value)) +
      theme_light(base_size = 12) +
      ggtitle(variable.labels[v]) +
      geom_boxplot() +
      stat_summary(fun.y="mean", color="blue")+
      facet_wrap(~ stat, scales = "free", ncol = ceiling(num_stats/2)) +
      geom_point(data = filter(hist_stats_aavg, variable == variables[v]),
                 color = "red") +
      labs(x = "Months", y = variable.units[v])


    ggsave(paste0(out.path,"monthly_stats_", variables[v],".png" ),
           height = 4 * ceiling(num_var_combs/2),
           width = 4 * ifelse(num_var_combs > 1, 2, 1))

    #### DAILY STATISTICS
    stats_df1v <- bind_rows(hist_stats_mon, sim_stats_mon_median) %>%
      filter(variable == variables[v]) %>%
      spread(key = type, value = value)

    p <- ggplot(stats_df1v, aes(x = Observed, y = Simulated)) +
      theme_light(base_size = 12) +
      ggtitle(variable.labels[v]) +
      geom_point(alpha = 0.4, size = 1.5) +
      geom_abline(color = "blue", size = 1) +
      labs(x = "Observed", y = "Simulated") +
      facet_wrap(stat ~ ., scales = "free", ncol = ceiling(num_stats/2)) +
      scale_x_continuous(breaks = scales::pretty_breaks(5), limits = c(0, NA)) +
      scale_y_continuous(breaks = scales::pretty_breaks(5), limits = c(0, NA)) +
      theme(plot.margin = unit(c(0.5,0.2,0.2,0.2), "cm"))

    ggsave(paste0(out.path,"daily_stats_", variables[v],".png" ),
           height = 4 * ceiling(num_var_combs/2),
           width = 4 * ifelse(num_var_combs > 1, 2, 1))

  }



}







