
#' Resample from WARM outputs
#'
#' @param series.sim  A numeric matrix, with simulated time-series.
#' @param series.obs  A numeric vector of observd time-series values.
#' @param save.plots A logical, to save the plots to file.
#' @param power.obs A numeric vector of power spectra of observed time-series.
#' @param power.sim A numeric matrix of power spectrum of simulated time-series.
#' @param power.period A time-series of power periods calculated.
#' @param power.signif A time-series of power significance.
#' @param sample.num A numeric value to define the final sample size.
#' @param seed A numeric value to define a seed for resampling.
#' @param save.series A logical to write the results to csv files.
#' @param verbose A logical to decide if further information to be displayed on the screen.
#' @param output.path Output folder path
#' @param padding placeholder
#' @param bounds placeholder
#'
#' @return
#' @export
#' @import dplyr
#' @import ggplot2
waveletARSubset <- function(
  series.obs = NULL,
  series.sim = NULL,
  power.obs = NULL,
  power.sim = NULL,
  power.period =  NULL,
  power.signif =  NULL,
  sample.num = 5,
  seed = NULL,
  save.plots = TRUE,
  save.series = TRUE,
  verbose = FALSE,
  output.path = NULL,
  padding = TRUE,
  bounds = list(
    mean = c(0.95,1.05),
    sd = c(0.90,1.10),
    min = c(0.90,1.10),
    max = c(0.90,1.10),
    power = c(0.60,2.25),
    nonsignif.threshold = 0.90)
  )

{

  # Workaround for rlang warning
  sim <- value <- yind <- par <- type <- variable <- y <- x <- 0


  sim.year.num = nrow(series.sim)

  # Statistics for simulated realizations
  stats_sim <- series.sim %>%
    as_tibble(.name_repair = ~as.character(1:ncol(series.sim))) %>%
    mutate(yind = 1:n()) %>%
    gather(key = sim, value = value, -yind) %>%
    mutate(sim = as.numeric(sim)) %>%
    group_by(sim) %>%
    summarize(mean = mean(value), sd = stats::sd(value),
              max = max(value), min = min(value))

  # Statistics for observed weather series
  stats_obs <- tibble(value = series.obs) %>%
    summarize(mean = mean(value), sd = stats::sd(value),
              max = max(value), min = min(value))

  # Significant periods
  periods_sig <- which(power.obs > power.signif)
  periods_sig <- periods_sig[periods_sig %in% 1:length(series.obs)]

  if(isTRUE(padding)) {
    periods_sig <- sort(intersect(unique(c(periods_sig-1, periods_sig, periods_sig+1)), 1:length(power.signif)))
  }
  periods_nonsig <- setdiff(1:length(power.signif),periods_sig)
  periods_nonsig <- periods_nonsig[periods_nonsig %in% 1:dim(power.sim)[1]]


  if (!is.null(bounds$power)) {

    # Filter scenarios have significant signals
    sub_power1 <- which(sapply(1:ncol(power.sim), function(x)
      any(power.sim[periods_sig,x] > power.signif[periods_sig])))

    # Signals within the bounds
    sub_power2 <- which(sapply(1:ncol(power.sim), function(x)
        all((power.sim[periods_sig,x] > power.obs[periods_sig] * bounds$power[1]) &
            (power.sim[periods_sig,x] < power.obs[periods_sig] * bounds$power[2]))))

    # Nonsignificant below threshold
    sub_power3 <- which(sapply(1:ncol(power.sim), function(x)
        all((power.sim[periods_nonsig,x] < power.signif[periods_nonsig]*bounds$nonsignif.threshold))))

    sub_power <- base::intersect(base::intersect(sub_power1, sub_power2), sub_power3)
  } else {
    sub_power  <- 1:ncol(series.sim)
  }


  if (!is.null(bounds$mean)) {
    sub_mean  <- which((stats_sim$mean > stats_obs$mean * bounds$mean[1]) &
                       (stats_sim$mean < stats_obs$mean * bounds$mean[2]))
  } else {
    sub_mean <- 1:ncol(series.sim)
  }


  if (!is.null(bounds$sd)) {
    sub_sd  <- which((stats_sim$sd > stats_obs$sd * bounds$sd[1]) &
                       (stats_sim$sd < stats_obs$sd * bounds$sd[2]))
  } else {
    sub_sd <- 1:ncol(series.sim)
  }


  if (!is.null(bounds$min)) {
    sub_min  <- which((stats_sim$min > stats_obs$min * bounds$min[1]) &
                      (stats_sim$min < stats_obs$min * bounds$min[2]))
  } else {
    sub_min  <- 1:ncol(series.sim)
  }

  if (!is.null(bounds$max)) {
    sub_max  <- which((stats_sim$max > stats_obs$max * bounds$max[1]) &
                      (stats_sim$max < stats_obs$max * bounds$max[2]))
  } else {
    sub_max  <- 1:ncol(series.sim)
  }

  #Select intersection
  sub_clim <- Reduce(base::intersect, list(sub_mean, sub_sd, sub_power,
                                     sub_min, sub_max))

  # Stochastically select from the initial dataset
  if(!is.null(seed)) set.seed(seed)
  sub_sample <- sample(sub_clim, min(sample.num, length(sub_clim)))

  if(isTRUE(verbose)) {
    print(tribble(
        ~"criteria", ~"# traces",
        "mean",  paste0(length(sub_mean), " out of ", ncol(series.sim)),
        "stdev", paste0(length(sub_sd), " out of ", ncol(series.sim)),
        "power", paste0(length(sub_power), " out of ", ncol(series.sim)),
        "min",   paste0(length(sub_min), " out of ", ncol(series.sim)),
        "max",   paste0(length(sub_max), " out of ", ncol(series.sim)),
        "final", paste0(length(sub_clim), " out of ", ncol(series.sim))))
  }

  if(length(sub_sample) < sample.num) {
    stop('subsetted traces less than the desired amount.
         Please relax the decision criteria and repeat.')
  }

  if(isTRUE(save.plots)) {

    ### Global Wavelet Spectral Plot
    pl <- min(length(power.period),dim(power.sim)[1])
    plr <- 5 * ceiling(power.period[pl] / 5)

    # For all matching realizations
    p <- waveletPlot(power.period = power.period[1:pl],
                   power.signif = power.signif[1:pl],
                   power.obs = power.obs[1:pl],
                   power.sim = power.sim[1:pl,sub_clim])  +
          scale_x_continuous(breaks=seq(5,plr,5), limits=c(0,plr), expand=c(0,0))

    ggsave(paste0(output.path, "warm_sim_matching_spectral.png"), width=8, height=6)

    # For subsetted realizations only
    p <- waveletPlot(power.period = power.period[1:pl],
                   power.signif = power.signif[1:pl],
                   power.obs = power.obs[1:pl],
                   power.sim = power.sim[1:pl,sub_sample])  +
          scale_x_continuous(breaks=seq(5,plr,5), limits=c(0,plr), expand=c(0,0))

    ggsave(paste0(output.path, "warm_sim_sample_spectral.png"), width=8, height=6)


    # Boxplots of all stats
    stats_obs_gg <- stats_obs %>% mutate(sim=1) %>%
      gather(key = par, value = value, -sim) %>% mutate(type = "Observed")
    stats_sim_gg <- stats_sim %>%
      gather(key = par, value = value, -sim) %>% mutate(type = "Simulated")

    stats_all <- bind_rows(stats_obs_gg, stats_sim_gg %>%
      filter(sim %in% sub_clim)) %>%
      mutate(type = factor(type, levels = c("Simulated", "Observed"))) %>%
      arrange(type)

    p <- ggplot(mapping = aes(x = par, y = value)) +
      theme_light() +
      facet_wrap(~par, scales = "free", drop = TRUE, nrow = 1) +
      geom_boxplot(data = stats_sim_gg, color = "gray60", outlier.shape = NA) +
      geom_violin(data = stats_sim_gg, color = "gray60") +
      #geom_jitter(data = stats_sim_gg, alpha = 0.1, color = "gray60") +
      geom_point(aes(fill = type), data = stats_all, size = 5, color = "white",
                 shape = 21) +
      scale_fill_manual(values = c("Simulated"="black","Observed"="blue")) +
      labs(x="", y = "", color = "", fill="") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())

      ggsave(paste0(output.path, "warm_sim_sample_stats.png"), width=8, height=6)

    # Plot simulated warm-series
    sub_clim_plot <- sub_clim[1:min(length(sub_clim), 50)]

    df1 <- series.sim[,sub_clim_plot] %>%
      as_tibble(.name_repair = ~paste0("rlz",1:length(sub_clim_plot))) %>%
      mutate(x = 1:sim.year.num) %>%
      gather(key = variable, value=y, -x) %>%
      mutate(y = y * 365)

    df2 <- tibble(x=1:length(series.obs), y = series.obs*365)

    p <- ggplot(df1, aes(x = x, y = y)) +
      theme_light(base_size = 12) +
      geom_line(aes(y = y, group = variable), color = "gray60", alpha = 0.6) +
      geom_line(aes(y=y), data = df2, color = "black", size = 1) +
      scale_x_continuous(limits = c(0,sim.year.num), breaks = seq(0,sim.year.num, 5)) +
      guides(color = "none") +
      labs(y = "Precipitation (mm/year)", x = "Year index")

    ggsave(paste0(output.path, "warm_sim_timeseries.png"), height = 5, width = 10)


  }

  if(isTRUE(save.series)) {

    utils::write.csv(x = series.sim[,sub_clim], row.names = FALSE,
        file = paste0(output.path, "warm_sim_matching.csv"))

    utils::write.csv(x = series.sim[,sub_sample], row.names = FALSE,
        file = paste0(output.path, "warm_sim_sample.csv"))

    }

  return(list(subsetted = series.sim[,sub_clim],
    sampled = series.sim[,sub_sample, drop=FALSE]))

}

