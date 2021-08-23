


#' A function to resample from WARM model results
#'
#' @param mean.bounds A numeric vector, defining the minimum and maximum limits for the mean value of the time-series when resampling.
#' @param sdev.bounds A numeric vector, defining the minimum and maximum limits for the standard deviation of the time-series when resampling.
#' @param max.bounds  A numeric vector, defining the minimum and maximum limits for the maximum value of the time-series when resampling.
#' @param min.bounds  A numeric vector, defining the minimum and maximum limits for the minimum value of the time-series when resampling.
#' @param series.sim  A numeric matrix, with simulated time-series.
#' @param series.obs  A numeric vector of observd time-series values.
#' @param seed A numeric value to define a seed for resampling.
#' @param save.plots A logical, to save the plots to file.
#' @param power.obs A numeric vector of power spectra of observed time-series.
#' @param power.sim A numeric matrix of power spectrum of simulated time-series.
#' @param power.period A time-series of power periods calculated.
#' @param power.signif A time-series of power significance.
#' @param power.bounds A numeric vector, defining the minimum and maximum limits for power spectra.
#' @param nonsig.threshold A numeric vector to define a resampling threshold for sampling.
#' @param nmax A numeric value to define the final sample size.
#' @param save.series A logical to write the results to csv files.
#' @param verbose A logical to decide if further information to be displayed on the screen.
#' @param out.path Output folder path
#'
#' @return
#' @export
#' @import dplyr
#' @importFrom utils write.csv
waveletARSubset <- function(
  series.obs = NULL,
  series.sim = NULL,
  power.obs = NULL,
  power.sim = NULL,
  power.period =  NULL,
  power.signif =  NULL,
  nmax = NULL,
  seed = NULL,
  save.plots = TRUE,
  save.series = TRUE,
  verbose = FALSE,
  out.path = out_path,
  mean.bounds = NULL,
  sdev.bounds = NULL,
  max.bounds  = NULL,
  min.bounds  = NULL,
  power.bounds = NULL,
  nonsig.threshold = NULL)

{

  # Statistics for simulated realizations
  stats_sim <- series.sim %>%
    as_tibble(.name_repair = ~as.character(1:ncol(series.sim))) %>%
    mutate(yind = 1:n()) %>%
    gather(key = sim, value = value, -yind) %>%
    mutate(sim = as.numeric(sim)) %>%
    group_by(sim) %>%
    summarize(mean = mean(value), sdev = sd(value),
              max = max(value), min = min(value))

  # Statistics for observed weather series
  stats_obs <- tibble(value = series.obs) %>%
    summarize(mean = mean(value), sdev = sd(value),
              max = max(value), min = min(value))

  # Significant periods
  periods_sig  <- which(power.obs > power.signif)

  if (!is.null(power.bounds)) {

    # Filter scenarios have significant signals
    sub_power1 <- which(sapply(1:ncol(power.sim), function(x)
      any(power.sim[periods_sig,x] > power.signif[periods_sig])))

    sub_power2 <- which(sapply(1:ncol(power.sim), function(x)
        all((power.sim[periods_sig,x] > power.obs[periods_sig] * power.bounds[1]) &
            (power.sim[periods_sig,x] < power.obs[periods_sig] * power.bounds[2]))))

    sub_power3 <- which(sapply(1:ncol(power.sim), function(x)
        all((power.sim[-periods_sig,x] < power.signif[-periods_sig]*nonsig.threshold))))

    sub_power <- intersect(intersect(sub_power1, sub_power2),sub_power3)
  } else {
    sub_power  <- 1:ncol(series.sim)
  }


  if (!is.null(mean.bounds)) {
    sub_mean  <- which((stats_sim$mean > stats_obs$mean * mean.bounds[1]) &
                       (stats_sim$mean < stats_obs$mean * mean.bounds[2]))
  } else {
    sub_mean <- 1:ncol(series.sim)
  }


  if (!is.null(sdev.bounds)) {
    sub_sdev  <- which((stats_sim$sdev > stats_obs$sdev * sdev.bounds[1]) &
                       (stats_sim$sdev < stats_obs$sdev * sdev.bounds[2]))
  } else {
    sub_sdev <- 1:ncol(series.sim)
  }


  if (!is.null(min.bounds)) {
    sub_min  <- which((stats_sim$min > stats_obs$min * min.bounds[1]) &
                      (stats_sim$min < stats_obs$min * min.bounds[2]))
  } else {
    sub_min  <- 1:ncol(series.sim)
  }

  if (!is.null(max.bounds)) {
    sub_max  <- which((stats_sim$max > stats_obs$max * max.bounds[1]) &
                      (stats_sim$max < stats_obs$max * max.bounds[2]))
  } else {
    sub_max  <- 1:ncol(series.sim)
  }

  #Select intersection
  sub_clim <- Reduce(intersect, list(sub_mean, sub_sdev, sub_power,
                                     sub_min, sub_max))

  # Stochastically select from the initial dataset
  if(!is.null(seed)) set.seed(seed)

  sub_sample <- sample(sub_clim, min(nmax, length(sub_clim)))

  if(isTRUE(verbose)) {
    print(tribble(
        ~"criteria", ~"# traces",
        "mean",  paste0(length(sub_mean), " out of ", ncol(series.sim)),
        "stdev", paste0(length(sub_sdev), " out of ", ncol(series.sim)),
        "power", paste0(length(sub_power), " out of ", ncol(series.sim)),
        "min",   paste0(length(sub_min), " out of ", ncol(series.sim)),
        "max",   paste0(length(sub_max), " out of ", ncol(series.sim)),
        "final", paste0(length(sub_clim), " out of ", ncol(series.sim))))
  }

  if(length(sub_sample) < nmax) {
    stop('subsetted traces less than the desired amount.
         Please relax the decision criteria and repeat.')
  }

  if(isTRUE(save.plots)) {

    ### Global Wavelet Spectral Plot
    pl <- length(power.period)
    plr <- 5 * ceiling(power.period[pl] / 5)

    p <- waveletPlot(power.period = power.period[1:pl],
                   power.signif = power.signif[1:pl],
                   power.obs = power.obs[1:pl],
                   power.sim = power.sim[1:pl,sub_clim])  +
          scale_x_continuous(breaks=seq(5,plr,5), limits=c(0,plr), expand=c(0,0))

    ggsave(paste0(out.path, "warm_subset_wavelet_spectra.png"), width=8, height=6)

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
      geom_jitter(data = stats_sim_gg, alpha = 0.1, color = "gray60") +
      geom_point(aes(fill = type), data = stats_all, size = 5, color = "white",
                 shape = 21) +
      scale_fill_manual(values = c("Simulated"="black","Observed"="blue")) +
      labs(x="", y = "", color = "", fill="") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())

      ggsave(paste0(out.path, "warm_subset_annual_stats.png"), width=8, height=6)

  }

  if(isTRUE(save.series)) {

    utils::write.csv(x = series.sim[,sub_clim], row.names = FALSE,
        file = paste0(out.path, "warm_sim_annual_set.csv"))

    utils::write.csv(x = series.sim[,sub_sample], row.names = FALSE,
        file = paste0(out.path, "warm_sim_annual_sample.csv"))

    }

  return(list(subsetted = series.sim[,sub_clim], sampled = series.sim[,sub_sample]))

}

