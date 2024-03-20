
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
    sd = c(0.80,1.20),
    min = c(0.80,1.20),
    max = c(0.80,1.20),
    power = c(0.40,3.00),
    nonsignif.threshold = 0.60))

{

  # Workaround for rlang warning
  sim <- value <- yind <- par <- type <- variable <- y <- x <- 0
  sim.year.num <- nrow(series.sim)

  # If no seed provided, sample a value
  if(is.null(seed)) seed <- sample.int(1e5,1)

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

  filterMatchingTS <- function(bounds = NA, stats_sim = NA, stats_obs = NA, series.sim = NA) {

    # Filter based on power spectra
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

    } else {sub_power  <- 1:ncol(series.sim)}

    # Filter based on means
    if (!is.null(bounds$mean)) {
      sub_mean  <- which((stats_sim$mean > stats_obs$mean * bounds$mean[1]) &
                         (stats_sim$mean < stats_obs$mean * bounds$mean[2]))
    } else {sub_mean <- 1:ncol(series.sim)}

    # Filter based on standard deviation
    if (!is.null(bounds$sd)) {
      sub_sd  <- which((stats_sim$sd > stats_obs$sd * bounds$sd[1]) &
                         (stats_sim$sd < stats_obs$sd * bounds$sd[2]))
    } else {sub_sd <- 1:ncol(series.sim)}

    # Filter based on minimum
    if (!is.null(bounds$min)) {
      sub_min  <- which((stats_sim$min > stats_obs$min * bounds$min[1]) &
                          (stats_sim$min < stats_obs$min * bounds$min[2]))
    } else {sub_min  <- 1:ncol(series.sim)}

    # Filter based on maximum
    if (!is.null(bounds$max)) {
      sub_max  <- which((stats_sim$max > stats_obs$max * bounds$max[1]) &
                          (stats_sim$max < stats_obs$max * bounds$max[2]))
    } else {sub_max  <- 1:ncol(series.sim)}

    #Select intersection of filtered
    output <- Reduce(base::intersect,
        list(sub_mean, sub_sd, sub_power, sub_min, sub_max))

    return(output)
  }

  ### Filter the values
  sub_clim <- filterMatchingTS(bounds = bounds, stats_sim = stats_sim,
                           stats_obs = stats_obs, series.sim = series.sim)

  set.seed(seed)
  sub_sample <- sample(sub_clim, min(sample.num, length(sub_clim)))

  if(length(sub_sample) < sample.num) {
    message('Not enough traces meeting criteria. Relaxing criteria and repeating subsetting')

    ### Filter the values
    bounds_rev = list(mean = c(0.90,1.10), sd = c(0.80,1.20), min = c(0.70,1.30),
                      max = c(0.70,1.30), power = NULL, nonsignif.threshold = NULL)

    sub_clim <- filterMatchingTS(bounds = bounds_rev, stats_sim = stats_sim,
                                 stats_obs = stats_obs, series.sim = series.sim)

    set.seed(seed)
    sub_sample <- sample(sub_clim, min(sample.num, length(sub_clim)))

  }

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

    ggsave(file.path(output.path, "warm_spectral_matching.png"), width=8, height=6)

    # For subsetted realizations only
    p <- waveletPlot(power.period = power.period[1:pl],
                   power.signif = power.signif[1:pl],
                   power.obs = power.obs[1:pl],
                   power.sim = power.sim[1:pl,sub_sample, drop = FALSE])  +
          scale_x_continuous(breaks=seq(5,plr,5), limits=c(0,plr), expand=c(0,0))

    ggsave(file.path(output.path, "warm_spectral_sampled.png"), width=8, height=6)

    # Boxplots of all stats
    stats_obs_gg <- stats_obs %>% mutate(sim=1) %>%
      gather(key = par, value = value, -sim) %>% mutate(type = "Observed") %>%
      mutate(type = factor(type, levels = c("Sampled", "Observed"))) %>%
      mutate(par = factor(par, levels = c("mean","sd", "min","max"),
        labels = c("Mean", "St. Deviation", "Minimum", "Maximum"))) %>%
      arrange(type)

    stats_sim_gg <- stats_sim %>%
      gather(key = par, value = value, -sim) %>% mutate(type = "Sampled") %>%
      mutate(type = factor(type, levels = c("Sampled", "Observed"))) %>%
      mutate(par = factor(par, levels = c("mean","sd", "min","max"),
        labels = c("Mean", "St. Deviation", "Minimum", "Maximum"))) %>%
      arrange(type)

    # Plot subsetted series statistics
    p <- ggplot(mapping = aes(x = par, y = value)) +
      theme_bw() +
      facet_wrap(~par, scales = "free", drop = TRUE, nrow = 1) +
      geom_violin(data = stats_sim_gg, color = "gray60") +
      geom_point(data = filter(stats_sim_gg, sim %in% sub_sample),
                 size = 3, color = "white", fill = "black", shape = 21) +
      geom_point(data = stats_obs_gg, size = 4, color = "white", fill = "blue", shape = 21) +
      labs(x="", y = "", color = "", fill="") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())

      ggsave(file.path(output.path, "warm_stats_sampled.png"), width = 8, height = 5)

    # Plot simulated warm-series
    sub_clim_plot <- sub_clim[1:min(length(sub_clim), 50)]

    df1 <- series.sim[,sub_sample] %>%
      as_tibble(.name_repair = ~paste0("rlz",1:length(sub_sample))) %>%
      mutate(x = 1:sim.year.num) %>%
      gather(key = variable, value=y, -x) %>%
      mutate(y = y * 365)

    df2 <- tibble(x=1:length(series.obs), y = series.obs*365)

    # Subsetted annual series
    p <- ggplot(df1, aes(x = x, y = y)) +
      theme_bw(base_size = 12) +
      geom_line(aes(y = y, group = variable, color = variable), alpha = 0.6) +
      geom_line(aes(y=y), data = df2, color = "black", linewidth = 1) +
      guides(color = "none") +
      labs(y = "Precipitation (mm/year)", x = "Year index")

    ggsave(file.path(output.path, "warm_annual_series.png"), height = 5, width = 8)


  }

  if(isTRUE(save.series)) {

    utils::write.csv(x = series.sim[,sub_clim], row.names = FALSE,
        file = file.path(output.path, "warm_output_matching.csv"))

    utils::write.csv(x = series.sim[,sub_sample], row.names = FALSE,
        file = file.path(output.path, "warm_output_sampled.csv"))

    }

  return(list(subsetted = series.sim[,sub_clim],
              sampled = series.sim[,sub_sample, drop=FALSE]))

}

