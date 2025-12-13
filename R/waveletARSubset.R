#' Subset and Sample Synthetic Time Series Based on Climate and Wavelet Criteria
#'
#' Selects and samples synthetic (simulated) time series that match statistical
#' and spectral (wavelet) criteria relative to observed data. Optionally saves
#' plots and selected series to disk.
#'
#' @param series.obs Numeric vector. Observed annual time series (length = years).
#' @param series.sim Numeric matrix. Simulated time series (nrow = years, ncol = realizations).
#' @param power.obs Numeric vector. Observed global wavelet power spectrum.
#' @param power.sim Numeric matrix. Simulated global wavelet spectra (nrow = periods, ncol = realizations).
#' @param power.period Numeric vector. Periods (years) for power spectra.
#' @param power.signif Numeric vector. Significance thresholds for power spectra.
#' @param sample.num Integer. Number of simulations to sample for final output.
#' @param seed Integer or NULL. Seed for reproducible sampling.
#' @param save.plots Logical. Whether to save plots (default: TRUE).
#' @param save.series Logical. Whether to write selected series to .csv (default: TRUE).
#' @param output.path String. Directory to write plots and .csv files to.
#' @param padding Logical. Whether to pad significant periods (default: TRUE).
#' @param bounds List. Bounds for filtering on mean, sd, min, max, sig.thr, nsig.thr.
#'
#' @return A list with two elements:
#'   \item{subsetted}{A numeric matrix containing all synthetic traces that pass the filtering criteria.}
#'   \item{sampled}{A numeric matrix of the sampled synthetic traces (number set by `sample.num`).}
#' @export
#' @import dplyr
#' @import ggplot2
waveletARSubset <- function(
    series.obs = NULL,
    series.sim = NULL,
    power.obs = NULL,
    power.sim = NULL,
    power.period = NULL,
    power.signif = NULL,
    sample.num = 5,
    seed = NULL,
    save.plots = TRUE,
    save.series = TRUE,
    output.path = NULL,
    padding = TRUE,
    bounds = list(
      mean = 0.1, sd = 0.2, min = 0.2, max = 0.2,
      sig.thr = 0.5, nsig.thr = 1.5
    )) {


  # Input checks
  if (is.null(series.obs) || !is.numeric(series.obs)) stop("'series.obs' must be a numeric vector.")
  if (is.null(series.sim) || !is.numeric(series.sim)) stop("'series.sim' must be a numeric matrix.")
  if (is.null(power.obs) || !is.numeric(power.obs)) stop("'power.obs' must be a numeric vector.")
  if (is.null(power.sim) || !is.numeric(power.sim)) stop("'power.sim' must be a numeric matrix.")
  if (is.null(power.period) || !is.numeric(power.period)) stop("'power.period' must be a numeric vector.")
  if (is.null(power.signif) || !is.numeric(power.signif)) stop("'power.signif' must be a numeric vector.")
  if (!is.numeric(sample.num) || sample.num < 1) stop("'sample.num' must be a positive integer.")

  # Workaround for rlang warning
  sim <- value <- yind <- par <- type <- variable <- y <- x <- obs <- 0

  sim.year.num <- nrow(series.sim)

  # Seed handling (local only)
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit(
      {
        .Random.seed <<- old_seed
      },
      add = TRUE
    )
    set.seed(seed)
  }

  # Statistics for observed series
  stats_obs <- tibble(value = series.obs) %>%
    summarize(
      mean = mean(value), sd = stats::sd(value),
      max = max(value), min = min(value)
    ) %>%
    gather(key = par, value = obs)

  # Statistics for synthetic series
  stats_sim <- tibble(sim = 1:ncol(series.sim)) %>%
    mutate(
      mean = apply(series.sim, 2, mean),
      sd = apply(series.sim, 2, sd),
      min = apply(series.sim, 2, min),
      max = apply(series.sim, 2, max)
    ) %>%
    group_by(sim) %>%
    gather(key = par, value = value, -sim) %>%
    left_join(stats_obs, by = "par") %>%
    mutate(value = 1 - (obs - value) / obs) %>%
    select(-obs) %>%
    pivot_wider(names_from = par, values_from = value)

  # Set significant and non-significant periods
  periods_sig <- which(power.obs > power.signif)
  periods_sig <- periods_sig[periods_sig %in% 1:length(series.obs)]

  if (isTRUE(padding)) {
    periods_sig <- sort(intersect(unique(c(periods_sig - 1, periods_sig, periods_sig + 1)), 1:length(power.signif)))
  }
  periods_nonsig <- setdiff(1:length(power.signif), periods_sig)
  periods_nonsig <- periods_nonsig[periods_nonsig %in% 1:dim(power.sim)[1]]

  ###################################################
  power_signif_max <- 10

  # Filter based on power spectra
  if (!is.null(bounds$sig.thr)) {
    # Filter scenarios have significant signals
    sub_power1 <- which(sapply(1:ncol(power.sim), function(x) {
      any(power.sim[periods_sig, x] > power.signif[periods_sig])
    }))

    # Significant signals within the bound
    sub_power2 <- which(sapply(1:ncol(power.sim), function(x) {
      all((power.sim[periods_sig, x] > power.obs[periods_sig] * bounds$sig.thr) &
        (power.sim[periods_sig, x] < power.obs[periods_sig] * power_signif_max))
    }))

    # Non-significant below threshold
    sub_power3 <- which(sapply(1:ncol(power.sim), function(x) {
      all((power.sim[periods_nonsig, x] < power.signif[periods_nonsig] * bounds$nsig.thr))
    }))

    sub_power <- base::intersect(base::intersect(sub_power1, sub_power2), sub_power3)
  } else {
    sub_power <- 1:ncol(series.sim)
  }

  # Filter based on means
  if (!is.null(bounds$mean)) {
    sub_mean <- which((stats_sim$mean > 1 - bounds$mean) & (stats_sim$mean < 1 + bounds$mean))
  } else {
    sub_mean <- 1:ncol(series.sim)
  }

  # Filter based on standard deviation
  if (!is.null(bounds$sd)) {
    sub_sd <- which((stats_sim$sd > 1 - bounds$sd) & (stats_sim$sd < 1 + bounds$sd))
  } else {
    sub_sd <- 1:ncol(series.sim)
  }

  # Filter based on minimum
  if (!is.null(bounds$min)) {
    sub_min <- which((stats_sim$min > 1 - bounds$min) & (stats_sim$min < 1 + bounds$min))
  } else {
    sub_min <- 1:ncol(series.sim)
  }

  # Filter based on maximum
  if (!is.null(bounds$max)) {
    sub_max <- which((stats_sim$max > 1 - bounds$max) & (stats_sim$max < 1 + bounds$max))
  } else {
    sub_max <- 1:ncol(series.sim)
  }

  # Select intersection of filtered
  sub_clim <- Reduce(
    base::intersect,
    list(sub_mean, sub_sd, sub_power, sub_min, sub_max)
  )

  set.seed(seed)
  sub_sample <- sample(sub_clim, min(sample.num, length(sub_clim)))

  if (length(sub_sample) < sample.num) {
    logger::log_info("\n[WARM] Not enough traces meeting criteria. Bypassing power constraint")

    # Select intersection of filtered
    sub_clim <- Reduce(base::intersect, list(sub_mean, sub_sd))

    set.seed(seed)
    sub_sample <- sample(sub_clim, min(sample.num, length(sub_clim)))
  }

  if (isTRUE(save.plots)) {

    suppressWarnings({

      ### Global Wavelet Spectral Plot
      pl <- min(length(power.period), dim(power.sim)[1])
      plr <- 5 * ceiling(power.period[pl] / 5)

      # For all matching realizations
      p <- waveletPlot(
        power.period = power.period[1:pl],
        power.signif = power.signif[1:pl],
        power.obs = power.obs[1:pl],
        power.sim = power.sim[1:pl, sub_clim]
      ) +
        scale_x_continuous(breaks = seq(5, plr, 5), limits = c(0, plr), expand = c(0, 0))

      ggsave(file.path(output.path, "warm_spectral_matching.png"), width = 8, height = 6)

      # For subsetted realizations only
      p <- waveletPlot(
        power.period = power.period[1:pl],
        power.signif = power.signif[1:pl],
        power.obs = power.obs[1:pl],
        power.sim = power.sim[1:pl, sub_sample, drop = FALSE]
      ) +
        scale_x_continuous(breaks = seq(5, plr, 5), limits = c(0, plr), expand = c(0, 0))

      ggsave(file.path(output.path, "warm_spectral_sampled.png"), width = 8, height = 6)

      # Boxplots of all stats
      par_labels <- c(`mean` = "Mean", `sd` = "StDev", `min` = "Minimum", `max` = "Maximum")

      stats_sim_gg <- stats_sim %>%
        pivot_longer(!sim, names_to = "par", values_to = "value") %>%
        mutate(par = factor(par, levels = c("mean", "sd", "min", "max"))) %>%
        mutate(value = value * 100 - 100)

      # Plot subsetted series statistics
      p <- ggplot(stats_sim_gg, aes(x = par, y = value)) +
        theme_bw() +
        geom_violin(color = "black", fill = "gray90") +
        facet_wrap(~par,
          scales = "free", drop = TRUE, nrow = 1,
          labeller = as_labeller(par_labels)
        ) +
        geom_hline(yintercept = 0, linewidth = 1, color = "blue") +
        geom_point(
          data = filter(stats_sim_gg, sim %in% sub_sample),
          size = 3, color = "white", fill = "black", shape = 21
        ) +
        scale_y_continuous(limits = c(-50, 50), breaks = seq(-50, 50, 25)) +
        labs(x = "", y = "Change (%)") +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )

      ggsave(file.path(output.path, "warm_stats_sampled.png"), width = 8, height = 6)

      # Plot simulated warm-series
      sub_clim_plot <- sub_clim[1:min(length(sub_clim), 50)]

      df1 <- series.sim[, sub_sample] %>%
        as_tibble(.name_repair = ~ paste0("rlz", 1:length(sub_sample))) %>%
        mutate(x = 1:sim.year.num) %>%
        gather(key = variable, value = y, -x) %>%
        mutate(y = y * 365)

      df2 <- tibble(x = 1:length(series.obs), y = series.obs * 365)

      # Subsetted annual series
      p <- ggplot(df1, aes(x = x, y = y)) +
        theme_bw(base_size = 12) +
        geom_line(aes(y = y, group = variable), color = "gray50", alpha = 0.5) +
        geom_line(aes(y = y), data = df2, color = "blue", linewidth = 1) +
        scale_x_continuous(expand = c(0, 0)) +
        guides(color = "none") +
        labs(y = "mm/year", x = "Serial year")

      ggsave(file.path(output.path, "warm_annual_series.png"), width = 8, height = 6)
    }) # close suppress warnings
  }

  if (isTRUE(save.series)) {
    utils::write.csv(
      x = series.sim[, sub_clim], row.names = FALSE,
      file = file.path(output.path, "warm_output_matching.csv")
    )

    utils::write.csv(
      x = series.sim[, sub_sample], row.names = FALSE,
      file = file.path(output.path, "warm_output_sampled.csv")
    )
  }

  return(list(
    subsetted = series.sim[, sub_clim],
    sampled = series.sim[, sub_sample, drop = FALSE]
  ))
}

