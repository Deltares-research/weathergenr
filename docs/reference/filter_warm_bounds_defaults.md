# WARM filtering default bounds

Returns internal default bounds for
[`filter_warm_pool()`](https://deltares-research.github.io/weathergenr/reference/filter_warm_pool.md).

The defaults are designed to be moderately selective for annual records
on the order of 50 to 100 years. Typical usage is to override only a
small subset of entries via `filter_bounds = list(...)` while keeping
the remaining defaults unchanged.

## Usage

``` r
filter_warm_bounds_defaults()
```

## Value

Named list of default bounds using snake case keys.

## Details

The returned list contains thresholds and controls used by the filtering
and relaxation logic.

Distributional tolerances:

- mean:

  Maximum absolute relative difference in the mean between a simulated
  realization and the observed series. Default 0.03.

- sd:

  Maximum absolute relative difference in standard deviation between a
  simulated realization and the observed series. Default 0.03.

Tail mass behavior:

- tail_low_p:

  Lower quantile used to define the low tail threshold in the observed
  series. Default 0.20.

- tail_high_p:

  Upper quantile used to define the high tail threshold in the observed
  series. Default 0.80.

- tail_tol_log:

  Maximum absolute log difference between simulated and observed tail
  mass metrics. Default log(1.03).

- tail_eps:

  Positive constant used for numerical stability in log transforms.
  Default 1e-5.

Spectral similarity:

- spectral_cor_min:

  Minimum correlation between log transformed observed and simulated
  global wavelet spectra. Default 0.60.

- spectral_eps:

  Positive constant used for numerical stability in spectral log
  transforms. Default 1e-10.

Peak matching:

- n_sig_peaks_max:

  Maximum number of significant observed spectral peaks to enforce.
  Default 2.

- peak_period_tol:

  Tolerance for matching peak periods in log2 period space. Default
  0.50.

- peak_mag_tol_log:

  Tolerance for matching peak magnitudes as the absolute log ratio
  between simulated and observed peak power. Default log(1.5).

- peak_match_frac_min:

  Minimum fraction of significant observed peaks that must be matched.
  Default 1.0.

Plot controls:

- plot_wavelet_q:

  Two probabilities used to summarize simulated spectra in diagnostic
  plots. Default c(0.50, 0.95).

Relaxation controls:

- relax_mult:

  Multiplicative factor applied when relaxing some bounds. Default 1.25.

- relax_mean_max:

  Maximum allowed mean tolerance during relaxation. Default 0.25.

- relax_sd_max:

  Maximum allowed sd tolerance during relaxation. Default 0.25.

- relax_tail_tol_log_max:

  Maximum tail tolerance during relaxation. Default log(2.0).

- relax_tail_p_step:

  Step size for relaxing tail quantiles. Default 0.02.

- relax_tail_p_low_max:

  Maximum lower tail quantile during relaxation. Default 0.40.

- relax_tail_p_high_min:

  Minimum upper tail quantile during relaxation. Default 0.40.

- relax_spectral_cor_step:

  Decrement applied to spectral_cor_min during wavelet relaxation.
  Default 0.05.

- relax_spectral_cor_min:

  Minimum spectral_cor_min allowed during relaxation. Default 0.30.

- relax_peak_match_frac_step:

  Decrement applied to peak_match_frac_min during wavelet relaxation.
  Default 0.10.

- relax_peak_match_frac_min:

  Minimum peak_match_frac_min allowed during relaxation. Default 0.00.

- relax_max_iter:

  Maximum number of relaxation iterations. Default 20.
