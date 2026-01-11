# Compute pass vectors for all filters

Evaluates which realizations pass each filter criterion.

## Usage

``` r
compute_pass_vectors(
  rel.diff.mean,
  rel.diff.sd,
  tail_metrics,
  P_sim_reg,
  P_sim_bg,
  presence_rpad,
  wavelet_diag,
  bounds,
  wavelet_active,
  wavelet_pars,
  n_realizations
)
```

## Arguments

- rel.diff.mean:

  Relative differences in mean

- rel.diff.sd:

  Relative differences in SD

- tail_metrics:

  List from compute_tailmass_metrics()

- P_sim_reg:

  Matrix of regional powers (n_realizations x n_regions)

- P_sim_bg:

  Vector of background powers (n_realizations)

- presence_rpad:

  Logical vector of presence indicators

- wavelet_diag:

  List of wavelet diagnostics

- bounds:

  Bounds environment or list

- wavelet_active:

  Logical for wavelet filter status

- wavelet_pars:

  List of wavelet parameters

- n_realizations:

  Number of realizations

## Value

Named list of logical vectors (pass\$mean, pass\$sd, pass\$tail_low,
etc.)
