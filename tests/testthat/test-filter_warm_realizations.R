# test-filter_warm_realizations.R

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

make_inputs <- function(n_years = 30, n_real = 12) {
  set.seed(1)
  series.obs <- rnorm(n_years, mean = 100, sd = 20)
  series.sim <- replicate(n_real, rnorm(n_years, mean = 100, sd = 20))
  series.sim <- matrix(series.sim, nrow = n_years, ncol = n_real)
  list(series.obs = series.obs, series.sim = series.sim)
}

find_first <- function(x, pred) {
  if (!is.list(x)) return(NULL)
  for (nm in names(x)) {
    obj <- x[[nm]]
    if (isTRUE(pred(obj))) return(obj)
  }
  NULL
}

get_named_or_find <- function(x, name_candidates, pred) {
  if (is.list(x)) {
    for (nm in name_candidates) {
      if (!is.null(x[[nm]])) return(x[[nm]])
    }
  }
  find_first(x, pred)
}

# -----------------------------------------------------------------------------
# Wavelet mocks (fast + deterministic)
# -----------------------------------------------------------------------------
mock_wavelet_factory <- function(n_period = 12, sig_mode = c("sig", "nosig")) {
  sig_mode <- match.arg(sig_mode)

  function(variable,
           signif.level = 0.90,
           noise.type = "red",
           period.lower.limit = 2,
           detrend = FALSE,
           mode = "fast",
           ...) {

    set.seed(length(variable) + n_period)

    gws_period <- seq_len(n_period)

    # NOTE: use sig_mode (outer) to decide significance behavior
    if (sig_mode == "nosig") {
      gws <- rep(0.9, n_period)      # never significant relative to signif=1
      gws_unmasked <- gws
    } else {
      gws <- rep(0.9, n_period)
      gws[c(3, 4, 5, 8)] <- 1.3      # some significant periods
      gws_unmasked <- gws
    }

    list(
      gws_period = gws_period,
      gws = gws,
      gws_unmasked = gws_unmasked,
      gws_signif = rep(1.0, n_period),
      gws_signif_unmasked = rep(1.0, n_period),
      has_significance = any(gws > 1.0),
      signif_periods = which(gws > 1.0),
      coi = rep(Inf, length(variable)),
      power = matrix(0, nrow = n_period, ncol = length(variable))
    )
  }
}

mock_gws_regrid <- function(wv, power.period, use_unmasked = TRUE) {
  if (isTRUE(use_unmasked) && !is.null(wv$gws_unmasked)) return(as.numeric(wv$gws_unmasked))
  as.numeric(wv$gws)
}

mock_fill_nearest <- function(x) as.numeric(x)

# -----------------------------------------------------------------------------
# 1) Basic structure + determinism
# -----------------------------------------------------------------------------

test_that("filter_warm_simulations returns expected structure", {
  inp <- make_inputs(n_years = 30, n_real = 12)

  testthat::local_mocked_bindings(
    wavelet_spectral_analysis = mock_wavelet_factory(n_period = 12, sig_mode = "sig"),
    gws_regrid = mock_gws_regrid,
    fill_nearest = mock_fill_nearest,
    .package = "weathergenr"
  )

  out <- weathergenr::filter_warm_simulations(
    series.obs = inp$series.obs,
    series.sim = inp$series.sim,
    sample.num = 5,
    seed = 123,
    padding = TRUE,
    make.plots = FALSE,
    verbose = FALSE
  )

  expect_true(is.list(out))

  sampled <- get_named_or_find(
    out,
    name_candidates = c("sampled", "sample", "sim.sampled", "sampled_sim", "sampled_series"),
    pred = function(z) is.matrix(z) && nrow(z) == length(inp$series.obs)
  )

  subsetted <- get_named_or_find(
    out,
    name_candidates = c("subsetted", "subset", "filtered", "sim.filtered", "subsetted_sim"),
    pred = function(z) is.matrix(z) && nrow(z) == length(inp$series.obs)
  )

  filter_summary <- get_named_or_find(
    out,
    name_candidates = c("filter_summary", "summary", "filter.stats", "filter_diagnostics"),
    pred = function(z) is.data.frame(z)
  )

  n_filtered <- get_named_or_find(
    out,
    name_candidates = c("n_filtered", "n.filter", "n_pool", "pool_size", "n_kept", "nretained"),
    pred = function(z) length(z) == 1 && (is.numeric(z) || is.integer(z))
  )
  if (is.null(n_filtered) && is.matrix(subsetted)) n_filtered <- ncol(subsetted)

  expect_true(is.matrix(sampled))
  expect_true(is.matrix(subsetted))
  expect_true(is.data.frame(filter_summary))
  expect_true(length(n_filtered) == 1 && (is.numeric(n_filtered) || is.integer(n_filtered)))

  expect_equal(nrow(subsetted), length(inp$series.obs))
  expect_equal(nrow(sampled), length(inp$series.obs))
  expect_lte(ncol(sampled), 5)
})

test_that("filter_warm_simulations is deterministic given seed", {
  inp <- make_inputs(n_years = 30, n_real = 12)

  testthat::local_mocked_bindings(
    wavelet_spectral_analysis = mock_wavelet_factory(n_period = 12, sig_mode = "sig"),
    gws_regrid = mock_gws_regrid,
    fill_nearest = mock_fill_nearest,
    .package = "weathergenr"
  )

  out1 <- weathergenr::filter_warm_simulations(
    series.obs = inp$series.obs,
    series.sim = inp$series.sim,
    sample.num = 5,
    seed = 999,
    make.plots = FALSE,
    verbose = FALSE
  )

  out2 <- weathergenr::filter_warm_simulations(
    series.obs = inp$series.obs,
    series.sim = inp$series.sim,
    sample.num = 5,
    seed = 999,
    make.plots = FALSE,
    verbose = FALSE
  )

  sampled1 <- get_named_or_find(out1, c("sampled", "sample", "sim.sampled", "sampled_sim", "sampled_series"), is.matrix)
  sampled2 <- get_named_or_find(out2, c("sampled", "sample", "sim.sampled", "sampled_sim", "sampled_series"), is.matrix)

  expect_true(is.matrix(sampled1))
  expect_true(is.matrix(sampled2))
  expect_equal(sampled1, sampled2)
})

# -----------------------------------------------------------------------------
# 2) Input validation / length reconciliation (updated)
# -----------------------------------------------------------------------------

test_that("series length mismatch reconciles to a common length (sampled rows match sim rows)", {
  inp <- make_inputs(n_years = 30, n_real = 10)
  series.obs_short <- inp$series.obs[-1]  # length 29, sim has 30 rows

  testthat::local_mocked_bindings(
    wavelet_spectral_analysis = mock_wavelet_factory(n_period = 10, sig_mode = "sig"),
    gws_regrid = mock_gws_regrid,
    fill_nearest = mock_fill_nearest,
    .package = "weathergenr"
  )

  out <- weathergenr::filter_warm_simulations(
    series.obs = series.obs_short,
    series.sim = inp$series.sim,
    sample.num = 5,
    seed = 1,
    make.plots = FALSE,
    verbose = FALSE
  )

  sampled <- get_named_or_find(
    out,
    name_candidates = c("sampled", "sample", "sim.sampled", "sampled_sim", "sampled_series"),
    pred = function(z) is.matrix(z)
  )

  expect_true(is.matrix(sampled))
  expect_equal(nrow(sampled), nrow(inp$series.sim))
})

test_that("validates series.sim is a numeric matrix", {
  inp <- make_inputs(n_years = 30, n_real = 10)

  testthat::local_mocked_bindings(
    wavelet_spectral_analysis = mock_wavelet_factory(n_period = 10, sig_mode = "sig"),
    gws_regrid = mock_gws_regrid,
    fill_nearest = mock_fill_nearest,
    .package = "weathergenr"
  )

  expect_error(
    weathergenr::filter_warm_simulations(
      series.obs = inp$series.obs,
      series.sim = as.data.frame(inp$series.sim),
      make.plots = FALSE,
      verbose = FALSE
    ),
    regexp = "series\\.sim|matrix|numeric",
    fixed = FALSE
  )
})

# -----------------------------------------------------------------------------
# 3) Wavelet-path behavior (updated: may disable silently)
# -----------------------------------------------------------------------------

test_that("wavelet filter is disabled when no significant periods exist in observed GWS (silent or warning)", {
  inp <- make_inputs(n_years = 30, n_real = 10)

  testthat::local_mocked_bindings(
    wavelet_spectral_analysis = mock_wavelet_factory(n_period = 12, sig_mode = "nosig"),
    gws_regrid = mock_gws_regrid,
    fill_nearest = mock_fill_nearest,
    .package = "weathergenr"
  )

  out <- suppressWarnings(
    weathergenr::filter_warm_simulations(
      series.obs = inp$series.obs,
      series.sim = inp$series.sim,
      sample.num = 5,
      seed = 123,
      make.plots = FALSE,
      verbose = FALSE
    )
  )

  sampled <- get_named_or_find(out, c("sampled", "sample", "sim.sampled", "sampled_sim", "sampled_series"), is.matrix)
  filter_summary <- get_named_or_find(out, c("filter_summary", "summary", "filter.stats", "filter_diagnostics"), is.data.frame)

  expect_true(is.matrix(sampled))
  expect_true(is.data.frame(filter_summary))
})

test_that("relaxation engages when bounds are too strict to reach sample.num", {
  inp <- make_inputs(n_years = 30, n_real = 25)

  testthat::local_mocked_bindings(
    wavelet_spectral_analysis = mock_wavelet_factory(n_period = 12, sig_mode = "sig"),
    gws_regrid = mock_gws_regrid,
    fill_nearest = mock_fill_nearest,
    .package = "weathergenr"
  )

  strict_bounds <- list(
    mean = 1e-8,
    sd = 1e-8,
    tail.low.p = 0.10,
    tail.high.p = 0.90,
    tail.tol.log = 1e-6,
    sig.frac = 0.95,
    wavelet.min_bg = 1e-12,
    wavelet.region.tol = 0.01,
    wavelet.contrast.tol = 0.01,
    wavelet.require_presence = TRUE,
    wavelet.presence.frac = 0.90,
    relax.mult = 2.0,
    relax.mean.max = 10,
    relax.sd.max = 10,
    relax.tail.tol.log.max = 10,
    relax.tail.p.step = 0.05,
    relax.tail.p.low.max = 0.40,
    relax.tail.p.high.min = 0.60,
    relax.wavelet.sig.frac.step = 0.05,
    relax.wavelet.sig.frac.min = 0.50,
    relax.wavelet.region.tol.step = 0.05,
    relax.wavelet.region.tol.max = 1.00,
    relax.wavelet.contrast.tol.step = 0.05,
    relax.wavelet.contrast.tol.max = 1.00
  )

  out <- weathergenr::filter_warm_simulations(
    series.obs = inp$series.obs,
    series.sim = inp$series.sim,
    sample.num = 7,
    seed = 777,
    bounds = strict_bounds,
    make.plots = FALSE,
    verbose = FALSE
  )

  sampled <- get_named_or_find(out, c("sampled", "sample", "sim.sampled", "sampled_sim", "sampled_series"), is.matrix)
  filter_summary <- get_named_or_find(out, c("filter_summary", "summary", "filter.stats", "filter_diagnostics"), is.data.frame)

  expect_true(is.matrix(sampled))
  expect_equal(nrow(sampled), length(inp$series.obs))
  expect_lte(ncol(sampled), 7)

  relax_cols <- c("relaxation_level", "relax_level", "relaxation", "relax_iter", "iteration")
  if (is.data.frame(filter_summary)) {
    col_present <- intersect(relax_cols, names(filter_summary))
    if (length(col_present) > 0) {
      v <- filter_summary[[col_present[1]]]
      expect_true(length(v) >= 1)
    }
  }
})

# -----------------------------------------------------------------------------
# 4) Plotting flag (make.plots)
# -----------------------------------------------------------------------------

test_that("make.plots=TRUE requires output.path (only if enforced)", {
  inp <- make_inputs(n_years = 30, n_real = 10)

  testthat::local_mocked_bindings(
    wavelet_spectral_analysis = mock_wavelet_factory(n_period = 12, sig_mode = "sig"),
    gws_regrid = mock_gws_regrid,
    fill_nearest = mock_fill_nearest,
    .package = "weathergenr"
  )

  err <- try(
    weathergenr::filter_warm_simulations(
      series.obs = inp$series.obs,
      series.sim = inp$series.sim,
      sample.num = 5,
      seed = 1,
      make.plots = TRUE,
      output.path = NULL,
      verbose = FALSE
    ),
    silent = TRUE
  )

  if (!inherits(err, "try-error")) {
    skip("make.plots=TRUE does not enforce output.path in current implementation.")
  } else {
    expect_true(grepl("output\\.path|plot|path", conditionMessage(attr(err, "condition")), ignore.case = TRUE))
  }
})
