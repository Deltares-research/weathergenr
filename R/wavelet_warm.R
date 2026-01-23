# ==============================================================================
# Wavelet Autoregressive Modeling (WARM) Simulation
# ==============================================================================
#
# This script implements a fast ARMA based simulator for WARM workflows.
#
# Core choices
# - Model fitting uses stats::arima over a small grid of ARMA(p,q) candidates and
#   selects the minimum AIC model.
# - Simulation uses a custom ARMA recursion that generates an entire
#   n by n_sim matrix per fitted model. This avoids S3 dispatch and scales well
#   when n_sim is large.
# - When ARMA fitting is not viable for a series (for example due to numerical
#   issues), the code falls back to a block bootstrap.
#
# Notes
# - In component mode, each component is modeled independently then summed.
# - MODWT MRA components are additive but not orthogonal; independent modeling
#   does not preserve cross component covariance unless you implement coupling
#   externally.
#

#' Simulate synthetic series with Wavelet Autoregressive Modeling (WARM)
#'
#' @description
#' Generates n_sim synthetic realizations of length n.
#'
#' Two operating modes are supported.
#'
#' Component mode (n >= bypass_n)
#' - Each wavelet component is centered.
#' - A stationary non seasonal ARMA(p,q) model is fit once per component using a
#'   grid search over p in 0..max_p and q in 0..max_q with AIC selection.
#' - Simulations are generated from the selected ARMA model using a fixed
#'   simulator (no S3 simulate dispatch).
#' - Component means are re added and components are summed.
#' - If fitting fails for a component, a block bootstrap fallback is used for
#'   that component only.
#'
#' Bypass mode (n < bypass_n)
#' - Component simulation is skipped.
#' - A stationary non seasonal ARMA(p,q) model is fit to the original series and
#'   simulated directly.
#' - If fitting fails, a block bootstrap fallback is used.
#'
#' @param components Matrix, data.frame, or list of numeric vectors. Wavelet
#'   components. In component mode, each component must have length n. In bypass
#'   mode, components are only used to reconstruct the original series if
#'   series_obs is not provided.
#' @param n Integer. Length of each simulated series.
#' @param n_sim Integer. Number of realizations to generate. Default is 1000.
#' @param seed Optional integer. Base RNG seed for reproducibility.
#' @param series_obs Optional numeric vector. Observed original series. Provide
#'   this in bypass mode to avoid reconstructing the series from components.
#' @param bypass_n Integer. Threshold for bypass mode. Default is 25.
#' @param match_variance Logical. If TRUE, rescale simulated outputs to match
#'   observed standard deviation when relative mismatch exceeds var_tol.
#' @param var_tol Numeric in [0, 1]. Relative standard deviation tolerance.
#' @param check_diagnostics Logical. If TRUE, runs a Ljung Box test on ARMA
#'   residuals and warns when residual autocorrelation remains.
#' @param verbose Logical. If TRUE, prints informational messages.
#'
#' @return Numeric matrix with dimension n by n_sim. Each column is a simulated
#'   realization.
#'
#' @details
#' ARMA configuration
#' - Non seasonal: seasonal terms are not considered.
#' - Stationary: fits enforce stationarity constraints via parameter transforms.
#' - No drift and no mean: include.mean is FALSE. The code centers series before
#'   fitting and re adds the observed mean after simulation.
#'
#' Fallback logic
#' - Fallback is applied per component in component mode. Components that fit
#'   successfully use ARMA simulation; failed components use block bootstrap.
#' - In bypass mode, fallback is applied to the whole series.
#'
#' @importFrom stats sd Box.test
#' @export
simulate_warm <- function(
    components = NULL,
    n = NULL,
    n_sim = 1000,
    seed = NULL,
    series_obs = NULL,
    bypass_n = 15L,
    match_variance = TRUE,
    var_tol = 0.1,
    check_diagnostics = FALSE,
    verbose = TRUE
) {

  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  if (is.null(n)) stop("Input 'n' must be specified.", call. = FALSE)
  if (!.is_int_scalar(n) || n < 1L) stop("'n' must be a positive integer.", call. = FALSE)
  n <- as.integer(n)

  if (!.is_int_scalar(n_sim) || n_sim < 1L) stop("'n_sim' must be a positive integer.", call. = FALSE)
  n_sim <- as.integer(n_sim)

  if (!.is_int_scalar(bypass_n) || bypass_n < 5L) stop("'bypass_n' must be an integer >= 5.", call. = FALSE)
  bypass_n <- as.integer(bypass_n)

  if (!is.logical(verbose) || length(verbose) != 1L) stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  verbose <- isTRUE(verbose)

  if (!is.logical(match_variance) || length(match_variance) != 1L) {
    stop("'match_variance' must be TRUE or FALSE.", call. = FALSE)
  }
  match_variance <- isTRUE(match_variance)

  if (!is.numeric(var_tol) || length(var_tol) != 1L || !is.finite(var_tol) || var_tol < 0 || var_tol > 1) {
    stop("'var_tol' must be between 0 and 1.", call. = FALSE)
  }

  if (!is.logical(check_diagnostics) || length(check_diagnostics) != 1L) {
    stop("'check_diagnostics' must be TRUE or FALSE.", call. = FALSE)
  }
  check_diagnostics <- isTRUE(check_diagnostics)

  # --------------------------------------------------------------------------
  # RNG management
  # --------------------------------------------------------------------------
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be NULL or a single finite number.", call. = FALSE)
    }
    seed <- as.integer(seed)

    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- .Random.seed
      has_seed <- TRUE
    } else {
      has_seed <- FALSE
    }
    on.exit({
      if (has_seed) .Random.seed <<- old_seed
    }, add = TRUE)
  }

  # --------------------------------------------------------------------------
  # Determine observed (fit) length
  # --------------------------------------------------------------------------
  n_fit <- NA_integer_

  if (!is.null(series_obs)) {
    if (!is.numeric(series_obs) || anyNA(series_obs)) {
      stop("'series_obs' must be numeric with no missing values.", call. = FALSE)
    }
    n_fit <- length(series_obs)
  } else if (!is.null(components)) {
    if (is.matrix(components) || is.data.frame(components)) {
      n_fit <- nrow(components)
    } else if (is.list(components)) {
      n_fit <- length(components[[1L]])
    }
  }

  if (!is.finite(n_fit) || n_fit < 1L) {
    stop("Cannot determine observed series length. Provide 'series_obs' or 'components'.", call. = FALSE)
  }
  n_fit <- as.integer(n_fit)













  # --------------------------------------------------------------------------
  # BYPASS MODE
  # --------------------------------------------------------------------------
  if (n_fit < bypass_n) {

    if (!is.null(series_obs)) {
      if (!is.numeric(series_obs) || anyNA(series_obs)) {
        stop("'series_obs' must be numeric with no missing values.", call. = FALSE)
      }
      series_obs <- as.numeric(series_obs)
      if (length(series_obs) != n_fit) stop("'series_obs' length must equal 'n'.", call. = FALSE)

    } else {

      if (is.null(components)) {
        stop("For n < bypass_n, provide 'series_obs' or 'components'.", call. = FALSE)
      }
      series_obs <- .reconstruct_series_from_components(components, n = n_fit)
    }

    if (!is.null(seed)) set.seed(seed)

    obs_mean <- mean(series_obs)
    x <- series_obs - obs_mean

    fit <- .fit_warm_arima_forecast(
      x = x,
      max_p = 2L,
      max_q = 2L
    )

    if (is.null(fit)) {
      bl <- .default_block_len(n)
      if (verbose) .log(paste0("[WARM] Bypass mode: ARMA fit failed, using block bootstrap (block_len=", bl, ")"))
      out <- matrix(NA_real_, nrow = n, ncol = n_sim)
      for (j in seq_len(n_sim)) {
        out[, j] <- .block_bootstrap(series_obs, n = n, block_len = bl)
      }
      return(out)
    }

    if (check_diagnostics && !is.null(fit$residuals) && length(fit$residuals) > 1L) {
      lj <- tryCatch(
        stats::Box.test(fit$residuals, type = "Ljung-Box", lag = min(10L, length(fit$residuals) - 1L)),
        error = function(e) NULL
      )
      if (!is.null(lj) && is.finite(lj$p.value) && lj$p.value < 0.05) {
        warning(
          "Bypass ARMA: residuals show significant autocorrelation (Ljung-Box p = ",
          round(lj$p.value, 4), "). Model may be inadequate.",
          call. = FALSE
        )
      }
    }

    sim_centered <- .warm_simulate_from_fit(fit, n = n, n_sim = n_sim)

    sim <- sim_centered + obs_mean

    if (match_variance) {
      target_sd <- stats::sd(series_obs)
      sim <- .variance_match_matrix(sim, target_sd = target_sd, tol = var_tol, target_mean = obs_mean)
    }

    if (verbose) .log(paste0("[WARM] Bypass mode used (n=", n, " < ", bypass_n, ")"))
    return(sim)
  }

  # --------------------------------------------------------------------------
  # COMPONENT MODE
  # --------------------------------------------------------------------------
  if (is.null(components)) stop("Input 'components' must not be NULL for component mode.", call. = FALSE)

  comp_names <- NULL

  if (is.matrix(components)) {
    if (!is.numeric(components)) stop("Matrix 'components' must be numeric.", call. = FALSE)
    if (nrow(components) != n_fit) stop("Matrix 'components' must have n rows.", call. = FALSE)
    ncomp <- ncol(components)
    comp_list <- vector("list", ncomp)
    for (k in seq_len(ncomp)) comp_list[[k]] <- components[, k]
    comp_names <- colnames(components)

  } else if (is.data.frame(components)) {
    col_types <- vapply(components, is.numeric, logical(1))
    if (!all(col_types)) stop("All columns of 'components' must be numeric.", call. = FALSE)
    if (nrow(components) != n_fit) stop("Data.frame 'components' must have n rows.", call. = FALSE)
    ncomp <- ncol(components)
    comp_list <- vector("list", ncomp)
    for (k in seq_len(ncomp)) comp_list[[k]] <- components[[k]]
    comp_names <- colnames(components)

  } else if (is.list(components)) {
    elem_types <- vapply(components, is.numeric, logical(1))
    if (!all(elem_types)) stop("All elements of 'components' must be numeric vectors.", call. = FALSE)
    ncomp <- length(components)
    comp_list <- components
    comp_names <- names(components)
    lens <- vapply(comp_list, length, integer(1))
    if (any(lens != n)) stop("All component vectors must have length 'n'.", call. = FALSE)

  } else {
    stop("'components' must be a matrix, data.frame, or list of numeric vectors.", call. = FALSE)
  }

  na_comp <- vapply(comp_list, function(x) anyNA(x), logical(1))
  if (any(na_comp)) {
    bad <- which(na_comp)
    stop(
      "Missing values detected in component(s): ",
      paste(bad, collapse = ", "),
      ". Remove or impute NAs before calling simulate_warm().",
      call. = FALSE
    )
  }

  output <- matrix(0, nrow = n, ncol = n_sim)
  variance_corrections <- FALSE

  for (k in seq_along(comp_list)) {

    component <- as.numeric(comp_list[[k]])
    comp_sd <- stats::sd(component)
    if (!is.finite(comp_sd)) stop("Component ", k, " has non-finite standard deviation.", call. = FALSE)

    nm <- ""
    if (!is.null(comp_names) && length(comp_names) >= k && !is.na(comp_names[k])) nm <- comp_names[k]
    is_smooth <- nzchar(nm) && grepl("^S", nm)

    # Constant component: carry through deterministically
    if (comp_sd < 1e-10) {
      output <- output + mean(component)
      next
    }

    if (!is.null(seed)) set.seed(seed + k * 1000L)

    comp_mean <- mean(component)
    centered <- component - comp_mean

    max_p_use <- if (is_smooth) 1L else 2L
    max_q_use <- if (is_smooth) 0L else 2L

    is_viable <- .warm_arima_viable(
      centered,
      min_n = max(8L, max_p_use + max_q_use + 2L),
      sd_eps = 1e-8,
      min_unique = 3L
    )

    fit <- NULL
    if (is_viable) {
      fit <- .warm_fit_arima_safe(
        x = centered,
        max_p = max_p_use,
        max_q = max_q_use,
        stationary = TRUE,
        include_mean = FALSE,
        allow_drift = FALSE
      )
    }

    if (is.null(fit)) {
      bl <- .default_block_len(n)
      if (verbose) {
        .log(paste0(
          "[WARM] Component ", k,
          if (nzchar(nm)) paste0(" (", nm, ")") else "",
          ": ARMA fit failed, using block bootstrap (block_len=", bl, ")"
        ))
      }
      sim_comp <- matrix(NA_real_, nrow = n, ncol = n_sim)
      for (j in seq_len(n_sim)) {
        sim_comp[, j] <- .warm_block_bootstrap(component, n = n)
      }
      output <- output + sim_comp
      next
    }

    if (check_diagnostics && !is.null(fit$residuals) && length(fit$residuals) > 1L) {
      lj <- tryCatch(
        stats::Box.test(fit$residuals, type = "Ljung-Box", lag = min(10L, length(fit$residuals) - 1L)),
        error = function(e) NULL
      )
      if (!is.null(lj) && is.finite(lj$p.value) && lj$p.value < 0.05) {
        warning(
          "Component ", k,
          if (nzchar(nm)) paste0(" (", nm, ")") else "",
          ": residuals show significant autocorrelation (Ljung-Box p = ",
          round(lj$p.value, 4), "). Model may be inadequate.",
          call. = FALSE
        )
      }
    }

    sim_centered <- .warm_simulate_from_fit(fit, n = n, n_sim = n_sim)

    sim_comp <- sim_centered + comp_mean

    if (match_variance) {
      sim_comp2 <- .variance_match_matrix(sim_comp, target_sd = comp_sd, tol = var_tol, target_mean = comp_mean)
      variance_corrections <- variance_corrections || !identical(sim_comp2, sim_comp)
      sim_comp <- sim_comp2
    }

    output <- output + sim_comp
  }

  if (verbose && match_variance && isTRUE(variance_corrections)) {
    .log("[WARM] Variance correction applied")
  }

  output
}


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

#' Check integer scalar
#'
#' @param x Object.
#' @return Logical scalar.
#' @keywords internal
.is_int_scalar <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x == as.integer(x)
}

#' Reconstruct original series from wavelet components
#'
#' @param components Matrix, data.frame, or list of numeric vectors.
#' @param n Integer. Expected length.
#' @return Numeric vector of length n.
#' @keywords internal
.reconstruct_series_from_components <- function(components, n) {

  if (is.matrix(components) || is.data.frame(components)) {
    if (nrow(components) != n) stop("Cannot reconstruct series: 'components' must have n rows.", call. = FALSE)
    return(as.numeric(rowSums(as.matrix(components))))

  } else if (is.list(components)) {
    lens <- vapply(components, length, integer(1))
    if (any(lens != n)) stop("Cannot reconstruct series: all component vectors must have length n.", call. = FALSE)
    m <- do.call(cbind, lapply(components, as.numeric))
    return(as.numeric(rowSums(m)))
  }

  stop("Cannot reconstruct series: unsupported 'components' type.", call. = FALSE)
}

#' Check ARIMA viability for WARM components
#'
#' @param x Numeric vector.
#' @param min_n Integer. Minimum length.
#' @param sd_eps Numeric. Minimum standard deviation threshold.
#' @param min_unique Integer. Minimum unique values.
#' @return Logical scalar.
#' @keywords internal
.warm_arima_viable <- function(x, min_n = 8L, sd_eps = 1e-8, min_unique = 3L) {
  if (!is.numeric(x) || anyNA(x)) return(FALSE)
  if (length(x) < as.integer(min_n)) return(FALSE)
  if (!is.finite(stats::sd(x)) || stats::sd(x) < sd_eps) return(FALSE)
  if (length(unique(x)) < as.integer(min_unique)) return(FALSE)
  TRUE
}

#' Fit ARIMA for WARM with safe fallback
#'
#' @param x Numeric vector (centered).
#' @param max_p Integer. Maximum AR order.
#' @param max_q Integer. Maximum MA order.
#' @param stationary Logical. Whether to enforce stationarity.
#' @param include_mean Logical. Whether to include mean.
#' @param allow_drift Logical. Whether to allow drift.
#' @return List with model or parameter entries, or NULL.
#' @keywords internal
.warm_fit_arima_safe <- function(x, max_p, max_q, stationary = TRUE, include_mean = FALSE, allow_drift = FALSE) {
  fit <- .fit_arma_aic(x, max_p = max_p, max_q = max_q)
  if (is.null(fit)) return(NULL)
  list(
    model = NULL,
    ar = fit$ar,
    ma = fit$ma,
    sigma2 = fit$sigma2,
    residuals = fit$residuals,
    order = fit$order,
    stationary = isTRUE(stationary),
    include_mean = isTRUE(include_mean),
    allow_drift = isTRUE(allow_drift)
  )
}

#' Fit ARIMA for bypass mode
#'
#' @param x Numeric vector (centered).
#' @param max_p Integer. Maximum AR order.
#' @param max_q Integer. Maximum MA order.
#' @param stationary Logical. Whether to enforce stationarity.
#' @param include_mean Logical. Whether to include mean.
#' @param allow_drift Logical. Whether to allow drift.
#' @return List with model or parameter entries, or NULL.
#' @keywords internal
.fit_warm_arima_forecast <- function(x, max_p, max_q, stationary = TRUE, include_mean = FALSE, allow_drift = FALSE) {
  .warm_fit_arima_safe(
    x = x,
    max_p = max_p,
    max_q = max_q,
    stationary = stationary,
    include_mean = include_mean,
    allow_drift = allow_drift
  )
}

#' Simulate from ARIMA fit for WARM
#'
#' @param fit List. Output from .warm_fit_arima_safe or .fit_warm_arima_forecast.
#' @param n Integer. Simulation length.
#' @param n_sim Integer. Number of realizations.
#' @return Numeric matrix n by n_sim.
#' @keywords internal
.warm_simulate_from_fit <- function(fit, n, n_sim) {
  if (!is.null(fit$model)) {
    out <- matrix(NA_real_, nrow = n, ncol = n_sim)
    for (j in seq_len(n_sim)) {
      out[, j] <- as.numeric(stats::simulate(fit$model, nsim = n))
    }
    return(out)
  }

  .simulate_arma_matrix(
    ar = fit$ar,
    ma = fit$ma,
    sd = sqrt(fit$sigma2),
    n = n,
    n_sim = n_sim,
    burnin = 100L
  )
}

#' Default block length for block bootstrap
#'
#' @param n Integer. Series length.
#' @return Integer block length.
#' @keywords internal
.default_block_len <- function(n) {
  n <- as.integer(n)
  bl <- floor(sqrt(n))
  bl <- max(3L, min(bl, n))
  bl
}

#' Block bootstrap for one realization
#'
#' @param x Numeric vector.
#' @param n Integer. Output length.
#' @param block_len Integer. Block length.
#' @return Numeric vector length n.
#' @keywords internal
.block_bootstrap <- function(x, n, block_len) {
  out <- .block_bootstrap_matrix(x, n = n, n_sim = 1L, block_len = block_len)
  as.numeric(out[, 1L])
}

#' Block bootstrap for WARM components
#'
#' @param x Numeric vector.
#' @param n Integer. Output length.
#' @param block_len Integer. Block length.
#' @return Numeric vector length n.
#' @keywords internal
.warm_block_bootstrap <- function(x, n, block_len = NULL) {
  if (is.null(block_len)) block_len <- .default_block_len(n)
  .block_bootstrap(x, n = n, block_len = block_len)
}
#' Block bootstrap for many realizations
#'
#' @description
#' Generates n_sim bootstrap realizations using contiguous blocks sampled with
#' replacement. The implementation is vectorized over realizations within each
#' block.
#'
#' @param x Numeric vector.
#' @param n Integer. Output length.
#' @param n_sim Integer. Number of realizations.
#' @param block_len Integer. Block length.
#' @return Numeric matrix n by n_sim.
#' @keywords internal
.block_bootstrap_matrix <- function(x, n, n_sim, block_len) {
  x <- as.numeric(x)
  n <- as.integer(n)
  n_sim <- as.integer(n_sim)
  block_len <- as.integer(block_len)

  m <- length(x)
  if (m < 1L) stop("Block bootstrap: input series is empty.", call. = FALSE)
  if (block_len < 1L) block_len <- 1L
  block_len <- min(block_len, m)

  n_blocks <- as.integer(ceiling(n / block_len))
  out <- matrix(NA_real_, nrow = n, ncol = n_sim)

  max_start <- m - block_len + 1L
  if (max_start < 1L) {
    # Degenerate case: block_len == m
    for (j in seq_len(n_sim)) out[, j] <- rep(x, length.out = n)
    return(out)
  }

  r0 <- 1L
  for (b in seq_len(n_blocks)) {
    r1 <- min(n, r0 + block_len - 1L)
    bl_now <- r1 - r0 + 1L

    starts <- sample.int(max_start, size = n_sim, replace = TRUE)
    idx <- outer(starts, seq_len(bl_now) - 1L, "+")
    seg <- x[idx]               # n_sim by bl_now
    out[r0:r1, ] <- t(seg)      # bl_now by n_sim

    r0 <- r1 + 1L
    if (r0 > n) break
  }

  out
}

#' Fit stationary non seasonal ARMA model by AIC grid search
#'
#' @description
#' Fits ARMA(p,q) models using stats::arima with include.mean = FALSE and
#' selects the model with minimum AIC.
#'
#' @param x Numeric vector. Should be centered.
#' @param max_p Integer. Maximum AR order.
#' @param max_q Integer. Maximum MA order.
#' @return List with ar, ma, sigma2, residuals, and order. Returns NULL if no
#'   candidate fit succeeds.
#' @keywords internal
.fit_arma_aic <- function(x, max_p = 2L, max_q = 2L) {
  x <- as.numeric(x)
  max_p <- as.integer(max_p)
  max_q <- as.integer(max_q)

  if (length(x) < 3L) return(NULL)
  if (!all(is.finite(x))) return(NULL)

  best_aic <- Inf
  best <- NULL

  for (p in 0:max_p) {
    for (q in 0:max_q) {

      fit <- tryCatch(
        stats::arima(
          x,
          order = c(p, 0L, q),
          seasonal = list(order = c(0L, 0L, 0L)),
          include.mean = FALSE,
          method = "ML",
          transform.pars = TRUE
        ),
        error = function(e) NULL
      )

      if (is.null(fit)) next
      aic <- tryCatch(as.numeric(stats::AIC(fit)), error = function(e) NA_real_)
      if (!is.finite(aic)) next

      if (aic < best_aic) {
        co <- fit$coef
        nm <- names(co)

        ar <- numeric(0)
        ma <- numeric(0)
        if (!is.null(nm)) {
          ar <- co[grepl("^ar", nm)]
          ma <- co[grepl("^ma", nm)]
        }

        best_aic <- aic
        best <- list(
          ar = as.numeric(ar),
          ma = as.numeric(ma),
          sigma2 = as.numeric(fit$sigma2),
          residuals = as.numeric(fit$residuals),
          order = c(p = as.integer(p), q = as.integer(q))
        )
      }
    }
  }

  best
}

#' Fast ARMA simulation for many realizations
#'
#' @description
#' Simulates a stationary ARMA process using direct recursion with Gaussian
#' innovations. This is intended for small p and q (0 to 2) and large n_sim.
#'
#' @param ar Numeric vector. AR coefficients.
#' @param ma Numeric vector. MA coefficients.
#' @param sd Numeric scalar. Innovation standard deviation.
#' @param n Integer. Output length.
#' @param n_sim Integer. Number of realizations.
#' @param burnin Integer. Burn in length.
#' @return Numeric matrix n by n_sim.
#' @keywords internal
.simulate_arma_matrix <- function(ar, ma, sd, n, n_sim, burnin = 100L) {
  ar <- as.numeric(ar)
  ma <- as.numeric(ma)
  sd <- as.numeric(sd)
  n <- as.integer(n)
  n_sim <- as.integer(n_sim)
  burnin <- as.integer(burnin)

  if (!is.finite(sd) || sd <= 0) sd <- 1
  if (n < 1L || n_sim < 1L) stop("ARMA simulation requires n and n_sim >= 1.", call. = FALSE)

  p <- length(ar)
  q <- length(ma)

  n_tot <- n + burnin
  max_lag <- max(p, q, 1L)

  # Innovations must be (n_tot + max_lag) by n_sim. Using only n_tot + max_lag
  # values would recycle and can cause subscript out of bounds errors.
  eps <- matrix(
    stats::rnorm((n_tot + max_lag) * n_sim, mean = 0, sd = sd),
    nrow = n_tot + max_lag,
    ncol = n_sim
  )
  x <- matrix(0, nrow = n_tot + max_lag, ncol = n_sim)

  for (t in (max_lag + 1L):(n_tot + max_lag)) {
    xt <- eps[t, ]

    if (q >= 1L) xt <- xt + ma[1L] * eps[t - 1L, ]
    if (q >= 2L) xt <- xt + ma[2L] * eps[t - 2L, ]

    if (p >= 1L) xt <- xt + ar[1L] * x[t - 1L, ]
    if (p >= 2L) xt <- xt + ar[2L] * x[t - 2L, ]

    if (p > 2L) {
      for (i in 3L:p) xt <- xt + ar[i] * x[t - i, ]
    }
    if (q > 2L) {
      for (i in 3L:q) xt <- xt + ma[i] * eps[t - i, ]
    }

    x[t, ] <- xt
  }

  x_out <- x[(max_lag + burnin + 1L):(max_lag + burnin + n), , drop = FALSE]
  x_out
}

#' Fast column standard deviations
#'
#' @param x Numeric matrix.
#' @return Numeric vector of column standard deviations.
#' @keywords internal
.fast_col_sd <- function(x) {
  m1 <- colMeans(x)
  m2 <- colMeans(x * x)
  v <- pmax(m2 - m1 * m1, 0)
  sqrt(v)
}

#' Vectorized variance matching for simulation matrices
#'
#' @description
#' If the relative standard deviation mismatch exceeds tol, rescales each column
#' to match target_sd while keeping the column mean at target_mean.
#'
#' @param sim Numeric matrix.
#' @param target_sd Numeric scalar.
#' @param tol Numeric scalar in [0, 1].
#' @param target_mean Numeric scalar.
#' @return Numeric matrix with the same dimension as sim.
#' @keywords internal
.variance_match_matrix <- function(sim, target_sd, tol, target_mean) {
  if (!is.finite(target_sd) || target_sd <= 0) return(sim)

  sim_sd <- .fast_col_sd(sim)
  sim_mean <- colMeans(sim)

  ok <- is.finite(sim_sd) & (sim_sd > 0)
  rel <- rep(0, length(sim_sd))
  rel[ok] <- abs(sim_sd[ok] - target_sd) / target_sd

  do_fix <- ok & (rel > tol)
  if (!any(do_fix)) return(sim)

  scale <- rep(1, length(sim_sd))
  scale[do_fix] <- target_sd / sim_sd[do_fix]

  sim2 <- sim
  sim2[, do_fix] <- sweep(sim2[, do_fix, drop = FALSE], 2, sim_mean[do_fix], "-")
  sim2[, do_fix] <- sweep(sim2[, do_fix, drop = FALSE], 2, scale[do_fix], "*")
  sim2[, do_fix] <- sweep(sim2[, do_fix, drop = FALSE], 2, target_mean, "+")

  sim2
}
