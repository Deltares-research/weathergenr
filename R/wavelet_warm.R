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
#'   residuals and warns when residual autocorrelation remains. In component
#'   mode, also warns when pre-correction aggregate variance mismatch is large.
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
#' Cross-component covariance correction (component mode only)
#' MODWT MRA components are additive but not orthogonal; fitting independent
#' ARMAs and summing inflates ensemble variance when cross-component covariances
#' are non-zero. Two tiers of correction are applied.
#' - Tier 2 (Cholesky): when >= 2 variable components are ARMA-viable and the
#'   component covariance matrix is well-conditioned (condition number < 1e6),
#'   components are decorrelated via X %*% L^{-1} (L = chol(cov(X))), ARMA is
#'   fit to the whitened series, and simulations are re-correlated by
#'   multiplying back by L before summing.
#' - Tier 1 (aggregate): after summing, the total simulated series is rescaled
#'   to match the observed total standard deviation when the relative mismatch
#'   exceeds var_tol. Applied in both the Tier 2 path (as a safety check) and
#'   the fallback path (as the primary correction).
#'
#' Fallback logic
#' - If any variable component fails the ARMA viability pre-scan, Tier 2 is
#'   skipped and per-component variance matching is applied instead.
#' - If the covariance matrix is ill-conditioned, Tier 2 is skipped; Tier 1
#'   remains active.
#' - In bypass mode, fallback is applied to the whole series (block bootstrap).
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

    bypass_max_p <- if (n_fit < 20L) 1L else 2L
    bypass_max_q <- if (n_fit < 20L) 1L else 2L

    fit <- .fit_warm_arima_forecast(
      x = x,
      max_p = bypass_max_p,
      max_q = bypass_max_q
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

  # Store original component means; used when re-adding means after re-correlation.
  comp_means_all <- vapply(comp_list, mean, numeric(1))

  # Identify variable components (non-constant).
  var_idx <- which(vapply(comp_list, function(x) stats::sd(x) >= 1e-10, logical(1)))
  nvar    <- length(var_idx)

  # --------------------------------------------------------------------------
  # Tier 2: Cholesky decorrelation setup
  # --------------------------------------------------------------------------
  # Attempt decorrelation only when >= 2 variable components exist and all of
  # them pass the ARMA viability check. The viability pre-scan lets us skip
  # decorrelation upfront rather than discovering a bootstrap fallback is
  # needed mid-loop (which would leave comp_sims in a mixed decorr/original
  # space that cannot be cleanly re-correlated).
  use_decorrelation <- FALSE
  L_chol            <- NULL

  if (nvar >= 2L) {
    # Pre-scan: all variable components must be ARMA-viable.
    all_viable <- TRUE
    for (ki in seq_along(var_idx)) {
      k    <- var_idx[ki]
      nm_k <- if (!is.null(comp_names) && length(comp_names) >= k && !is.na(comp_names[k])) comp_names[k] else ""
      is_smooth_k <- nzchar(nm_k) && grepl("^S", nm_k)
      max_p_k <- if (is_smooth_k) 1L else 2L
      max_q_k <- if (is_smooth_k) 0L else 2L
      ctr_k <- comp_list[[k]] - comp_means_all[k]
      if (!.warm_arima_viable(ctr_k, min_n = max(8L, max_p_k + max_q_k + 2L), sd_eps = 1e-8, min_unique = 3L)) {
        all_viable <- FALSE
        break
      }
    }

    if (all_viable) {
      comp_mat_var <- do.call(cbind, lapply(var_idx, function(k) comp_list[[k]] - comp_means_all[k]))
      comp_cov     <- stats::cov(comp_mat_var)
      cond_num     <- tryCatch(kappa(comp_cov, exact = FALSE), error = function(e) Inf)

      if (is.finite(cond_num) && cond_num < 1e6) {
        L_chol <- tryCatch(chol(comp_cov), error = function(e) NULL)  # upper: t(L) %*% L = cov
        if (!is.null(L_chol)) {
          # Z = X %*% L^{-1}  =>  cov(Z) ~ I
          L_inv    <- backsolve(L_chol, diag(nvar))
          Z_mat    <- comp_mat_var %*% L_inv  # n_fit x nvar, decorrelated, centered
          for (ki in seq_along(var_idx)) comp_list[[var_idx[ki]]] <- Z_mat[, ki]
          use_decorrelation <- TRUE
          if (verbose) .log(paste0(
            "[WARM] Cholesky decorrelation active (cond=", round(cond_num, 1), ", nvar=", nvar, ")"
          ))
        }
      }

      if (!use_decorrelation && verbose) .log(paste0(
        "[WARM] Cholesky decorrelation skipped (cond=",
        if (is.finite(cond_num)) round(cond_num, 1) else "Inf",
        "); falling back to per-component variance correction"
      ))
    } else if (verbose) {
      .log("[WARM] Cholesky decorrelation skipped: one or more variable components not ARMA-viable")
    }
  }

  # comp_sims[[k]] stores the centered simulation for component k (n x n_sim).
  # Means are accounted for separately and added back after the loop.
  comp_sims          <- vector("list", length(comp_list))
  variance_corrections <- FALSE

  for (k in seq_along(comp_list)) {

    component <- as.numeric(comp_list[[k]])
    # When decorrelation is active, var_idx components hold Z columns (centered).
    # Constant components and non-variable components are untouched.
    comp_sd <- stats::sd(component)
    if (!is.finite(comp_sd)) stop("Component ", k, " has non-finite standard deviation.", call. = FALSE)

    nm <- ""
    if (!is.null(comp_names) && length(comp_names) >= k && !is.na(comp_names[k])) nm <- comp_names[k]
    is_smooth <- nzchar(nm) && grepl("^S", nm)

    # Constant component: store a zero matrix; mean is added back post-loop.
    if (comp_sd < 1e-10) {
      comp_sims[[k]] <- matrix(0, nrow = n, ncol = n_sim)
      next
    }

    if (!is.null(seed)) set.seed(seed + k * 1000L)

    # Decorrelated components are already centered (mean ~ 0); others need centering.
    is_decorr_comp <- use_decorrelation && (k %in% var_idx)
    comp_mean  <- if (is_decorr_comp) 0 else comp_means_all[k]
    centered   <- component - comp_mean

    # After Cholesky whitening, a decorrelated component is a linear combination
    # of all original components and is no longer smooth regardless of its name.
    # Restrict max_p/max_q only for non-decorrelated smooth components.
    is_smooth_effective <- is_smooth && !is_decorr_comp
    max_p_use <- if (is_smooth_effective) 1L else 2L
    max_q_use <- if (is_smooth_effective) 0L else 2L

    is_viable <- .warm_arima_viable(
      centered,
      min_n      = max(8L, max_p_use + max_q_use + 2L),
      sd_eps     = 1e-8,
      min_unique = 3L
    )

    fit <- NULL
    if (is_viable) {
      fit <- .warm_fit_arima_safe(
        x            = centered,
        max_p        = max_p_use,
        max_q        = max_q_use,
        stationary   = TRUE,
        include_mean = FALSE,
        allow_drift  = FALSE
      )
    }

    if (is.null(fit)) {
      # Bootstrap fallback. use_decorrelation cannot be TRUE here because the
      # pre-scan ensures all variable components are ARMA-viable. This branch
      # is therefore only reached in non-decorrelation mode.
      bl <- .default_block_len(n)
      if (verbose) {
        .log(paste0(
          "[WARM] Component ", k,
          if (nzchar(nm)) paste0(" (", nm, ")") else "",
          ": ARMA fit failed, using block bootstrap (block_len=", bl, ")"
        ))
      }
      sim_comp <- matrix(NA_real_, nrow = n, ncol = n_sim)
      for (j in seq_len(n_sim)) sim_comp[, j] <- .warm_block_bootstrap(component, n = n)
      # Bootstrap output already includes the component mean; zero out the stored mean.
      comp_means_all[k] <- 0
      comp_sims[[k]] <- sim_comp
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

    # In non-decorrelation mode, apply per-component variance correction now.
    # In decorrelation mode, variance is corrected after re-correlation (Tier 1).
    if (!use_decorrelation && match_variance) {
      sim_comp2 <- .variance_match_matrix(
        sim_centered + comp_mean,
        target_sd   = stats::sd(comp_list[[k]]),
        tol         = var_tol,
        target_mean = comp_mean
      )
      variance_corrections <- variance_corrections || !identical(sim_comp2, sim_centered + comp_mean)
      # Store mean-subtracted version so post-loop mean addition is consistent.
      comp_sims[[k]] <- sim_comp2 - comp_mean
    } else {
      comp_sims[[k]] <- sim_centered  # centered simulation; mean added post-loop
    }
  }

  # --------------------------------------------------------------------------
  # Post-loop: aggregate and apply variance correction
  # --------------------------------------------------------------------------

  if (use_decorrelation) {
    # Tier 2: re-correlate.
    # X_sim (sum over components) = sum_k  sum_j Z_sim_j * L_chol[j, k]
    #                              = sum_j Z_sim_j * rowSums(L_chol)[j]
    L_row_sums <- rowSums(L_chol)
    output <- matrix(0, nrow = n, ncol = n_sim)
    for (ki in seq_along(var_idx)) {
      output <- output + comp_sims[[var_idx[ki]]] * L_row_sums[ki]
    }
    # Add constant-component contributions.
    const_idx <- setdiff(seq_along(comp_sims), var_idx)
    for (k in const_idx) output <- output + comp_sims[[k]]
    # Re-add original means (sum of all components' means = observed total mean).
    output <- output + sum(comp_means_all)
    if (verbose) .log("[WARM] Cholesky re-correlation applied")

  } else {
    # No decorrelation: simple sum; per-component means are stored separately.
    output <- matrix(0, nrow = n, ncol = n_sim)
    for (k in seq_along(comp_sims)) output <- output + comp_sims[[k]] + comp_means_all[k]
    if (verbose && match_variance && isTRUE(variance_corrections)) {
      .log("[WARM] Per-component variance correction applied")
    }
  }

  # Tier 1: aggregate variance correction applied in both paths.
  # In decorrelation mode this is a safety check; in fallback mode it is the
  # primary correction for cross-component covariance leakage.
  if (match_variance) {
    obs_total      <- .reconstruct_series_from_components(components, n = n_fit)
    obs_total_sd   <- stats::sd(obs_total)
    obs_total_mean <- mean(obs_total)
    sim_total_sd   <- stats::sd(output)
    variance_identity_error <- NA_real_
    if (is.finite(obs_total_sd) && obs_total_sd > 0 && is.finite(sim_total_sd)) {
      variance_identity_error <- abs(sim_total_sd - obs_total_sd) / obs_total_sd
    }
    if (check_diagnostics &&
        is.finite(variance_identity_error) &&
        variance_identity_error > max(var_tol, 0.1)) {
      warning(
        "WARM pre-correction variance mismatch is high (relative error = ",
        round(variance_identity_error, 4),
        "). This can indicate non-negligible cross-component covariance.",
        call. = FALSE
      )
    }
    output2 <- .variance_match_matrix(
      output,
      target_sd   = obs_total_sd,
      tol         = var_tol,
      target_mean = obs_total_mean
    )
    if (!identical(output2, output)) {
      if (verbose) .log("[WARM] Aggregate variance correction applied")
    }
    output <- output2
  }

  output
}


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

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
  n_rows <- n_tot + max_lag

  # Chunk columns to avoid allocating a full innovations matrix for very large n_sim.
  chunk_size <- min(1000L, n_sim)
  n_chunks <- as.integer(ceiling(n_sim / chunk_size))
  x_out <- matrix(0, nrow = n, ncol = n_sim)

  for (ch in seq_len(n_chunks)) {
    j0 <- (ch - 1L) * chunk_size + 1L
    j1 <- min(ch * chunk_size, n_sim)
    idx <- seq.int(j0, j1)
    n_col <- length(idx)

    eps <- matrix(
      stats::rnorm(n_rows * n_col, mean = 0, sd = sd),
      nrow = n_rows,
      ncol = n_col
    )
    x <- matrix(0, nrow = n_rows, ncol = n_col)

    for (t in (max_lag + 1L):n_rows) {
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

    x_out[, idx] <- x[(max_lag + burnin + 1L):(max_lag + burnin + n), , drop = FALSE]
  }

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
  sim_fix <- sim2[, do_fix, drop = FALSE]
  sim_fix <- (sim_fix - rep(sim_mean[do_fix], each = nrow(sim_fix))) *
    rep(scale[do_fix], each = nrow(sim_fix)) + target_mean
  sim2[, do_fix] <- sim_fix

  sim2
}
