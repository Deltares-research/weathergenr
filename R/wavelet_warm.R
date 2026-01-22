#' Simulate synthetic series with Wavelet Autoregressive Modeling (WARM)
#'
#' @description
#' Generates \code{n_sim} synthetic realizations of length \code{n} based on an observed
#' time series represented by wavelet components.
#'
#' Two operating modes are supported:
#'
#' - Component mode (\code{n >= bypass_n}): each wavelet component is fit with an ARIMA model,
#'   simulated forward, and then components are summed to reconstruct synthetic series.
#'
#' - Bypass mode (\code{n < bypass_n}): component-level simulation is skipped. A stationary
#'   AR/ARMA model is fit to the original series and simulated directly. This is intended
#'   to avoid unstable low-frequency behavior on short records.
#'
#' @param components Matrix, data.frame, or list of numeric vectors. Wavelet components
#'   (for example, MODWT-MRA detail and smooth components). In component mode, each component
#'   must have length \code{n}. In bypass mode, components are only used to reconstruct the
#'   original series if \code{series_obs} is not provided.
#' @param n Integer. Length of each simulated series.
#' @param n_sim Integer. Number of realizations to generate. Default is \code{1000}.
#' @param seed Optional integer. Base RNG seed for reproducibility. If provided, seeds are set
#'   deterministically per component in component mode.
#' @param series_obs Optional numeric vector. Observed original series. Recommended in bypass
#'   mode so the function does not depend on reconstructing the original series from
#'   \code{components}.
#' @param bypass_n Integer. Threshold for bypass mode. If \code{n < bypass_n}, fit AR/ARMA
#'   directly to the original series. Default is \code{25}.
#' @param match_variance Logical. If \code{TRUE} (default), variance matching may be applied
#'   to simulated outputs when relative SD mismatch exceeds \code{var_tol}.
#' @param var_tol Numeric in [0, 1]. Relative tolerance for SD mismatch that triggers variance
#'   rescaling. Default is \code{0.1} (10\%).
#' @param check_diagnostics Logical. If \code{TRUE}, runs a simple residual autocorrelation
#'   check (Ljung-Box) on fitted ARIMA models and warns if residual structure remains. Default
#'   is \code{FALSE}.
#' @param verbose Logical. If \code{TRUE} (default), prints informational messages.
#'
#' @return A numeric matrix of dimension \code{n x n_sim}. Each column is a simulated realization.
#'
#' @details
#' Component mode (\code{n >= bypass_n})
#'
#' - Each component is centered prior to ARIMA fitting; the observed component mean is
#'   re-added after simulation.
#' - Components are simulated independently and then summed. Because MODWT-MRA components are
#'   additive but not orthogonal, this does not preserve cross-component covariance unless an
#'   explicit coupling step is implemented externally.
#' - Smooth components (names starting with \code{"S"}) are constrained to AR(1)
#'   (\code{max.p = 1, max.q = 0}) to reduce spurious low-frequency oscillations.
#' - Mean and drift terms are disabled during fitting (\code{include.mean = FALSE},
#'   \code{allowdrift = FALSE}) because components are explicitly centered.
#' - If ARIMA fitting fails for a component, the function falls back to bootstrap resampling
#'   (with replacement) from the observed component values.
#' - If \code{match_variance = TRUE}, simulated component variance may be rescaled to match the
#'   observed component SD when relative mismatch exceeds \code{var_tol}.
#'
#' Bypass mode (\code{n < bypass_n})
#'
#' - Fits a stationary AR/ARMA model (\code{max.p = 2, max.q = 2}) to the original series and
#'   simulates forward using \code{stats::simulate()}.
#' - If model fitting fails, the function falls back to bootstrap resampling from the observed
#'   series.
#' - Optional variance matching is applied to the simulated series when enabled.
#'
#' @section When to use bypass mode:
#' Use bypass mode for short records where component-level ARIMA fits are unstable and may
#' introduce artificial low-frequency resonance. Provide \code{series_obs} explicitly to avoid
#' ambiguity from reconstructing the original series from \code{components}.
#'
#' @importFrom stats sd simulate Box.test var
#' @importFrom forecast auto.arima
#' @export
simulate_warm <- function(
    components = NULL,
    n = NULL,
    n_sim = 1000,
    seed = NULL,
    series_obs = NULL,
    bypass_n = 25L,
    match_variance = TRUE,
    var_tol = 0.1,
    check_diagnostics = FALSE,
    verbose = TRUE
) {

  if (is.null(n)) stop("Input 'n' must be specified.")
  if (!.is_int_scalar(n) || n < 1L) stop("'n' must be a positive integer.", call. = FALSE)
  n <- as.integer(n)

  if (!.is_int_scalar(n_sim) || n_sim < 1L) stop("'n_sim' must be a positive integer.", call. = FALSE)
  n_sim <- as.integer(n_sim)

  if (!.is_int_scalar(bypass_n) || bypass_n < 5L) stop("'bypass_n' must be an integer >= 5.", call. = FALSE)
  bypass_n <- as.integer(bypass_n)

  if (!is.logical(verbose) || length(verbose) != 1L) stop("'verbose' must be TRUE/FALSE.", call. = FALSE)

  if (!is.logical(match_variance) || length(match_variance) != 1L) {
    stop("'match_variance' must be TRUE/FALSE.", call. = FALSE)
  }

  if (!is.numeric(var_tol) || length(var_tol) != 1L || !is.finite(var_tol) || var_tol < 0 || var_tol > 1) {
    stop("'var_tol' must be between 0 and 1.", call. = FALSE)
  }

  if (!is.logical(check_diagnostics) || length(check_diagnostics) != 1L) {
    stop("'check_diagnostics' must be TRUE/FALSE.", call. = FALSE)
  }

  # RNG state management (restore old seed on exit)
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
  # BYPASS MODE: n < bypass_n -> fit AR/ARMA on original series directly
  # --------------------------------------------------------------------------
  if (n < bypass_n) {

    # Prefer explicit observed series; otherwise reconstruct from components.
    if (!is.null(series_obs)) {
      if (!is.numeric(series_obs) || anyNA(series_obs)) stop("'series_obs' must be numeric with no missing values.", call. = FALSE)
      series_obs <- as.numeric(series_obs)
      if (length(series_obs) != n) stop("'series_obs' length must equal 'n'.", call. = FALSE)

    } else {

      if (is.null(components)) stop("For n < bypass_n, provide 'series_obs' or 'components' to reconstruct the series.", call. = FALSE)
      series_obs <- .reconstruct_series_from_components(components, n = n)
    }

    if (!is.null(seed)) set.seed(seed)

    # Fit a simple stationary AR/ARMA model.
    model <- tryCatch(
      {
        suppressPackageStartupMessages({
          forecast::auto.arima(
            series_obs,
            max.p = 2,
            max.q = 2,
            max.P = 0,
            max.Q = 0,
            stationary = TRUE,
            seasonal = FALSE,
            allowdrift = FALSE
          )
        })
      },
      error = function(e) {
        warning(
          "Bypass AR/ARMA fitting failed (n=", n, "): ", e$message,
          "\nFalling back to resampling with replacement from the observed series.",
          call. = FALSE
        )
        NULL
      }
    )

    out <- matrix(0, nrow = n, ncol = n_sim)

    if (is.null(model)) {
      for (j in seq_len(n_sim)) out[, j] <- sample(series_obs, n, replace = TRUE)
      return(out)
    }

    # Optional diagnostics
    if (check_diagnostics && length(model$residuals) > 1) {
      lj <- tryCatch(
        stats::Box.test(model$residuals, type = "Ljung-Box",
                        lag = min(10, length(model$residuals) - 1)),
        error = function(e) NULL
      )
      if (!is.null(lj) && is.finite(lj$p.value) && lj$p.value < 0.05) {
        warning(
          "Bypass AR/ARMA: residuals show significant autocorrelation ",
          "(Ljung-Box p = ", round(lj$p.value, 4), "). Model may be inadequate.",
          call. = FALSE
        )
      }
    }

    target_sd <- stats::sd(series_obs)

    for (j in seq_len(n_sim)) {
      sim <- as.numeric(stats::simulate(model, nsim = n))

      if (isTRUE(match_variance) && is.finite(target_sd) && target_sd > 0) {
        sim_sd <- stats::sd(sim)
        if (is.finite(sim_sd) && sim_sd > 0) {
          rel <- abs(sim_sd - target_sd) / target_sd
          if (rel > var_tol) {
            sim_mean <- mean(sim)
            obs_mean <- mean(series_obs)
            sim <- obs_mean + (sim - sim_mean) * (target_sd / sim_sd)
          }
        }
      }

      out[, j] <- sim
    }

    if (verbose) .log(paste0("[WARM] Bypass mode used (n=", n, " < ", bypass_n, "): AR/ARMA on original series"))
    return(out)
  }

  # --------------------------------------------------------------------------
  # COMPONENT MODE: n >= bypass_n -> simulate per component and sum
  # --------------------------------------------------------------------------
  if (is.null(components)) stop("Input 'components' must not be NULL for component mode.", call. = FALSE)

  # Efficient component list conversion + capture component names if available
  comp_names <- NULL

  if (is.matrix(components)) {
    if (!is.numeric(components)) stop("Matrix 'components' must be numeric.", call. = FALSE)
    if (nrow(components) != n) stop("Matrix 'components' must have n rows (same as 'n').", call. = FALSE)
    ncomp <- ncol(components)
    comp_list <- vector("list", ncomp)
    for (k in seq_len(ncomp)) comp_list[[k]] <- components[, k]
    comp_names <- colnames(components)

  } else if (is.data.frame(components)) {
    col_types <- vapply(components, is.numeric, logical(1))
    if (!all(col_types)) stop("All columns of 'components' must be numeric.", call. = FALSE)
    if (nrow(components) != n) stop("Data.frame 'components' must have n rows (same as 'n').", call. = FALSE)
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
      ". Remove/impute NAs before calling simulate_warm().",
      call. = FALSE
    )
  }

  output <- matrix(0, nrow = n, ncol = n_sim)
  variance_corrections <- logical(ncomp)

  for (k in seq_len(ncomp)) {

    component <- as.numeric(comp_list[[k]])
    comp_sd <- stats::sd(component)
    if (!is.finite(comp_sd)) stop("Component ", k, " has non-finite standard deviation.", call. = FALSE)

    nm <- if (!is.null(comp_names) && length(comp_names) >= k && !is.na(comp_names[k])) comp_names[k] else ""
    is_smooth <- nzchar(nm) && grepl("^S", nm)

    if (comp_sd < 1e-10) {
      warning(
        "Component ", k, if (nzchar(nm)) paste0(" (", nm, ")") else "",
        " is essentially constant. Returning constant values without ARIMA modeling.",
        call. = FALSE
      )
      output <- output + mean(component)
      next
    }

    if (!is.null(seed)) set.seed(seed + k * 1000L)

    # Center for fitting, re-add observed mean after simulation.
    comp_mean <- mean(component)
    centered <- component - comp_mean
    target_sd <- comp_sd

    # Constrain smooth components to AR(1) (avoid low-frequency resonance)
    max_p_use <- if (is_smooth) 1 else 2
    max_q_use <- if (is_smooth) 0 else 2

    MODEL <- tryCatch(
      {
        suppressPackageStartupMessages({
          forecast::auto.arima(
            centered,
            max.p = max_p_use,
            max.q = max_q_use,
            max.P = 0,
            max.Q = 0,
            stationary = TRUE,
            seasonal = FALSE,
            include.mean = FALSE,
            allowdrift = FALSE
          )
        })
      },
      error = function(e) {
        warning(
          "ARIMA fitting failed for component ", k,
          if (nzchar(nm)) paste0(" (", nm, ")") else "",
          ": ", e$message,
          "\nFalling back to resampling with replacement.",
          call. = FALSE
        )
        NULL
      }
    )

    if (is.null(MODEL)) {
      for (j in seq_len(n_sim)) {
        output[, j] <- output[, j] + sample(component, n, replace = TRUE)
      }
      next
    }

    if (check_diagnostics && length(MODEL$residuals) > 1) {
      ljung_test <- tryCatch(
        stats::Box.test(
          MODEL$residuals,
          type = "Ljung-Box",
          lag = min(10, length(MODEL$residuals) - 1)
        ),
        error = function(e) NULL
      )

      if (!is.null(ljung_test) && is.finite(ljung_test$p.value) && ljung_test$p.value < 0.05) {
        warning(
          "Component ", k,
          if (nzchar(nm)) paste0(" (", nm, ")") else "",
          ": residuals show significant autocorrelation ",
          "(Ljung-Box p = ", round(ljung_test$p.value, 4), "). Model may be inadequate.",
          call. = FALSE
        )
      }
    }

    needs_variance_check <- isTRUE(match_variance)

    for (j in seq_len(n_sim)) {
      sim_centered <- as.numeric(stats::simulate(MODEL, nsim = n))
      simulated <- sim_centered + comp_mean

      if (needs_variance_check) {
        sim_sd <- stats::sd(simulated)
        if (is.finite(sim_sd) && sim_sd > 0) {
          relative_diff <- abs(sim_sd - target_sd) / target_sd
          if (relative_diff > var_tol) {
            sim_mean <- mean(simulated)
            simulated <- comp_mean + (simulated - sim_mean) * (target_sd / sim_sd)
            if (!variance_corrections[k]) variance_corrections[k] <- TRUE
          }
        }
      }

      output[, j] <- output[, j] + simulated
    }
  }

  if (verbose && match_variance && any(variance_corrections)) {
    .log("[WARM] Variance correction applied")
  }

  output
}


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

#' @keywords internal
.is_int_scalar <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x == as.integer(x)
}

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
