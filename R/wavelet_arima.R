#' Wavelet Autoregressive Modeling (WARM)
#'
#' Simulates synthetic time series by modeling each wavelet component (signal or noise)
#' with an ARIMA model, and then summing the simulated components. Optionally enforces
#' variance matching to preserve statistical properties of the original components.
#'
#' @param wavelet.components A list or matrix where each column (or list element) is a
#'   numeric vector corresponding to a wavelet component (low-frequency signal or noise).
#' @param sim.year.num Integer. Desired length (number of years or timesteps) of each
#'   simulated series.
#' @param sim.num Integer. Number of synthetic series to produce. Default: 1000.
#' @param seed Optional. Integer random seed for reproducibility.
#' @param match.variance Logical. If TRUE, rescale simulated components to match the
#'   variance of the original components (default: TRUE). Recommended for preserving
#'   statistical properties.
#' @param variance.tolerance Numeric. Relative tolerance for variance matching
#'   (default: 0.1 = 10\%). Only applies if match.variance = TRUE.
#' @param check.diagnostics Logical. If TRUE, perform basic ARIMA model diagnostics
#'   and issue warnings if models appear inadequate (default: FALSE).
#'
#' @return A matrix of dimension \code{sim.year.num} x \code{sim.num}, where each
#'   column is a synthetic time series realization.
#'
#' @details
#' The function implements the Wavelet Autoregressive Modeling (WARM) approach:
#'
#' 1. Each wavelet component is centered (mean removed)
#'
#' 2. An ARIMA(p,0,q) model is fitted to the centered component
#'
#' 3. Multiple realizations are simulated from the fitted model
#'
#' 4. Mean is restored to each simulation
#'
#' 5. Optionally, variance is corrected to match the original component
#'
#' 6. Components are summed to produce final synthetic series
#'
#' Variance matching ensures that ARIMA models that slightly underfit do not produce
#' under-dispersed simulations. This is particularly important for maintaining the
#' correct amplitude of low-frequency oscillations.
#'
#' @references
#' Steinschneider, S. and Lall, U. (2015). A hierarchical Bayesian regional model
#' for nonstationary precipitation extremes in Northern California conditioned on
#' tropical moisture exports. \emph{Water Resources Research}, 51(3), 1472-1492.
#'
#' @examples
#' \dontrun{
#' # Simulate 10 synthetic annual series (length 30) from two wavelet components
#' set.seed(42)
#' wavelet.components <- list(
#'   signal = sin(2 * pi * 1:30 / 8) + rnorm(30, 0, 0.2),
#'   noise  = rnorm(30, 0, 0.5)
#' )
#'
#' # With variance matching (recommended)
#' result <- waveletARIMA(
#'   wavelet.components,
#'   sim.year.num = 30,
#'   sim.num = 10,
#'   seed = 123,
#'   match.variance = TRUE
#' )
#'
#' # Check variance preservation
#' original_var <- var(rowSums(sapply(wavelet.components, identity)))
#' simulated_var <- apply(result, 2, var)
#' cat("Original variance:", round(original_var, 3), "\n")
#' cat("Simulated variance range:",
#'     round(range(simulated_var), 3), "\n")
#'
#' # Visualize
#' matplot(result, type = "l", lty = 1,
#'         col = adjustcolor("black", alpha = 0.3),
#'         main = "Synthetic Wavelet-ARIMA Realizations")
#' lines(rowSums(sapply(wavelet.components, identity)),
#'       col = "red", lwd = 2)
#' legend("topright", legend = c("Observed", "Simulated"),
#'        col = c("red", "black"), lwd = c(2, 1))
#' }
#'
#' @importFrom stats sd simulate
#' @importFrom forecast auto.arima
#' @export
wavelet_arima <- function(wavelet.components = NULL,
                          sim.year.num = NULL,
                          sim.num = 1000,
                          seed = NULL,
                          match.variance = TRUE,
                          variance.tolerance = 0.1,
                          check.diagnostics = FALSE) {

  # ============================================================================
  # Input Validation
  # ============================================================================

  if (is.null(wavelet.components)) {
    stop("Input 'wavelet.components' must not be NULL.")
  }

  if (is.null(sim.year.num)) {
    stop("Input 'sim.year.num' must be specified.")
  }

  # --- FIX 1: enforce integer-ish sim.year.num and sim.num with clean messages
  .is_int_scalar <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && (x %% 1 == 0)
  }

  if (!.is_int_scalar(sim.year.num) || sim.year.num < 1L) {
    stop("'sim.year.num' must be a positive integer.")
  }
  sim.year.num <- as.integer(sim.year.num)

  if (!.is_int_scalar(sim.num) || sim.num < 1L) {
    stop("'sim.num' must be a positive integer.")
  }
  sim.num <- as.integer(sim.num)

  if (!is.logical(match.variance) || length(match.variance) != 1L) {
    stop("'match.variance' must be logical (TRUE/FALSE).")
  }

  if (!is.numeric(variance.tolerance) || length(variance.tolerance) != 1L ||
      !is.finite(variance.tolerance) || variance.tolerance < 0 || variance.tolerance > 1) {
    stop("'variance.tolerance' must be between 0 and 1.")
  }

  if (!is.logical(check.diagnostics) || length(check.diagnostics) != 1L) {
    stop("'check.diagnostics' must be logical (TRUE/FALSE).")
  }

  # ============================================================================
  # Efficient component list conversion
  # ============================================================================

  if (is.matrix(wavelet.components)) {
    if (!is.numeric(wavelet.components)) {
      stop("Matrix 'wavelet.components' must be numeric.")
    }
    ncomp <- ncol(wavelet.components)
    comp_list <- vector("list", ncomp)
    for (k in seq_len(ncomp)) comp_list[[k]] <- wavelet.components[, k]

  } else if (is.data.frame(wavelet.components)) {
    col_types <- vapply(wavelet.components, is.numeric, logical(1))
    if (!all(col_types)) stop("All columns of 'wavelet.components' must be numeric.")
    ncomp <- ncol(wavelet.components)
    comp_list <- vector("list", ncomp)
    for (k in seq_len(ncomp)) comp_list[[k]] <- wavelet.components[[k]]

  } else if (is.list(wavelet.components)) {
    elem_types <- vapply(wavelet.components, is.numeric, logical(1))
    if (!all(elem_types)) stop("All elements of 'wavelet.components' must be numeric vectors.")
    ncomp <- length(wavelet.components)
    comp_list <- wavelet.components

  } else {
    stop("'wavelet.components' must be a matrix, data.frame, or list of numeric vectors.")
  }

  # --- FIX 2: explicit NA check inside each component vector (before sd/mean)
  na_comp <- vapply(comp_list, function(x) anyNA(x), logical(1))
  if (any(na_comp)) {
    bad <- which(na_comp)
    stop(
      "Missing values detected in wavelet component(s): ",
      paste(bad, collapse = ", "),
      ". Remove/impute NAs before calling wavelet_arima().",
      call. = FALSE
    )
  }

  # ============================================================================
  # RNG State Management
  # ============================================================================

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be NULL or a single finite number.", call. = FALSE)
    }

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

  # ============================================================================
  # Pre-allocate output matrix
  # ============================================================================

  output <- matrix(0, nrow = sim.year.num, ncol = sim.num)
  variance_corrections <- logical(ncomp)

  # ============================================================================
  # Component Simulation Loop
  # ============================================================================

  for (k in seq_len(ncomp)) {

    component <- as.numeric(comp_list[[k]])
    n_obs <- length(component)

    # constant component check is now safe because NA is excluded above
    comp_sd <- stats::sd(component)
    if (!is.finite(comp_sd)) {
      stop("Component ", k, " has non-finite standard deviation.", call. = FALSE)
    }

    if (comp_sd < 1e-10) {
      warning(
        "Component ", k, " is essentially constant. ",
        "Returning constant values without ARIMA modeling.",
        call. = FALSE
      )
      output <- output + mean(component)
      next
    }

    if (n_obs < 10) {
      warning(
        "Component ", k, " has only ", n_obs,
        " observations. ARIMA modeling may be unreliable.",
        call. = FALSE
      )
    }

    comp_mean <- mean(component)
    centered <- component - comp_mean
    target_sd <- comp_sd

    if (!is.null(seed)) set.seed(as.integer(seed) + k * 1000L)

    MODEL <- tryCatch(
      {
        suppressPackageStartupMessages({
          forecast::auto.arima(
            centered,
            max.p = 2,
            max.q = 2,
            max.P = 0,
            max.Q = 0,
            stationary = TRUE,
            seasonal = FALSE
          )
        })
      },
      error = function(e) {
        warning(
          "ARIMA fitting failed for component ", k, ": ", e$message,
          "\nFalling back to resampling with replacement.",
          call. = FALSE
        )
        NULL
      }
    )

    if (is.null(MODEL)) {
      for (j in seq_len(sim.num)) {
        output[, j] <- output[, j] + sample(component, sim.year.num, replace = TRUE)
      }
      next
    }

    if (check.diagnostics) {
      if (!is.null(MODEL$convergence) && MODEL$convergence != 0) {
        warning("ARIMA model for component ", k, " did not converge properly.", call. = FALSE)
      }

      if (length(MODEL$residuals) > 1) {
        ljung_test <- tryCatch(
          stats::Box.test(MODEL$residuals, type = "Ljung-Box",
                          lag = min(10, length(MODEL$residuals) - 1)),
          error = function(e) NULL
        )

        if (!is.null(ljung_test) && is.finite(ljung_test$p.value) && ljung_test$p.value < 0.05) {
          warning(
            "Component ", k, ": Residuals show significant autocorrelation ",
            "(Ljung-Box p = ", round(ljung_test$p.value, 4), "). ",
            "Model may be inadequate.",
            call. = FALSE
          )
        }
      }
    }

    intercept <- 0
    if ("intercept" %in% names(MODEL$coef)) {
      intercept <- MODEL$coef[["intercept"]]
    } else if ("drift" %in% names(MODEL$coef)) {
      intercept <- MODEL$coef[["drift"]]
    }

    model_sd <- sqrt(MODEL$sigma2)
    needs_variance_check <- isTRUE(match.variance)

    for (j in seq_len(sim.num)) {

      simulated <- as.numeric(stats::simulate(MODEL, sim.year.num, sd = model_sd)) +
        intercept + comp_mean

      if (needs_variance_check) {
        sim_sd <- stats::sd(simulated)

        # If sim_sd is 0/non-finite, rescaling is impossible; skip correction safely
        if (is.finite(sim_sd) && sim_sd > 0) {
          relative_diff <- abs(sim_sd - target_sd) / target_sd
          if (relative_diff > variance.tolerance) {
            sim_mean <- mean(simulated)
            simulated <- comp_mean + (simulated - sim_mean) * (target_sd / sim_sd)
            if (!variance_corrections[k]) variance_corrections[k] <- TRUE
          }
        }
      }

      output[, j] <- output[, j] + simulated
    }
  }

  if (match.variance && any(variance_corrections) && requireNamespace("logger", quietly = TRUE)) {
    corrected_comps <- which(variance_corrections)
    logger::log_info(
      "[WARM] Variance correction applied to component(s): {paste(corrected_comps, collapse = ', ')}"
    )
  }

  output
}
