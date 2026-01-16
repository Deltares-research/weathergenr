test_that("simulate_warm validates required inputs", {
  suppressMessages(suppressWarnings({

    skip_if_not(exists("simulate_warm", where = asNamespace("weathergenr")))

    expect_error(
      weathergenr::simulate_warm(components = NULL, n = 10),
      "must not be NULL",
      fixed = FALSE
    )

    expect_error(
      weathergenr::simulate_warm(
        components = list(a = rnorm(20)),
        n = NULL,
        verbose = FALSE
      ),
      "must be specified",
      fixed = FALSE
    )

  }))
})

test_that("simulate_warm validates scalar numeric arguments", {
  suppressMessages(suppressWarnings({

    comps <- list(a = rnorm(20))

    expect_error(
      weathergenr::simulate_warm(comps, n = "10", verbose = FALSE),
      "n.*positive integer",
      fixed = FALSE
    )

    expect_error(
      weathergenr::simulate_warm(comps, n = 0, verbose = FALSE),
      "n.*positive integer",
      fixed = FALSE
    )

    expect_error(
      weathergenr::simulate_warm(comps, n = 10, n_sim = 0, verbose = FALSE),
      "n_sim.*positive integer",
      fixed = FALSE
    )

    expect_error(
      weathergenr::simulate_warm(comps, n = 10, match_variance = 1, verbose = FALSE),
      "match_variance.*TRUE/FALSE",
      fixed = FALSE
    )

    expect_error(
      weathergenr::simulate_warm(comps, n = 10, var_tol = -0.1, verbose = FALSE),
      "var_tol.*between 0 and 1",
      fixed = FALSE
    )

    expect_error(
      weathergenr::simulate_warm(comps, n = 10, var_tol = 1.1, verbose = FALSE),
      "var_tol.*between 0 and 1",
      fixed = FALSE
    )

  }))
})

test_that("simulate_warm rejects unsupported components types", {
  suppressMessages(suppressWarnings({

    expect_error(
      weathergenr::simulate_warm(components = 1:10, n = 10, verbose = FALSE),
      "must be a matrix, data\\.frame, or list",
      fixed = FALSE
    )

  }))
})

test_that("simulate_warm accepts matrix, data.frame, and list inputs and returns correct dimensions", {
  suppressMessages(suppressWarnings({

    sim_year <- 25
    sim_num  <- 7

    comps_list <- list(
      signal = sin(2 * pi * (1:30) / 8) + rnorm(30, 0, 0.1),
      noise  = rnorm(30, 0, 0.5)
    )

    comps_mat <- cbind(comps_list$signal, comps_list$noise)
    comps_df  <- data.frame(signal = comps_list$signal, noise = comps_list$noise)

    out_list <- weathergenr::simulate_warm(
      comps_list, n = sim_year, n_sim = sim_num, seed = 1, verbose = FALSE
    )
    out_mat  <- weathergenr::simulate_warm(
      comps_mat,  n = sim_year, n_sim = sim_num, seed = 1, verbose = FALSE
    )
    out_df   <- weathergenr::simulate_warm(
      comps_df,   n = sim_year, n_sim = sim_num, seed = 1, verbose = FALSE
    )

    expect_true(is.matrix(out_list))
    expect_true(is.numeric(out_list))
    expect_equal(dim(out_list), c(sim_year, sim_num))

    expect_true(is.matrix(out_mat))
    expect_equal(dim(out_mat), c(sim_year, sim_num))

    expect_true(is.matrix(out_df))
    expect_equal(dim(out_df), c(sim_year, sim_num))

  }))
})

test_that("simulate_warm is deterministic given seed", {
  suppressMessages(suppressWarnings({

    comps <- list(
      signal = sin(2 * pi * (1:40) / 10) + rnorm(40, 0, 0.1),
      noise  = rnorm(40, 0, 0.3)
    )

    out1 <- weathergenr::simulate_warm(comps, n = 30, n_sim = 5, seed = 123, verbose = FALSE)
    out2 <- weathergenr::simulate_warm(comps, n = 30, n_sim = 5, seed = 123, verbose = FALSE)
    out3 <- weathergenr::simulate_warm(comps, n = 30, n_sim = 5, seed = 124, verbose = FALSE)

    expect_identical(out1, out2)
    expect_false(isTRUE(all.equal(out1, out3)))

  }))
})

test_that("simulate_warm restores .Random.seed when seed is provided", {
  suppressMessages(suppressWarnings({

    comps <- list(a = rnorm(30), b = rnorm(30))

    set.seed(999)
    seed_before <- .Random.seed

    invisible(
      weathergenr::simulate_warm(comps, n = 20, n_sim = 3, seed = 1, verbose = FALSE)
    )

    expect_identical(.Random.seed, seed_before)

    set.seed(999)
    x1 <- runif(5)

    set.seed(999)
    invisible(
      weathergenr::simulate_warm(comps, n = 20, n_sim = 3, seed = 1, verbose = FALSE)
    )
    x2 <- runif(5)

    expect_identical(x1, x2)

  }))
})

test_that("constant components skip ARIMA and add constant mean with warning", {

  comps <- list(constant = rep(5, 30))

  expect_warning(
    out <- weathergenr::simulate_warm(comps, n = 10, n_sim = 4, seed = 7, verbose = FALSE),
    "essentially constant",
    fixed = FALSE
  )

  expect_equal(dim(out), c(10, 4))
  expect_true(all(out == 5))
})

test_that("short components (<10 obs) emit a warning", {

  comps <- list(short = rnorm(9))

  expect_warning(
    out <- weathergenr::simulate_warm(comps, n = 12, n_sim = 3, seed = 11, verbose = FALSE),
    "only 9 observations",
    fixed = TRUE
  )

  expect_equal(dim(out), c(12, 3))
})

test_that("NA inside a component causes an error (regression guard)", {
  suppressMessages(suppressWarnings({

    expect_error(
      weathergenr::simulate_warm(
        list(bad = c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)),
        n = 10,
        n_sim = 2,
        seed = 1, verbose = FALSE
      ),
      "Missing values detected in component\\(s\\)",
      fixed = FALSE
    )

  }))
})

test_that("non-integer n fails validation (regression guard)", {
  suppressMessages(suppressWarnings({

    expect_error(
      weathergenr::simulate_warm(
        list(a = rnorm(30)),
        n = 10.5,
        n_sim = 2,
        seed = 1, verbose = FALSE
      ),
      "n.*positive integer",
      fixed = FALSE
    )

  }))
})
