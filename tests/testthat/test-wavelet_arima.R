test_that("wavelet_arima validates required inputs", {
  suppressMessages(suppressWarnings({

    skip_if_not(exists("wavelet_arima", where = asNamespace("weathergenr")))

    expect_error(
      weathergenr::wavelet_arima(wavelet.components = NULL, sim.year.num = 10),
      "must not be NULL",
      fixed = FALSE
    )

    expect_error(
      weathergenr::wavelet_arima(
        wavelet.components = list(a = rnorm(20)),
        sim.year.num = NULL,
        verbose = FALSE
      ),
      "must be specified",
      fixed = FALSE
    )

  }))
})

test_that("wavelet_arima validates scalar numeric arguments", {
  suppressMessages(suppressWarnings({

    comps <- list(a = rnorm(20))

    expect_error(
      weathergenr::wavelet_arima(comps, sim.year.num = "10", verbose = FALSE),
      "sim.year.num.*positive integer",
      fixed = FALSE
    )

    expect_error(
      weathergenr::wavelet_arima(comps, sim.year.num = 0, verbose = FALSE),
      "sim.year.num.*positive integer",
      fixed = FALSE
    )

    expect_error(
      weathergenr::wavelet_arima(comps, sim.year.num = 10, sim.num = 0, verbose = FALSE),
      "sim.num.*positive integer",
      fixed = FALSE
    )

    expect_error(
      weathergenr::wavelet_arima(comps, sim.year.num = 10, match.variance = 1, verbose = FALSE),
      "match\\.variance.*logical",
      fixed = FALSE
    )

    expect_error(
      weathergenr::wavelet_arima(comps, sim.year.num = 10, variance.tolerance = -0.1, verbose = FALSE),
      "variance\\.tolerance.*between 0 and 1",
      fixed = FALSE
    )

    expect_error(
      weathergenr::wavelet_arima(comps, sim.year.num = 10, variance.tolerance = 1.1, verbose = FALSE),
      "variance\\.tolerance.*between 0 and 1",
      fixed = FALSE
    )

  }))
})

test_that("wavelet_arima rejects unsupported wavelet.components types", {
  suppressMessages(suppressWarnings({

    expect_error(
      weathergenr::wavelet_arima(wavelet.components = 1:10, sim.year.num = 10, verbose = FALSE),
      "must be a matrix, data\\.frame, or list",
      fixed = FALSE
    )

  }))
})

test_that("wavelet_arima accepts matrix, data.frame, and list inputs and returns correct dimensions", {
  suppressMessages(suppressWarnings({

    sim_year <- 25
    sim_num  <- 7

    comps_list <- list(
      signal = sin(2 * pi * (1:30) / 8) + rnorm(30, 0, 0.1),
      noise  = rnorm(30, 0, 0.5)
    )

    comps_mat <- cbind(comps_list$signal, comps_list$noise)
    comps_df  <- data.frame(signal = comps_list$signal, noise = comps_list$noise)

    out_list <- weathergenr::wavelet_arima(
      comps_list, sim.year.num = sim_year, sim.num = sim_num, seed = 1, verbose = FALSE
    )
    out_mat  <- weathergenr::wavelet_arima(
      comps_mat,  sim.year.num = sim_year, sim.num = sim_num, seed = 1, verbose = FALSE
    )
    out_df   <- weathergenr::wavelet_arima(
      comps_df,   sim.year.num = sim_year, sim.num = sim_num, seed = 1, verbose = FALSE
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

test_that("wavelet_arima is deterministic given seed", {
  suppressMessages(suppressWarnings({

    comps <- list(
      signal = sin(2 * pi * (1:40) / 10) + rnorm(40, 0, 0.1),
      noise  = rnorm(40, 0, 0.3)
    )

    out1 <- weathergenr::wavelet_arima(comps, sim.year.num = 30, sim.num = 5, seed = 123, verbose = FALSE)
    out2 <- weathergenr::wavelet_arima(comps, sim.year.num = 30, sim.num = 5, seed = 123, verbose = FALSE)
    out3 <- weathergenr::wavelet_arima(comps, sim.year.num = 30, sim.num = 5, seed = 124, verbose = FALSE)

    expect_identical(out1, out2)
    expect_false(isTRUE(all.equal(out1, out3)))

  }))
})

test_that("wavelet_arima restores .Random.seed when seed is provided", {
  suppressMessages(suppressWarnings({

    comps <- list(a = rnorm(30), b = rnorm(30))

    set.seed(999)
    seed_before <- .Random.seed

    invisible(
      weathergenr::wavelet_arima(comps, sim.year.num = 20, sim.num = 3, seed = 1, verbose = FALSE)
    )

    expect_identical(.Random.seed, seed_before)

    set.seed(999)
    x1 <- runif(5)

    set.seed(999)
    invisible(
      weathergenr::wavelet_arima(comps, sim.year.num = 20, sim.num = 3, seed = 1, verbose = FALSE)
    )
    x2 <- runif(5)

    expect_identical(x1, x2)

  }))
})

test_that("constant components skip ARIMA and add constant mean with warning", {

  comps <- list(constant = rep(5, 30))

  expect_warning(
    out <- weathergenr::wavelet_arima(comps, sim.year.num = 10, sim.num = 4, seed = 7, verbose = FALSE),
    "essentially constant",
    fixed = FALSE
  )

  expect_equal(dim(out), c(10, 4))
  expect_true(all(out == 5))
})

test_that("short components (<10 obs) emit a warning", {

  comps <- list(short = rnorm(9))

  expect_warning(
    out <- weathergenr::wavelet_arima(comps, sim.year.num = 12, sim.num = 3, seed = 11, verbose = FALSE),
    "only 9 observations",
    fixed = TRUE
  )

  expect_equal(dim(out), c(12, 3))
})

test_that("NA inside a component currently causes a logical NA failure (regression guard)", {
  suppressMessages(suppressWarnings({

    expect_error(
      weathergenr::wavelet_arima(
        list(bad = c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)),
        sim.year.num = 10,
        sim.num = 2,
        seed = 1, verbose = FALSE
      ),
      "Missing values detected in wavelet component\\(s\\)",
      fixed = FALSE
    )

  }))
})

test_that("non-integer sim.year.num currently fails downstream (regression guard)", {
  suppressMessages(suppressWarnings({

    expect_error(
      weathergenr::wavelet_arima(
        list(a = rnorm(30)),
        sim.year.num = 10.5,
        sim.num = 2,
        seed = 1, verbose = FALSE
      ),
      "sim\\.year\\.num.*positive integer",
      fixed = FALSE
    )

  }))
})
