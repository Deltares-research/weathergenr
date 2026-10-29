test_that("R planner agrees with the language-agnostic measured-width fixtures", {
  fixtures <- utils::read.csv(test_path("fixtures", "panel-layout.csv"))
  for (i in seq_len(nrow(fixtures))) {
    row <- fixtures[i, ]
    call <- function() .choose_panel_grid(row$panel_count, c(3, 2),
      c(row$width_3_mm, row$width_2_mm), row$max_width_mm)
    if (row$expected_columns == 0) {
      expect_error(call(), "do not fit")
    } else {
      expect_equal(unname(call()), c(row$expected_columns, row$expected_rows),
                   info = row$case)
    }
  }
  expect_equal(.choose_panel_grid(1, 1, 90, 180), c(columns = 1, rows = 1))
  expect_equal(.choose_panel_grid(2, c(3, 2), c(250, 170), 180),
               c(columns = 2, rows = 1))
})

test_that("physical panel layouts preserve size, type and input plots", {
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point() + theme_weathergenr(base_size = 7)
  for (count in c(1, 2, 3, 4, 9)) {
    width <- if (count == 9) 270 else 180
    result <- arrange_weather_panels(rep(list(p), count), max_width_mm = width)
    layout <- attr(result, "weathergenr_panel_layout")
    expect_equal(layout$columns, if (count == 1) 1 else if (count == 9) 3 else 2)
    expect_equal(layout$rows, ceiling(count / layout$columns))
    expect_lte(layout$canvas_mm[["width"]], width + 1e-6)
    expect_equal(layout$panel_width_mm, 70)
    expect_equal(layout$panel_height_mm, 70)
    expect_equal(p$theme$text$size, 7)
  }
  expect_error(arrange_weather_panels(list(p, p), max_width_mm = 100), "do not fit")
  expect_error(arrange_weather_panels(list()), "nonempty")
  expect_error(arrange_weather_panels(list(p), preferred_columns = 1), "at least two")
})

test_that("compact style is opt-in and preserves annual values, RNG and defaults", {
  sim <- expand.grid(year = 1:6, rlz = 1:2)
  sim$stat <- "mean"
  sim$variable <- "precip"
  sim$Simulated <- seq_len(nrow(sim)) / 2
  obs <- data.frame(year = 1:6, stat = "mean", variable = "precip", Observed = 1:6)
  original <- .create_annual_precip_plot(sim, obs, list(), FALSE, FALSE, NULL)
  before <- ggplot2::ggplot_build(original)
  aes_before <- lapply(original$layers, function(layer) layer$aes_params)
  .set_seed_fixed_kind(42)
  seed <- .Random.seed
  result <- style_weather_plot(original, "annual_precip")
  expect_identical(.Random.seed, seed)
  expect_identical(lapply(original$layers, function(layer) layer$aes_params), aes_before)
  after <- ggplot2::ggplot_build(result)
  for (i in seq_along(before$data)) {
    expect_identical(before$data[[i]][c("x", "y", "group")],
                     after$data[[i]][c("x", "y", "group")])
  }
  expect_equal(result$layers[[1]]$aes_params$alpha, 0.6)
  expect_equal(result$layers[[2]]$aes_params$linewidth, 0.45)
  expect_equal(result$labels$x, "Serial year")
  expect_error(style_weather_plot(original, "unknown"), "stable plot ID")
  short <- .create_annual_precip_plot(sim[sim$year <= 2, ], obs[obs$year <= 2, ],
    list(), FALSE, FALSE, NULL)
  expect_no_error(style_weather_plot(short, "annual_precip"))
  empty <- .create_annual_precip_plot(sim[sim$year == 1, ], obs[obs$year == 1, ],
    list(), FALSE, FALSE, NULL)
  expect_no_error(style_weather_plot(empty, "annual_precip"))
})

test_that("compact facets handle four variables and keep singleton units", {
  daily <- expand.grid(variable = c("precip", "temp", "temp_min", "temp_max"),
    mon = 1:3, rlz = 1:2)
  daily$stat <- "mean"
  daily$Observed <- daily$mon + as.integer(daily$variable)
  daily$Simulated <- daily$Observed + daily$rlz / 10
  original <- .create_daily_mean_plot(daily, list(alpha = 0.3), FALSE, FALSE, NULL)
  result <- style_weather_plot(original, "daily_mean")
  expect_equal(attr(result, "weathergenr_panel_layout")$rows, 2)
  expect_equal(attr(result, "weathergenr_panel_layout")$columns, 2)
  expect_equal(.diagnostic_panel_letter(27), "aa")
  singleton <- daily[daily$variable == "temp_min", ]
  p <- .create_daily_mean_plot(singleton, list(alpha = 0.3), FALSE, FALSE, NULL)
  styled <- style_weather_plot(p, "daily_mean", plot_style = list(
    theme = ggplot2::theme(strip.text = ggplot2::element_text(size = 9))))
  expect_match(styled[[1]]$labels$x, "Observed daily mean")
  expect_match(styled[[1]]$labels$x, "\u00b0C", fixed = TRUE)
  expect_s3_class(styled[[1]]$theme$strip.text, "element_blank")
  monthly <- .create_monthly_cycle_plot(singleton, list(), FALSE, FALSE, NULL)
  styled <- style_weather_plot(monthly, "monthly_cycle")
  expect_equal(styled[[1]]$labels$y, "Minimum temperature (\u00b0C)")
  simulated <- expand.grid(variable = "temp_min", rlz = 1:2,
    year = 2001:2004, mon = 1:12)
  simulated$stat <- "mean"
  simulated$Simulated <- simulated$mon + simulated$rlz + simulated$year / 1000
  observed <- simulated[simulated$rlz == 1, ]
  observed$rlz <- NULL
  observed$Observed <- observed$Simulated - 0.2
  observed$Simulated <- NULL
  p <- .create_monthly_pattern_plot(simulated, observed, "temp_min",
    list(colors = c(Observed = "blue", Simulated = "black")),
    FALSE, FALSE, NULL)
  styled <- style_weather_plot(p, "annual_pattern_temp_min")
  expect_equal(styled[[1]]$labels$y, "Mean (\u00b0C)")
})
