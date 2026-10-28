test_that("plot styles reject unknown fields, roles and complete themes", {
  expect_error(.validate_plot_style(list(colours = "red")), "named list")
  expect_error(.validate_plot_style(list(colors = c(other = "red"))), "unique names")
  expect_error(.validate_plot_style(list(colors = c(observed = "nonsense"))), "invalid color")
  expect_error(.validate_plot_style(list(power_palette = NA_character_)), "color strings")
  expect_error(.validate_plot_style(list(theme = ggplot2::theme_bw())), "partial")
  expect_silent(.validate_plot_style(list()))
  expect_identical(.plot_color(list(colors = c(observed = "red")), "simulated", "black"), "black")
})

test_that("public export uses millimetres, device, background and stored geometry", {
  calls <- list()
  local_mocked_bindings(ggsave = function(...) {
    calls[[length(calls) + 1L]] <<- list(...)
  })
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point() + theme_weathergenr()
  p <- .export_figure(p, "unused.png", NULL, family = "wide")
  export_weather_plot(p, "figure.png", width = 254, height = 127,
                      dpi = 96, device = grDevices::png, bg = "transparent")
  expect_equal(c(calls[[1]]$width, calls[[1]]$height), c(10, 5))
  expect_equal(calls[[1]]$dpi, 96)
  expect_identical(calls[[1]]$device, grDevices::png)
  expect_identical(calls[[1]]$bg, "transparent")
  export_weather_plot(p, "figure.png")
  expect_equal(c(calls[[2]]$width, calls[[2]]$height), c(8, 4.5))
  expect_error(export_weather_plot(p, "figure.png", width = NA_real_), "positive")
  expect_error(export_weather_plot(NULL, "figure.png"), "ggplot")
})

test_that("partial themes preserve local layout and explicit font sizes at export", {
  local_mocked_bindings(ggsave = function(...) invisible(NULL))
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point() + theme_weathergenr() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  p <- .style_plot(p, list(theme = ggplot2::theme(
    text = ggplot2::element_text(size = 8, family = "sans"),
    plot.title = ggplot2::element_text(size = 11)
  )))
  out <- export_weather_plot(p, "figure.png", width = 300, height = 200)
  expect_equal(out$theme$text$size, 8)
  expect_equal(out$theme$plot.title$size, 11)
  expect_s3_class(out$theme$axis.text.x, "element_blank")
  changed <- p + ggplot2::theme(text = ggplot2::element_text(size = 7))
  out <- export_weather_plot(changed, "figure.png", width = 100, height = 100)
  expect_equal(out$theme$text$size, 7)
  bare <- ggplot2::ggplot() + ggplot2::theme(text = ggplot2::element_text(size = 6))
  out <- export_weather_plot(bare, "figure.png", width = 100, height = 100)
  expect_equal(out$theme$text$size, 6)
})

test_that("wavelet style reaches both panels and semantic layers", {
  x <- seq_len(32)
  series <- 100 + sin(x / 2)
  w <- analyze_wavelet_spectrum(series, mode = "complete", seed = 1)
  args <- list(series = series, period = w$period, power = w$power,
               gws = w$gws, gws_signif = w$gws_signif,
               coi = w$coi, signif_mask = w$sigm)
  style <- list(theme = ggplot2::theme(text = ggplot2::element_text(size = 9)),
                colors = c(observed = "orange", simulated = "purple", reference = "green"),
                power_palette = c("white", "navy"))
  default <- do.call(plot_wavelet_power, args)
  styled <- suppressMessages(do.call(plot_wavelet_power, c(args, list(plot_style = style))))
  for (i in 1:2) {
    expect_equal(styled[[i]]$theme$text$size, 9)
    expect_equal(styled[[i]]$data, default[[i]]$data)
  }
  expect_identical(styled[[2]]$layers[[1]]$aes_params$colour, "orange")
  expect_identical(styled[[2]]$layers[[2]]$aes_params$colour, "green")
  expect_identical(styled[[1]]$layers[[2]]$aes_params$colour, "green")
  expect_no_error(ggplot2::ggplot_build(styled[[1]]))
  p <- plot_wavelet_global_spectrum(w$period, w$gws_signif, w$gws,
                                    cbind(w$gws, w$gws), plot_style = style)
  expect_identical(p$layers[[3]]$aes_params$colour, "purple")
})
