# Run weathergenr end-to-end (generate and evaluate)

Top-level convenience wrapper that runs:
[`generate_weather()`](https://deltares-research.github.io/weathergenr/reference/generate_weather.md)
[`prepare_evaluation_data()`](https://deltares-research.github.io/weathergenr/reference/prepare_evaluation_data.md)
[`evaluate_weather_generator()`](https://deltares-research.github.io/weathergenr/reference/evaluate_weather_generator.md)

All execution and evaluation settings are taken from `config`. Optional
logging writes console output to a timestamped file in `out_dir` while
continuing to display output on the console.

## Usage

``` r
run_weather_generator(
  obs_data,
  obs_grid,
  obs_dates,
  out_dir,
  config,
  eval_max_grids = 25L,
  log_messages = FALSE
)
```

## Arguments

- obs_data:

  Observed data (e.g. `ncdata$data`).

- obs_grid:

  Observed grid metadata (e.g. `ncdata$grid`).

- obs_dates:

  Observed dates (e.g. `ncdata$date`).

- out_dir:

  Character. Output directory.

- config:

  List. Full simulation/evaluation configuration.

- eval_max_grids:

  Integer. Maximum number of grids to evaluate.

- log_messages:

  Logical. If TRUE, save console output to `log_YYYYMMDD_HHMMSS.txt` in
  `out_dir`.

## Value

A list with components:

- `gen_output`: output of
  [`generate_weather()`](https://deltares-research.github.io/weathergenr/reference/generate_weather.md)

- `evaluation`: output of
  [`evaluate_weather_generator()`](https://deltares-research.github.io/weathergenr/reference/evaluate_weather_generator.md)

- `log_path`: path to the log file (or NULL if `log_messages=FALSE`)
