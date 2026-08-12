# ==============================================================================
# End-to-end numeric baseline
# ==============================================================================
#
# Guards the numeric output of the full generate -> evaluate pipeline against a
# recorded baseline, so a refactor that is supposed to be behaviour-preserving
# can be shown to be behaviour-preserving rather than assumed to be.
#
# This is an OPT-IN LOCAL GATE, not a suite test. It costs roughly 40 s and does
# not run under devtools::test(), R CMD check, or CI. Run it deliberately, both
# before and after a change that touches numeric behaviour:
#
#   WEATHERGENR_BASELINE=1 Rscript -e 'devtools::test(filter = "baseline-e2e")'
#
# A green CI therefore does NOT mean the baseline was checked.
#
# When it fails, read the reported keys before re-recording. A failure is a
# question ("did I mean to change this?"), not a chore. Re-record only after
# deciding the new numbers are correct:
#
#   Rscript tools/record_baseline.R --dry-run   # inspect the delta first
#   Rscript tools/record_baseline.R --force
#
# Scenario definitions and the fingerprint live in helper-baseline.R, which the
# recorder sources too -- there is exactly one definition of both.

skip_baseline_unless_requested <- function() {
  testthat::skip_on_cran()
  # skip_on_ci() keys on the CI env var and is NOT overridable by
  # WEATHERGENR_BASELINE. The CI door is closed deliberately: the exact-hash tier
  # has only ever been verified on one platform (R 4.6.0, x86_64-w64-mingw32), so
  # running it across the 3-OS x 3-R matrix could fail on floating-point
  # differences rather than on a real regression. Remove this line only after
  # confirming the hashes reproduce on another platform.
  testthat::skip_on_ci()
  if (!nzchar(Sys.getenv("WEATHERGENR_BASELINE"))) {
    testthat::skip(paste(
      "End-to-end baseline is opt-in (~40 s).",
      "Set WEATHERGENR_BASELINE=1 to run it."
    ))
  }
  # nolint start: object_usage_linter.
  # baseline_path() lives in helper-baseline.R; lintr does not resolve symbols
  # across testthat helper files, so it reads as undefined here.
  if (!file.exists(baseline_path())) {
    testthat::skip(paste(
      "No baseline recorded at", baseline_path(),
      "-- run: Rscript tools/record_baseline.R"
    ))
  }
  # nolint end
}

# Compare one scenario's fingerprint against the stored rows and build a message
# that names what moved, so a failure localises itself.
compare_fingerprint <- function(current, stored, scenario) {

  merged <- merge(
    stored[stored$scenario == scenario, c("key", "value"), drop = FALSE],
    current[current$scenario == scenario, c("key", "value"), drop = FALSE],
    by = "key", all = TRUE, suffixes = c(".stored", ".current")
  )

  bad <- merged[is.na(merged$value.stored) | is.na(merged$value.current) |
                  merged$value.stored != merged$value.current, , drop = FALSE]

  if (!nrow(bad)) return(NULL)

  # Put the exact-tier hashes first: they say *whether* it moved, the rounded
  # statistics below them say *how*.
  bad <- bad[order(!grepl("hash", bad$key), bad$key), , drop = FALSE]

  lines <- sprintf("  %-42s stored=%s current=%s",
                   bad$key,
                   ifelse(is.na(bad$value.stored), "<absent>", bad$value.stored),
                   ifelse(is.na(bad$value.current), "<absent>", bad$value.current))

  paste0(
    "Scenario '", scenario, "': ", nrow(bad), " of ", nrow(merged),
    " baseline keys differ.\n",
    paste(lines, collapse = "\n"),
    "\n\nIf this change is intended, re-record with:",
    "\n  Rscript tools/record_baseline.R --dry-run",
    "\n  Rscript tools/record_baseline.R --force"
  )
}


test_that("the recorded baseline matches the fixture it was recorded from", {

  skip_baseline_unless_requested()

  stored <- baseline_read()
  prov <- baseline_split(stored)$provenance
  recorded_md5 <- prov$value[prov$key == "fixture_md5"]

  skip_if(!length(recorded_md5) || is.na(recorded_md5),
          "Baseline predates fixture-hash provenance; re-record it.")

  # Checked separately from the numbers: if the fixture itself moved, every
  # downstream key would differ and the real cause would be buried.
  expect_identical(
    unname(tools::md5sum(baseline_fixture_path())),
    recorded_md5,
    info = paste(
      "The fixture NetCDF differs from the one the baseline was recorded",
      "against, so the numeric comparison below is meaningless. Re-record."
    )
  )
})


for (scenario_name in names(baseline_scenarios())) {

  # Force the value so each generated test closes over its own scenario name.
  local({
    nm <- scenario_name

    test_that(paste0("end-to-end output is unchanged (", nm, ")"), {

      skip_baseline_unless_requested()

      stored <- baseline_split(baseline_read())$fingerprint
      skip_if(!any(stored$scenario == nm),
              paste0("Baseline has no rows for scenario '", nm,
                     "' -- re-record with tools/record_baseline.R"))

      cfg <- baseline_scenarios()[[nm]]
      ncdata <- read_netcdf(nc_path = baseline_fixture_path())
      res <- baseline_run_scenario(cfg, ncdata)
      current <- baseline_fingerprint(nm, res, cfg = cfg)

      failure <- compare_fingerprint(current, stored, nm)
      expect_null(failure, info = failure)
    })
  })
}


test_that("parallel execution reproduces the recorded sequential baseline", {

  skip_baseline_unless_requested()

  # Asserted against the SEQUENTIAL water_year rows on purpose. That pins two
  # properties with one comparison: parallel agrees with sequential, and both
  # still agree with the recorded numbers. A separately recorded parallel
  # fingerprint would only ever restate the sequential one.
  nm <- "water_year"
  stored <- baseline_split(baseline_read())$fingerprint
  skip_if(!any(stored$scenario == nm),
          paste0("Baseline has no rows for '", nm, "'."))

  skip_if_not(
    baseline_workers_match_master(2L),
    paste(
      "PSOCK workers are running a different build from the master.",
      "devtools::load_all() does not reach workers -- they load the INSTALLED",
      "package. Run devtools::install() first, or rely on R CMD check, which",
      "installs before testing."
    )
  )

  cfg <- baseline_scenarios()[[nm]]
  ncdata <- read_netcdf(nc_path = baseline_fixture_path())
  res <- baseline_run_scenario(cfg, ncdata, parallel = TRUE, n_cores = 2L)
  current <- baseline_fingerprint(nm, res, cfg = cfg)

  failure <- compare_fingerprint(current, stored, nm)
  expect_null(failure, info = failure)
})


test_that("baseline fingerprinting is itself deterministic", {

  # Cheap guard on the harness rather than the package: if the fingerprint
  # depended on run order, locale, or a timestamp, every future baseline failure
  # would be untrustworthy. Uses a stub, so it costs nothing and always runs.
  fake <- list(
    gen = list(
      resampled = data.frame(
        rlz_1 = as.Date(c("2001-01-01", "2001-01-02", NA)),
        rlz_2 = as.Date(c("2002-03-04", "2002-03-05", "2002-03-06"))
      ),
      dates = as.Date(c("2020-10-01", "2020-10-02", "2020-10-03"))
    ),
    fit_summary = data.frame(
      rlz = c(1L, 2L),
      rank = c(1L, 2L),
      overall_score = c(0.5, 0.5),
      mae_mean_precip = c(2.766108123456, 2.840301987654)
    )
  )

  a <- baseline_fingerprint("stub", fake)
  b <- baseline_fingerprint("stub", fake)

  expect_identical(a, b)
  expect_true(all(c("gen.hash_resampled", "eval.hash_fit_summary") %in% a$key))

  # The NA in rlz_1 must be visible in the fingerprint, not silently dropped:
  # unfilled simulated days are exactly the kind of regression this guards.
  expect_identical(a$value[a$key == "gen.rlz_1.n_na"], "1")

  # Rounding must actually bite, otherwise the "platform-tolerant" claim in
  # helper-baseline.R is false.
  fake_jittered <- fake
  fake_jittered$fit_summary$mae_mean_precip <-
    fake$fit_summary$mae_mean_precip * (1 + 1e-14)
  expect_identical(
    baseline_fingerprint("stub", fake_jittered)$value,
    a$value
  )
})
