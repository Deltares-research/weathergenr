# ==============================================================================
# Record the end-to-end numeric baseline
# ==============================================================================
#
# Repo-only operator script. Recording a baseline BLESSES whatever the package
# currently produces, so it is a deliberate act, never a chore: run it only
# after you have decided that a numeric change is correct and intended.
#
# Usage (from the package root):
#
#   Rscript tools/record_baseline.R              # record, refusing to overwrite
#   Rscript tools/record_baseline.R --force      # re-record over an existing file
#   Rscript tools/record_baseline.R --dry-run    # show the diff, write nothing
#
# Checking a baseline is the other half and lives in the test suite:
#
#   WEATHERGENR_BASELINE=1 Rscript -e 'devtools::test(filter = "baseline-e2e")'
#
# Runtime is roughly 40 s: two scenarios, each a full generate + evaluate pass
# over inst/extdata/ntoum_era5_data.nc.

args     <- commandArgs(trailingOnly = TRUE)
force    <- "--force"   %in% args
dry_run  <- "--dry-run" %in% args

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# The helper is the single definition of scenarios and fingerprinting; the
# recorder must never grow its own copy.
source(file.path("tests", "testthat", "helper-baseline.R"))

# helper-baseline.R reaches for testthat::test_path() to locate the artifact,
# which resolves differently outside a test run. Pin it explicitly here.
out_path <- file.path("tests", "testthat", "baseline", "e2e-baseline.csv")

if (file.exists(out_path) && !force && !dry_run) {
  stop("Baseline already exists at ", out_path, ".\n",
       "Re-record only when a numeric change is intended and reviewed, ",
       "then pass --force.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Provenance
# ------------------------------------------------------------------------------

git_val <- function(cmd) {
  out <- suppressWarnings(
    tryCatch(system2("git", cmd, stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
  )
  if (!length(out) || is.na(out[1])) NA_character_ else out[1]
}

fixture <- baseline_fixture_path()

provenance <- data.frame(
  scenario = "provenance",
  key = c("recorded_at", "package_version", "r_version", "platform",
          "git_commit", "git_branch", "git_dirty",
          "fixture_file", "fixture_md5", "signif_digits"),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    as.character(utils::packageVersion("weathergenr")),
    R.version.string,
    R.version$platform,
    git_val(c("rev-parse", "HEAD")),
    git_val(c("rev-parse", "--abbrev-ref", "HEAD")),
    as.character(nzchar(paste(
      suppressWarnings(system2("git", c("status", "--porcelain"),
                               stdout = TRUE, stderr = FALSE)),
      collapse = ""))),
    basename(fixture),
    unname(tools::md5sum(fixture)),
    as.character(BASELINE_SIGNIF)
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Run every scenario
# ------------------------------------------------------------------------------

ncdata <- read_netcdf(nc_path = fixture)
scenarios <- baseline_scenarios()

fingerprints <- list()
for (nm in names(scenarios)) {
  message("Recording scenario: ", nm, " ...")
  t0 <- Sys.time()
  res <- baseline_run_scenario(scenarios[[nm]], ncdata)
  fingerprints[[nm]] <- baseline_fingerprint(nm, res, cfg = scenarios[[nm]])
  message(sprintf("  done in %.1f s (%d keys)",
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  nrow(fingerprints[[nm]])))
}

# Scenario configuration is emitted by baseline_fingerprint() itself, so the
# recorder and the checker cannot disagree about it.
new_baseline <- rbind(provenance, do.call(rbind, fingerprints))
rownames(new_baseline) <- NULL

# ------------------------------------------------------------------------------
# Report the delta before writing
# ------------------------------------------------------------------------------

old <- baseline_read(out_path)

if (is.null(old)) {
  message("\nNo existing baseline -- this will be the first record.")
} else {
  old_fp <- baseline_split(old)$fingerprint
  new_fp <- baseline_split(new_baseline)$fingerprint

  merged <- merge(old_fp, new_fp, by = c("scenario", "key"),
                  all = TRUE, suffixes = c(".old", ".new"))
  changed <- merged[is.na(merged$value.old) | is.na(merged$value.new) |
                      merged$value.old != merged$value.new, , drop = FALSE]

  if (!nrow(changed)) {
    message("\nNo change: the current output matches the stored baseline.")
  } else {
    message("\n", nrow(changed), " of ", nrow(merged), " fingerprint keys differ:")
    for (i in seq_len(min(nrow(changed), 40L))) {
      message(sprintf("  %-16s %-38s %s -> %s",
                      changed$scenario[i], changed$key[i],
                      changed$value.old[i], changed$value.new[i]))
    }
    if (nrow(changed) > 40L) message("  ... and ", nrow(changed) - 40L, " more")
  }
}

if (dry_run) {
  message("\n--dry-run: nothing written.")
} else {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(new_baseline, out_path, row.names = FALSE)
  message("\nWrote ", out_path, " (", nrow(new_baseline), " rows).")
  message("Review the diff with `git diff -- ", out_path, "` before committing.")
}
