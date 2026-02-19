#!/usr/bin/env Rscript

# =============================================================================
# Parse Preprocessed Brain Volumes
# =============================================================================
# Aggregates brain volume data from preprocessed segmentation files
#
# Inputs:
#   - data/raw/pp_vols/*.csv: Individual participant volume files
#
# Outputs:
#   - data/raw/ukb-lng_icc-scale.csv: Aggregated ICC volumes and scale factors
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("01_parse_volumes.R")

# Load configuration
config <- load_config()

# ----- Input/Output Paths -----
preproc_dir <- get_data_path("raw", "pp_vols_dir")
log_info("Checking for preprocessed volumes directory")

if (!dir.exists(preproc_dir)) {
  log_error("Preprocessed volumes directory not found: %s", preproc_dir)
  stop("Missing required directory: ", preproc_dir, call. = FALSE)
}

log_info("Found preprocessed volumes directory")

# ----- Main Processing -----
output <- get_data_path("processed", "icc_scale")

if (!check_files_exist(output, stop_on_missing = FALSE)) {
  log_section("Aggregating volume data")

  columns.v <- c("SubjectID", "VisitID", "ICC_vol", "ScaleFactor")

  csv_files <- list.files(preproc_dir, pattern = "\\.csv$", full.names = TRUE)
  log_info("Found %d volume files to process", length(csv_files))

  if (length(csv_files) == 0) {
    log_error("No CSV files found in %s", preproc_dir)
    stop("No volume files to process", call. = FALSE)
  }

  log_time("Reading and aggregating volume files", {
    volumes.dt <- csv_files |>
      lapply(fread, select = columns.v, encoding = "UTF-8") |>  # BP7: explicit encoding
      rbindlist()
  })

  # Validate aggregated data
  log_section("Validating aggregated data")
  validate_columns(volumes.dt, columns.v, "volumes")
  validate_not_empty(volumes.dt, "volumes")
  log_info("Aggregated data has %d rows", nrow(volumes.dt))

  # Write output
  write_csv_safe(volumes.dt, output, "Aggregated ICC volumes")

  rm(volumes.dt, columns.v, csv_files)
} else {
  log_info("Output already exists, skipping: %s", basename(output))
}

# ----- Cleanup -----
rm(preproc_dir, output)
log_script_end("01_parse_volumes.R", success = TRUE)
