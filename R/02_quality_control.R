#!/usr/bin/env Rscript

# =============================================================================
# Quality Control Filtering
# =============================================================================
# Combines DARQ and AssemblyNet quality control metrics
#
# Inputs:
#   - data/lists/darq_out_20240408_lng_r152.csv: DARQ QC scores
#   - data/fst/ukb_assemblynet.fst: AssemblyNet QC
#
# Outputs:
#   - data/derivatives/qc_darq-assemblynet.rds: Combined QC metrics
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
library(fst)
library(stringr)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("02_quality_control.R")

# Load configuration
config <- load_config()

# ----- Input Files -----
log_section("Checking requirements")
darq_file <- get_data_path("raw", "darq_qc")
assemblynet_file <- get_data_path("raw", "assemblynet_fst")

check_files_exist(c(darq_file, assemblynet_file))

# ----- DARQ Quality Control -----
log_section("Processing DARQ QC scores")
darq.dt <- read_csv_safe(darq_file, description = "DARQ QC scores")

log_info("Extracting EID and instance from filenames")
darq.dt[, EID := as.integer(stringr::str_extract(fn, "(?<=-)\\d{7}"))]
darq.dt[, INSTANCE := str_extract(fn, "ses-\\d")]
darq.dt[, fn := NULL]

darq.dt |>
  setcolorder(2:3) |>
  setnames("darq", "DARQ")

# Validate DARQ data
validate_columns(darq.dt, c("EID", "INSTANCE", "DARQ"), "DARQ")
validate_not_empty(darq.dt, "DARQ")
validate_num_range(
  darq.dt, "DARQ",
  min_val = 0, max_val = 1, na_allowed = TRUE
)
log_info(
  "DARQ data: %d rows, mean DARQ = %.3f",
  nrow(darq.dt),
  mean(darq.dt$DARQ, na.rm = TRUE)
)

# ----- AssemblyNet Quality Control -----
log_section("Processing AssemblyNet QC scores")
asnet.dt <- read_fst_safe(
  assemblynet_file,
  as_data_table = TRUE, description = "AssemblyNet covariates"
) |>
  (\(DT) DT[, c(1:2, 8)])() |>
  setnames(c("EID", "INSTANCE", "ASBLYNET"))

log_info("Cleaning EID format")
asnet.dt[, EID := as.integer(str_remove(EID, "sub-"))]

# Validate AssemblyNet data
validate_columns(asnet.dt, c("EID", "INSTANCE", "ASBLYNET"), "AssemblyNet")
validate_not_empty(asnet.dt, "AssemblyNet")
log_info("AssemblyNet data: %d rows", nrow(asnet.dt))

# Log QC grade distribution
if ("ASBLYNET" %in% names(asnet.dt)) {
  qc_counts <- asnet.dt[, .N, by = ASBLYNET]
  log_info("AssemblyNet QC distribution:")
  for (i in seq_len(nrow(qc_counts))) {
    log_info("  Grade %s: %d", qc_counts[i, ASBLYNET], qc_counts[i, N])
  }
}

# ----- Combine QC Metrics -----
log_section("Merging QC metrics")
qc.dt <- asnet.dt[darq.dt, on = .(EID, INSTANCE)]

validate_not_empty(qc.dt, "combined QC")
log_info("Combined QC data: %d rows", nrow(qc.dt))

# Check for missing data
n_missing_darq <- qc.dt[is.na(DARQ), .N]
n_missing_asnet <- qc.dt[is.na(ASBLYNET), .N]
if (n_missing_darq > 0) log_warn("Missing DARQ scores: %d rows", n_missing_darq)
if (n_missing_asnet > 0) {
  log_warn(
    "Missing AssemblyNet scores: %d rows", n_missing_asnet
  )
}

rm(darq.dt, asnet.dt, darq_file, assemblynet_file)

# ----- Output -----
log_section("Saving combined QC data")
output <- get_data_path("processed", "qc_combined")
write_rds_safe(qc.dt, output, "Combined QC metrics")

rm(qc.dt, output)
log_script_end("02_quality_control.R", success = TRUE)
