#!/usr/bin/env Rscript

# =============================================================================
# Sex-Balance Audit of the UKB Normative Pipeline
# =============================================================================
# Reports per-stage sample composition by sex along the upstream pipeline
# (R/02_quality_control.R + R/05_adjust_headsize.R), and explains why the
# final GAMLSS normative sample is 56.2% female.
#
# Stages audited:
#   1. Raw UKB imaging cohort (any ses-2 or ses-3 MRI date)
#   2. DARQ scan-level QC failures (DARQ < darq_threshold)
#   3. Manually-curated bad-EID exclusion list
#   4. Volume-based outlier removal (|Z| > outlier_sd_threshold on
#      residuals of lm(VAL ~ SEX + AGE + ICC))
#   5. AssemblyNet QC failures (grade != assemblynet_rating)
#   6-9. ICD-10 exclusions (F, G, Q0, primary union)
#   10-12. Final GAMLSS sample (all / train / test splits)
#
# Inputs (all from pipeline_config.yaml):
#   - data/fst/ukb_covars.fst
#   - data/derivatives/qc_darq-assemblynet.rds
#   - data/raw/ukb-lng_icc-scale.csv
#   - data/derivatives/hclvag_segmentations.csv
#   - data/lists/exclude_id.txt
#   - models/fits/gamlss.rds
#
# Output:
#   - outputs/audit/sex_balance.rds (per-stage counts + QC failure rates)
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
library(stringr)
library(lubridate)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("audit_sex_balance.R")

config <- load_config()

# ----- Constants -----
DARQ_THRESHOLD    <- get_parameter("qc", "darq_threshold",       default = 0.25)
ASNET_PASS        <- get_parameter("qc", "assemblynet_rating",   default = "A")
OUTLIER_THRESHOLD <- get_parameter("qc", "outlier_sd_threshold", default = 3.0)
ICD_PRIMARY       <- get_parameter("exclusions", "icd10_patterns",
                                   default = "F|G|Q0")

log_info("DARQ threshold: %.2f", DARQ_THRESHOLD)
log_info("AssemblyNet pass grade: %s", ASNET_PASS)
log_info("Outlier threshold: |Z| > %.1f", OUTLIER_THRESHOLD)
log_info("ICD-10 primary exclusion pattern: %s", ICD_PRIMARY)

# ----- Input Files -----
log_section("Checking inputs")
covars_file        <- get_data_path("processed", "covars_fst")
qc_file            <- get_data_path("processed", "qc_combined")
icc_scale_file     <- get_data_path("processed", "icc_scale")
segmentations_file <- get_data_path("processed", "hclvag_segmentations")
exclude_ids_file   <- get_data_path("raw",       "exclude_ids")
gamlss_file        <- get_data_path("processed", "gamlss_fits")

required_files <- c(covars_file, qc_file, icc_scale_file,
                    segmentations_file, exclude_ids_file, gamlss_file)
check_files_exist(required_files, stop_on_missing = TRUE) |> invisible()

# ----- Load Data -----
log_section("Loading data")

covars.dt <- read_fst_safe(covars_file, description = "covariates")
validate_columns(
  covars.dt,
  c("EID", "SEX_r", "DATE_mri_ses2", "DATE_mri_ses3", "ICD_10"),
  "covariates"
)

qc.dt <- read_rds_safe(qc_file, description = "combined QC")
validate_columns(qc.dt, c("EID", "INSTANCE", "DARQ", "ASBLYNET"), "QC")
# Align column name with cnn/asblynet datasets so anti-joins use both keys
setnames(qc.dt, "INSTANCE", "INST")

icc_scl.dt <- read_csv_safe(
  icc_scale_file,
  col.names = c("EID", "INST", "ICC", "SCALE"),
  description = "ICC and scale factors"
)
icc_scl.dt[, EID := as.integer(str_remove(EID, "sub-"))]
icc_scl.dt[, ICC := ICC / 1000]
setkey(icc_scl.dt, EID, INST)

vols.dt <- read_csv_safe(segmentations_file, description = "CNN segmentations")

exclude_ids.v <- as.integer(readLines(exclude_ids_file))
log_info("Bad-EID exclusion list: %d IDs", length(exclude_ids.v))

gamlss.lst <- read_rds_safe(gamlss_file, description = "GAMLSS fits")
final.dt   <- as.data.table(gamlss.lst$DATA$CRS)
validate_columns(final.dt, c("EID", "SEX", "SPLIT"), "GAMLSS sample")
rm(gamlss.lst)

# ----- Prepare CNN Volumes (mirrors R/05_adjust_headsize.R) -----
log_section("Preparing CNN volumes for outlier detection")

vols.dt[, c("l_amy", "r_amy") := NULL]
old_names.v <- names(vols.dt)[-1]
new_names.v <- old_names.v |>
  str_replace("^(.)_(..)_(.)$", "\\2_\\1_\\3") |>
  str_replace_all("vc", "lv") |>
  str_to_upper() |>
  str_replace_all(c("T$" = "tail", "B$" = "body", "H$" = "head"))
setnames(vols.dt, old_names.v, new_names.v)
vols.dt[, (new_names.v) := lapply(.SD, \(v) v / 1000), .SDcols = -1]

# Bilateral totals — these are the columns the outlier filter operates on
vols.dt[, HC   := rowSums(.SD), .SDcols = patterns("^HC")]
vols.dt[, HC_L := rowSums(.SD), .SDcols = patterns("^HC_L")]
vols.dt[, HC_R := rowSums(.SD), .SDcols = patterns("^HC_R")]
vols.dt[, LV   := rowSums(.SD), .SDcols = patterns("^LV")]
vols.dt[, LV_L := rowSums(.SD), .SDcols = patterns("^LV_L")]
vols.dt[, LV_R := rowSums(.SD), .SDcols = patterns("^LV_R")]

setnames(vols.dt, "id", "EID")
vols.dt[, let(
  INST = str_extract(EID, "ses-\\d"),
  EID  = as.integer(str_extract(EID, "(?<=sub-)\\d*"))
)]
setcolorder(vols.dt, "INST", after = "EID")
setkey(vols.dt, EID, INST)
rm(old_names.v, new_names.v)

# Merge ICC/Scale and bring CNN back to subject space
cnn.dt <- icc_scl.dt[vols.dt, nomatch = NULL]
roi_cols.v <- grep("^(HC|LV)", names(cnn.dt), value = TRUE)
cnn.dt[, (roi_cols.v) := .SD / SCALE, .SDcols = roi_cols.v]
rm(vols.dt, icc_scl.dt, roi_cols.v)

# Apply DARQ + bad-EID filters (R/05_adjust_headsize.R line 168)
setkey(qc.dt, EID, INST)
cnn.dt <- cnn.dt[!qc.dt[DARQ < DARQ_THRESHOLD]][!EID %in% exclude_ids.v]
log_info("Post-DARQ + post-bad-EID CNN observations: %d", nrow(cnn.dt))

# ----- Age & Sex (mirrors R/05_adjust_headsize.R lines 173-196) -----
age_sex.dt <- covars.dt[
  !is.na(DATE_mri_ses2) | !is.na(DATE_mri_ses3),
  .(BIRTH_m, BIRTH_y, DATE_mri_ses2, DATE_mri_ses3),
  by = .(EID, SEX = SEX_r)
][
  , BIRTH_my := my(sprintf("%i-%i", match(BIRTH_m, month.name), BIRTH_y))
] |>
  melt(measure = patterns("DATE"), variable = "INST") |>
  (\(dt) {
    dt[, INST := fifelse(INST %like% 2, "ses-2", "ses-3")]
    dt[, value := str_extract(value, "^[^ ]+")]
    dt[, AGE := as.duration(ymd(value) - BIRTH_my) / dyears(1)]
    dt[, c("BIRTH_m", "BIRTH_y", "BIRTH_my", "value") := NULL]
  })() |>
  na.omit() |>
  setkey(EID, INST)

# ----- Stage Counts -----
log_section("Computing per-stage counts by sex")

sex.dt <- covars.dt[, .(EID, SEX = SEX_r)]
setkey(sex.dt, EID)

# Stage 1: raw imaging cohort
imaging.dt <- covars.dt[
  !is.na(DATE_mri_ses2) | !is.na(DATE_mri_ses3),
  .(EID, SEX = SEX_r)
]
imaging_eids.v   <- imaging.dt$EID
stage_imaging.dt <- imaging.dt[, .N, keyby = SEX]
log_info("Raw imaging cohort: %d (%.1f%% female)",
         sum(stage_imaging.dt$N),
         100 * stage_imaging.dt[SEX == "Female", N] / sum(stage_imaging.dt$N))

# Stage 2: DARQ QC failures
qc_sex.dt    <- merge(qc.dt, sex.dt, by = "EID")
darq_fail.dt <- qc_sex.dt[DARQ < DARQ_THRESHOLD, .N, keyby = SEX]

# Stage 3: bad-EID exclusion list
bad_eid.dt <- sex.dt[EID %in% exclude_ids.v, .N, keyby = SEX]

# Stage 4: outlier removal (mirrors R/05_adjust_headsize.R lines 199-219)
cnn_long.dt <- merge(
  melt(cnn.dt, measure = patterns("HC|LV")),
  age_sex.dt,
  by = c("EID", "INST")
)
cnn_long.dt[
  !variable %like% "tail|body|head",
  Z_res := lm(value ~ SEX + AGE + ICC) |> residuals() |> scale(),
  by = variable
]
outlier_ids.v   <- cnn_long.dt[abs(Z_res) > OUTLIER_THRESHOLD, unique(EID)]
outlier_fail.dt <- sex.dt[EID %in% outlier_ids.v, .N, keyby = SEX]
log_info("Volume outliers (|Z| > %.1f): %d", OUTLIER_THRESHOLD,
         length(outlier_ids.v))
rm(cnn_long.dt, cnn.dt, age_sex.dt)

# Stage 5: AssemblyNet QC failures
asnet_fail.dt <- qc_sex.dt[ASBLYNET != ASNET_PASS, .N, keyby = SEX]

# Within-sex QC failure rates
qc_rates.dt <- qc_sex.dt[, .(
  total          = .N,
  darq_fail_pct  = round(100 * mean(DARQ < DARQ_THRESHOLD,  na.rm = TRUE), 2),
  asnet_fail_pct = round(100 * mean(ASBLYNET != ASNET_PASS, na.rm = TRUE), 2)
), keyby = SEX]

# Stages 6-9: ICD-10 exclusions within imaging cohort
icd_dt.fn <- function(pattern.s) {
  covars.dt[
    EID %in% imaging_eids.v & ICD_10 %like% pattern.s,
    .N, keyby = .(SEX = SEX_r)
  ]
}
icd_f.dt   <- icd_dt.fn("^F")
icd_g.dt   <- icd_dt.fn("^G")
icd_q0.dt  <- icd_dt.fn("^Q0")
icd_pri.dt <- icd_dt.fn(ICD_PRIMARY)

# Stages 10-12: final GAMLSS sample
final_all.dt   <- unique(final.dt[, .(EID, SEX)])[, .N, keyby = SEX]
final_train.dt <- unique(final.dt[SPLIT == "train", .(EID, SEX)])[, .N, keyby = SEX]
final_test.dt  <- unique(final.dt[SPLIT == "test",  .(EID, SEX)])[, .N, keyby = SEX]

# ----- Assemble Summary -----
# ----- Cumulative Survivor Counts (anchors the PRISMA boxes) -----
log_section("Computing cumulative survivor counts by sex")

get_n.fn <- function(dt, sex.s) {
  if (nrow(dt) == 0) return(0L)
  v <- dt[SEX == sex.s, N]
  if (length(v) == 0) 0L else as.integer(v)
}

# Build the per-EID worst-instance flags that map onto PRISMA stages.
# Each EID is considered to "fail QC" if any of its sessions fail DARQ or
# AssemblyNet, or if it appears on the bad-EID list, or if its post-DARQ
# residuals flag it as an outlier (using the upstream definition).
cnn_eids.v        <- unique(qc.dt$EID)
asnet_pass_eids.v <- unique(qc.dt[ASBLYNET == ASNET_PASS, EID])
darq_pass_eids.v  <- unique(qc.dt[DARQ >= DARQ_THRESHOLD, EID])

# Anchor stage A: CNN-segmented imaging cohort (matches PRISMA: 47,398).
# Defined as EIDs with QC records and ICC/SCALE data — i.e. those whose
# CNN segmentations exist (the pipeline starting point in
# R/05_adjust_headsize.R after the icc_scl <-> vols inner join).
seg_eids.v <- intersect(cnn_eids.v, imaging_eids.v)
stageA.dt  <- sex.dt[EID %in% seg_eids.v, .N, keyby = SEX]

# Stage B: passed all QC (DARQ + bad-EID + outlier + AssemblyNet),
# matching PRISMA "Passed all processing" (38,567).
qc_pass_eids.v <- Reduce(intersect, list(
  seg_eids.v,
  darq_pass_eids.v,
  asnet_pass_eids.v,
  setdiff(seg_eids.v, exclude_ids.v),
  setdiff(seg_eids.v, outlier_ids.v)
))
stageB.dt <- sex.dt[EID %in% qc_pass_eids.v, .N, keyby = SEX]

# Stage C: sensitivity sample (post-QC minus ICD G|Q0, matching PRISMA 30,032)
icd_sens_eids.v   <- covars.dt[ICD_10 %like% "G|Q0", EID]
sens_eids.v       <- setdiff(qc_pass_eids.v, icd_sens_eids.v)
stageC.dt         <- sex.dt[EID %in% sens_eids.v, .N, keyby = SEX]

# Stage D: primary sample (sensitivity minus ICD F, matching PRISMA 27,680)
icd_f_only_eids.v <- covars.dt[ICD_10 %like% "F", EID]
primary_eids.v    <- setdiff(sens_eids.v, icd_f_only_eids.v)
stageD.dt         <- sex.dt[EID %in% primary_eids.v, .N, keyby = SEX]

# Stage E: GAMLSS final sample (already loaded from gamlss.rds)
stageE.dt <- final_all.dt

survivors.dt <- data.table(
  stage = c(
    "A. CNN-segmented imaging cohort",
    "B. Passed all QC (DARQ + bad-EID + outlier + AssemblyNet)",
    "C. Sensitivity sample (post-QC minus ICD G|Q0)",
    "D. Primary sample (post-QC minus ICD F|G|Q0)",
    "E. Final GAMLSS sample (primary minus missing education/site)"
  ),
  female = vapply(
    list(stageA.dt, stageB.dt, stageC.dt, stageD.dt, stageE.dt),
    get_n.fn, integer(1), sex.s = "Female"
  ),
  male = vapply(
    list(stageA.dt, stageB.dt, stageC.dt, stageD.dt, stageE.dt),
    get_n.fn, integer(1), sex.s = "Male"
  )
)
survivors.dt[, total      := female + male]
survivors.dt[, pct_female := round(100 * female / pmax(total, 1L), 1)]

message("\n=== Cumulative survivors by sex (PRISMA anchors) ===")
print(survivors.dt, row.names = FALSE)

log_section("Assembling exclusions summary table")

stages.v <- c(
  "1. Raw imaging cohort (ses-2 or ses-3 MRI)",
  "2. DARQ QC failures (DARQ < 0.25)",
  "3. Manually-curated bad-EID exclusion",
  "4. Volume outlier removal (|Z| > 3)",
  "5. AssemblyNet QC failures (grade != A)",
  "6. ICD-10 F-code (any)",
  "7. ICD-10 G-code (any)",
  "8. ICD-10 Q0-code (any)",
  "9. ICD-10 primary exclusion (F|G|Q0)",
  "10. Final GAMLSS sample (all splits)",
  "11. Final GAMLSS sample (train only)",
  "12. Final GAMLSS sample (test only)"
)
stage_dts.lst <- list(
  stage_imaging.dt,
  darq_fail.dt,
  bad_eid.dt,
  outlier_fail.dt,
  asnet_fail.dt,
  icd_f.dt,
  icd_g.dt,
  icd_q0.dt,
  icd_pri.dt,
  final_all.dt,
  final_train.dt,
  final_test.dt
)

summary.dt <- data.table(
  stage  = stages.v,
  female = vapply(stage_dts.lst, get_n.fn, integer(1), sex.s = "Female"),
  male   = vapply(stage_dts.lst, get_n.fn, integer(1), sex.s = "Male")
)
summary.dt[, total := female + male]
summary.dt[, pct_female := round(100 * female / pmax(total, 1L), 1)]

# Console summary (kept for interactive runs; logger captures the rest)
message("\n=== UKB sex-balance audit ===")
print(summary.dt, row.names = FALSE)
message("\n=== QC failure rates within sex (%) ===")
print(qc_rates.dt, row.names = FALSE)

# Anchor the final-sample claim in the log
final_total <- final_all.dt[, sum(N)]
final_pct_f <- round(100 * final_all.dt[SEX == "Female", N] / final_total, 1)
log_info("Final GAMLSS sample: %d (%.1f%% female)", final_total, final_pct_f)

# ----- Save -----
log_section("Saving audit output")
out_path <- here("outputs/audit/sex_balance.rds")
write_rds_safe(
  list(
    stages    = summary.dt,
    survivors = survivors.dt,
    qc_rates  = qc_rates.dt
  ),
  out_path,
  description = "sex-balance audit"
)

log_script_end("audit_sex_balance.R", success = TRUE)
