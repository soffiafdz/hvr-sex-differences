#!/usr/bin/env Rscript

# =============================================================================
# Head-size Adjustment and Matching
# =============================================================================
# Applies head-size adjustments and performs sex/age/ICC matching
#
# Inputs:
#   - data/raw/ukb-lng_icc-scale.csv: ICC and scale factors
#   - data/derivatives/hclvag_segmentations.csv: Segmentation volumes
#   - data/fst/ukb_assemblynet.fst: AssemblyNet segmentations
#   - data/derivatives/qc_darq-assemblynet.rds: QC data
#   - data/fst/ukb_covars.fst: Covariates
#
# Outputs:
#   - data/rds/sex-age-icc_matching.rds: Matched samples
#   - data/derivatives/hc-hvr_adj.rds: Adjusted HC/HVR values
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
library(stringr)
library(lubridate)
library(progress)
library(MatchIt)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("05_adjust_headsize.R")

# Load configuration
config <- load_config()

# Set random seed
set_seed()

# ----- Constants -----
REMATCH <- get_script_setting("adjust_headsize", "rematch_samples")
log_info("Rematch samples: %s", REMATCH)

# ----- Input Files -----
assemblynet_file <- get_data_path("raw", "assemblynet_fst")
exclude_ids_file <- get_data_path("raw", "exclude_ids")
icc_scale_file <- get_data_path("processed", "icc_scale")
segmentations_file <- get_data_path("processed", "hclvag_segmentations")
qc_file <- get_data_path("processed", "qc_combined")
covars_file <- get_data_path("processed", "covars_fst")

# Check required files exist
required_files <- c(icc_scale_file, segmentations_file, assemblynet_file)
check_files_exist(required_files, stop_on_missing = TRUE) |> invisible()

# ----- ICC Volume and Scale Factors -----
log_info("Loading ICC and scale factors")
icc_scl.dt <- read_csv_safe(
  icc_scale_file,
  col.names = c("EID", "INST", "ICC", "SCALE"),
  description = "ICC scale factors"
)

# ----- Segmentations -----
log_info("Loading segmentations")
# CNN
vols.dt <- read_csv_safe(segmentations_file, description = "CNN segmentations")

# AssemblyNet
asblynet.dt <- read_fst_safe(
  assemblynet_file,
  as_data_table = TRUE, description = "AssemblyNet segmentations"
)

# ----- QC Data -----
if (!check_files_exist(qc_file, stop_on_missing = FALSE)) {
  log_warn("QC file not found, running 02_quality_control.R")
  source(here("R/02_quality_control.R"))
}
qc.dt <- read_rds_safe(qc_file, description = "QC data")
validate_not_empty(qc.dt, "QC data")
validate_columns(qc.dt, c("DARQ", "ASBLYNET"), "QC data")

# ----- Covariates -----
if (!check_files_exist(covars_file, stop_on_missing = FALSE)) {
  log_warn("Covariates file not found, running 03_parse_covariates.R")
  source(here("R/03_parse_covariates.R"))
}
covars.dt <- read_fst_safe(
  covars_file,
  as_data_table = TRUE, description = "Covariates"
)
validate_not_empty(covars.dt, "covariates")
validate_columns(covars.dt, c("EID", "SEX_r", "DATE_mri_ses2", "ICD_10"), "covariates")

# ----- Data Cleaning -----
# AssemblyNet
asbly_cols.v <- asblynet.dt |>
  names() |>
  grep(pattern = "sub|vis|IC|Hipp|Inf. lat", value = TRUE) |>
  grep(pattern = "%|total|asymm", invert = TRUE, value = TRUE)
asblynet.dt <- asblynet.dt[, ..asbly_cols.v]
setnames(
  asblynet.dt, asbly_cols.v,
  c("EID", "INST", "ICC", "HC_L", "HC_R", "LV_L", "LV_R")
)

invisible({
  asblynet.dt[, HC := rowSums(.SD), .SDcols = patterns("^HC")]
  asblynet.dt[, LV := rowSums(.SD), .SDcols = patterns("^LV")]
  asblynet.dt[, EID := as.integer(str_remove(EID, "sub-"))]
  asblynet.dt[, SEGM := "ABLYNET"]
})

setkey(asblynet.dt, EID, INST)
rm(asbly_cols.v)

# CNN volumes
vols.dt[, c("l_amy", "r_amy") := NULL] |> invisible()
old_names.v <- names(vols.dt)[-1]
new_names.v <- old_names.v |>
  str_replace("^(.)_(..)_(.)$", "\\2_\\1_\\3") |>
  str_replace_all("vc", "lv") |>
  str_to_upper() |>
  str_replace_all(c("T$" = "tail", "B$" = "body", "H$" = "head"))
setnames(vols.dt, old_names.v, new_names.v)

invisible({
  vols.dt[, (new_names.v) := lapply(.SD, \(v) v / 1000), .SDcols = -1]
  vols.dt[, HC := rowSums(.SD), .SDcols = patterns("^HC")]
  vols.dt[, HC_L := rowSums(.SD), .SDcols = patterns("^HC_L")]
  vols.dt[, HC_R := rowSums(.SD), .SDcols = patterns("^HC_R")]
  vols.dt[, LV := rowSums(.SD), .SDcols = patterns("^LV")]
  vols.dt[, LV_L := rowSums(.SD), .SDcols = patterns("^LV_L")]
  vols.dt[, LV_R := rowSums(.SD), .SDcols = patterns("^LV_R")]
})

setnames(vols.dt, "id", "EID")
vols.dt[, let(
  EID   = as.integer(str_extract(EID, "(?<=sub-)\\d*")),
  INST  = str_extract(EID, "ses-\\d")
)] |> invisible()
setcolorder(vols.dt, "INST", after = "EID")
setkey(vols.dt, EID, INST)
rm(old_names.v, new_names.v)

# Merge ICC/Scale & volumes; bring back to subject space
invisible({
  icc_scl.dt[, EID := as.integer(str_remove(EID, "sub-"))]
  icc_scl.dt[, ICC := ICC / 1000] |> invisible()
  setkey(icc_scl.dt, EID, INST)
  cnn.dt <- icc_scl.dt[vols.dt, nomatch = NULL]
  roi_cols.v <- cnn.dt |>
    names() |>
    grep(patt = "^(HC|LV)", val = TRUE)
  cnn.dt[, (roi_cols.v) := .SD / SCALE, .SDcols = roi_cols.v]
  cnn.dt[, SEGM := "LPP_CNN"]
})
rm(roi_cols.v, vols.dt, icc_scl.dt)

# Quality control
log_info("Applying quality control filters")
exc_id <- readLines(exclude_ids_file) |> as.numeric()
cnn.dt <- cnn.dt[!qc.dt[DARQ < .25]][!EID %in% exc_id]
asblynet.dt <- asblynet.dt[!qc.dt[ASBLYNET != "A"]]
rm(qc.dt)

# Exclude outliers:
# Age & Sex
age_sex.dt <-
  covars.dt[
    !is.na(DATE_mri_ses2) | !is.na(DATE_mri_ses3),
    .(BIRTH_m, BIRTH_y, DATE_mri_ses2, DATE_mri_ses3),
    by = .(EID, SEX = SEX_r)
  ][
    # Convert month name to number (match returns 1 for "January", 2 for "February", etc.)
    , BIRTH_my := my(sprintf("%i-%i", match(BIRTH_m, month.name), BIRTH_y))
  ] |>
  melt(measure = patterns("DATE"), variable = "INST") |>
  (\(data_table) {
    data_table[
      , INST := fifelse(INST %like% 2, "ses-2", "ses-3")
    ][
      , value := str_extract(value, "^[^ ]+")
    ][
      , AGE := as.duration(ymd(value) - BIRTH_my) / dyears(1)
    ][
      , c("BIRTH_m", "BIRTH_y", "BIRTH_my", "value") := NULL
    ]
  })() |>
  na.omit() |>
  setkey(EID, INST)

# Total volumes (bilateral)
cnn_long.dt <- cnn.dt |>
  melt(measure = patterns("HC|LV")) |>
  merge(age_sex.dt)

# Z-scores adjusted for age & head-size
invisible(
  cnn_long.dt[
    !variable %like% "tail|body|head",
    Z_res := lm(value ~ SEX + AGE + ICC) |> residuals() |> scale(),
    by = variable
  ]
)
# Cutoff: |Z| > threshold from config
outlier_threshold <- get_parameter("qc", "outlier_sd_threshold", default = 3.0)
outliers_ids.v <- cnn_long.dt[abs(Z_res) > outlier_threshold, unique(EID)]
log_info("Outlier threshold: |Z| > %.1f (excluded %d participants)",
         outlier_threshold, length(outliers_ids.v))

cnn.dt <- cnn_long.dt[!EID %in% outliers_ids.v, -"Z_res"] |>
  dcast(... ~ variable, value.var = "value") |>
  setkey(EID, INST)
rm(cnn_long.dt, outliers_ids.v)

# Clinical history exclusions
# Primary sample: exclude F (psychiatric), G (neurological), Q0 (CNS malformations)
# Sensitivity sample: exclude only G and Q0 (include psychiatric F-codes)
icd10_primary <- get_parameter("exclusions", "icd10_patterns")
icd10_sens <- get_parameter("exclusions", "icd10_patterns_sensitivity")

log_info("Primary exclusion pattern: %s", icd10_primary)
log_info("Sensitivity exclusion pattern: %s", icd10_sens)

# IDs to exclude for each sample
icd_primary.dt <- covars.dt[ICD_10 %like% icd10_primary, .(EID)]
icd_sens.dt <- covars.dt[ICD_10 %like% icd10_sens, .(EID)]

log_info("Primary exclusions (F|G|Q0): %d participants", nrow(icd_primary.dt))
log_info("Sensitivity exclusions (G|Q0 only): %d participants", nrow(icd_sens.dt))
log_info("F-code only participants (included in sensitivity): %d",
         nrow(icd_primary.dt) - nrow(icd_sens.dt))

# Create copies for sensitivity sample before primary exclusion
cnn_sens.dt <- copy(cnn.dt)[!icd_sens.dt, on = "EID"]
asblynet_sens.dt <- copy(asblynet.dt)[!icd_sens.dt, on = "EID"]

# Apply primary exclusions
cnn.dt <- cnn.dt[!icd_primary.dt, on = "EID"]
asblynet.dt <- asblynet.dt[!icd_primary.dt, on = "EID"]

log_info("Primary sample CNN: %d observations", nrow(cnn.dt))
log_info("Sensitivity sample CNN: %d observations", nrow(cnn_sens.dt))

rm(covars.dt, icd_primary.dt, icd_sens.dt)

# Keep union of segmentation methods - PRIMARY sample
hclv.dt <- rbind(
  cnn.dt[asblynet.dt[, 1:2], nomatch = NULL],
  asblynet.dt[cnn.dt[, c(1:2, 4)], nomatch = NULL],
  use.names = TRUE, fill = TRUE
) |>
  setcolorder("SEGM", after = "INST") |>
  setkey(EID, INST)

# Keep union of segmentation methods - SENSITIVITY sample
hclv_sens.dt <- rbind(
  cnn_sens.dt[asblynet_sens.dt[, 1:2], nomatch = NULL],
  asblynet_sens.dt[cnn_sens.dt[, c(1:2, 4)], nomatch = NULL],
  use.names = TRUE, fill = TRUE
) |>
  setcolorder("SEGM", after = "INST") |>
  setkey(EID, INST)

rm(cnn.dt, asblynet.dt, cnn_sens.dt, asblynet_sens.dt)

# ----- Subcohorts -----
# Which ICC measure to use.
# LPP_CNN ~ assemblynet < .9
icc_measure <- "LPP_CNN"

# Cross-sectional cohort - PRIMARY
crs.dt <- hclv.dt[
  icc_measure,
  on = "SEGM",
  .(EID, INST, ICC)
][
  !duplicated(EID)
][
  age_sex.dt,
  nomatch = NULL
] |>
  setcolorder("SEX", after = "EID")

# Cross-sectional cohort - SENSITIVITY (includes F-codes)
crs_sens.dt <- hclv_sens.dt[
  icc_measure,
  on = "SEGM",
  .(EID, INST, ICC)
][
  !duplicated(EID)
][
  age_sex.dt,
  nomatch = NULL
] |>
  setcolorder("SEX", after = "EID")

# Longitudinal cohort - PRIMARY
lng.dt <- hclv.dt[
  icc_measure,
  on = "SEGM",
  .(EID, INST, ICC)
][
  ,
  if (.N == 2) .SD, EID
][
  age_sex.dt,
  on = .(EID, INST),
  nomatch = NULL
] |>
  setcolorder("SEX", after = "EID") |>
  setkey(EID, INST)

# Longitudinal cohort - SENSITIVITY (includes F-codes)
lng_sens.dt <- hclv_sens.dt[
  icc_measure,
  on = "SEGM",
  .(EID, INST, ICC)
][
  ,
  if (.N == 2) .SD, EID
][
  age_sex.dt,
  on = .(EID, INST),
  nomatch = NULL
] |>
  setcolorder("SEX", after = "EID") |>
  setkey(EID, INST)

log_info("Primary CRS: %d subjects, Sensitivity CRS: %d subjects",
         nrow(crs.dt), nrow(crs_sens.dt))
log_info("Primary LNG: %d subjects, Sensitivity LNG: %d subjects",
         lng.dt[, uniqueN(EID)], lng_sens.dt[, uniqueN(EID)])

rm(age_sex.dt, icc_measure)

# ----- Matching Algorithm -----
# Using MatchIt package with optimal matching (minimizes total Mahalanobis distance)
# This replaces the previous greedy algorithm for more statistically rigorous pairing.
#
# Method: Optimal 1:1 matching without replacement
# Distance: Mahalanobis distance on AGE and ICV
# Calipers: Applied to enforce tolerance bounds from config

log_section("Sex/Age/ICC Matching (MatchIt Optimal)")
matching_output <- get_data_path("processed", "matching_pairs")

if (all(!REMATCH, check_files_exist(matching_output, stop_on_missing = FALSE))) {
  log_info("Loading existing matching results")
  match.lst <- read_rds_safe(matching_output, description = "Matched samples")
  validate_not_empty(match.lst$CRS, "cross-sectional matched pairs")
  validate_not_empty(match.lst$LNG, "longitudinal matched pairs")
  log_info(
    "Loaded %d cross-sectional and %d longitudinal matched pairs",
    nrow(match.lst$CRS), nrow(match.lst$LNG)
  )
} else {
  log_info("Running MatchIt optimal matching algorithm")

  # Matching tolerances from config (used as calipers)
  tol <- get_script_setting("adjust_headsize", "matching_tolerance")
  log_info(
    "Matching tolerances - Age: %.1f years, ICC: %.1f cc, Offset: %.1f years",
    tol$age, tol$icc, tol$offset
  )

  # Helper function to run MatchIt and extract pairs
  run_matchit <- function(dt, match_vars, caliper_list, label) {
    # Prepare data for MatchIt (SEX must be 0/1 treatment indicator)
    match_data <- copy(dt)
    match_data[, TREAT := fifelse(SEX %like% "M", 1L, 0L)]

    # Build formula
    formula_str <- paste("TREAT ~", paste(match_vars, collapse = " + "))
    match_formula <- as.formula(formula_str)

    # Compute SD of each variable for caliper (MatchIt uses SD units by default)
    # We want raw units, so convert: caliper_sd = caliper_raw / sd(var)
    caliper_sd <- sapply(names(caliper_list), function(v) {
      caliper_list[[v]] / sd(match_data[[v]], na.rm = TRUE)
    })

    log_info("  %s: Running optimal matching on %d subjects", label, nrow(match_data))
    log_info("  Variables: %s", paste(match_vars, collapse = ", "))
    log_info("  Calipers (SD units): %s",
             paste(sprintf("%s=%.3f", names(caliper_sd), caliper_sd), collapse = ", "))

    # Run MatchIt with optimal matching
    m_out <- tryCatch({
      matchit(
        match_formula,
        data = match_data,
        method = "optimal",
        distance = "mahalanobis",
        caliper = caliper_sd,
        std.caliper = TRUE,  # caliper_sd is already in SD units
        ratio = 1
      )
    }, error = function(e) {
      log_warn("  MatchIt optimal failed: %s", e$message)
      log_info("  Falling back to nearest-neighbor matching")
      matchit(
        match_formula,
        data = match_data,
        method = "nearest",
        distance = "mahalanobis",
        caliper = caliper_sd,
        std.caliper = TRUE,
        ratio = 1,
        replace = FALSE
      )
    })

    # Extract matched pairs
    matched_data <- match.data(m_out)
    n_matched <- sum(matched_data$TREAT == 1)
    log_info("  %s: Matched %d pairs", label, n_matched)

    if (n_matched == 0) {
      return(NULL)
    }

    # Get match IDs (subclass indicates pair membership)
    matched_data <- as.data.table(matched_data)
    matched_data[, subclass := as.integer(as.character(subclass))]

    # Split into male/female and merge to get pairs
    males <- matched_data[TREAT == 1, .(EID_m = EID, subclass,
                                         AGE_m = AGE, ICC_m = ICC)]
    females <- matched_data[TREAT == 0, .(EID_f = EID, subclass,
                                           AGE_f = AGE, ICC_f = ICC)]

    # Handle optional OFF column for longitudinal
    if ("OFF" %in% names(matched_data)) {
      males[, OFF_m := matched_data[TREAT == 1, OFF]]
      females[, OFF_f := matched_data[TREAT == 0, OFF]]
    }

    pairs <- merge(males, females, by = "subclass")
    pairs[, ID := .I]
    pairs[, DIFF_age := AGE_m - AGE_f]
    pairs[, DIFF_icc := ICC_m - ICC_f]

    if ("OFF_m" %in% names(pairs)) {
      pairs[, DIFF_off := OFF_m - OFF_f]
      return(pairs[, .(ID, EID_m, EID_f, DIFF_off, DIFF_age, DIFF_icc)])
    } else {
      return(pairs[, .(ID, EID_m, EID_f, DIFF_age, DIFF_icc)])
    }
  }

  # ----- Cross-sectional matching -----
  log_info("Cross-sectional matching:")
  crs_match <- run_matchit(
    dt = crs.dt[, .(EID, SEX, AGE, ICC)],
    match_vars = c("AGE", "ICC"),
    caliper_list = list(AGE = tol$age, ICC = tol$icc),
    label = "CRS"
  )

  # ----- Longitudinal matching -----
  log_info("Longitudinal matching:")
  # Prepare longitudinal data with offset
  lng_b.dt <- lng.dt[, .(OFF = diff(AGE)), by = EID][lng.dt["ses-2", on = "INST"]]

  lng_match <- run_matchit(
    dt = lng_b.dt[, .(EID, SEX, AGE, ICC, OFF)],
    match_vars = c("AGE", "ICC", "OFF"),
    caliper_list = list(AGE = tol$age, ICC = tol$icc, OFF = tol$offset),
    label = "LNG"
  )

  # Compile results

  match.lst <- list(
    CRS = if (!is.null(crs_match)) crs_match else data.table(
      ID = integer(), EID_m = integer(), EID_f = integer(),
      DIFF_age = numeric(), DIFF_icc = numeric()
    ),
    LNG = if (!is.null(lng_match)) lng_match else data.table(
      ID = integer(), EID_m = integer(), EID_f = integer(),
      DIFF_off = numeric(), DIFF_age = numeric(), DIFF_icc = numeric()
    )
  )

  rm(lng_b.dt)

  # Save matching results
  write_rds_safe(match.lst, matching_output, description = "Matched samples")
  log_info(
    "Matched %d cross-sectional and %d longitudinal pairs (optimal matching)",
    nrow(match.lst$CRS), nrow(match.lst$LNG)
  )
}

# ----- Head-size Adjustment -----
log_section("Head-size Adjustment")
log_info("Applying adjustment methods: NON, PRP, STX, RES")

# Helper function to melt volume data to long format
melt_hclv <- function(dt) {
  dt |>
    melt(measure = patterns("HC|LV"), value = "CC") |>
    (\(data_table) {
      data_table[
        , c("ROI", "SIDE", "SUBFIELD") := tstrsplit(variable, split = "_")
      ][
        is.na(SUBFIELD), SUBFIELD := "total"
      ][
        is.na(SIDE), SIDE := "LR"
      ][
        , c("SEX", "AGE", "variable") := NULL
      ]
      data_table[!is.na(CC)]
    })() |>
    setcolorder("CC", after = "SIDE") |>
    setkey(EID, INST)
}

# Make long-format for easier adjustment - PRIMARY and SENSITIVITY
hclv.dt <- melt_hclv(hclv.dt)
hclv_sens.dt <- melt_hclv(hclv_sens.dt)

# Reconcile subcohorts in a nested list
hclv.lst <- list(LPP_CNN = list(), ABLYNET = list())
for (i in seq_along(hclv.lst)) {
  hclv.lst[[i]] <- list(CRS = NULL, LNG = NULL)
  for (j in seq_along(hclv.lst[[i]])) {
    # PRIMARY: ALL sample (excludes F|G|Q0)
    hclv.lst[[i]][[j]]$ALL <- merge(
      list(crs.dt, lng.dt)[[j]][, -"ICC"],
      hclv.dt[names(hclv.lst)[i], on = "SEGM", -"SEGM"],
      all = FALSE
    )
    # PRIMARY: MTCH sample (matched pairs from primary)
    hclv.lst[[i]][[j]]$MTCH <- match.lst[[j]][, .(MATCH = ID, EID_m, EID_f)] |>
      melt(id = "MATCH", value = "EID") |>
      merge(hclv.lst[[i]][[j]]$ALL, by = "EID") |>
      setkey(EID, INST)
    hclv.lst[[i]][[j]]$MTCH[, variable := NULL]
    # SENSITIVITY: SENS sample (excludes only G|Q0, includes F-codes)
    hclv.lst[[i]][[j]]$SENS <- merge(
      list(crs_sens.dt, lng_sens.dt)[[j]][, -"ICC"],
      hclv_sens.dt[names(hclv.lst)[i], on = "SEGM", -"SEGM"],
      all = FALSE
    )
    if (i == 2) lapply(hclv.lst[[i]][[j]], \(DT) DT[, SUBFIELD := NULL])
    rm(j)
  }
  rm(i)
}
rm(crs.dt, lng.dt, crs_sens.dt, lng_sens.dt, hclv.dt, hclv_sens.dt, REMATCH)

# Apply adjustment methods
hc_hvr.lst <- hclv.lst |>
  lapply(\(segm_lvl) {
    lapply(
      segm_lvl, \(dsgn_lvl) {
        lapply(
          dsgn_lvl, \(match_lvl) {
            cols_by.v <- match_lvl |>
              names() |>
              grep(pattern = "ROI|SIDE|SUBFIELD", value = TRUE)
            match_lvl[, ICC_mean := mean(ICC)]
            # Extract ICC coefficient using named indexing (more robust than numeric indices)
            match_lvl[, B := coef(lm(CC ~ ICC))["ICC"], cols_by.v]
            match_lvl[, let(
              NON_cc = CC,
              PRP_cc = CC / ICC,
              STX_cc = CC * SCALE,
              RES_cc = CC - B * (ICC - ICC_mean)
            )]
            match_lvl[, c("CC", "B", "ICC_mean") := NULL]
            match_lvl |>
              melt(measure = patterns("_cc$"), variable = "ADJ") |>
              dcast(... ~ ROI, value = "value") |>
              (\(data_table) {
                data_table[, ADJ := str_remove(ADJ, "_cc$")]
                # HVR is calculated from RAW (NON) HC and LV values ONLY
                # HVR = HC_raw / (HC_raw + LV_raw) - this is the original definition
                # HVR is self-normalizing and does NOT have adjustment methods
                # Only store HVR for ADJ == "NON" to avoid confusion and reduce storage
                hvr_vals <- data_table[ADJ == "NON", .(EID, INST, SIDE, HVR = HC / (HC + LV))]
                if ("SUBFIELD" %in% cols_by.v) {
                  hvr_vals <- data_table[ADJ == "NON" & SUBFIELD == "total",
                                         .(EID, INST, SIDE, SUBFIELD, HVR = HC / (HC + LV))]
                  data_table <- merge(data_table, hvr_vals,
                                      by = c("EID", "INST", "SIDE", "SUBFIELD"), all.x = TRUE)
                } else {
                  data_table <- merge(data_table, hvr_vals,
                                      by = c("EID", "INST", "SIDE"), all.x = TRUE)
                }
                # Set HVR to NA for non-NON rows (HVR only valid for unadjusted values)
                data_table[ADJ != "NON", HVR := NA_real_]
                setkey(data_table, EID, INST)
                data_table
              })()
          }
        )
      }
    )
  })

rm(hclv.lst)

# ----- Outputs -----
log_section("Saving outputs")

# Validate output structure before saving
for (segm in names(hc_hvr.lst)) {
  for (cohort in names(hc_hvr.lst[[segm]])) {
    for (match_type in names(hc_hvr.lst[[segm]][[cohort]])) {
      validate_not_empty(
        hc_hvr.lst[[segm]][[cohort]][[match_type]],
        sprintf("hc_hvr %s/%s/%s", segm, cohort, match_type)
      )
    }
  }
}

log_info("Output validation passed for all segmentation/cohort/match combinations")
log_info("LPP_CNN CRS ALL: %d subjects", hc_hvr.lst$LPP_CNN$CRS$ALL[!duplicated(EID), .N])
log_info("LPP_CNN CRS SENS: %d subjects", hc_hvr.lst$LPP_CNN$CRS$SENS[!duplicated(EID), .N])
log_info("LPP_CNN CRS MTCH: %d pairs", hc_hvr.lst$LPP_CNN$CRS$MTCH[!duplicated(MATCH), .N])
log_info("LPP_CNN LNG ALL: %d subjects", hc_hvr.lst$LPP_CNN$LNG$ALL[!duplicated(EID), .N])
log_info("LPP_CNN LNG SENS: %d subjects", hc_hvr.lst$LPP_CNN$LNG$SENS[!duplicated(EID), .N])

# Note: matching_pairs was saved during matching algorithm if computed
# Save head-size adjusted data
write_rds_safe(
  hc_hvr.lst,
  get_data_path("processed", "hc_hvr_adjusted"),
  description = "Head-size adjusted HC/HVR data"
)

log_script_end("05_adjust_headsize.R", success = TRUE)
