#!/usr/bin/env Rscript

# =============================================================================
# Pre-compute Manuscript Objects
# =============================================================================
# Assembles all analysis outputs into a single environment for the Quarto
# manuscript, so that .qmd files load one RDS instead of re-running analyses.
#
# Inputs:
#   - All processed data from steps 05-11 (demographics, sex differences,
#     SEM fits, GAMLSS models, HVR comparisons, cognitive reliability)
#
# Outputs:
#   - data/derivatives/manuscript_env.rds  (list with $scalars, $data, $tables)
#
# Dependencies: 07, 08, 09, 10, 11
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/validation.R"))
source(here("R/utils/formatting.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/tables_core.R"))
source(here("R/utils/tables_gt.R"))
source(here("R/utils/tables_data.R"))
source(here("R/utils/tables_normative.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("12_manuscript_objects.R")

config <- load_config()

log_info("Pre-computing manuscript objects")

# =============================================================================
# Section A: Data Loading
# =============================================================================

log_info("Loading analysis data...")
analysis_data <- load_analysis_data()
demog_data    <- analysis_data$demog_data
sex_diff      <- analysis_data$sex_diff
sem_fit       <- analysis_data$sem_fit
sem_params    <- analysis_data$sem_params
cfa_minv      <- analysis_data$cfa_minv

hvr_comparison    <- load_hvr_comparison()
brain_data        <- load_brain_volumes()
matched_brain_data <- load_matched_brain_volumes(demog_data)
norm_tables       <- load_norm_tables()
gamlss_calib      <- load_gamlss_calibration()
cog_reliability   <- load_cog_reliability()

# =============================================================================
# Section B: Derived Scalars
# =============================================================================

log_info("Computing derived scalars...")

# --- Sample sizes ---
all_demog <- demog_data$CRS$ALL
setDT(all_demog)
n_primary  <- nrow(all_demog)
n_female   <- nrow(all_demog[SEX == "Female"])
n_male     <- nrow(all_demog[SEX == "Male"])
pct_female <- round(n_female / n_primary * 100, 1)

# Age stats
age_min    <- round(min(all_demog$AGE), 1)
age_max    <- round(max(all_demog$AGE), 1)
age_mean_f <- round(mean(all_demog[SEX == "Female"]$AGE), 1)
age_mean_m <- round(mean(all_demog[SEX == "Male"]$AGE), 1)
age_sd_f   <- round(sd(all_demog[SEX == "Female"]$AGE), 1)
age_sd_m   <- round(sd(all_demog[SEX == "Male"]$AGE), 1)

# England subsample
eng_demog <- demog_data$CRS$ENG
if (!is.null(eng_demog)) {
  setDT(eng_demog)
  n_england <- nrow(eng_demog)
} else {
  n_england <- NA
}

# Matched sample
mtch_demog <- demog_data$CRS$MTCH
setDT(mtch_demog)
n_matched <- nrow(mtch_demog)
n_initial_imaging <- 47398

sample_sizes    <- extract_sample_sizes(demog_data)
n_matched_pairs <- sample_sizes$n_matched_pairs

# --- Effect size rows ---
bilateral <- sex_diff$OVERALL[SIDE == "LR"]
setDT(bilateral)

hc_non  <- bilateral[ROI == "HC" & ADJ == "NON"]
hc_prp  <- bilateral[ROI == "HC" & ADJ == "PRP"]
hc_stx  <- bilateral[ROI == "HC" & ADJ == "STX"]
hc_res  <- bilateral[ROI == "HC" & ADJ == "RES"]
lv_non  <- bilateral[ROI == "LV" & ADJ == "NON"]
lv_prp  <- bilateral[ROI == "LV" & ADJ == "PRP"]
lv_stx  <- bilateral[ROI == "LV" & ADJ == "STX"]
lv_res  <- bilateral[ROI == "LV" & ADJ == "RES"]
hvr_non <- bilateral[ROI == "HVR" & ADJ == "NON"]

matched_hvr    <- sex_diff$COMPARISON[ROI == "HVR" & ADJ == "NON"]
matched_hc_res <- sex_diff$COMPARISON[ROI == "HC" & ADJ == "RES"]
matched_lv_res <- sex_diff$COMPARISON[ROI == "LV" & ADJ == "RES"]
matched_hc_non <- sex_diff$COMPARISON[ROI == "HC" & ADJ == "NON"]
matched_lv_non <- sex_diff$COMPARISON[ROI == "LV" & ADJ == "NON"]

# --- Derived percentages ---
hvr_icv <- get_hvr_icv_validation(sex_diff)

r_hvr_icv <- hvr_icv[VARIABLE == "HVR", CORRELATION]
r_hc_icv  <- hvr_icv[VARIABLE == "HC", CORRELATION]
r_lv_icv  <- hvr_icv[VARIABLE == "LV", CORRELATION]
icv_reduction_pct <- round((1 - abs(r_hvr_icv) / abs(r_hc_icv)) * 100, 0)

hvr_d_full    <- hvr_non$ESTIMATE
hvr_d_matched <- matched_hvr$ESTIMATE_MTCH
hvr_reduction_pct <- round((1 - abs(hvr_d_matched) / abs(hvr_d_full)) * 100, 0)

hc_d_range <- max(bilateral[ROI == "HC"]$ESTIMATE) - min(bilateral[ROI == "HC"]$ESTIMATE)

# --- SEM paths (pooled) ---
hc_g     <- get_sem_path(sem_params, "HC_COG_POOLED", "g", "HC")
hvr_g    <- get_sem_path(sem_params, "HVR_COG_POOLED", "g", "HVR")
hc_res_g <- get_sem_path(sem_params, "HC_RES_COG_POOLED", "g", "HC_RES")

# SEM fit measures
sem_fit_hc  <- sem_fit$HC_COG
sem_fit_hvr <- sem_fit$HVR_COG

# --- Age interactions ---
age_int  <- hvr_comparison$age_sex_interactions
hvr_int  <- age_int[VARIABLE == "HVR"]
hc_int   <- age_int[VARIABLE == "HC"]
lv_int   <- age_int[VARIABLE == "LV"]

# Age-stratified effect sizes
hvr_age_strat <- sex_diff$AGE_STRATIFIED[ROI == "HVR" & ADJ == "NON" & SIDE == "LR"]
setDT(hvr_age_strat)
hvr_d_young <- hvr_age_strat[AGE_BIN == "[45,55]"]
hvr_d_old   <- hvr_age_strat[AGE_BIN == "(75,85]"]

# =============================================================================
# Section C: Heavy Computations
# =============================================================================

# --- C1: Temporal Stability (ICC) ---
log_info("Computing temporal stability metrics...")

calc_temporal_stability <- function() {
  gamlss_site_path <- get_data_path("models", "fit", "gamlss_site")
  gamlss_site <- tryCatch(
    read_rds_safe(gamlss_site_path, "gamlss_site"),
    error = function(e) NULL
  )

  if (is.null(gamlss_site) || is.null(gamlss_site$DATA$LNG)) {
    return(NULL)
  }

  lng <- gamlss_site$DATA$LNG
  setDT(lng)

  calc_retest <- function(data, roi_name, side = "LR") {
    roi_lng <- data[ROI == roi_name & ADJ == "NON" & SIDE == side]
    if (nrow(roi_lng) == 0) return(list(n = NA, r = NA, icc = NA))

    roi_wide <- dcast(roi_lng, EID + SEX ~ INST, value.var = c("AGE", "VAL"))
    roi_wide <- roi_wide[complete.cases(roi_wide)]

    if (nrow(roi_wide) < 10) return(list(n = nrow(roi_wide), r = NA, icc = NA))

    r <- cor(roi_wide$`VAL_ses-2`, roi_wide$`VAL_ses-3`)

    vals <- as.matrix(roi_wide[, .(`VAL_ses-2`, `VAL_ses-3`)])
    n <- nrow(vals)
    k <- 2

    grand_mean <- mean(vals)
    row_means  <- rowMeans(vals)
    col_means  <- colMeans(vals)

    SS_total    <- sum((vals - grand_mean)^2)
    SS_rows     <- k * sum((row_means - grand_mean)^2)
    SS_cols     <- n * sum((col_means - grand_mean)^2)
    SS_residual <- SS_total - SS_rows - SS_cols

    MS_rows     <- SS_rows / (n - 1)
    MS_residual <- SS_residual / ((n - 1) * (k - 1))
    MS_cols     <- SS_cols / (k - 1)

    icc <- (MS_rows - MS_residual) / (MS_rows + (k - 1) * MS_residual + k * (MS_cols - MS_residual) / n)

    list(n = n, r = r, icc = icc)
  }

  hvr_lr <- calc_retest(lng, "HVR", "LR")
  hc_lr  <- calc_retest(lng, "HC", "LR")
  hc_l   <- calc_retest(lng, "HC", "L")
  hc_r   <- calc_retest(lng, "HC", "R")
  lv_lr  <- calc_retest(lng, "LV", "LR")
  lv_l   <- calc_retest(lng, "LV", "L")
  lv_r   <- calc_retest(lng, "LV", "R")

  results <- data.table(
    ROI = c("Hippocampal-to-Ventricle Ratio",
            rep("Hippocampus", 3),
            rep("Lateral Ventricles", 3)),
    Laterality = c("Bilateral",
                   "Bilateral", "Left", "Right",
                   "Bilateral", "Left", "Right"),
    N = c(hvr_lr$n,
          hc_lr$n, hc_l$n, hc_r$n,
          lv_lr$n, lv_l$n, lv_r$n),
    `Pearson r` = c(hvr_lr$r,
                    hc_lr$r, hc_l$r, hc_r$r,
                    lv_lr$r, lv_l$r, lv_r$r),
    `ICC` = c(hvr_lr$icc,
              hc_lr$icc, hc_l$icc, hc_r$icc,
              lv_lr$icc, lv_l$icc, lv_r$icc)
  )

  results[!is.na(N)]
}

stability_dt <- calc_temporal_stability()

# Extract inline scalars from stability
if (!is.null(stability_dt)) {
  hvr_r   <- stability_dt[ROI == "Hippocampal-to-Ventricle Ratio" & Laterality == "Bilateral", `Pearson r`]
  hvr_icc <- stability_dt[ROI == "Hippocampal-to-Ventricle Ratio" & Laterality == "Bilateral", ICC]
  hc_r    <- stability_dt[ROI == "Hippocampus" & Laterality == "Bilateral", `Pearson r`]
  hc_icc  <- stability_dt[ROI == "Hippocampus" & Laterality == "Bilateral", ICC]
  lv_r    <- stability_dt[ROI == "Lateral Ventricles" & Laterality == "Bilateral", `Pearson r`]
  lv_icc  <- stability_dt[ROI == "Lateral Ventricles" & Laterality == "Bilateral", ICC]
} else {
  hvr_r <- hvr_icc <- hc_r <- hc_icc <- lv_r <- lv_icc <- NA
}

# --- C2: GAMLSS Model Coefficients ---
log_info("Extracting GAMLSS model coefficients...")

# Need gamlss library loaded for coefficient extraction
library(gamlss)

invisible(capture.output({
  hvr_coefs <- get_gamlss_model_summary(roi = "HVR", adj = "NON", side = "LR")
}, type = "output"))

invisible(capture.output({
  hc_coefs <- get_gamlss_model_summary(roi = "HC", adj = "NON", side = "LR")
}, type = "output"))

# --- C3: SEM Parameter Extractions ---
log_info("Extracting SEM parameters...")

# Speed paths (pooled)
hc_speed     <- get_sem_path(sem_params, "HC_COG_POOLED", "PRSP_s", "HC")
hvr_speed    <- get_sem_path(sem_params, "HVR_COG_POOLED", "PRSP_s", "HVR")
hc_res_speed <- get_sem_path(sem_params, "HC_RES_COG_POOLED", "PRSP_s", "HC_RES")

# Sex-stratified paths
hc_g_f       <- get_sem_sex_path(sem_params, "HC_COG", "g", "HC", 1)
hc_g_m       <- get_sem_sex_path(sem_params, "HC_COG", "g", "HC", 2)
hvr_g_f      <- get_sem_sex_path(sem_params, "HVR_COG", "g", "HVR", 1)
hvr_g_m      <- get_sem_sex_path(sem_params, "HVR_COG", "g", "HVR", 2)
hc_res_g_f   <- get_sem_sex_path(sem_params, "HC_RES_COG", "g", "HC_RES", 1)
hc_res_g_m   <- get_sem_sex_path(sem_params, "HC_RES_COG", "g", "HC_RES", 2)
hc_speed_f   <- get_sem_sex_path(sem_params, "HC_COG", "PRSP_s", "HC", 1)
hc_speed_m   <- get_sem_sex_path(sem_params, "HC_COG", "PRSP_s", "HC", 2)
hc_res_speed_m <- get_sem_sex_path(sem_params, "HC_RES_COG", "PRSP_s", "HC_RES", 2)

# SEM fit table
sem_fit_dt <- load_sem_fit()

# SEM full params + derived tables
sem_full_params <- load_sem_full_params()
loadings_dt     <- extract_sem_factor_loadings(sem_full_params, model_filter = "HVR_COG")
covar_dt        <- extract_sem_covariate_effects(sem_full_params, model_filter = "HVR_COG")

# Age-squared coefficients
hvr_age_sq      <- get_sem_age_sq_coef(sem_params, "HVR_COG_POOLED", "HVR")
g_age_sq_hvr    <- get_sem_age_sq_coef(sem_params, "HVR_COG_POOLED", "g")
speed_age_sq_hvr <- get_sem_age_sq_coef(sem_params, "HVR_COG_POOLED", "PRSP_s")
hc_res_age_sq   <- get_sem_age_sq_coef(sem_params, "HC_RES_COG_POOLED", "HC_RES")
g_age_sq_hc     <- get_sem_age_sq_coef(sem_params, "HC_RES_COG_POOLED", "g")
speed_age_sq_hc <- get_sem_age_sq_coef(sem_params, "HC_RES_COG_POOLED", "PRSP_s")

# =============================================================================
# Section D: Assemble and Save
# =============================================================================

log_info("Assembling manuscript environment...")

ms <- list(
  scalars = list(
    # Sample sizes
    n_primary = n_primary, n_female = n_female, n_male = n_male,
    pct_female = pct_female, n_england = n_england,
    n_matched = n_matched, n_matched_pairs = n_matched_pairs,
    n_initial_imaging = n_initial_imaging,
    # Age stats
    age_min = age_min, age_max = age_max,
    age_mean_f = age_mean_f, age_mean_m = age_mean_m,
    age_sd_f = age_sd_f, age_sd_m = age_sd_m,
    # Effect size rows
    hc_non = hc_non, hc_prp = hc_prp, hc_stx = hc_stx, hc_res = hc_res,
    lv_non = lv_non, lv_prp = lv_prp, lv_stx = lv_stx, lv_res = lv_res,
    hvr_non = hvr_non,
    matched_hvr = matched_hvr, matched_hc_res = matched_hc_res,
    matched_lv_res = matched_lv_res,
    matched_hc_non = matched_hc_non, matched_lv_non = matched_lv_non,
    # Derived percentages
    icv_reduction_pct = icv_reduction_pct, hvr_reduction_pct = hvr_reduction_pct,
    hc_d_range = hc_d_range,
    # HVR-ICV validation
    hvr_icv = hvr_icv,
    r_hvr_icv = r_hvr_icv, r_hc_icv = r_hc_icv, r_lv_icv = r_lv_icv,
    # SEM paths (pooled)
    hc_g = hc_g, hvr_g = hvr_g, hc_res_g = hc_res_g,
    # SEM fit
    sem_fit_hc = sem_fit_hc, sem_fit_hvr = sem_fit_hvr,
    # Age interactions
    hvr_int = hvr_int, hc_int = hc_int, lv_int = lv_int,
    hvr_age_strat = hvr_age_strat,
    hvr_d_young = hvr_d_young, hvr_d_old = hvr_d_old,
    # Speed paths
    hc_speed = hc_speed, hvr_speed = hvr_speed, hc_res_speed = hc_res_speed,
    # Sex-stratified paths
    hc_g_f = hc_g_f, hc_g_m = hc_g_m,
    hvr_g_f = hvr_g_f, hvr_g_m = hvr_g_m,
    hc_res_g_f = hc_res_g_f, hc_res_g_m = hc_res_g_m,
    hc_speed_f = hc_speed_f, hc_speed_m = hc_speed_m,
    hc_res_speed_m = hc_res_speed_m,
    # Temporal stability inline values
    hvr_r = hvr_r, hvr_icc = hvr_icc,
    hc_r = hc_r, hc_icc = hc_icc,
    lv_r = lv_r, lv_icc = lv_icc,
    # Age-squared coefficients
    hvr_age_sq = hvr_age_sq, g_age_sq_hvr = g_age_sq_hvr,
    speed_age_sq_hvr = speed_age_sq_hvr,
    hc_res_age_sq = hc_res_age_sq, g_age_sq_hc = g_age_sq_hc,
    speed_age_sq_hc = speed_age_sq_hc
  ),
  data = list(
    # Demographics (for tbl-demographics)
    demog_data = demog_data,
    all_demog = all_demog,
    eng_demog = eng_demog,
    mtch_demog = mtch_demog,
    # Sex differences (for figures and tables)
    sex_diff = sex_diff,
    bilateral = bilateral,
    # SEM params (for tbl-sem-paths, fig-sem-forest, format_sem_paths_by_sex_gt)
    sem_params = sem_params,
    sem_fit = sem_fit,
    cfa_minv = cfa_minv,
    # Brain data (for figures)
    brain_data = brain_data,
    matched_brain_data = matched_brain_data,
    # Normative tables (for centile tables + figures)
    norm_tables = norm_tables,
    # GAMLSS calibration (for calibration figure)
    gamlss_calib = gamlss_calib,
    # HVR comparison (for age interaction text)
    hvr_comparison = hvr_comparison,
    # Cognitive reliability
    cog_reliability = cog_reliability,
    # SEM full params (for loadings/covariates tables)
    sem_full_params = sem_full_params
  ),
  tables = list(
    stability_dt = stability_dt,
    hvr_coefs = hvr_coefs,
    hc_coefs = hc_coefs,
    sem_fit_dt = sem_fit_dt,
    loadings_dt = loadings_dt,
    covar_dt = covar_dt
  )
)

out_path <- here("data/derivatives/manuscript_env.rds")
write_rds_safe(ms, out_path, "manuscript environment")

log_script_end("12_manuscript_objects.R", success = TRUE)
