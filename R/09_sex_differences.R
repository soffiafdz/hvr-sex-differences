#!/usr/bin/env Rscript

# =============================================================================
# Sex Differences Analysis
# =============================================================================
# Analyzes sex differences in brain volumes across head-size adjustment methods
#
# Key analyses:
# - Effect sizes (Cohen's d) by ROI × Adjustment × Hemisphere
# - Age-stratified effect sizes
# - Matched vs unmatched sample comparisons
# - Longitudinal stability of sex differences
#
# Inputs:
#   - data/derivatives/hc-hvr_adj.rds: Head-size adjusted volumes
#   - data/derivatives/z_scores.rds: Normative z-scores (optional)
#
# Outputs:
#   - data/derivatives/sex_differences.rds: Effect sizes and statistics
#   - outputs/tables/sex_differences/*.{html,tex}: Publication tables
#   - outputs/figures/sex_differences/*.png: Effect size plots
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
library(effsize)
library(lme4)
library(ggplot2)
library(GGally)
library(ggtext)
library(gt)
library(progress)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))
source(here("R/utils/tables_core.R"))
source(here("R/utils/plotting_core.R"))
source(here("R/utils/plotting_pipeline.R"))
source(here("R/utils/export.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("09_sex_differences.R")

# Load configuration
config <- load_config()
set_seed()

# ----- Constants -----
FORCE_REGENERATE <- get_script_setting("force_regenerate", "sex_differences", default = FALSE)
REDO_TABLES <- get_script_setting("sex_differences", "redo_tables", default = TRUE)
REDO_PLOTS <- get_script_setting("sex_differences", "redo_plots", default = TRUE)
NETWORK <- get_script_setting("sex_differences", "network", default = "LPP_CNN")
RUN_SENSITIVITY <- get_script_setting("sex_differences", "sensitivity_f_codes", default = FALSE)

ROIS <- get_parameter("rois")
ADJS <- get_parameter("adjustment_methods")
ICD10_PATTERN_PRIMARY <- get_parameter("exclusions", "icd10_patterns")
ICD10_PATTERN_SENSITIVITY <- get_parameter("exclusions", "icd10_patterns_sensitivity")

log_info("Network: %s", NETWORK)
log_info("Force regenerate: %s, Redo tables: %s, Redo plots: %s", FORCE_REGENERATE, REDO_TABLES, REDO_PLOTS)
log_info("Run sensitivity (include F-codes): %s", RUN_SENSITIVITY)

# ===========================================================================
# Load Data
# ===========================================================================
log_section("Loading input data")

# Head-size adjusted volumes
hc_hvr_path <- get_data_path("processed", "hc_hvr_adjusted")
if (!check_files_exist(hc_hvr_path, stop_on_missing = FALSE)) {
  log_warn("Head-size adjusted data not found, running 05_adjust_headsize.R")
  source(here("R", "05_adjust_headsize.R"))
}
hc_hvr.lst <- read_rds_safe(hc_hvr_path, "Head-size adjusted data")
validate_not_empty(hc_hvr.lst, "Head-size adjusted data")

# Load covariates for scanner site (Field 54 - Assessment Centre)
covars_path <- get_data_path("processed", "covars_fst")
if (!check_files_exist(covars_path, stop_on_missing = FALSE)) {
  log_error("Covariates file not found: %s", covars_path)
  stop("Covariates required for site adjustment", call. = FALSE)
}
covars.dt <- read_fst_safe(covars_path, as_data_table = TRUE, description = "Covariates")

# Extract scanner site (Assessment Centre, Field 54) for random effect
# SITE_ses2 is the assessment centre for imaging session 2
site.dt <- covars.dt[, .(SITE = as.factor(SITE_ses2)), EID]
site.dt <- site.dt[!is.na(SITE)]
setkey(site.dt, EID)
log_info("Scanner site data: %d subjects with valid SITE (Field 54)", nrow(site.dt))
log_info("  Sites: %s", paste(levels(site.dt$SITE), collapse = ", "))
rm(covars.dt)

# ===========================================================================
# Data Preparation
# ===========================================================================
log_section("Preparing data for sex differences analysis")

# Extract cross-sectional data
crs_data_wide.dt <- hc_hvr.lst[[NETWORK]]$CRS$ALL

# Validate
validate_columns(crs_data_wide.dt, c("EID", "SEX", "AGE", "INST", "ADJ", "SIDE", "HC", "LV", "HVR"),
                "Cross-sectional data")
validate_categorical(crs_data_wide.dt, "SEX", c("Female", "Male"), "sex")

# Reshape to long format for ROI analysis
crs_data.dt <- melt(
  crs_data_wide.dt,
  id.vars = c("EID", "INST", "SEX", "AGE", "ICC", "SCALE", "SIDE", "SUBFIELD", "ADJ"),
  measure.vars = c("HC", "LV", "HVR"),
  variable.name = "ROI",
  value.name = "VAL",
  na.rm = TRUE
)

# Filter to SUBFIELD == "total" only
crs_data.dt <- crs_data.dt[SUBFIELD == "total", .(EID, INST, SEX, AGE, ICC, ADJ, ROI, SIDE, VAL)]

# HVR is self-normalizing and does NOT have adjustment methods
# It's calculated from raw HC/LV and is identical across all ADJ rows
# Keep only one HVR row per observation (use NON as placeholder)
crs_data.dt <- crs_data.dt[!(ROI == "HVR" & ADJ != "NON")]

setkey(crs_data.dt, EID, INST, ROI, ADJ, SIDE)

# Merge scanner site (Field 54) for site-adjusted effect sizes
crs_data.dt <- merge(crs_data.dt, site.dt, by = "EID", all.x = TRUE)
n_with_site <- crs_data.dt[!is.na(SITE), uniqueN(EID)]
log_info("Merged scanner site: %d/%d subjects have SITE",
         n_with_site, crs_data.dt[, uniqueN(EID)])

log_info("Data loaded: %d observations from %d subjects",
         nrow(crs_data.dt), crs_data.dt[, uniqueN(EID)])
log_info("  Females: %d, Males: %d",
         crs_data.dt[SEX == "Female", uniqueN(EID)],
         crs_data.dt[SEX == "Male", uniqueN(EID)])
log_info("  Sides: %s", paste(unique(crs_data.dt$SIDE), collapse = ", "))

# ===========================================================================
# HVR Self-Normalizing Verification
# ===========================================================================
# Verify that HVR is independent of head size (ICV) as claimed.
# Self-normalization requires that HC and LV scale similarly with ICV
# such that their ratio cancels out head-size effects.
#
# We compute correlations of HVR, HC (all adjustment methods), and LV with ICV:
# - r(HVR, ICV) should be comparable to r(HC_RES, ICV) ≈ 0
# - This validates the "self-normalizing" claim in the manuscript

log_section("HVR self-normalizing verification")

# Helper function to compute correlation with confidence intervals
compute_cor_ci <- function(x, y, conf.level = 0.95) {
  complete <- complete.cases(x, y)
  x <- x[complete]
  y <- y[complete]
  test <- cor.test(x, y, conf.level = conf.level)
  list(
    r = test$estimate,
    ci_lower = test$conf.int[1],
    ci_upper = test$conf.int[2],
    p = test$p.value,
    n = length(x)
  )
}

# Get bilateral values for correlation analysis
# HC: all adjustment methods (NON, PRP, STX, RES)
# LV: unadjusted only (NON)
# HVR: self-normalizing (no adjustment)

# Extract HC for each adjustment method
hc_by_adj.dt <- crs_data.dt[SIDE == "LR" & ROI == "HC", .(EID, ADJ, HC = VAL, ICC)]
hc_wide.dt <- dcast(hc_by_adj.dt, EID + ICC ~ ADJ, value.var = "HC")
setnames(hc_wide.dt, c("NON", "PRP", "STX", "RES"), c("HC", "HC_PRP", "HC_STX", "HC_RES"))

# Extract LV unadjusted
lv_data.dt <- crs_data.dt[SIDE == "LR" & ADJ == "NON" & ROI == "LV", .(EID, LV = VAL)]

# Extract HVR
hvr_data.dt <- crs_data.dt[SIDE == "LR" & ROI == "HVR", .(EID, HVR = VAL)]

# Merge all
hvr_icv_data.dt <- hc_wide.dt[lv_data.dt, on = "EID"][hvr_data.dt, on = "EID"]

# Compute correlations with ICV for all variables
variables_to_test <- c("HC", "HC_PRP", "HC_STX", "HC_RES", "LV", "HVR")

hvr_icv_results.dt <- rbindlist(lapply(variables_to_test, function(var) {
  test <- compute_cor_ci(hvr_icv_data.dt[[var]], hvr_icv_data.dt$ICC)
  data.table(
    SAMPLE = "Primary",
    VARIABLE = var,
    CORRELATION = test$r,
    CI_LOWER = test$ci_lower,
    CI_UPPER = test$ci_upper,
    P_VALUE = test$p,
    N = test$n
  )
}))

# Compute variance explained
hvr_icv_results.dt[, R_SQUARED := CORRELATION^2]
hvr_icv_results.dt[, VARIANCE_EXPLAINED_PCT := R_SQUARED * 100]

# Compute reduction in ICV dependence relative to unadjusted HC (within Primary sample)
r_hc <- hvr_icv_results.dt[SAMPLE == "Primary" & VARIABLE == "HC", CORRELATION]
hvr_icv_results.dt[SAMPLE == "Primary", ICV_REDUCTION_PCT := (1 - abs(CORRELATION) / abs(r_hc)) * 100]
hvr_icv_results.dt[SAMPLE == "Primary" & VARIABLE == "HC", ICV_REDUCTION_PCT := NA_real_]  # Reference category

# Log results (Primary sample)
log_info("ICV correlation verification - Primary sample (comparing adjustment methods):")
for (var in variables_to_test) {
  row <- hvr_icv_results.dt[SAMPLE == "Primary" & VARIABLE == var]
  log_info("  r(%s, ICV) = %.3f [%.3f, %.3f], p = %.2e (%.1f%% variance)",
           var, row$CORRELATION, row$CI_LOWER, row$CI_UPPER, row$P_VALUE, row$VARIANCE_EXPLAINED_PCT)
}

r_hvr <- hvr_icv_results.dt[SAMPLE == "Primary" & VARIABLE == "HVR", CORRELATION]
r_hc_res <- hvr_icv_results.dt[SAMPLE == "Primary" & VARIABLE == "HC_RES", CORRELATION]
log_info("  HVR r(ICV) = %.3f vs HC_RES r(ICV) = %.3f", r_hvr, r_hc_res)

# Interpretation
if (abs(r_hvr) < 0.1) {
  log_info("  INTERPRETATION: HVR is effectively independent of head size (|r| < 0.1)")
} else if (abs(r_hvr) <= abs(r_hc_res) * 1.5) {
  log_info("  INTERPRETATION: HVR achieves ICV independence comparable to residualization")
} else {
  log_warn("  INTERPRETATION: HVR retains more ICV dependence than residualization")
}

rm(hc_by_adj.dt, hc_wide.dt, lv_data.dt, hvr_data.dt, hvr_icv_data.dt)  # Clean up

# ===========================================================================
# Calculate Effect Sizes
# ===========================================================================
# EFFECT SIZE CONVENTION (CL5):
# Cohen's d from effsize::cohen.d(VAL ~ SEX) computes (Mean_Female - Mean_Male) / SD_pooled
# when SEX factor levels are c("Female", "Male") - i.e., first level minus second level.
# Positive d = Females have LARGER values than Males
# Negative d = Males have LARGER values than Females
#
# For paired comparisons: d = (Female - Male) / SD_differences
# Positive d = Females have LARGER values within matched pairs
#
# Site-adjusted effect sizes (from lme4) are computed from SEXMale coefficient
# and NEGATED to match the cohen.d convention (positive = Females larger)

log_section("Calculating sex difference effect sizes")

sex_diff_path <- get_data_path("processed", "sex_differences")

if (FORCE_REGENERATE || !file.exists(sex_diff_path)) {
  log_info("Computing Cohen's d effect sizes%s",
           if (FORCE_REGENERATE) " (force_regenerate = TRUE)" else "")

  # Calculate effect sizes for each ROI × ADJ × SIDE
  effs.dt <- crs_data.dt[
    !is.na(VAL),
    {
      res <- cohen.d(VAL ~ SEX, na.rm = TRUE)
      list(
        ESTIMATE = res$estimate,
        MAGNITUDE = res$magnitude,
        CI_LOWER = res$conf.int[1],
        CI_UPPER = res$conf.int[2],
        N_FEMALE = sum(SEX == "Female"),
        N_MALE = sum(SEX == "Male")
      )
    },
    by = .(ROI, ADJ, SIDE)
  ]

  # Add labels
  effs.dt[, LABEL := sprintf("d = %.2f [%.2f, %.2f]", ESTIMATE, CI_LOWER, CI_UPPER)]

  # Add ROI and ADJ full names (ensure atomic character vectors)
  effs.dt[, ROI_LABEL := as.character(ROIS[ROI])]
  effs.dt[, ADJ_LABEL := as.character(ADJS[ADJ])]
  # HVR is self-normalizing, not "Unadjusted" - fix the label
  effs.dt[ROI == "HVR", ADJ_LABEL := "HVR (Self-Normalizing)"]

  log_info("Effect sizes calculated: %d combinations", nrow(effs.dt))

  # -------------------------------------------------------------------------
  # Site-Adjusted Effect Sizes (Mixed-Effects Model)
  # -------------------------------------------------------------------------
  # Scanner site (Field 54 Assessment Centre) can introduce systematic variance
  # Use lme4 to account for site as random intercept
  # This is consistent with the GAMLSS normative models which use random(SITE)
  log_info("Computing site-adjusted effect sizes")

  # Compute site-adjusted Cohen's d using mixed-effects model
  # Model: VAL ~ SEX + (1|SITE)
  # Cohen's d approximated as: beta_SEX / residual_SD
  effs_site_adj.dt <- crs_data.dt[
    !is.na(VAL) & !is.na(SITE),
    {
      # Fit mixed-effects model with site as random intercept
      tryCatch({
        mod <- lmer(VAL ~ SEX + (1|SITE), data = .SD, REML = TRUE)

        # Extract fixed effect for SEX (Male - Female, since Female is reference)
        beta_sex <- fixef(mod)["SEXMale"]

        # Get residual SD for effect size calculation
        sigma_resid <- sigma(mod)

        # Approximate Cohen's d = beta / sigma_residual
        # Negate to match cohen.d convention: positive = Females larger
        d_estimate <- -beta_sex / sigma_resid

        # Confidence interval via profile likelihood or Wald approximation
        se_beta <- sqrt(vcov(mod)["SEXMale", "SEXMale"])
        d_se <- se_beta / sigma_resid
        ci_lower <- d_estimate - 1.96 * d_se
        ci_upper <- d_estimate + 1.96 * d_se

        # Site variance proportion (ICC for site)
        var_comp <- as.data.frame(VarCorr(mod))
        var_site <- var_comp[var_comp$grp == "SITE", "vcov"]
        var_resid <- sigma_resid^2
        icc_site <- var_site / (var_site + var_resid)

        list(
          ESTIMATE = d_estimate,
          CI_LOWER = ci_lower,
          CI_UPPER = ci_upper,
          BETA_SEX = beta_sex,
          SE_BETA = se_beta,
          SIGMA_RESID = sigma_resid,
          ICC_SITE = icc_site,
          N_SITES = uniqueN(SITE),
          N_FEMALE = sum(SEX == "Female"),
          N_MALE = sum(SEX == "Male"),
          CONVERGED = TRUE
        )
      }, error = function(e) {
        list(
          ESTIMATE = NA_real_,
          CI_LOWER = NA_real_,
          CI_UPPER = NA_real_,
          BETA_SEX = NA_real_,
          SE_BETA = NA_real_,
          SIGMA_RESID = NA_real_,
          ICC_SITE = NA_real_,
          N_SITES = uniqueN(SITE),
          N_FEMALE = sum(SEX == "Female"),
          N_MALE = sum(SEX == "Male"),
          CONVERGED = FALSE
        )
      })
    },
    by = .(ROI, ADJ, SIDE)
  ]

  effs_site_adj.dt[, LABEL := sprintf("d = %.2f [%.2f, %.2f]", ESTIMATE, CI_LOWER, CI_UPPER)]
  effs_site_adj.dt[, ROI_LABEL := as.character(ROIS[ROI])]
  effs_site_adj.dt[, ADJ_LABEL := as.character(ADJS[ADJ])]
  effs_site_adj.dt[ROI == "HVR", ADJ_LABEL := "HVR (Self-Normalizing)"]

  log_info("Site-adjusted effect sizes calculated: %d combinations", nrow(effs_site_adj.dt))
  log_info("  Models converged: %d/%d", sum(effs_site_adj.dt$CONVERGED), nrow(effs_site_adj.dt))
  log_info("  Mean site ICC: %.3f", mean(effs_site_adj.dt$ICC_SITE, na.rm = TRUE))

  # Age-stratified effect sizes
  log_info("Computing age-stratified effect sizes")

  # Define age bins (BP1: configurable age stratification)
  age_breaks <- get_parameter("age_stratification", "breaks", default = c(45, 55, 65, 75, 85))
  log_info("Age bins: %s", paste(age_breaks, collapse = ", "))
  crs_data.dt[, AGE_BIN := cut(AGE, breaks = age_breaks, include.lowest = TRUE)]

  effs_age.dt <- crs_data.dt[
    !is.na(VAL) & !is.na(AGE_BIN),
    {
      res <- cohen.d(VAL ~ SEX, na.rm = TRUE)
      list(
        ESTIMATE = res$estimate,
        CI_LOWER = res$conf.int[1],
        CI_UPPER = res$conf.int[2],
        N_FEMALE = sum(SEX == "Female"),
        N_MALE = sum(SEX == "Male")
      )
    },
    by = .(ROI, ADJ, SIDE, AGE_BIN)
  ]

  effs_age.dt[, ROI_LABEL := as.character(ROIS[ROI])]
  effs_age.dt[, ADJ_LABEL := as.character(ADJS[ADJ])]
  effs_age.dt[ROI == "HVR", ADJ_LABEL := "HVR (Self-Normalizing)"]

  log_info("Age-stratified effect sizes calculated: %d combinations", nrow(effs_age.dt))

  # -------------------------------------------------------------------------
  # Matched-Sample Effect Sizes
  # -------------------------------------------------------------------------
  log_info("Computing matched-sample effect sizes")

  # Extract matched data (Sex/Age/ICV matched pairs)
  mtch_data_wide.dt <- hc_hvr.lst[[NETWORK]]$CRS$MTCH

  if (!is.null(mtch_data_wide.dt) && nrow(mtch_data_wide.dt) > 0) {
    # Validate matched data
    validate_columns(mtch_data_wide.dt, c("EID", "SEX", "AGE", "MATCH", "ADJ", "SIDE", "HC", "LV", "HVR"),
                    "Matched data")

    # Reshape to long format
    mtch_data.dt <- melt(
      mtch_data_wide.dt,
      id.vars = c("EID", "INST", "SEX", "AGE", "ICC", "SCALE", "SIDE", "SUBFIELD", "ADJ", "MATCH"),
      measure.vars = c("HC", "LV", "HVR"),
      variable.name = "ROI",
      value.name = "VAL",
      na.rm = TRUE
    )

    # Filter to SUBFIELD == "total" only
    mtch_data.dt <- mtch_data.dt[SUBFIELD == "total",
                                  .(EID, INST, SEX, AGE, ICC, ADJ, ROI, SIDE, VAL, MATCH)]

    # HVR has no adjustment method - deduplicate (keep NON as placeholder)
    mtch_data.dt <- mtch_data.dt[!(ROI == "HVR" & ADJ != "NON")]

    setkey(mtch_data.dt, EID, INST, ROI, ADJ, SIDE)

    log_info("Matched data: %d observations from %d pairs",
             nrow(mtch_data.dt), mtch_data.dt[, uniqueN(MATCH)])
    log_info("  Females: %d, Males: %d",
             mtch_data.dt[SEX == "Female", uniqueN(EID)],
             mtch_data.dt[SEX == "Male", uniqueN(EID)])

    # Calculate effect sizes for matched sample
    effs_mtch.dt <- mtch_data.dt[
      !is.na(VAL),
      {
        res <- cohen.d(VAL ~ SEX, na.rm = TRUE)
        list(
          ESTIMATE = res$estimate,
          MAGNITUDE = res$magnitude,
          CI_LOWER = res$conf.int[1],
          CI_UPPER = res$conf.int[2],
          N_FEMALE = sum(SEX == "Female"),
          N_MALE = sum(SEX == "Male"),
          N_PAIRS = uniqueN(MATCH)
        )
      },
      by = .(ROI, ADJ, SIDE)
    ]

    effs_mtch.dt[, LABEL := sprintf("d = %.2f [%.2f, %.2f]", ESTIMATE, CI_LOWER, CI_UPPER)]
    effs_mtch.dt[, ROI_LABEL := as.character(ROIS[ROI])]
    effs_mtch.dt[, ADJ_LABEL := as.character(ADJS[ADJ])]
    effs_mtch.dt[ROI == "HVR", ADJ_LABEL := "HVR (Self-Normalizing)"]

    log_info("Matched sample effect sizes calculated: %d combinations", nrow(effs_mtch.dt))

    # Paired t-test for matched pairs (within-pair differences)
    log_info("Computing paired effect sizes (within-pair differences)")

    # Create wide format with female and male values side by side
    mtch_paired.dt <- mtch_data.dt[SIDE == "LR", .(MATCH, ROI, ADJ, SEX, VAL)]
    mtch_paired_wide.dt <- dcast(mtch_paired.dt, MATCH + ROI + ADJ ~ SEX, value.var = "VAL")

    effs_paired.dt <- mtch_paired_wide.dt[
      !is.na(Female) & !is.na(Male),
      {
        diff <- Female - Male
        mean_diff <- mean(diff)
        sd_diff <- sd(diff)
        n <- .N

        # Protection against division by zero (CL2)
        # If SD is zero (all differences identical) or n < 2, effect size is undefined
        if (n < 2 || sd_diff < .Machine$double.eps) {
          list(
            ESTIMATE_PAIRED = NA_real_,
            CI_LOWER_PAIRED = NA_real_,
            CI_UPPER_PAIRED = NA_real_,
            MEAN_DIFF = mean_diff,
            SD_DIFF = sd_diff,
            T_STAT = NA_real_,
            P_VALUE = NA_real_,
            N_PAIRS = n
          )
        } else {
          se <- sd_diff / sqrt(n)
          # Cohen's d for paired samples (using SD of differences)
          d_paired <- mean_diff / sd_diff
          # 95% CI for paired d
          ci_lower <- d_paired - 1.96 * sqrt(1/n + d_paired^2/(2*n))
          ci_upper <- d_paired + 1.96 * sqrt(1/n + d_paired^2/(2*n))
          # Paired t-test
          t_stat <- mean_diff / se
          p_val <- 2 * pt(-abs(t_stat), df = n - 1)
          list(
            ESTIMATE_PAIRED = d_paired,
            CI_LOWER_PAIRED = ci_lower,
            CI_UPPER_PAIRED = ci_upper,
            MEAN_DIFF = mean_diff,
            SD_DIFF = sd_diff,
            T_STAT = t_stat,
            P_VALUE = p_val,
            N_PAIRS = n
          )
        }
      },
      by = .(ROI, ADJ)
    ]

    effs_paired.dt[, ROI_LABEL := as.character(ROIS[ROI])]
    effs_paired.dt[, ADJ_LABEL := as.character(ADJS[ADJ])]
    effs_paired.dt[ROI == "HVR", ADJ_LABEL := "HVR (Self-Normalizing)"]

    log_info("Paired effect sizes calculated: %d combinations", nrow(effs_paired.dt))

    # Compare matched vs unmatched effect sizes
    log_info("Computing matched vs unmatched comparison")
    comparison.dt <- merge(
      effs.dt[SIDE == "LR", .(ROI, ADJ, ESTIMATE_ALL = ESTIMATE, CI_LOWER_ALL = CI_LOWER, CI_UPPER_ALL = CI_UPPER, N_ALL = N_FEMALE + N_MALE)],
      effs_mtch.dt[SIDE == "LR", .(ROI, ADJ, ESTIMATE_MTCH = ESTIMATE, CI_LOWER_MTCH = CI_LOWER, CI_UPPER_MTCH = CI_UPPER, N_MTCH = N_FEMALE + N_MALE)],
      by = c("ROI", "ADJ")
    )
    comparison.dt[, DIFF := ESTIMATE_MTCH - ESTIMATE_ALL]
    comparison.dt[, ROI_LABEL := as.character(ROIS[ROI])]
    comparison.dt[, ADJ_LABEL := as.character(ADJS[ADJ])]
    comparison.dt[ROI == "HVR", ADJ_LABEL := "HVR (Self-Normalizing)"]

    log_info("Comparison computed: %d combinations", nrow(comparison.dt))

  } else {
    log_warn("No matched sample data found - skipping matched-sample analyses")
    effs_mtch.dt <- NULL
    effs_paired.dt <- NULL
    comparison.dt <- NULL
  }

  # -------------------------------------------------------------------------
  # F-Code Sensitivity Analysis (Include Psychiatric Disorders)
  # -------------------------------------------------------------------------
  effs_sens.dt <- NULL
  sens_comparison.dt <- NULL

  if (RUN_SENSITIVITY) {
    log_info("Running F-code sensitivity analysis")

    # Check if sensitivity sample exists
    sens_data_wide.dt <- hc_hvr.lst[[NETWORK]]$CRS$SENS

    if (!is.null(sens_data_wide.dt) && nrow(sens_data_wide.dt) > 0) {
      log_info("Sensitivity sample found in data - using pre-computed sample")

      # Reshape to long format
      sens_data.dt <- melt(
        sens_data_wide.dt,
        id.vars = c("EID", "INST", "SEX", "AGE", "ICC", "SCALE", "SIDE", "SUBFIELD", "ADJ"),
        measure.vars = c("HC", "LV", "HVR"),
        variable.name = "ROI",
        value.name = "VAL",
        na.rm = TRUE
      )

      setkey(sens_data.dt, EID, INST, ROI, ADJ, SIDE)

      # Filter to SUBFIELD == "total" only
      # The source data already has L, R, and LR (sum of sides) - no recomputation needed
      sens_data.dt <- sens_data.dt[SUBFIELD == "total",
                                    .(EID, INST, SEX, AGE, ICC, ADJ, ROI, SIDE, VAL)]

      setkey(sens_data.dt, EID, INST, ROI, ADJ, SIDE)

    } else {
      log_info("Sensitivity sample not found - creating from covariates")

      # Load covariates to identify F-code participants
      covars_path <- get_data_path("processed", "covars_fst")
      if (file.exists(covars_path)) {
        covars.dt <- fst::read_fst(covars_path, as.data.table = TRUE)

        # Identify participants excluded ONLY due to F-codes (not G or Q0)
        f_only_eids.dt <- covars.dt[
          grepl("F", ICD_10, ignore.case = FALSE) &
          !grepl("G|Q0", ICD_10, ignore.case = FALSE),
          .(EID)
        ]

        log_info("Found %d participants with F-codes only", nrow(f_only_eids.dt))

        # Get volume data for F-code participants from segmentations
        segs_path <- get_data_path("processed", "hclvag_segmentations")
        if (file.exists(segs_path)) {
          segs.dt <- read_csv_safe(segs_path, description = "Segmentations for F-code analysis")

          # Merge with covariates for F-code participants
          f_data.dt <- merge(segs.dt, f_only_eids.dt, by = "EID")
          f_data.dt <- merge(f_data.dt, covars.dt[, .(EID, SEX, AGE, ICC, SCALE)], by = "EID")

          if (nrow(f_data.dt) > 0) {
            log_info("F-code sensitivity data: %d observations", nrow(f_data.dt))

            # Combine with primary sample for sensitivity analysis
            sens_combined.dt <- rbind(crs_data.dt[, .(EID, INST, SEX, AGE, ICC, ADJ, ROI, SIDE, VAL)],
                                       f_data.dt[, .(EID, INST, SEX, AGE, ICC, ADJ, ROI, SIDE, VAL)],
                                       fill = TRUE)
            sens_data.dt <- sens_combined.dt
          } else {
            log_warn("No F-code participant data found in segmentations")
            sens_data.dt <- NULL
          }
        } else {
          log_warn("Segmentations file not found: %s", segs_path)
          sens_data.dt <- NULL
        }
      } else {
        log_warn("Covariates file not found: %s", covars_path)
        sens_data.dt <- NULL
      }
    }

    # Calculate effect sizes for sensitivity sample if available
    if (!is.null(sens_data.dt) && nrow(sens_data.dt) > 0) {
      log_info("Computing sensitivity sample effect sizes")
      log_info("Sensitivity data: %d observations from %d subjects",
               nrow(sens_data.dt), sens_data.dt[, uniqueN(EID)])
      log_info("  Females: %d, Males: %d",
               sens_data.dt[SEX == "Female", uniqueN(EID)],
               sens_data.dt[SEX == "Male", uniqueN(EID)])

      effs_sens.dt <- sens_data.dt[
        !is.na(VAL),
        {
          res <- cohen.d(VAL ~ SEX, na.rm = TRUE)
          list(
            ESTIMATE = res$estimate,
            MAGNITUDE = res$magnitude,
            CI_LOWER = res$conf.int[1],
            CI_UPPER = res$conf.int[2],
            N_FEMALE = sum(SEX == "Female"),
            N_MALE = sum(SEX == "Male")
          )
        },
        by = .(ROI, ADJ, SIDE)
      ]

      effs_sens.dt[, LABEL := sprintf("d = %.2f [%.2f, %.2f]", ESTIMATE, CI_LOWER, CI_UPPER)]
      effs_sens.dt[, ROI_LABEL := as.character(ROIS[ROI])]
      effs_sens.dt[, ADJ_LABEL := as.character(ADJS[ADJ])]
      effs_sens.dt[ROI == "HVR", ADJ_LABEL := "HVR (Self-Normalizing)"]

      log_info("Sensitivity effect sizes calculated: %d combinations", nrow(effs_sens.dt))

      # Compare primary vs sensitivity
      sens_comparison.dt <- merge(
        effs.dt[SIDE == "LR", .(ROI, ADJ, ESTIMATE_PRIMARY = ESTIMATE, CI_LOWER_PRIMARY = CI_LOWER, CI_UPPER_PRIMARY = CI_UPPER, N_PRIMARY = N_FEMALE + N_MALE)],
        effs_sens.dt[SIDE == "LR", .(ROI, ADJ, ESTIMATE_SENS = ESTIMATE, CI_LOWER_SENS = CI_LOWER, CI_UPPER_SENS = CI_UPPER, N_SENS = N_FEMALE + N_MALE)],
        by = c("ROI", "ADJ")
      )
      sens_comparison.dt[, DIFF := ESTIMATE_SENS - ESTIMATE_PRIMARY]
      sens_comparison.dt[, ROI_LABEL := as.character(ROIS[ROI])]
      sens_comparison.dt[, ADJ_LABEL := as.character(ADJS[ADJ])]
      sens_comparison.dt[ROI == "HVR", ADJ_LABEL := "HVR (Self-Normalizing)"]

      log_info("Sensitivity comparison computed")

      # Compute HVR-ICV validation for sensitivity sample
      log_info("Computing HVR-ICV validation for sensitivity sample")

      # Extract data for ICV correlation (same approach as primary sample)
      sens_hc_by_adj.dt <- sens_data.dt[SIDE == "LR" & ROI == "HC", .(EID, ADJ, HC = VAL, ICC)]
      sens_hc_wide.dt <- dcast(sens_hc_by_adj.dt, EID + ICC ~ ADJ, value.var = "HC")
      if (all(c("NON", "PRP", "STX", "RES") %in% names(sens_hc_wide.dt))) {
        setnames(sens_hc_wide.dt, c("NON", "PRP", "STX", "RES"), c("HC", "HC_PRP", "HC_STX", "HC_RES"))

        sens_lv_data.dt <- sens_data.dt[SIDE == "LR" & ADJ == "NON" & ROI == "LV", .(EID, LV = VAL)]
        sens_hvr_data.dt <- sens_data.dt[SIDE == "LR" & ROI == "HVR", .(EID, HVR = VAL)]

        sens_hvr_icv_data.dt <- sens_hc_wide.dt[sens_lv_data.dt, on = "EID"][sens_hvr_data.dt, on = "EID"]

        # Compute correlations for sensitivity sample
        sens_hvr_icv_results.dt <- rbindlist(lapply(variables_to_test, function(var) {
          if (var %in% names(sens_hvr_icv_data.dt)) {
            test <- compute_cor_ci(sens_hvr_icv_data.dt[[var]], sens_hvr_icv_data.dt$ICC)
            data.table(
              SAMPLE = "Sensitivity",
              VARIABLE = var,
              CORRELATION = test$r,
              CI_LOWER = test$ci_lower,
              CI_UPPER = test$ci_upper,
              P_VALUE = test$p,
              N = test$n
            )
          }
        }))

        # Compute derived metrics for sensitivity sample
        sens_hvr_icv_results.dt[, R_SQUARED := CORRELATION^2]
        sens_hvr_icv_results.dt[, VARIANCE_EXPLAINED_PCT := R_SQUARED * 100]
        r_hc_sens <- sens_hvr_icv_results.dt[VARIABLE == "HC", CORRELATION]
        sens_hvr_icv_results.dt[, ICV_REDUCTION_PCT := (1 - abs(CORRELATION) / abs(r_hc_sens)) * 100]
        sens_hvr_icv_results.dt[VARIABLE == "HC", ICV_REDUCTION_PCT := NA_real_]

        # Append to primary results
        hvr_icv_results.dt <- rbind(hvr_icv_results.dt, sens_hvr_icv_results.dt)

        log_info("HVR-ICV validation for sensitivity sample computed (N = %d)",
                 sens_hvr_icv_data.dt[, .N])

        rm(sens_hc_by_adj.dt, sens_hc_wide.dt, sens_lv_data.dt, sens_hvr_data.dt, sens_hvr_icv_data.dt)
      } else {
        log_warn("Sensitivity sample missing required columns for HVR-ICV validation")
      }
    } else {
      log_warn("Sensitivity data not available - skipping F-code sensitivity analysis")
    }
  }

  # Store results
  sex_diffs.lst <- list(
    OVERALL = effs.dt,
    SITE_ADJUSTED = effs_site_adj.dt,
    AGE_STRATIFIED = effs_age.dt,
    MATCHED = effs_mtch.dt,
    PAIRED = effs_paired.dt,
    COMPARISON = comparison.dt,
    SENSITIVITY = effs_sens.dt,
    SENS_COMPARISON = sens_comparison.dt,
    HVR_ICV_VALIDATION = hvr_icv_results.dt
  )

  # Save
  ensure_directory(dirname(sex_diff_path))
  write_rds_safe(sex_diffs.lst, sex_diff_path, "Sex difference effect sizes")

} else {
  log_info("Loading existing sex difference effect sizes")
  sex_diffs.lst <- read_rds_safe(sex_diff_path, "Sex difference effect sizes")
}

effs.dt <- sex_diffs.lst$OVERALL
effs_site_adj.dt <- sex_diffs.lst$SITE_ADJUSTED
effs_age.dt <- sex_diffs.lst$AGE_STRATIFIED
effs_mtch.dt <- sex_diffs.lst$MATCHED
effs_paired.dt <- sex_diffs.lst$PAIRED
comparison.dt <- sex_diffs.lst$COMPARISON
effs_sens.dt <- sex_diffs.lst$SENSITIVITY
sens_comparison.dt <- sex_diffs.lst$SENS_COMPARISON
hvr_icv_results.dt <- sex_diffs.lst$HVR_ICV_VALIDATION

# ===========================================================================
# Generate Summary Tables
# ===========================================================================
if (REDO_TABLES) {
  log_section("Generating sex difference tables")

  tables_dir <- get_output_path("tables")
  ensure_directory(tables_dir)

  # ----- Table 1: Overall Effect Sizes (Bilateral) -----
  log_info("Creating overall effect sizes table (bilateral)")

  overall_bilat.dt <- effs.dt[SIDE == "LR", .(ROI, ADJ, ROI_LABEL, ADJ_LABEL,
                                               ESTIMATE, CI_LOWER, CI_UPPER,
                                               MAGNITUDE, N_FEMALE, N_MALE)]

  # Create gt table
  overall_gt <- overall_bilat.dt |>
    gt(groupname_col = "ROI_LABEL") |>
    tab_header(
      title = "Sex Differences in Brain Volumes",
      subtitle = "Effect sizes (Cohen's d) across adjustment methods"
    ) |>
    cols_hide(columns = c("ROI", "ADJ")) |>
    cols_label(
      ADJ_LABEL = "Adjustment",
      ESTIMATE = "Cohen's d",
      CI_LOWER = "95% CI Lower",
      CI_UPPER = "95% CI Upper",
      MAGNITUDE = "Magnitude",
      N_FEMALE = "N Female",
      N_MALE = "N Male"
    ) |>
    fmt_number(
      columns = c(ESTIMATE, CI_LOWER, CI_UPPER),
      decimals = 3
    ) |>
    fmt_number(
      columns = c(N_FEMALE, N_MALE),
      decimals = 0
    ) |>
    tab_source_note("Positive values indicate larger volumes in females.") |>
    tab_source_note("Effect size interpretation: |d| < 0.2 = negligible, 0.2-0.5 = small, 0.5-0.8 = medium, > 0.8 = large") |>
    opt_table_font(font = "Times New Roman")

  # Save
  gt::gtsave(overall_gt, file.path(tables_dir, "effect_sizes_overall.html"))
  gt::gtsave(overall_gt, file.path(tables_dir, "effect_sizes_overall.tex"))

  log_info("Overall effect sizes table saved")

  # ----- Table 1b: Site-Adjusted Effect Sizes (Bilateral) -----
  log_info("Creating site-adjusted effect sizes table (bilateral)")

  if (!is.null(effs_site_adj.dt) && nrow(effs_site_adj.dt) > 0) {
    site_adj_bilat.dt <- effs_site_adj.dt[SIDE == "LR", .(ROI, ADJ, ROI_LABEL, ADJ_LABEL,
                                                          ESTIMATE, CI_LOWER, CI_UPPER,
                                                          ICC_SITE, N_SITES, N_FEMALE, N_MALE)]

    # Create gt table
    site_adj_gt <- site_adj_bilat.dt |>
      gt(groupname_col = "ROI_LABEL") |>
      tab_header(
        title = "Site-Adjusted Sex Differences in Brain Volumes",
        subtitle = "Effect sizes (Cohen's d) with scanner site as random intercept"
      ) |>
      cols_hide(columns = c("ROI", "ADJ")) |>
      cols_label(
        ADJ_LABEL = "Adjustment",
        ESTIMATE = "Cohen's d",
        CI_LOWER = "95% CI Lower",
        CI_UPPER = "95% CI Upper",
        ICC_SITE = "Site ICC",
        N_SITES = "N Sites",
        N_FEMALE = "N Female",
        N_MALE = "N Male"
      ) |>
      fmt_number(
        columns = c(ESTIMATE, CI_LOWER, CI_UPPER, ICC_SITE),
        decimals = 3
      ) |>
      fmt_number(
        columns = c(N_SITES, N_FEMALE, N_MALE),
        decimals = 0
      ) |>
      tab_source_note("Effect sizes from mixed-effects models: VAL ~ SEX + (1|SITE)") |>
      tab_source_note("Site ICC = proportion of variance attributable to scanner site") |>
      opt_table_font(font = "Times New Roman")

    # Save
    gt::gtsave(site_adj_gt, file.path(tables_dir, "effect_sizes_site_adjusted.html"))
    gt::gtsave(site_adj_gt, file.path(tables_dir, "effect_sizes_site_adjusted.tex"))

    log_info("Site-adjusted effect sizes table saved")

    # ----- Comparison: Unadjusted vs Site-Adjusted -----
    log_info("Creating site adjustment comparison table")

    site_comparison.dt <- merge(
      effs.dt[SIDE == "LR", .(ROI, ADJ, d_UNADJ = ESTIMATE, CI_L_UNADJ = CI_LOWER, CI_U_UNADJ = CI_UPPER)],
      effs_site_adj.dt[SIDE == "LR", .(ROI, ADJ, d_SITE = ESTIMATE, CI_L_SITE = CI_LOWER, CI_U_SITE = CI_UPPER, ICC_SITE)],
      by = c("ROI", "ADJ")
    )
    site_comparison.dt[, DIFF := d_SITE - d_UNADJ]
    site_comparison.dt[, PCT_CHANGE := 100 * (d_SITE - d_UNADJ) / abs(d_UNADJ)]
    site_comparison.dt[, ROI_LABEL := as.character(ROIS[ROI])]
    site_comparison.dt[, ADJ_LABEL := as.character(ADJS[ADJ])]
    site_comparison.dt[ROI == "HVR", ADJ_LABEL := "HVR (Self-Normalizing)"]

    site_comp_gt <- site_comparison.dt |>
      gt(groupname_col = "ROI_LABEL") |>
      tab_header(
        title = "Impact of Site Adjustment on Sex Difference Effect Sizes",
        subtitle = "Comparison of unadjusted vs site-adjusted Cohen's d"
      ) |>
      cols_hide(columns = c("ROI", "ADJ", "CI_L_UNADJ", "CI_U_UNADJ", "CI_L_SITE", "CI_U_SITE")) |>
      cols_label(
        ADJ_LABEL = "Adjustment",
        d_UNADJ = "d (Unadjusted)",
        d_SITE = "d (Site-Adj)",
        DIFF = "Difference",
        PCT_CHANGE = "% Change",
        ICC_SITE = "Site ICC"
      ) |>
      fmt_number(
        columns = c(d_UNADJ, d_SITE, DIFF, ICC_SITE),
        decimals = 3
      ) |>
      fmt_number(
        columns = PCT_CHANGE,
        decimals = 1
      ) |>
      tab_source_note("Negative % change indicates attenuation after site adjustment") |>
      opt_table_font(font = "Times New Roman")

    gt::gtsave(site_comp_gt, file.path(tables_dir, "effect_sizes_site_comparison.html"))
    gt::gtsave(site_comp_gt, file.path(tables_dir, "effect_sizes_site_comparison.tex"))

    log_info("Site comparison table saved")
  } else {
    log_warn("Site-adjusted effect sizes not available - skipping site tables")
  }

  # ----- Table 2: Age-Stratified Effect Sizes -----
  log_info("Creating age-stratified effect sizes table")

  age_strat_bilat.dt <- effs_age.dt[SIDE == "LR" & ROI != "ICC"]

  # Format age bins for display
  age_strat_bilat.dt[, AGE_BIN_LABEL := gsub("[\\[\\(\\]\\)]", "", as.character(AGE_BIN))]
  age_strat_bilat.dt[, AGE_BIN_LABEL := gsub(",", "-", AGE_BIN_LABEL)]

  age_strat_gt <- age_strat_bilat.dt |>
    gt(groupname_col = "ROI_LABEL") |>
    tab_header(
      title = "Age-Stratified Sex Differences",
      subtitle = "Effect sizes by age group and adjustment method"
    ) |>
    cols_hide(columns = c("ROI", "ADJ", "SIDE", "AGE_BIN")) |>
    cols_label(
      ADJ_LABEL = "Adjustment",
      AGE_BIN_LABEL = "Age Range",
      ESTIMATE = "Cohen's d",
      CI_LOWER = "95% CI Lower",
      CI_UPPER = "95% CI Upper",
      N_FEMALE = "N Female",
      N_MALE = "N Male"
    ) |>
    fmt_number(
      columns = c(ESTIMATE, CI_LOWER, CI_UPPER),
      decimals = 3
    ) |>
    fmt_number(
      columns = c(N_FEMALE, N_MALE),
      decimals = 0
    ) |>
    tab_source_note("Age bins: 5-year intervals from 45 to 85 years") |>
    opt_table_font(font = "Times New Roman")

  # Save
  gt::gtsave(age_strat_gt, file.path(tables_dir, "effect_sizes_age_stratified.html"))
  gt::gtsave(age_strat_gt, file.path(tables_dir, "effect_sizes_age_stratified.tex"))

  log_info("Age-stratified effect sizes table saved")

  # ----- Table 3: Hemisphere Comparison -----
  log_info("Creating hemisphere comparison table")

  hemi_comp.dt <- effs.dt[ROI %in% c("HC", "LV") & SIDE %in% c("L", "R")]

  hemi_gt <- hemi_comp.dt |>
    gt(groupname_col = "ROI_LABEL") |>
    tab_header(
      title = "Sex Differences by Hemisphere",
      subtitle = "Left vs Right hemisphere comparison"
    ) |>
    cols_hide(columns = c("ROI", "ADJ")) |>
    cols_label(
      ADJ_LABEL = "Adjustment",
      SIDE = "Hemisphere",
      ESTIMATE = "Cohen's d",
      CI_LOWER = "95% CI Lower",
      CI_UPPER = "95% CI Upper",
      MAGNITUDE = "Magnitude"
    ) |>
    fmt_number(
      columns = c(ESTIMATE, CI_LOWER, CI_UPPER),
      decimals = 3
    ) |>
    opt_table_font(font = "Times New Roman")

  # Save
  gt::gtsave(hemi_gt, file.path(tables_dir, "effect_sizes_hemisphere.html"))
  gt::gtsave(hemi_gt, file.path(tables_dir, "effect_sizes_hemisphere.tex"))

  log_info("Hemisphere comparison table saved")

  # ----- Table 4: Matched Sample Effect Sizes -----
  if (!is.null(effs_mtch.dt) && nrow(effs_mtch.dt) > 0) {
    log_info("Creating matched sample effect sizes table")

    mtch_bilat.dt <- effs_mtch.dt[SIDE == "LR", .(ROI, ADJ, ROI_LABEL, ADJ_LABEL,
                                                   ESTIMATE, CI_LOWER, CI_UPPER,
                                                   MAGNITUDE, N_FEMALE, N_MALE, N_PAIRS)]

    mtch_gt <- mtch_bilat.dt |>
      gt(groupname_col = "ROI_LABEL") |>
      tab_header(
        title = "Sex Differences in Matched Sample",
        subtitle = "Effect sizes after matching on Age and ICV"
      ) |>
      cols_hide(columns = c("ROI", "ADJ")) |>
      cols_label(
        ADJ_LABEL = "Adjustment",
        ESTIMATE = "Cohen's d",
        CI_LOWER = "95% CI Lower",
        CI_UPPER = "95% CI Upper",
        MAGNITUDE = "Magnitude",
        N_FEMALE = "N Female",
        N_MALE = "N Male",
        N_PAIRS = "N Pairs"
      ) |>
      fmt_number(
        columns = c(ESTIMATE, CI_LOWER, CI_UPPER),
        decimals = 3
      ) |>
      fmt_number(
        columns = c(N_FEMALE, N_MALE, N_PAIRS),
        decimals = 0
      ) |>
      tab_source_note("Matched sample: females and males matched on age (±1 year) and ICV (±25 cc).") |>
      tab_source_note("Positive values indicate larger volumes in females.") |>
      opt_table_font(font = "Times New Roman")

    gt::gtsave(mtch_gt, file.path(tables_dir, "effect_sizes_matched.html"))
    gt::gtsave(mtch_gt, file.path(tables_dir, "effect_sizes_matched.tex"))

    log_info("Matched sample effect sizes table saved")
  }

  # ----- Table 5: Paired Effect Sizes -----
  if (!is.null(effs_paired.dt) && nrow(effs_paired.dt) > 0) {
    log_info("Creating paired effect sizes table")

    paired_gt <- effs_paired.dt |>
      gt(groupname_col = "ROI_LABEL") |>
      tab_header(
        title = "Paired Sex Differences (Within-Pair)",
        subtitle = "Effect sizes based on female-male differences within matched pairs"
      ) |>
      cols_hide(columns = c("ROI", "ADJ")) |>
      cols_label(
        ADJ_LABEL = "Adjustment",
        ESTIMATE_PAIRED = "Cohen's d",
        CI_LOWER_PAIRED = "95% CI Lower",
        CI_UPPER_PAIRED = "95% CI Upper",
        MEAN_DIFF = "Mean Diff (F-M)",
        SD_DIFF = "SD Diff",
        T_STAT = "t",
        P_VALUE = "p",
        N_PAIRS = "N Pairs"
      ) |>
      fmt_number(
        columns = c(ESTIMATE_PAIRED, CI_LOWER_PAIRED, CI_UPPER_PAIRED),
        decimals = 3
      ) |>
      fmt_number(
        columns = c(MEAN_DIFF, SD_DIFF, T_STAT),
        decimals = 2
      ) |>
      fmt_scientific(
        columns = P_VALUE,
        decimals = 2
      ) |>
      fmt_number(
        columns = N_PAIRS,
        decimals = 0
      ) |>
      tab_source_note("Paired analysis: Cohen's d computed from within-pair differences.") |>
      tab_source_note("Positive values indicate larger volumes in females.") |>
      opt_table_font(font = "Times New Roman")

    gt::gtsave(paired_gt, file.path(tables_dir, "effect_sizes_paired.html"))
    gt::gtsave(paired_gt, file.path(tables_dir, "effect_sizes_paired.tex"))

    log_info("Paired effect sizes table saved")
  }

  # ----- Table 6: Matched vs Unmatched Comparison -----
  if (!is.null(comparison.dt) && nrow(comparison.dt) > 0) {
    log_info("Creating matched vs unmatched comparison table")

    comp_gt <- comparison.dt |>
      gt(groupname_col = "ROI_LABEL") |>
      tab_header(
        title = "Comparison: Matched vs Full Sample",
        subtitle = "Effect of Age/ICV matching on sex difference estimates"
      ) |>
      cols_hide(columns = c("ROI", "ADJ")) |>
      cols_label(
        ADJ_LABEL = "Adjustment",
        ESTIMATE_ALL = "d (Full)",
        CI_LOWER_ALL = "CI Low",
        CI_UPPER_ALL = "CI High",
        ESTIMATE_MTCH = "d (Matched)",
        CI_LOWER_MTCH = "CI Low",
        CI_UPPER_MTCH = "CI High",
        DIFF = "Difference",
        N_ALL = "N Full",
        N_MTCH = "N Matched"
      ) |>
      fmt_number(
        columns = c(ESTIMATE_ALL, CI_LOWER_ALL, CI_UPPER_ALL,
                   ESTIMATE_MTCH, CI_LOWER_MTCH, CI_UPPER_MTCH, DIFF),
        decimals = 3
      ) |>
      fmt_number(
        columns = c(N_ALL, N_MTCH),
        decimals = 0
      ) |>
      tab_spanner(
        label = "Full Sample",
        columns = c(ESTIMATE_ALL, CI_LOWER_ALL, CI_UPPER_ALL, N_ALL)
      ) |>
      tab_spanner(
        label = "Matched Sample",
        columns = c(ESTIMATE_MTCH, CI_LOWER_MTCH, CI_UPPER_MTCH, N_MTCH)
      ) |>
      tab_source_note("Difference = d(Matched) - d(Full). Positive difference indicates larger effect after matching.") |>
      opt_table_font(font = "Times New Roman")

    gt::gtsave(comp_gt, file.path(tables_dir, "effect_sizes_comparison.html"))
    gt::gtsave(comp_gt, file.path(tables_dir, "effect_sizes_comparison.tex"))

    log_info("Matched vs unmatched comparison table saved")
  }

  # ----- Table 7: F-Code Sensitivity Effect Sizes -----
  if (!is.null(effs_sens.dt) && nrow(effs_sens.dt) > 0) {
    log_info("Creating F-code sensitivity effect sizes table")

    sens_bilat.dt <- effs_sens.dt[SIDE == "LR", .(ROI, ADJ, ROI_LABEL, ADJ_LABEL,
                                                   ESTIMATE, CI_LOWER, CI_UPPER,
                                                   MAGNITUDE, N_FEMALE, N_MALE)]

    sens_gt <- sens_bilat.dt |>
      gt(groupname_col = "ROI_LABEL") |>
      tab_header(
        title = "Sex Differences: Sensitivity Analysis (Including Psychiatric Disorders)",
        subtitle = "Effect sizes with F-code participants included"
      ) |>
      cols_hide(columns = c("ROI", "ADJ")) |>
      cols_label(
        ADJ_LABEL = "Adjustment",
        ESTIMATE = "Cohen's d",
        CI_LOWER = "95% CI Lower",
        CI_UPPER = "95% CI Upper",
        MAGNITUDE = "Magnitude",
        N_FEMALE = "N Female",
        N_MALE = "N Male"
      ) |>
      fmt_number(
        columns = c(ESTIMATE, CI_LOWER, CI_UPPER),
        decimals = 3
      ) |>
      fmt_number(
        columns = c(N_FEMALE, N_MALE),
        decimals = 0
      ) |>
      tab_source_note("Sensitivity analysis: includes participants with ICD-10 F-codes (psychiatric disorders).") |>
      tab_source_note("Primary analysis excluded F-codes due to associations with hippocampal volume.") |>
      opt_table_font(font = "Times New Roman")

    gt::gtsave(sens_gt, file.path(tables_dir, "effect_sizes_sensitivity.html"))
    gt::gtsave(sens_gt, file.path(tables_dir, "effect_sizes_sensitivity.tex"))

    log_info("F-code sensitivity effect sizes table saved")
  }

  # ----- Table 8: Primary vs Sensitivity Comparison -----
  if (!is.null(sens_comparison.dt) && nrow(sens_comparison.dt) > 0) {
    log_info("Creating primary vs sensitivity comparison table")

    sens_comp_gt <- sens_comparison.dt |>
      gt(groupname_col = "ROI_LABEL") |>
      tab_header(
        title = "Comparison: Primary vs Sensitivity Analysis",
        subtitle = "Effect of including psychiatric disorder participants"
      ) |>
      cols_hide(columns = c("ROI", "ADJ")) |>
      cols_label(
        ADJ_LABEL = "Adjustment",
        ESTIMATE_PRIMARY = "d (Primary)",
        CI_LOWER_PRIMARY = "CI Low",
        CI_UPPER_PRIMARY = "CI High",
        ESTIMATE_SENS = "d (Sensitivity)",
        CI_LOWER_SENS = "CI Low",
        CI_UPPER_SENS = "CI High",
        DIFF = "Difference",
        N_PRIMARY = "N Primary",
        N_SENS = "N Sensitivity"
      ) |>
      fmt_number(
        columns = c(ESTIMATE_PRIMARY, CI_LOWER_PRIMARY, CI_UPPER_PRIMARY,
                   ESTIMATE_SENS, CI_LOWER_SENS, CI_UPPER_SENS, DIFF),
        decimals = 3
      ) |>
      fmt_number(
        columns = c(N_PRIMARY, N_SENS),
        decimals = 0
      ) |>
      tab_spanner(
        label = "Primary (F-codes excluded)",
        columns = c(ESTIMATE_PRIMARY, CI_LOWER_PRIMARY, CI_UPPER_PRIMARY, N_PRIMARY)
      ) |>
      tab_spanner(
        label = "Sensitivity (F-codes included)",
        columns = c(ESTIMATE_SENS, CI_LOWER_SENS, CI_UPPER_SENS, N_SENS)
      ) |>
      tab_source_note("Difference = d(Sensitivity) - d(Primary). Minimal differences suggest robust findings.") |>
      opt_table_font(font = "Times New Roman")

    gt::gtsave(sens_comp_gt, file.path(tables_dir, "effect_sizes_sens_comparison.html"))
    gt::gtsave(sens_comp_gt, file.path(tables_dir, "effect_sizes_sens_comparison.tex"))

    log_info("Primary vs sensitivity comparison table saved")
  }

  log_info("Tables saved to: %s", tables_dir)
}

# ===========================================================================
# Generate Plots
# ===========================================================================
if (REDO_PLOTS) {
  log_section("Generating sex difference plots")

  plots_dir <- get_output_path("figures")
  ensure_directory(plots_dir)

  # Color palettes
  sex_colors <- get_palette("sex")
  adj_colors <- get_palette("adjustment")  # Named colors for adjustment methods incl HVR

  # ----- Plot 1: Effect Sizes by Adjustment Method (Bilateral) -----
  log_info("Creating effect sizes by adjustment plot")

  plot_data <- effs.dt[SIDE == "LR" & ROI != "ICC"]

  # Ensure labels are atomic character vectors
  plot_data[, ROI_LABEL := as.character(ROI_LABEL)]
  plot_data[, ADJ_LABEL := as.character(ADJ_LABEL)]

  # Reorder adjustments: NON first, include HVR (Self-Normalizing) for HVR rows
  adj_order <- c("NON", setdiff(names(ADJS), "NON"))
  adj_label_order <- c(as.character(ADJS[adj_order]), "HVR (Self-Normalizing)")
  plot_data[, ADJ := factor(ADJ, levels = adj_order)]
  plot_data[, ADJ_LABEL := factor(ADJ_LABEL, levels = adj_label_order)]
  plot_data[, ROI_LABEL := factor(ROI_LABEL)]

  p1 <- ggplot(plot_data, aes(x = ADJ_LABEL, y = ESTIMATE, color = ROI_LABEL, group = ROI_LABEL)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = CI_LOWER, ymax = CI_UPPER),
      width = 0.2,
      linewidth = 0.6
    ) +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
    labs(
      title = "Sex Differences by Adjustment Method",
      subtitle = "Bilateral volumes (left + right average)",
      x = "Head-Size Adjustment Method",
      y = "Effect Size (Cohen's d)\nPositive = Females > Males",
      color = "Brain Region",
      caption = paste(
        "NON = Unadjusted; PRP = Proportions (V/ICV);",
        "STX = Stereotaxic (V/ICV^1.15); RES = Residuals"
      )
    ) +
    theme_publication(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.caption = element_text(size = 8, hjust = 0)
    )

  save_plot(p1, file.path(plots_dir, "effect_sizes_by_adjustment.png"),
            width = 8, height = 6)

  # ----- Plot 2: Effect Sizes by ROI -----
  log_info("Creating effect sizes by ROI plot")

  p2 <- ggplot(plot_data, aes(x = ROI_LABEL, y = ESTIMATE, color = ADJ_LABEL)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(
      aes(ymin = CI_LOWER, ymax = CI_UPPER),
      position = position_dodge(width = 0.5),
      width = 0.2,
      linewidth = 0.6
    ) +
    scale_color_manual(values = adj_colors) +
    labs(
      title = "Sex Differences Across Brain Regions",
      x = "Brain Region",
      y = "Effect Size (Cohen's d)\nPositive = Females > Males",
      color = "Adjustment Method"
    ) +
    theme_publication(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    coord_flip()

  save_plot(p2, file.path(plots_dir, "effect_sizes_by_roi.png"),
            width = 8, height = 6)

  # ----- Plot 3: Age Trajectory -----
  log_info("Creating age trajectory plot")

  age_plot_data <- effs_age.dt[SIDE == "LR" & ROI != "ICC"]
  age_plot_data[, AGE_MID := sapply(AGE_BIN, function(x) {
    mean(as.numeric(unlist(strsplit(gsub("[\\[\\(\\]\\)]", "", as.character(x)), ","))))
  })]

  # Ensure labels are atomic character vectors
  age_plot_data[, ROI_LABEL := as.character(ROI_LABEL)]
  age_plot_data[, ADJ_LABEL := as.character(ADJ_LABEL)]

  age_plot_data[, ADJ := factor(ADJ, levels = adj_order)]
  age_plot_data[, ADJ_LABEL := factor(ADJ_LABEL, levels = adj_label_order)]
  age_plot_data[, ROI_LABEL := factor(ROI_LABEL)]

  p3 <- ggplot(age_plot_data, aes(x = AGE_MID, y = ESTIMATE, color = ADJ_LABEL)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = CI_LOWER, ymax = CI_UPPER, fill = ADJ_LABEL),
                alpha = 0.2, color = NA) +
    facet_wrap(~ ROI_LABEL, scales = "free_y", ncol = 2) +
    scale_color_manual(values = adj_colors) +
    scale_fill_manual(values = adj_colors) +
    labs(
      title = "Age Trajectory of Sex Differences",
      subtitle = "Effect sizes across age groups",
      x = "Age (years)",
      y = "Effect Size (Cohen's d)",
      color = "Adjustment",
      fill = "Adjustment"
    ) +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )

  save_plot(p3, file.path(plots_dir, "age_trajectory.png"),
            width = 9, height = 6)

  # ----- Plot 4: Scatter Matrix (Unadjusted) -----
  log_info("Creating scatter matrix plot")

  # Get unadjusted bilateral data (only HC and LV have bilateral averages)
  scatter_data <- crs_data.dt[ADJ == "NON" & SIDE == "LR" & ROI %in% c("HC", "LV")]

  # Create wide format
  scatter_wide <- dcast(scatter_data, EID + SEX + AGE ~ ROI, value.var = "VAL")

  # Select ROIs that exist in bilateral data
  scatter_rois <- c("HC", "LV")
  scatter_labels <- ROIS[scatter_rois]

  # Create scatter plot matrix
  p4 <- ggpairs(
    scatter_wide,
    columns = scatter_rois,
    aes(color = SEX, alpha = 0.3),
    upper = list(continuous = wrap("cor", size = 3)),
    lower = list(continuous = wrap("points", alpha = 0.2, size = 0.3)),
    diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
    columnLabels = as.character(scatter_labels)
  ) +
    scale_color_manual(values = sex_colors) +
    scale_fill_manual(values = sex_colors) +
    theme_minimal(base_size = 9) +
    labs(title = "Brain Volume Distributions by Sex (Unadjusted)")

  save_plot(p4, file.path(plots_dir, "scatter_matrix_unadjusted.png"),
            width = 10, height = 10)

  # ----- Plot 5: Hemisphere Comparison -----
  log_info("Creating hemisphere comparison plot")

  hemi_data <- effs.dt[ROI %in% c("HC", "LV") & SIDE %in% c("L", "R")]

  # Ensure labels are atomic character vectors
  hemi_data[, ROI_LABEL := as.character(ROI_LABEL)]
  hemi_data[, ADJ_LABEL := as.character(ADJ_LABEL)]

  hemi_data[, ADJ := factor(ADJ, levels = adj_order)]
  hemi_data[, ADJ_LABEL := factor(ADJ_LABEL, levels = adj_label_order)]
  hemi_data[, ROI_LABEL := factor(ROI_LABEL)]

  p5 <- ggplot(hemi_data, aes(x = SIDE, y = ESTIMATE, fill = ADJ_LABEL)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = CI_LOWER, ymax = CI_UPPER),
      position = position_dodge(width = 0.8),
      width = 0.2
    ) +
    facet_wrap(~ ROI_LABEL, scales = "free_y") +
    scale_fill_manual(values = adj_colors) +
    scale_x_discrete(labels = c("L" = "Left", "R" = "Right")) +
    labs(
      title = "Sex Differences by Hemisphere",
      x = "Hemisphere",
      y = "Effect Size (Cohen's d)",
      fill = "Adjustment"
    ) +
    theme_publication(base_size = 11) +
    theme(
      legend.position = "bottom"
    )

  save_plot(p5, file.path(plots_dir, "hemisphere_comparison.png"),
            width = 8, height = 5)

  # ----- Plot 6: Matched vs Unmatched Comparison -----
  if (!is.null(comparison.dt) && nrow(comparison.dt) > 0) {
    log_info("Creating matched vs unmatched comparison plot")

    # Prepare data for comparison plot
    comp_plot.dt <- melt(
      comparison.dt,
      id.vars = c("ROI", "ADJ", "ROI_LABEL", "ADJ_LABEL"),
      measure.vars = c("ESTIMATE_ALL", "ESTIMATE_MTCH"),
      variable.name = "SAMPLE",
      value.name = "ESTIMATE"
    )
    comp_plot.dt[, SAMPLE_LABEL := fifelse(SAMPLE == "ESTIMATE_ALL", "Full Sample", "Matched Sample")]

    # Get CIs
    ci_all.dt <- comparison.dt[, .(ROI, ADJ, SAMPLE = "ESTIMATE_ALL", CI_LOWER = CI_LOWER_ALL, CI_UPPER = CI_UPPER_ALL)]
    ci_mtch.dt <- comparison.dt[, .(ROI, ADJ, SAMPLE = "ESTIMATE_MTCH", CI_LOWER = CI_LOWER_MTCH, CI_UPPER = CI_UPPER_MTCH)]
    ci.dt <- rbind(ci_all.dt, ci_mtch.dt)

    comp_plot.dt <- merge(comp_plot.dt, ci.dt, by = c("ROI", "ADJ", "SAMPLE"))
    comp_plot.dt[, ADJ := factor(ADJ, levels = adj_order)]
    comp_plot.dt[, ADJ_LABEL := factor(ADJ_LABEL, levels = adj_label_order)]
    comp_plot.dt[, ROI_LABEL := factor(ROI_LABEL)]

    p6 <- ggplot(comp_plot.dt, aes(x = ADJ_LABEL, y = ESTIMATE, color = SAMPLE_LABEL, group = SAMPLE_LABEL)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_line(linewidth = 0.7, alpha = 0.7) +
      geom_point(size = 3, position = position_dodge(width = 0.2)) +
      geom_errorbar(
        aes(ymin = CI_LOWER, ymax = CI_UPPER),
        width = 0.15,
        linewidth = 0.5,
        position = position_dodge(width = 0.2)
      ) +
      facet_wrap(~ ROI_LABEL, scales = "free_y", ncol = 2) +
      scale_color_manual(values = c("Full Sample" = "#1f77b4", "Matched Sample" = "#ff7f0e")) +
      labs(
        title = "Effect of Age/ICV Matching on Sex Differences",
        subtitle = "Comparison of full sample vs matched sample effect sizes",
        x = "Head-Size Adjustment Method",
        y = "Effect Size (Cohen's d)",
        color = "Sample"
      ) +
      theme_publication(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.text = element_text(face = "bold")
      )

    save_plot(p6, file.path(plots_dir, "matched_vs_unmatched.png"),
              width = 9, height = 6)

    log_info("Matched vs unmatched comparison plot saved")
  }

  # ----- Plot 7: Paired Effect Sizes (Forest Plot) -----
  if (!is.null(effs_paired.dt) && nrow(effs_paired.dt) > 0) {
    log_info("Creating paired effect sizes forest plot")

    paired_plot.dt <- copy(effs_paired.dt)
    paired_plot.dt[, ADJ := factor(ADJ, levels = adj_order)]
    paired_plot.dt[, ADJ_LABEL := factor(ADJ_LABEL, levels = adj_label_order)]
    paired_plot.dt[, ROI_LABEL := factor(ROI_LABEL)]

    # Create label combining ROI and ADJ
    paired_plot.dt[, LABEL := paste0(ROI_LABEL, " (", ADJ_LABEL, ")")]
    paired_plot.dt[, LABEL := factor(LABEL, levels = unique(LABEL[order(ROI, ADJ)]))]

    p7 <- ggplot(paired_plot.dt, aes(x = ESTIMATE_PAIRED, y = reorder(LABEL, -as.numeric(LABEL)))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_errorbarh(
        aes(xmin = CI_LOWER_PAIRED, xmax = CI_UPPER_PAIRED, color = ROI_LABEL),
        height = 0.2,
        linewidth = 0.6
      ) +
      geom_point(aes(color = ROI_LABEL), size = 3) +
      scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
      labs(
        title = "Paired Sex Differences (Within-Pair Analysis)",
        subtitle = "Effect sizes from matched female-male pairs",
        x = "Effect Size (Cohen's d)\nPositive = Females > Males",
        y = "",
        color = "Brain Region"
      ) +
      theme_publication(base_size = 10) +
      theme(
        legend.position = "bottom"
      )

    save_plot(p7, file.path(plots_dir, "paired_effect_sizes.png"),
              width = 8, height = 6)

    log_info("Paired effect sizes forest plot saved")
  }

  # ----- Plot 8: Primary vs Sensitivity Comparison -----
  if (!is.null(sens_comparison.dt) && nrow(sens_comparison.dt) > 0) {
    log_info("Creating primary vs sensitivity comparison plot")

    # Prepare data for sensitivity comparison plot
    sens_comp_plot.dt <- melt(
      sens_comparison.dt,
      id.vars = c("ROI", "ADJ", "ROI_LABEL", "ADJ_LABEL"),
      measure.vars = c("ESTIMATE_PRIMARY", "ESTIMATE_SENS"),
      variable.name = "SAMPLE",
      value.name = "ESTIMATE"
    )
    sens_comp_plot.dt[, SAMPLE_LABEL := fifelse(SAMPLE == "ESTIMATE_PRIMARY", "Primary (F-codes excluded)", "Sensitivity (F-codes included)")]

    # Get CIs
    ci_primary.dt <- sens_comparison.dt[, .(ROI, ADJ, SAMPLE = "ESTIMATE_PRIMARY", CI_LOWER = CI_LOWER_PRIMARY, CI_UPPER = CI_UPPER_PRIMARY)]
    ci_sens.dt <- sens_comparison.dt[, .(ROI, ADJ, SAMPLE = "ESTIMATE_SENS", CI_LOWER = CI_LOWER_SENS, CI_UPPER = CI_UPPER_SENS)]
    ci_sens.dt <- rbind(ci_primary.dt, ci_sens.dt)

    sens_comp_plot.dt <- merge(sens_comp_plot.dt, ci_sens.dt, by = c("ROI", "ADJ", "SAMPLE"))
    sens_comp_plot.dt[, ADJ := factor(ADJ, levels = adj_order)]
    sens_comp_plot.dt[, ADJ_LABEL := factor(ADJ_LABEL, levels = adj_label_order)]
    sens_comp_plot.dt[, ROI_LABEL := factor(ROI_LABEL)]

    p8 <- ggplot(sens_comp_plot.dt, aes(x = ADJ_LABEL, y = ESTIMATE, color = SAMPLE_LABEL, group = SAMPLE_LABEL)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_line(linewidth = 0.7, alpha = 0.7) +
      geom_point(size = 3, position = position_dodge(width = 0.2)) +
      geom_errorbar(
        aes(ymin = CI_LOWER, ymax = CI_UPPER),
        width = 0.15,
        linewidth = 0.5,
        position = position_dodge(width = 0.2)
      ) +
      facet_wrap(~ ROI_LABEL, scales = "free_y", ncol = 2) +
      scale_color_manual(values = c("Primary (F-codes excluded)" = "#2ca02c", "Sensitivity (F-codes included)" = "#d62728")) +
      labs(
        title = "Sensitivity Analysis: Effect of Including Psychiatric Disorders",
        subtitle = "Comparison of primary vs sensitivity analysis effect sizes",
        x = "Head-Size Adjustment Method",
        y = "Effect Size (Cohen's d)",
        color = "Sample"
      ) +
      theme_publication(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.text = element_text(face = "bold")
      )

    save_plot(p8, file.path(plots_dir, "sensitivity_comparison.png"),
              width = 9, height = 6)

    log_info("Primary vs sensitivity comparison plot saved")
  }

  # ----- Plot 9: Distribution Overlap (HVR vs HC Methods) -----
  log_info("Creating distribution overlap plot")

  # Use crs_data.dt which has the long-format data
  if (exists("crs_data.dt") && nrow(crs_data.dt) > 0) {
    # Filter to bilateral only for cleaner visualization
    dist_data <- crs_data.dt[SIDE == "LR"]

    # Get effect sizes for annotation
    dist_effects <- effs.dt[SIDE == "LR", .(ROI, ADJ, ESTIMATE)]

    # Create the HVR vs HC distribution comparison plot
    p9 <- plot_hvr_vs_hc_distributions(
      data = dist_data,
      effect_sizes = dist_effects,
      sex_colors = get_palette("sex")
    )

    save_plot(p9, file.path(plots_dir, "distribution_overlap_hc_hvr.png"),
              width = 14, height = 5)

    log_info("Distribution overlap plot saved")

    # Also create the full faceted version (all ROIs × all adjustments)
    p10 <- plot_distribution_overlap_faceted(
      data = dist_data,
      effect_sizes = dist_effects,
      title = "Sex Differences in Brain Volumes by Adjustment Method",
      sex_colors = get_palette("sex")
    )

    save_plot(p10, file.path(plots_dir, "distribution_overlap_faceted.png"),
              width = 12, height = 10)

    log_info("Faceted distribution overlap plot saved")
  } else {
    log_warn("Distribution data not available - skipping overlap plots")
  }

  log_info("Plots saved to: %s", plots_dir)
}

# ===========================================================================
# Summary Statistics
# ===========================================================================
log_section("Sex Differences Summary")

log_info("Overall effect sizes (bilateral):")
for (roi in names(ROIS)) {
  for (adj in names(ADJS)) {
    eff <- effs.dt[ROI == roi & ADJ == adj & SIDE == "LR"]
    if (nrow(eff) > 0) {
      log_info("  %s (%s): d = %.3f [%.3f, %.3f] - %s",
               ROIS[[roi]], ADJS[[adj]],
               eff$ESTIMATE, eff$CI_LOWER, eff$CI_UPPER,
               eff$MAGNITUDE)
    }
  }
}

log_info("\nAge effects on sex differences:")
log_info("  Checking for age × sex interactions...")
age_range <- range(effs_age.dt$ESTIMATE, na.rm = TRUE)
log_info("  Effect size range across age bins: [%.3f, %.3f]", age_range[1], age_range[2])

log_info("\nHemisphere asymmetry in sex differences:")
hemi_asym <- effs.dt[ROI %in% c("HC", "LV") & SIDE %in% c("L", "R")]
hemi_summary <- hemi_asym[, .(
  MEAN_L = mean(ESTIMATE[SIDE == "L"]),
  MEAN_R = mean(ESTIMATE[SIDE == "R"]),
  DIFF = mean(ESTIMATE[SIDE == "L"]) - mean(ESTIMATE[SIDE == "R"])
), by = ROI]
log_info("  Left-Right differences:")
for (i in seq_len(nrow(hemi_summary))) {
  log_info("    %s: L = %.3f, R = %.3f, Diff = %.3f",
           ROIS[[hemi_summary$ROI[i]]],
           hemi_summary$MEAN_L[i],
           hemi_summary$MEAN_R[i],
           hemi_summary$DIFF[i])
}

# Matched sample summary
if (!is.null(effs_mtch.dt) && nrow(effs_mtch.dt) > 0) {
  log_info("\nMatched sample effect sizes (bilateral):")
  for (roi in names(ROIS)) {
    for (adj in names(ADJS)) {
      eff <- effs_mtch.dt[ROI == roi & ADJ == adj & SIDE == "LR"]
      if (nrow(eff) > 0) {
        log_info("  %s (%s): d = %.3f [%.3f, %.3f] - %s (N pairs = %d)",
                 ROIS[[roi]], ADJS[[adj]],
                 eff$ESTIMATE, eff$CI_LOWER, eff$CI_UPPER,
                 eff$MAGNITUDE, eff$N_PAIRS)
      }
    }
  }
}

# Matched vs Full comparison
if (!is.null(comparison.dt) && nrow(comparison.dt) > 0) {
  log_info("\nMatched vs Full Sample comparison:")
  for (roi in names(ROIS)) {
    log_info("  %s:", ROIS[[roi]])
    for (adj in names(ADJS)) {
      comp <- comparison.dt[ROI == roi & ADJ == adj]
      if (nrow(comp) > 0) {
        log_info("    %s: Full d = %.3f, Matched d = %.3f, Diff = %.3f",
                 ADJS[[adj]],
                 comp$ESTIMATE_ALL, comp$ESTIMATE_MTCH, comp$DIFF)
      }
    }
  }
}

# Paired effect sizes
if (!is.null(effs_paired.dt) && nrow(effs_paired.dt) > 0) {
  log_info("\nPaired effect sizes (within-pair differences):")
  for (roi in names(ROIS)) {
    for (adj in names(ADJS)) {
      eff <- effs_paired.dt[ROI == roi & ADJ == adj]
      if (nrow(eff) > 0) {
        log_info("  %s (%s): d = %.3f [%.3f, %.3f], t = %.2f, p = %.2e (N = %d pairs)",
                 ROIS[[roi]], ADJS[[adj]],
                 eff$ESTIMATE_PAIRED, eff$CI_LOWER_PAIRED, eff$CI_UPPER_PAIRED,
                 eff$T_STAT, eff$P_VALUE, eff$N_PAIRS)
      }
    }
  }
}

# F-code sensitivity summary
if (!is.null(effs_sens.dt) && nrow(effs_sens.dt) > 0) {
  log_info("\nSensitivity analysis effect sizes (F-codes included, bilateral):")
  for (roi in names(ROIS)) {
    for (adj in names(ADJS)) {
      eff <- effs_sens.dt[ROI == roi & ADJ == adj & SIDE == "LR"]
      if (nrow(eff) > 0) {
        log_info("  %s (%s): d = %.3f [%.3f, %.3f] - %s",
                 ROIS[[roi]], ADJS[[adj]],
                 eff$ESTIMATE, eff$CI_LOWER, eff$CI_UPPER,
                 eff$MAGNITUDE)
      }
    }
  }
}

# Primary vs Sensitivity comparison
if (!is.null(sens_comparison.dt) && nrow(sens_comparison.dt) > 0) {
  log_info("\nPrimary vs Sensitivity comparison:")
  for (roi in names(ROIS)) {
    log_info("  %s:", ROIS[[roi]])
    for (adj in names(ADJS)) {
      comp <- sens_comparison.dt[ROI == roi & ADJ == adj]
      if (nrow(comp) > 0) {
        log_info("    %s: Primary d = %.3f, Sensitivity d = %.3f, Diff = %.3f",
                 ADJS[[adj]],
                 comp$ESTIMATE_PRIMARY, comp$ESTIMATE_SENS, comp$DIFF)
      }
    }
  }
}

# ===========================================================================
# Export Manuscript Assets
# ===========================================================================
log_section("Exporting manuscript figures")

# Prepare data environment for figure generation
# The _common.R functions need: sex_diff and brain_data
export_data <- list(
  sex_diff = sex_diffs.lst,
  brain_data = hc_hvr.lst[[NETWORK]]$CRS$ALL
)

# Export figures using _common.R functions with config-defined names
export_manuscript_figures("sex_differences", export_data)

log_info("Manuscript assets exported to: %s", get_output_path("figures"))

log_script_end("09_sex_differences.R", success = TRUE)
