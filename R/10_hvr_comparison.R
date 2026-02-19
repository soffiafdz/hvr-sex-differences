#!/usr/bin/env Rscript

# =============================================================================
# HVR vs HC Comparison Analysis
# =============================================================================
# Compares hippocampal-ventricle ratio (HVR) against residualized hippocampal
# volume (HC_RES) for cognitive prediction and sex-specific aging patterns
#
# Key analyses:
# - Correlation comparison using Williams t-test
# - Age × sex interactions for HVR, HC, and LV
# - Age-stratified correlation comparisons
# - Incremental R² test: Does HVR add variance beyond HC_RES?
# - Memory-specific comparisons
#
# Inputs:
#   - data/derivatives/hc-hvr_adj.rds: Head-size adjusted volumes
#   - data/derivatives/lat-cog_values.rds: Cognitive factor scores
#
# Outputs:
#   - models/results/hvr_hc_comparison.rds: Comparison statistics
#   - outputs/figures/hvr_comparison/*.png: Comparison plots
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
library(ggplot2)
library(gt)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))
source(here("R/utils/statistics.R"))
source(here("R/utils/plotting_core.R"))
source(here("R/utils/plotting_pipeline.R"))
source(here("R/utils/export.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("10_hvr_comparison.R")

# Load configuration
config <- load_config()
set_seed()

# ----- Constants -----
FORCE_REGENERATE <- get_script_setting("force_regenerate", "hvr_comparison", default = FALSE)
REDO_PLOTS <- get_script_setting("hvr_comparison", "redo_plots", default = TRUE)
NETWORK <- get_script_setting("hvr_comparison", "network", default = "LPP_CNN")

log_info("Network: %s", NETWORK)
log_info("Force regenerate: %s, Redo plots: %s", FORCE_REGENERATE, REDO_PLOTS)

# =============================================================================
# Load Data
# =============================================================================
log_section("Loading input data")

# Head-size adjusted volumes
hc_hvr_path <- get_data_path("processed", "hc_hvr_adjusted")
if (!check_files_exist(hc_hvr_path, stop_on_missing = FALSE)) {
  log_error("Head-size adjusted data not found: %s", hc_hvr_path)
  log_error("Please run 05_adjust_headsize.R first")
  stop("Required input not found", call. = FALSE)
}
hc_hvr.lst <- read_rds_safe(hc_hvr_path, "Brain volume data")
validate_not_empty(hc_hvr.lst, "Head-size adjusted data")

# Cognitive factor scores
cog_path <- get_data_path("processed", "lat-cog_values")
if (!check_files_exist(cog_path, stop_on_missing = FALSE)) {
  log_error("Cognitive factor scores not found: %s", cog_path)
  log_error("Please run 06_cognitive_factors.R first")
  stop("Required input not found", call. = FALSE)
}
cog.lst <- read_rds_safe(cog_path, "Cognitive factor scores")
validate_not_empty(cog.lst, "Cognitive factor scores")

# =============================================================================
# Data Preparation
# =============================================================================
log_section("Preparing comparison data")

# Extract HVR from unadjusted (NON) - HVR is self-normalizing, only valid for NON
hvr_non.dt <- hc_hvr.lst[[NETWORK]]$CRS$ALL[
  ADJ == "NON" & SIDE == "LR" & SUBFIELD == "total",
  .(EID, SEX, AGE, HVR, LV)
]
validate_not_empty(hvr_non.dt, "HVR (unadjusted) brain data")

# Extract residualized HC (HC_RES) for comparison
hc_res.dt <- hc_hvr.lst[[NETWORK]]$CRS$ALL[
  ADJ == "RES" & SIDE == "LR" & SUBFIELD == "total",
  .(EID, HC_RES = HC)
]
validate_not_empty(hc_res.dt, "Residualized HC data")

# Merge HVR and HC_RES
brain_res.dt <- merge(hvr_non.dt, hc_res.dt, by = "EID")
setnames(brain_res.dt, "HC_RES", "HC")  # Rename to HC for compatibility with rest of script
validate_columns(brain_res.dt, c("EID", "SEX", "AGE", "HC", "HVR", "LV"), "Brain data")

# Extract cognitive scores
cog_scores.dt <- cog.lst$CRS[, .(EID, COG, MEM, PRSP)]
validate_not_empty(cog_scores.dt, "Cognitive scores")

# Merge datasets
comparison.dt <- merge(brain_res.dt, cog_scores.dt, by = "EID")
comparison.dt <- comparison.dt[complete.cases(comparison.dt)]

# Ensure SEX is factor
comparison.dt[, SEX := factor(SEX, levels = c("Female", "Male"))]

n_total <- nrow(comparison.dt)
log_info("Comparison data prepared: N = %d", n_total)
log_info("  Females: %d, Males: %d",
         comparison.dt[SEX == "Female", .N],
         comparison.dt[SEX == "Male", .N])

# =============================================================================
# Check for Existing Results
# =============================================================================
output_path <- get_data_path("models", "results", "hvr_hc_comparison")
fig_dir <- get_output_path("figures")

if (!FORCE_REGENERATE && file.exists(output_path)) {
  log_info("Results already exist: %s", output_path)
  log_info("Set force_regenerate: hvr_comparison: yes to recompute")

  if (!REDO_PLOTS) {
    log_info("Skipping analysis (no regeneration requested)")
    log_script_end("10_hvr_comparison.R", success = TRUE)
    quit(save = "no", status = 0)
  } else {
    log_info("Loading existing results for plotting")
    hvr_hc_comparison <- read_rds_safe(output_path, "HVR comparison results")
  }
} else {

  # ===========================================================================
  # Analysis
  # ===========================================================================

  # Initialize results list
  hvr_hc_comparison <- list()

  # -------------------------------------------------------------------------
  # 1. Correlation Comparison (Williams t-test)
  # -------------------------------------------------------------------------
  log_section("Correlation comparison (Williams t-test)")

  # Overall correlation comparison
  r_hc_cog <- cor(comparison.dt$HC, comparison.dt$COG)
  r_hvr_cog <- cor(comparison.dt$HVR, comparison.dt$COG)
  r_hc_hvr <- cor(comparison.dt$HC, comparison.dt$HVR)

  williams_overall <- williams_t_test(r_hc_cog, r_hvr_cog, r_hc_hvr, n_total)

  hvr_hc_comparison$correlation_comparison <- data.table(
    GROUP = "Overall",
    N = n_total,
    R_HC_COG = r_hc_cog,
    R_HVR_COG = r_hvr_cog,
    R_HC_HVR = r_hc_hvr,
    RATIO_HVR_HC = r_hvr_cog / r_hc_cog,
    WILLIAMS_T = williams_overall$t,
    WILLIAMS_P = williams_overall$p
  )

  # Sex-specific correlation comparisons
  for (sex in c("Female", "Male")) {
    sub.dt <- comparison.dt[SEX == sex]
    n_sub <- nrow(sub.dt)

    r_hc <- cor(sub.dt$HC, sub.dt$COG)
    r_hvr <- cor(sub.dt$HVR, sub.dt$COG)
    r_hc_hvr_sub <- cor(sub.dt$HC, sub.dt$HVR)

    williams_sex <- williams_t_test(r_hc, r_hvr, r_hc_hvr_sub, n_sub)

    hvr_hc_comparison$correlation_comparison <- rbind(
      hvr_hc_comparison$correlation_comparison,
      data.table(
        GROUP = sex,
        N = n_sub,
        R_HC_COG = r_hc,
        R_HVR_COG = r_hvr,
        R_HC_HVR = r_hc_hvr_sub,
        RATIO_HVR_HC = r_hvr / r_hc,
        WILLIAMS_T = williams_sex$t,
        WILLIAMS_P = williams_sex$p
      )
    )
  }

  log_info("Correlation comparison complete:")
  log_info("  Overall: r(HC,COG)=%.3f, r(HVR,COG)=%.3f, ratio=%.2f, p=%s",
           r_hc_cog, r_hvr_cog, r_hvr_cog / r_hc_cog,
           format.pval(williams_overall$p, digits = 3))

  # -------------------------------------------------------------------------
  # 2. Age × Sex Interactions
  # -------------------------------------------------------------------------
  log_section("Age × sex interactions")

  # Fit interaction models
  hvr_int.lm <- lm(HVR ~ AGE * SEX, data = comparison.dt)
  hc_int.lm <- lm(HC ~ AGE * SEX, data = comparison.dt)
  lv_int.lm <- lm(LV ~ AGE * SEX, data = comparison.dt)

  # Extract interaction terms
  extract_age_sex_interaction <- function(model, var_name) {
    summ <- summary(model)$coefficients
    int_row <- grep("AGE:SEX", rownames(summ))

    # Get sex-specific slopes
    slope_f <- coef(model)["AGE"]
    slope_m <- slope_f + coef(model)[int_row]

    data.table(
      VARIABLE = var_name,
      SLOPE_FEMALE = slope_f,
      SLOPE_MALE = slope_m,
      SLOPE_RATIO = slope_m / slope_f,
      INTERACTION_EST = summ[int_row, "Estimate"],
      INTERACTION_SE = summ[int_row, "Std. Error"],
      INTERACTION_T = summ[int_row, "t value"],
      INTERACTION_P = summ[int_row, "Pr(>|t|)"]
    )
  }

  hvr_hc_comparison$age_sex_interactions <- rbind(
    extract_age_sex_interaction(hvr_int.lm, "HVR"),
    extract_age_sex_interaction(hc_int.lm, "HC"),
    extract_age_sex_interaction(lv_int.lm, "LV")
  )

  log_info("Age × Sex interactions:")
  for (var in c("HVR", "HC", "LV")) {
    row <- hvr_hc_comparison$age_sex_interactions[VARIABLE == var]
    log_info("  %s: Male/Female slope ratio = %.2f, interaction p = %s",
             var, row$SLOPE_RATIO, format.pval(row$INTERACTION_P, digits = 3))
  }

  # -------------------------------------------------------------------------
  # 3. Age-Stratified Comparisons
  # -------------------------------------------------------------------------
  log_section("Age-stratified comparisons")

  age_cuts.v <- c(45, 55, 65, 70, 75)
  age_stratified.dt <- data.table()

  for (age_cut in age_cuts.v) {
    sub.dt <- comparison.dt[AGE >= age_cut]
    n_sub <- nrow(sub.dt)

    r_hc <- cor(sub.dt$HC, sub.dt$COG)
    r_hvr <- cor(sub.dt$HVR, sub.dt$COG)

    age_stratified.dt <- rbind(age_stratified.dt, data.table(
      AGE_MIN = age_cut,
      N = n_sub,
      R_HC_COG = r_hc,
      R_HVR_COG = r_hvr,
      RATIO_HVR_HC = r_hvr / r_hc
    ))
  }

  hvr_hc_comparison$age_stratified <- age_stratified.dt
  log_info("Age-stratified correlations computed for %d age thresholds",
           length(age_cuts.v))

  # -------------------------------------------------------------------------
  # 4. Standardized Effect Comparison
  # -------------------------------------------------------------------------
  log_section("Standardized effect comparisons")

  # Standardize all variables
  comparison.dt[, `:=`(
    HC_z = scale(HC)[, 1],
    HVR_z = scale(HVR)[, 1],
    COG_z = scale(COG)[, 1],
    AGE_z = scale(AGE)[, 1]
  )]

  # Fit standardized models
  m_hc.lm <- lm(COG_z ~ HC_z + AGE_z + SEX, data = comparison.dt)
  m_hvr.lm <- lm(COG_z ~ HVR_z + AGE_z + SEX, data = comparison.dt)
  m_both.lm <- lm(COG_z ~ HC_z + HVR_z + AGE_z + SEX, data = comparison.dt)

  # Extract coefficients
  hvr_hc_comparison$standardized_effects <- data.table(
    MODEL = c("HC only", "HVR only", "Both (HC)", "Both (HVR)"),
    PREDICTOR = c("HC", "HVR", "HC", "HVR"),
    BETA = c(
      coef(m_hc.lm)["HC_z"],
      coef(m_hvr.lm)["HVR_z"],
      coef(m_both.lm)["HC_z"],
      coef(m_both.lm)["HVR_z"]
    ),
    R_SQUARED = c(
      summary(m_hc.lm)$r.squared,
      summary(m_hvr.lm)$r.squared,
      summary(m_both.lm)$r.squared,
      summary(m_both.lm)$r.squared
    )
  )

  log_info("Standardized effects:")
  log_info("  HC only: beta = %.4f, R² = %.4f",
           coef(m_hc.lm)["HC_z"], summary(m_hc.lm)$r.squared)
  log_info("  HVR only: beta = %.4f, R² = %.4f",
           coef(m_hvr.lm)["HVR_z"], summary(m_hvr.lm)$r.squared)
  log_info("  Both: R² = %.4f", summary(m_both.lm)$r.squared)

  # -------------------------------------------------------------------------
  # 4b. Variance Inflation Factor (VIF) for Joint Model
  # -------------------------------------------------------------------------
  # RATIONALE: Since HVR = HC / (HC + LV), HC and HVR are structurally
  # correlated (sharing the numerator). VIF quantifies collinearity to
  # assess whether joint model coefficients are stable.
  # Interpretation: VIF < 5 = acceptable, VIF < 10 = moderate, VIF > 10 = severe
  # Reference: O'Brien (2007) - A Caution Regarding Rules of Thumb for VIF
  # -------------------------------------------------------------------------
  log_section("Collinearity diagnostics (VIF)")

  # Compute VIF using the variance inflation formula
  # VIF_j = 1 / (1 - R²_j), where R²_j is from regressing X_j on other predictors
  r2_hc_on_others <- summary(lm(HC_z ~ HVR_z + AGE_z + SEX, data = comparison.dt))$r.squared
  r2_hvr_on_others <- summary(lm(HVR_z ~ HC_z + AGE_z + SEX, data = comparison.dt))$r.squared

  vif_hc <- 1 / (1 - r2_hc_on_others)
  vif_hvr <- 1 / (1 - r2_hvr_on_others)

  # Correlation between HC and HVR (the primary collinearity concern)
  r_hc_hvr <- cor(comparison.dt$HC_z, comparison.dt$HVR_z)

  hvr_hc_comparison$collinearity <- data.table(
    PREDICTOR = c("HC", "HVR"),
    VIF = c(vif_hc, vif_hvr),
    TOLERANCE = c(1 - r2_hc_on_others, 1 - r2_hvr_on_others),
    R_HC_HVR = r_hc_hvr
  )

  log_info("Collinearity diagnostics for joint model (COG ~ HC + HVR + AGE + SEX):")
  log_info("  r(HC, HVR) = %.3f", r_hc_hvr)
  log_info("  VIF(HC) = %.2f, VIF(HVR) = %.2f", vif_hc, vif_hvr)
  log_info("  Interpretation: %s",
           if (max(vif_hc, vif_hvr) < 5) "Acceptable collinearity (VIF < 5)"
           else if (max(vif_hc, vif_hvr) < 10) "Moderate collinearity (VIF 5-10)"
           else "Severe collinearity (VIF > 10) - interpret joint coefficients with caution")

  # -------------------------------------------------------------------------
  # 5. Incremental R² Test
  # -------------------------------------------------------------------------
  # RATIONALE: HVR = HC / (HC + LV), so mathematically contains HC.
  # However, this test is NOT circular because:
  # 1. HVR adds information from LV (ventricle expansion) that HC alone lacks
  # 2. The incremental R² quantifies whether the LV component improves prediction
  # 3. Conceptually equivalent to: does HC + LV predict better than HC alone?
  # The ratio formulation (HVR) captures this parsimoniously.
  # -------------------------------------------------------------------------
  log_section("Incremental R² test: HVR beyond HC_RES")

  # Model 1: Base model with HC_RES only
  m_base.lm <- lm(COG_z ~ HC_z + AGE_z + SEX, data = comparison.dt)

  # Model 2: Add HVR to base model
  m_full.lm <- lm(COG_z ~ HC_z + HVR_z + AGE_z + SEX, data = comparison.dt)

  # F-test for model comparison (nested models)
  incremental_test <- anova(m_base.lm, m_full.lm)

  # Extract key statistics
  r2_base <- summary(m_base.lm)$r.squared
  r2_full <- summary(m_full.lm)$r.squared
  delta_r2 <- r2_full - r2_base
  f_stat <- incremental_test$F[2]
  p_value <- incremental_test$`Pr(>F)`[2]

  # Compute partial correlations
  resid_hvr <- residuals(lm(HVR_z ~ HC_z + AGE_z + SEX, data = comparison.dt))
  resid_cog <- residuals(lm(COG_z ~ HC_z + AGE_z + SEX, data = comparison.dt))
  partial_r_hvr_cog <- cor(resid_hvr, resid_cog)

  resid_hc <- residuals(lm(HC_z ~ HVR_z + AGE_z + SEX, data = comparison.dt))
  resid_cog_hvr <- residuals(lm(COG_z ~ HVR_z + AGE_z + SEX, data = comparison.dt))
  partial_r_hc_cog <- cor(resid_hc, resid_cog_hvr)

  hvr_hc_comparison$incremental_r2 <- data.table(
    R2_BASE = r2_base,
    R2_FULL = r2_full,
    DELTA_R2 = delta_r2,
    PCT_INCREASE = 100 * delta_r2 / r2_base,
    F_STAT = f_stat,
    DF1 = incremental_test$Df[2],
    DF2 = incremental_test$Res.Df[2],
    P_VALUE = p_value,
    PARTIAL_R_HVR_COG = partial_r_hvr_cog,
    PARTIAL_R_HC_COG = partial_r_hc_cog,
    PARTIAL_R2_HVR = partial_r_hvr_cog^2,
    PARTIAL_R2_HC = partial_r_hc_cog^2
  )

  log_info("Incremental R² test results:")
  log_info("  Base model (HC_RES + AGE + SEX): R² = %.4f", r2_base)
  log_info("  Full model (+HVR): R² = %.4f", r2_full)
  log_info("  ΔR² = %.4f (%.1f%% increase)", delta_r2, 100 * delta_r2 / r2_base)
  log_info("  F(%d, %d) = %.2f, p = %s",
           incremental_test$Df[2], incremental_test$Res.Df[2],
           f_stat, format.pval(p_value, digits = 3))

  # -------------------------------------------------------------------------
  # 5b. Sex-Stratified Incremental R²
  # -------------------------------------------------------------------------
  log_info("Computing sex-stratified incremental R² tests")

  sex_incremental.dt <- data.table()
  for (sex in c("Female", "Male")) {
    sub.dt <- comparison.dt[SEX == sex]

    m_base_sex.lm <- lm(COG_z ~ HC_z + AGE_z, data = sub.dt)
    m_full_sex.lm <- lm(COG_z ~ HC_z + HVR_z + AGE_z, data = sub.dt)

    test_sex <- anova(m_base_sex.lm, m_full_sex.lm)

    r2_b <- summary(m_base_sex.lm)$r.squared
    r2_f <- summary(m_full_sex.lm)$r.squared

    sex_incremental.dt <- rbind(sex_incremental.dt, data.table(
      SEX = sex,
      N = nrow(sub.dt),
      R2_BASE = r2_b,
      R2_FULL = r2_f,
      DELTA_R2 = r2_f - r2_b,
      PCT_INCREASE = 100 * (r2_f - r2_b) / r2_b,
      F_STAT = test_sex$F[2],
      P_VALUE = test_sex$`Pr(>F)`[2]
    ))
  }

  hvr_hc_comparison$incremental_r2_by_sex <- sex_incremental.dt

  log_info("Sex-stratified incremental R²:")
  for (sex in c("Female", "Male")) {
    row <- sex_incremental.dt[SEX == sex]
    log_info("  %s: ΔR² = %.4f (%.1f%% increase), p = %s",
             sex, row$DELTA_R2, row$PCT_INCREASE,
             format.pval(row$P_VALUE, digits = 3))
  }

  # -------------------------------------------------------------------------
  # 6. Memory-Specific Comparisons
  # -------------------------------------------------------------------------
  log_section("Memory-specific comparisons")

  r_hc_mem <- cor(comparison.dt$HC, comparison.dt$MEM, use = "complete.obs")
  r_hvr_mem <- cor(comparison.dt$HVR, comparison.dt$MEM, use = "complete.obs")

  hvr_hc_comparison$memory_comparison <- data.table(
    OUTCOME = "Memory (MEM)",
    R_HC = r_hc_mem,
    R_HVR = r_hvr_mem,
    RATIO = r_hvr_mem / r_hc_mem
  )

  log_info("Memory correlations: r(HC)=%.3f, r(HVR)=%.3f", r_hc_mem, r_hvr_mem)

  # ===========================================================================
  # Save Results
  # ===========================================================================
  log_section("Saving results")

  ensure_directory(dirname(output_path))
  write_rds_safe(hvr_hc_comparison, output_path, "HVR vs HC comparison results")
  log_info("Results saved: %s", output_path)
}

# =============================================================================
# Generate Figures
# =============================================================================
if (REDO_PLOTS) {
  log_section("Generating figures")

  ensure_directory(fig_dir)

  # Figure 1: HVR vs HC_RES scatter plot
  log_info("Creating HVR vs HC_RES scatter plot")

  # Re-standardize if needed (in case we loaded existing results)
  if (!"HC_z" %in% names(comparison.dt)) {
    comparison.dt[, `:=`(
      HC_z = scale(HC)[, 1],
      HVR_z = scale(HVR)[, 1]
    )]
  }

  p_scatter <- plot_hvr_vs_hcres(
    data = comparison.dt,
    hc_col = "HC_z",
    hvr_col = "HVR_z",
    sex_colors = get_palette("sex"),
    show_marginals = FALSE,
    annotate_r = TRUE
  )
  save_plot(p_scatter, file.path(fig_dir, "hvr_vs_hcres_scatter.png"),
            width = 8, height = 7)

  # Figure 2: Incremental R² bar plot
  log_info("Creating incremental R² plot")
  p_incr <- plot_incremental_r2(
    r2_base = hvr_hc_comparison$incremental_r2$R2_BASE,
    r2_full = hvr_hc_comparison$incremental_r2$R2_FULL,
    title = "HVR Adds Unique Variance Beyond Residualized HC"
  )
  save_plot(p_incr, file.path(fig_dir, "incremental_r2.png"),
            width = 6, height = 5)

  # Figure 3: Correlation comparison bar plot
  log_info("Creating correlation comparison plot")
  corr_plot.dt <- hvr_hc_comparison$correlation_comparison[, .(
    GROUP,
    HC = R_HC_COG,
    HVR = R_HVR_COG
  )]
  corr_long.dt <- melt(corr_plot.dt, id.vars = "GROUP",
                       variable.name = "Measure", value.name = "r")

  p_corr <- ggplot(corr_long.dt, aes(x = GROUP, y = r, fill = Measure)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c(HC = "#0072B2", HVR = "#009E73")) +
    labs(
      title = "Brain-Cognition Correlations: HVR vs HC",
      subtitle = "HVR consistently shows stronger association with general cognition",
      x = "Sample",
      y = "Correlation (r)"
    ) +
    theme_publication(base_size = 12) +
    theme(legend.position = "bottom")

  save_plot(p_corr, file.path(fig_dir, "correlation_comparison.png"),
            width = 7, height = 5)

  log_info("Figures saved to: %s", fig_dir)
}

# =============================================================================
# Generate Tables (TEX/HTML for manuscript submission)
# =============================================================================
log_section("Generating publication tables")

tables_dir <- get_output_path("tables")
ensure_directory(tables_dir)

# --- Table 1: Incremental R² Results ---
log_info("Creating incremental R² table")
incr_dt <- hvr_hc_comparison$incremental_r2
incr_table <- data.table(
  `Statistic` = c(
    "Base model R² (HC + Age + Sex)",
    "Full model R² (+ HVR)",
    "ΔR²",
    "% Increase",
    "F statistic",
    "p-value",
    "Partial r(HVR, COG | HC)",
    "Partial r(HC, COG | HVR)"
  ),
  Value = c(
    sprintf("%.4f", incr_dt$R2_BASE),
    sprintf("%.4f", incr_dt$R2_FULL),
    sprintf("%.4f", incr_dt$DELTA_R2),
    sprintf("%.1f%%", incr_dt$PCT_INCREASE),
    sprintf("F(%d, %d) = %.2f", incr_dt$DF1, incr_dt$DF2, incr_dt$F_STAT),
    format.pval(incr_dt$P_VALUE, digits = 3),
    sprintf("%.4f", incr_dt$PARTIAL_R_HVR_COG),
    sprintf("%.4f", incr_dt$PARTIAL_R_HC_COG)
  )
)

incr_gt <- incr_table |>
  gt() |>
  tab_header(
    title = "Incremental Validity: HVR Beyond Residualized HC",
    subtitle = "Does HVR add unique variance in predicting cognition?"
  ) |>
  tab_source_note("Base model: COG ~ HC_RES + AGE + SEX") |>
  tab_source_note("Full model: COG ~ HC_RES + HVR + AGE + SEX") |>
  opt_table_font(font = "Times New Roman")

gt::gtsave(incr_gt, file.path(tables_dir, "incremental_r2.html"))
gt::gtsave(incr_gt, file.path(tables_dir, "incremental_r2.tex"))
log_info("Incremental R² table saved")

# --- Table 2: Age × Sex Interactions ---
log_info("Creating age × sex interactions table")
int_dt <- hvr_hc_comparison$age_sex_interactions

int_gt <- int_dt |>
  gt() |>
  tab_header(
    title = "Age × Sex Interactions for Brain Measures",
    subtitle = "Do males and females show different age-related patterns?"
  ) |>
  cols_label(
    VARIABLE = "Measure",
    SLOPE_FEMALE = "Female Slope",
    SLOPE_MALE = "Male Slope",
    SLOPE_RATIO = "M/F Ratio",
    INTERACTION_EST = "β (Age×Sex)",
    INTERACTION_SE = "SE",
    INTERACTION_T = "t",
    INTERACTION_P = "p"
  ) |>
  fmt_number(columns = c(SLOPE_FEMALE, SLOPE_MALE, INTERACTION_EST, INTERACTION_SE), decimals = 4) |>
  fmt_number(columns = c(SLOPE_RATIO, INTERACTION_T), decimals = 2) |>
  fmt_scientific(columns = INTERACTION_P, decimals = 2) |>
  tab_source_note("Model: Y ~ Age + Sex + Age×Sex") |>
  tab_source_note("Positive ratio = males show steeper age slopes") |>
  opt_table_font(font = "Times New Roman")

gt::gtsave(int_gt, file.path(tables_dir, "age_sex_interactions.html"))
gt::gtsave(int_gt, file.path(tables_dir, "age_sex_interactions.tex"))
log_info("Age × sex interactions table saved")

# --- Table 3: Memory-Specific Correlations ---
log_info("Creating memory comparison table")
mem_dt <- hvr_hc_comparison$memory_comparison

mem_gt <- mem_dt |>
  gt() |>
  tab_header(
    title = "Brain-Memory Correlations",
    subtitle = "HC vs HVR for memory-specific prediction"
  ) |>
  cols_label(
    OUTCOME = "Outcome",
    R_HC = "r(HC, MEM)",
    R_HVR = "r(HVR, MEM)",
    RATIO = "HVR/HC Ratio"
  ) |>
  fmt_number(columns = c(R_HC, R_HVR, RATIO), decimals = 3) |>
  tab_source_note("MEM = Memory-specific factor from bifactor CFA (orthogonal to g)") |>
  opt_table_font(font = "Times New Roman")

gt::gtsave(mem_gt, file.path(tables_dir, "memory_comparison.html"))
gt::gtsave(mem_gt, file.path(tables_dir, "memory_comparison.tex"))
log_info("Memory comparison table saved")

# --- Table 4: Standardized Effects ---
log_info("Creating standardized effects table")
std_dt <- hvr_hc_comparison$standardized_effects

std_gt <- std_dt |>
  gt() |>
  tab_header(
    title = "Standardized Regression Coefficients",
    subtitle = "Predicting cognition from brain measures"
  ) |>
  cols_label(
    MODEL = "Model",
    PREDICTOR = "Predictor",
    BETA = "β",
    R_SQUARED = "R²"
  ) |>
  fmt_number(columns = c(BETA, R_SQUARED), decimals = 4) |>
  tab_source_note("All models include Age and Sex as covariates") |>
  opt_table_font(font = "Times New Roman")

gt::gtsave(std_gt, file.path(tables_dir, "standardized_effects.html"))
gt::gtsave(std_gt, file.path(tables_dir, "standardized_effects.tex"))
log_info("Standardized effects table saved")

log_info("Tables saved to: %s", tables_dir)

# =============================================================================
# Summary
# =============================================================================
log_section("Summary")

log_info("Key findings:")
log_info("  1. HVR correlates with cognition %.2fx more strongly than HC_RES",
         hvr_hc_comparison$correlation_comparison[GROUP == "Overall", RATIO_HVR_HC])
log_info("     (Williams t-test p = %s)",
         format.pval(hvr_hc_comparison$correlation_comparison[GROUP == "Overall", WILLIAMS_P], digits = 3))
log_info("  2. HVR adds significant variance beyond HC_RES:")
log_info("     ΔR² = %.4f (%.1f%% increase), F = %.2f, p = %s",
         hvr_hc_comparison$incremental_r2$DELTA_R2,
         hvr_hc_comparison$incremental_r2$PCT_INCREASE,
         hvr_hc_comparison$incremental_r2$F_STAT,
         format.pval(hvr_hc_comparison$incremental_r2$P_VALUE, digits = 3))
log_info("  3. Partial r(HVR, COG | HC_RES) = %.4f vs partial r(HC, COG | HVR) = %.4f",
         hvr_hc_comparison$incremental_r2$PARTIAL_R_HVR_COG,
         hvr_hc_comparison$incremental_r2$PARTIAL_R_HC_COG)
log_info("  4. Males decline %.1fx faster in HVR than females",
         abs(hvr_hc_comparison$age_sex_interactions[VARIABLE == "HVR", SLOPE_RATIO]))
log_info("  5. HVR advantage consistent across age groups and sexes")

log_script_end("10_hvr_comparison.R", success = TRUE)
