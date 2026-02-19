# =============================================================================
# GT Table Formatting Functions
# =============================================================================
# Functions for formatting analysis outputs as gt tables for display.
# Part of the R/utils/ module system.
# =============================================================================

#' Apply manuscript-style formatting to gt tables
#'
#' @param gt_tbl A gt table object
#' @param table_type Type: "main" (larger), "supplementary" (smaller), "compact" (smallest)
#' @return Formatted gt table
#' @export
style_manuscript_table <- function(gt_tbl, table_type = "main") {
  # Font sizes based on table type (in pt for PDF compatibility)
  if (table_type == "main") {
    title_size <- 12
    subtitle_size <- 10
    body_size <- 9
    footnote_size <- 8
    source_size <- 8
  } else if (table_type == "supplementary") {
    title_size <- 11
    subtitle_size <- 9
    body_size <- 8
    footnote_size <- 7
    source_size <- 7
  } else if (table_type == "compact") {
    title_size <- 10
    subtitle_size <- 8
    body_size <- 7
    footnote_size <- 6
    source_size <- 6
  } else {
    stop("table_type must be 'main', 'supplementary', or 'compact'")
  }

  gt_tbl |>
    gt::tab_options(
      table.font.size = gt::px(body_size),
      quarto.disable_processing = TRUE,
      heading.title.font.size = gt::px(title_size),
      heading.title.font.weight = "bold",
      heading.subtitle.font.size = gt::px(subtitle_size),
      column_labels.font.size = gt::px(body_size),
      column_labels.font.weight = "bold",
      footnotes.font.size = gt::px(footnote_size),
      footnotes.multiline = FALSE,
      source_notes.font.size = gt::px(source_size),
      row_group.font.size = gt::px(body_size),
      row_group.font.weight = "bold",
      table.border.top.style = "solid",
      table.border.bottom.style = "solid",
      heading.border.bottom.style = "solid"
    ) |>
    # Right-align all columns by default (most are numeric or formatted numeric)
    gt::cols_align(align = "right", columns = everything()) |>
    # Left-align only the first column (row labels/identifiers)
    gt::cols_align(align = "left", columns = 1)
}

#' Alias for backward compatibility
#' @param gt_tbl A gt table object
#' @return Formatted gt table
#' @export
gt_pdf_style <- function(gt_tbl) {
  style_manuscript_table(gt_tbl, table_type = "main")
}

#' Format centile table for display with gt
#' @param centile_dt Centile data.table from get_centile_table()
#' @param title Table title
#' @param ages_subset Optional vector of ages to display (default: every 5 years)
#' @return gt table object
#' @export
format_centile_gt <- function(centile_dt, title = "Normative Centiles",
                               ages_subset = seq(50, 80, by = 5),
                               decimals = 3) {
  if (is.null(centile_dt) || nrow(centile_dt) == 0) {
    return(gt(data.table(Message = "Data not available")) |>
             tab_header(title = title))
  }

  # Make a copy to avoid modifying original
  dt <- copy(centile_dt)

  # Subset ages if specified
  if (!is.null(ages_subset)) {
    dt <- dt[AGE %in% ages_subset]
  }

  # Select key centiles for display (match actual column names in data)
  display_cols <- c("AGE", "p5", "p10", "p25", "p50", "p75", "p90", "p95")
  available_cols <- intersect(display_cols, names(dt))

  if (length(available_cols) <= 1) {
    return(gt(data.table(Message = "Centile columns not found")) |>
             tab_header(title = title))
  }

  # Build column labels dynamically based on available columns
  col_labels <- list(AGE = "Age")
  if ("p5" %in% available_cols) col_labels$p5 <- "5th"
  if ("p10" %in% available_cols) col_labels$p10 <- "10th"
  if ("p25" %in% available_cols) col_labels$p25 <- "25th"
  if ("p50" %in% available_cols) col_labels$p50 <- "50th"
  if ("p75" %in% available_cols) col_labels$p75 <- "75th"
  if ("p90" %in% available_cols) col_labels$p90 <- "90th"
  if ("p95" %in% available_cols) col_labels$p95 <- "95th"

  dt[, ..available_cols] |>
    gt() |>
    fmt_number(columns = -AGE, decimals = decimals) |>
    cols_label(.list = col_labels) |>
    tab_header(title = title) |>
    tab_spanner(label = "Percentile", columns = -AGE) |>
    gt_pdf_style()
}

#' Format hemisphere-specific centile table with Bilateral/L/R as row groups
#' @param norm_tables Normative tables from load_norm_tables()
#' @param roi ROI name (e.g., "HC", "HVR", "LV")
#' @param adj Adjustment method
#' @param sex Sex ("Female" or "Male")
#' @param title Table title
#' @param ages_subset Optional vector of ages to display
#' @param centiles Which centiles to show
#' @param decimals Number of decimal places
#' @return gt table object
#' @export
format_centile_hemisphere_gt <- function(norm_tables, roi, adj, sex,
                                          title = NULL,
                                          ages_subset = seq(45, 80, by = 5),
                                          centiles = c("p10", "p25", "p50", "p75", "p90"),
                                          decimals = 3) {
  if (is.null(norm_tables)) {
    return(gt(data.table(Message = "Data not available")) |>
             tab_header(title = title %||% paste(roi, "Centiles")))
  }

  # Get bilateral, left, and right data
  bilateral_dt <- tryCatch(copy(norm_tables[[sex]][[roi]][[adj]][["LR"]]), error = function(e) NULL)
  left_dt <- tryCatch(copy(norm_tables[[sex]][[roi]][[adj]][["L"]]), error = function(e) NULL)
  right_dt <- tryCatch(copy(norm_tables[[sex]][[roi]][[adj]][["R"]]), error = function(e) NULL)

  if (is.null(bilateral_dt) && is.null(left_dt) && is.null(right_dt)) {
    return(gt(data.table(Message = "Hemisphere data not available")) |>
             tab_header(title = title %||% paste(roi, "Centiles")))
  }

  # Find available centiles from first non-null table
  ref_dt <- bilateral_dt %||% left_dt %||% right_dt
  available_centiles <- intersect(centiles, names(ref_dt))
  if (length(available_centiles) == 0) {
    return(gt(data.table(Message = "Centile columns not found")) |>
             tab_header(title = title %||% paste(roi, "Centiles")))
  }

  # Helper to prepare each hemisphere's data
  prep_hemi <- function(dt, hemi_label) {
    if (is.null(dt)) return(NULL)
    dt <- copy(dt)
    if (!is.null(ages_subset)) dt <- dt[AGE %in% ages_subset]
    cols <- c("AGE", available_centiles)
    dt <- dt[, ..cols]
    dt[, Hemisphere := hemi_label]
    dt
  }

  # Stack all three
  combined <- rbindlist(list(
    prep_hemi(bilateral_dt, "Bilateral"),
    prep_hemi(left_dt, "Left"),
    prep_hemi(right_dt, "Right")
  ), fill = TRUE)

  combined[, Hemisphere := factor(Hemisphere, levels = c("Bilateral", "Left", "Right"))]
  setcolorder(combined, c("Hemisphere", "AGE", available_centiles))

  # Build column labels
  centile_labels <- c("p5" = "5th", "p10" = "10th", "p25" = "25th",
                      "p50" = "50th", "p75" = "75th", "p90" = "90th", "p95" = "95th")
  col_labels <- list(AGE = "Age")
  for (c in available_centiles) {
    col_labels[[c]] <- centile_labels[c]
  }

  default_title <- paste0(roi, " Normative Centiles: ", sex, "s (by Hemisphere)")

  combined |>
    gt(groupname_col = "Hemisphere") |>
    fmt_number(columns = all_of(available_centiles), decimals = decimals) |>
    cols_label(.list = col_labels) |>
    tab_header(title = title %||% default_title) |>
    tab_spanner(label = "Percentile", columns = all_of(available_centiles)) |>
    style_manuscript_table(table_type = "supplementary")
}

#' Format structural paths as gt table
#' @param paths_dt data.table from get_sem_structural_paths()
#' @return gt table object
#' @export
format_sem_paths_gt <- function(paths_dt) {
  if (is.null(paths_dt) || nrow(paths_dt) == 0) {
    return(gt(data.table(Message = "SEM paths not available")))
  }

  # Variable labels
  var_labels <- c(
    "g" = "General Cognition (g)",
    "MEM_s" = "Memory-Specific",
    "PRSP_s" = "Processing Speed-Specific",
    "HC" = "Hippocampus",
    "HC_RES" = "Hippocampus (Residualized)",
    "HVR" = "HVR",
    "AGE" = "Age",
    "ICV" = "ICV",
    "EDUC" = "Education",
    "IMDP" = "Deprivation"
  )

  dt <- copy(paths_dt)
  dt[, Outcome_Label := fifelse(Outcome %in% names(var_labels), var_labels[Outcome], Outcome)]
  dt[, Predictor_Label := fifelse(Predictor %in% names(var_labels), var_labels[Predictor], Predictor)]

  # Handle NA/Inf values gracefully
  dt[, `:=`(
    beta_fmt = fifelse(is.na(Std_Est) | is.infinite(Std_Est), "—", sprintf("%.3f", Std_Est)),
    ci_fmt = fifelse(is.na(CI_Lower) | is.na(CI_Upper) | is.infinite(CI_Lower) | is.infinite(CI_Upper),
                     "—", sprintf("[%.3f, %.3f]", CI_Lower, CI_Upper)),
    p_fmt = fifelse(is.na(P_Value) | is.infinite(P_Value), "—",
                    fifelse(P_Value < 0.001, sprintf("%.2e", P_Value), sprintf("%.3f", P_Value)))
  )]

  dt[, .(
    Outcome = Outcome_Label,
    Predictor = Predictor_Label,
    Sex,
    Beta_std = beta_fmt,
    CI_95 = ci_fmt,
    p = p_fmt
  )] |>
    gt() |>
    cols_label(Beta_std = "B (std)", CI_95 = "95% CI") |>
    tab_header(
      title = "Structural Path Estimates",
      subtitle = "Multi-group SEM with bootstrap CIs (2,000 resamples)"
    ) |>
    tab_spanner(
      label = "Parameter",
      columns = c(Outcome, Predictor, Sex)
    ) |>
    tab_spanner(
      label = "Estimate",
      columns = c(Beta_std, CI_95, p)
    ) |>
    gt_pdf_style()
}

#' Format SEM fit indices table for display
#' @param fit_dt data.table from load_sem_fit() or get_sem_fit_summary()
#' @return gt table object
#' @export
format_sem_fit_gt <- function(fit_dt) {
  if (is.null(fit_dt) || nrow(fit_dt) == 0) {
    return(gt(data.table(Message = "SEM fit data not available")))
  }

  # Create display-friendly model names and identify model type
  # POOLED models have suffix _POOLED, stratified models don't (multigroup)
  brain_labels <- c(
    "HC_COG" = "Hippocampus (Unadjusted)",
    "HC_RES_COG" = "Hippocampus (Residualized)",
    "HVR_COG" = "Hippocampal-to-Ventricle Ratio"
  )

  dt <- copy(fit_dt)

  # Determine if model is pooled or stratified based on MODEL name
  dt[, `:=`(
    IS_POOLED = grepl("_POOLED$", MODEL),
    BASE_MODEL = gsub("_POOLED$", "", MODEL)
  )]

  # Create labels
  dt[, Brain_Measure := fifelse(BASE_MODEL %in% names(brain_labels),
                                 brain_labels[BASE_MODEL], BASE_MODEL)]
  dt[, Analysis_Type := fifelse(IS_POOLED, "Pooled", "Sex-Stratified")]

  # Set factor order: HVR, HC (Residualized), HC (Unadjusted)
  dt[, Brain_Measure := factor(Brain_Measure,
                                levels = c("Hippocampal-to-Ventricle Ratio",
                                           "Hippocampus (Residualized)",
                                           "Hippocampus (Unadjusted)"))]
  setorder(dt, Analysis_Type, Brain_Measure)

  dt[, .(
    Analysis = Analysis_Type,
    `Brain Measure` = Brain_Measure,
    CFI = round(cfi, 3),
    TLI = round(tli, 3),
    RMSEA = sprintf("%.3f [%.3f, %.3f]", rmsea, `rmsea.ci.lower`, `rmsea.ci.upper`),
    SRMR = round(srmr, 3)
  )] |>
    gt(groupname_col = "Analysis") |>
    tab_header(title = "SEM Model Fit Indices") |>
    gt_pdf_style()
}

#' Format CFA fit indices table for display
#' @param fit_dt data.table from load_cfa_fit()
#' @return gt table object
#' @export
format_cfa_fit_gt <- function(fit_dt) {
  if (is.null(fit_dt) || nrow(fit_dt) == 0) {
    return(gt(data.table(Message = "CFA fit data not available")))
  }

  # Model labels - use N from data if available, otherwise just descriptive names
  base_labels <- c(
    "CRS" = "Cross-sectional",
    "LNG" = "Longitudinal"
  )

  dt <- copy(fit_dt)

  # Build label with N if available in data
  if ("N" %in% names(dt)) {
    dt[, Model := fifelse(
      MODEL %in% names(base_labels),
      paste0(base_labels[MODEL], " (N = ", format(N, big.mark = ","), ")"),
      MODEL
    )]
  } else {
    dt[, Model := fifelse(MODEL %in% names(base_labels), base_labels[MODEL], MODEL)]
  }

  dt[, .(
    Model,
    CFI = round(CFI, 3),
    TLI = round(TLI, 3),
    RMSEA = sprintf("%.3f [%.3f, %.3f]", RMSEA, RMSEA_CI_LOWER, RMSEA_CI_UPPER),
    SRMR = round(SRMR, 3)
  )] |>
    gt() |>
    tab_header(title = "Cognitive CFA Model Fit Indices") |>
    gt_pdf_style()
}

#' Format measurement invariance table for display
#' @param minv_dt data.table from load_measurement_invariance()
#' @return gt table object
#' @export
format_minv_gt <- function(minv_dt) {
  if (is.null(minv_dt) || nrow(minv_dt) == 0) {
    return(gt(data.table(Message = "Invariance data not available")))
  }

  # Model labels
  model_labels <- c(
    "CONFIG" = "Configural",
    "METRIC" = "Metric (equal loadings)",
    "SCALAR" = "Scalar (equal intercepts)",
    "STRICT" = "Strict (equal residuals)"
  )

  dt <- copy(minv_dt)
  dt[, Model := fifelse(MODEL %in% names(model_labels),
                        model_labels[MODEL], MODEL)]

  dt[, .(
    Model,
    CFI = round(CFI, 3),
    RMSEA = round(RMSEA, 3),
    Delta_CFI = round(DELTA_CFI, 3),
    Delta_RMSEA = round(DELTA_RMSEA, 3)
  )] |>
    gt() |>
    cols_label(Delta_CFI = "ΔCFI", Delta_RMSEA = "ΔRMSEA") |>
    tab_header(title = "Cognitive CFA Measurement Invariance Across Sex") |>
    sub_missing(missing_text = "—") |>
    gt_pdf_style()
}

#' Format cognitive reliability table for display
#' @param rel_dt data.table from load_cog_reliability()
#' @return gt table object
#' @export
format_cog_reliability_gt <- function(rel_dt) {
  if (is.null(rel_dt) || nrow(rel_dt) == 0) {
    return(gt(data.table(Message = "Reliability data not available")))
  }

  # Factor labels
  factor_labels <- c(
    "COG" = "General cognition (g)",
    "MEM" = "Memory-specific (MEM_s)",
    "PRSP" = "Processing speed-specific (PRSP_s)"
  )

  dt <- copy(rel_dt)
  dt[, Factor := factor_labels[FACTOR]]

  dt[, .(
    Factor,
    N,
    `r [95% CI]` = sprintf("%.3f [%.3f, %.3f]", R, R_CI_LOWER, R_CI_UPPER),
    `ICC [95% CI]` = sprintf("%.3f [%.3f, %.3f]", ICC, ICC_CI_LOWER, ICC_CI_UPPER),
    Interpretation = INTERPRETATION
  )] |>
    gt() |>
    tab_header(title = "Cognitive Factor Test-Retest Reliability") |>
    gt_pdf_style()
}

#' Format HVR vs HC SEM comparison table
#' @param comparison_dt data.table from get_hvr_hc_sem_comparison()
#' @return gt table object
#' @export
format_hvr_hc_sem_gt <- function(comparison_dt) {
  if (is.null(comparison_dt) || nrow(comparison_dt) == 0) {
    return(gt(data.table(Message = "HVR vs HC comparison data not available")))
  }

  # Parameter labels
  param_labels <- c(
    "diff_g" = "D(HVR - HC) -> g",
    "diff_mem" = "D(HVR - HC) -> MEM_s",
    "diff_prs" = "D(HVR - HC) -> PRSP_s",
    "path_hc_g" = "HC -> g",
    "path_hvr_g" = "HVR -> g",
    "path_hc_mem" = "HC -> MEM_s",
    "path_hvr_mem" = "HVR -> MEM_s",
    "path_hc_prs" = "HC -> PRSP_s",
    "path_hvr_prs" = "HVR -> PRSP_s"
  )

  dt <- copy(comparison_dt)
  dt[, Path := fifelse(Parameter %in% names(param_labels),
                       param_labels[Parameter], Parameter)]

  # Separate differences from individual paths
  diffs <- dt[Parameter %like% "^diff_"]

  # Format differences table
  diffs[, .(
    Path,
    Beta_std = round(Estimate_Std, 3),
    CI_95 = sprintf("[%.3f, %.3f]", CI_Lower, CI_Upper),
    p = format_p(P_Value)
  )] |>
    gt() |>
    cols_label(Beta_std = "β (std)", CI_95 = "95% CI") |>
    tab_header(title = "HVR vs HC: Path Coefficient Comparison") |>
    gt_pdf_style()
}

#' Format hemisphere comparison table
#' @param hemi_dt data.table from get_hemisphere_comparison()
#' @return gt table object
#' @export
format_hemisphere_gt <- function(hemi_dt) {
  if (is.null(hemi_dt) || nrow(hemi_dt) == 0) {
    return(gt(data.table(Message = "Hemisphere data not available")))
  }

  # Make a copy and apply ROI ordering
  dt <- copy(hemi_dt)
  if ("ROI" %in% names(dt)) {
    dt <- apply_roi_order(dt, "ROI")
    # Create ROI_LABEL if not present, using correct order
    if (!"ROI_LABEL" %in% names(dt)) {
      dt[, ROI_LABEL := fifelse(ROI %in% names(ROI_LABELS), ROI_LABELS[ROI], as.character(ROI))]
    }
  }

  # Ensure ROI_LABEL follows ROI order
  if ("ROI" %in% names(dt)) {
    dt[, ROI_LABEL := factor(ROI_LABEL, levels = unique(ROI_LABEL))]
  }

  # Reshape to wide format
  dt_wide <- dcast(dt, ROI_LABEL ~ SIDE, value.var = c("ESTIMATE", "CI_LOWER", "CI_UPPER"))

  dt_wide[, .(
    Region = ROI_LABEL,
    `Left d` = round(ESTIMATE_L, 2),
    `Left 95% CI` = sprintf("[%.2f, %.2f]", CI_LOWER_L, CI_UPPER_L),
    `Right d` = round(ESTIMATE_R, 2),
    `Right 95% CI` = sprintf("[%.2f, %.2f]", CI_LOWER_R, CI_UPPER_R),
    `L-R Diff` = round(ESTIMATE_L - ESTIMATE_R, 2)
  )] |>
    gt() |>
    tab_header(title = "Sex Differences by Hemisphere") |>
    gt_pdf_style()
}

#' Format age-stratified effect sizes table
#' @param age_dt data.table from get_age_stratified()
#' @return gt table object
#' @export
format_age_stratified_gt <- function(age_dt) {
  if (is.null(age_dt) || nrow(age_dt) == 0) {
    return(gt(data.table(Message = "Age-stratified data not available")))
  }

  age_dt[, .(
    Region = ROI_LABEL,
    `Age Group` = AGE_BIN,
    d = round(ESTIMATE, 2),
    `95% CI` = sprintf("[%.2f, %.2f]", CI_LOWER, CI_UPPER),
    N = N_FEMALE + N_MALE
  )] |>
    gt() |>
    tab_header(title = "Sex Differences by Age Group (Residualized)") |>
    gt_pdf_style()
}

#' Format sensitivity comparison table
#' @param sens_dt data.table from get_sensitivity_comparison()
#' @return gt table object
#' @export
format_sensitivity_gt <- function(sens_dt) {
  if (is.null(sens_dt) || nrow(sens_dt) == 0) {
    return(gt(data.table(Message = "Sensitivity data not available")))
  }

  # Make a copy to avoid modifying original
  dt <- copy(sens_dt)

  # Apply ROI ordering
  if ("ROI" %in% names(dt)) {
    dt <- apply_roi_order(dt, "ROI")
  }
  # Apply ADJ ordering
  if ("ADJ" %in% names(dt)) {
    dt <- apply_adj_order(dt, "ADJ")
  }

  # Handle column name variations (ESTIMATE_ALL vs ESTIMATE_PRIMARY)
  prim_col <- if ("ESTIMATE_ALL" %in% names(dt)) "ESTIMATE_ALL" else "ESTIMATE_PRIMARY"
  diff_col <- if ("DIFF" %in% names(dt)) "DIFF" else {
    # Calculate difference if not present
    dt[, DIFF_CALC := get(prim_col) - ESTIMATE_SENS]
    "DIFF_CALC"
  }

  # Check for required label columns and create display labels
  if ("ROI" %in% names(dt) && !"ROI_LABEL" %in% names(dt)) {
    dt[, ROI_LABEL := fifelse(ROI %in% names(ROI_LABELS), ROI_LABELS[ROI], as.character(ROI))]
  }
  if ("ADJ" %in% names(dt) && !"ADJ_LABEL" %in% names(dt)) {
    dt[, ADJ_LABEL := fifelse(ADJ %in% names(ADJ_LABELS), ADJ_LABELS[ADJ], as.character(ADJ))]
  }

  roi_col <- if ("ROI_LABEL" %in% names(dt)) "ROI_LABEL" else "ROI"
  adj_col <- if ("ADJ_LABEL" %in% names(dt)) "ADJ_LABEL" else "ADJ"

  dt[, .(
    Region = get(roi_col),
    Method = get(adj_col),
    `Primary d` = round(get(prim_col), 2),
    `Sensitivity d` = round(ESTIMATE_SENS, 2),
    `Difference` = round(get(diff_col), 3)
  )] |>
    gt() |>
    tab_header(title = "Primary vs Sensitivity Sample") |>
    gt_pdf_style()
}

#' Format HVR-ICV validation table for display
#' @param hvr_icv_dt data.table from get_hvr_icv_validation()
#' @return gt table object
#' @export
format_hvr_icv_validation_gt <- function(hvr_icv_dt) {
  if (is.null(hvr_icv_dt) || nrow(hvr_icv_dt) == 0) {
    return(gt(data.table(Message = "HVR-ICV validation data not available. Run R/09_sex_differences.R")))
  }

  # Variable labels
  var_labels <- c(
    "HVR" = "Hippocampal-Ventricle Ratio",
    "HC" = "Hippocampal Volume (raw)",
    "LV" = "Lateral Ventricle Volume (raw)"
  )

  dt <- copy(hvr_icv_dt)
  dt[, Variable := var_labels[VARIABLE]]

  dt[, .(
    Variable,
    r_ICV = round(CORRELATION, 3),
    CI_95 = sprintf("[%.3f, %.3f]", CI_LOWER, CI_UPPER),
    R2 = sprintf("%.1f%%", VARIANCE_EXPLAINED_PCT),
    N = format(N, big.mark = ",")
  )] |>
    gt() |>
    cols_label(r_ICV = "r(X, ICV)", CI_95 = "95% CI", R2 = "R²") |>
    tab_header(
      title = "HVR Self-Normalizing Verification",
      subtitle = "Correlations with intracranial volume (ICV)"
    ) |>
    gt_pdf_style()
}

#' Format cognitive paths table
#' @param paths_dt data.table from get_cognitive_paths()
#' @return gt table object
#' @export
format_cognitive_paths_gt <- function(paths_dt) {
  if (is.null(paths_dt) || nrow(paths_dt) == 0) {
    return(gt(data.table(Message = "Cognitive path data not available")))
  }

  # Outcome labels
  outcome_labels <- c(
    "g" = "General Cognition (g)",
    "MEM_s" = "Memory-Specific",
    "PRSP_s" = "Processing Speed-Specific"
  )

  paths_dt[, .(
    Outcome = fifelse(Outcome %in% names(outcome_labels),
                      outcome_labels[Outcome], Outcome),
    Sex,
    Beta = round(beta, 3),
    CI_95 = sprintf("[%.3f, %.3f]", ci_lower, ci_upper),
    p = fifelse(pvalue < 0.001, sprintf("%.2e", pvalue), sprintf("%.3f", pvalue))
  )] |>
    gt() |>
    cols_label(Beta = "β", CI_95 = "95% CI") |>
    tab_header(title = "Hippocampal Volume → Cognitive Domain Paths") |>
    gt_pdf_style()
}

# =============================================================================
# SEM Full Parameter Tables
# =============================================================================

#' Format SEM factor loadings table
#' @param loadings_dt data.table from extract_sem_factor_loadings()
#' @param model_name Model name for title
#' @return gt table object
#' @export
format_sem_loadings_gt <- function(loadings_dt, model_name = "SEM") {
  if (is.null(loadings_dt) || nrow(loadings_dt) == 0) {
    return(gt(data.table(Message = "Factor loadings not available")))
  }

  # Factor labels (publication-ready) - use markdown for cross-format compatibility
  factor_labels <- c(
    "g" = "General Cognition (<em>g</em>)",
    "MEM_s" = "Memory-Specific",
    "PRSP_s" = "Processing Speed-Specific"
  )

  # Cognitive test indicator labels (publication-ready)
  indicator_labels <- c(
    "PRS_mean_time" = "Pairs Matching: Time",
    "PRS_mean_inc" = "Pairs Matching: Errors",
    "NUM" = "Numeric Memory",
    "PRMEM_res_n" = "Prospective Memory",
    "FLINT" = "Fluid Intelligence",
    "MATS_corr" = "Matrix Pattern Completion",
    "TRLS_alnum_time" = "Trail Making: Alphanumeric",
    "TOWER_corr" = "Tower Rearranging",
    "REACT" = "Reaction Time",
    "TRLS_num_time" = "Trail Making: Numeric",
    "SYM_corr" = "Symbol-Digit: Correct",
    "SYM_try" = "Symbol-Digit: Attempted"
  )

  # Create display table with proper labels
  dt <- copy(loadings_dt)
  dt[, `:=`(
    Factor_Label = fifelse(FACTOR %in% names(factor_labels), factor_labels[FACTOR], FACTOR),
    Indicator_Label = fifelse(INDICATOR %in% names(indicator_labels),
                               indicator_labels[INDICATOR], INDICATOR),
    Sex = fifelse(is.na(GROUP) | GROUP == 0, "Pooled",
                  fifelse(GROUP == 1, "Female", "Male")),
    # Separate loading and SE columns
    Loading_fmt = fifelse(is.na(LOADING), "—", sprintf("%.3f", LOADING)),
    SE_fmt = fifelse(is.na(SE), "—", sprintf("%.3f", SE)),
    # P-value formatting
    p_fmt = fifelse(is.na(PVALUE), "—",
                    fifelse(PVALUE == 0, "< 2.2e-16",
                    fifelse(PVALUE < 0.001, sprintf("%.2e", PVALUE), sprintf("%.3f", PVALUE)))),
    # Flag significant
    is_sig = !is.na(PVALUE) & PVALUE < 0.05
  )]

  # Reshape to wide format (separate columns for Female/Male)
  dt_wide <- dcast(
    dt,
    Factor_Label + Indicator_Label ~ Sex,
    value.var = c("Loading_fmt", "SE_fmt", "p_fmt", "is_sig"),
    fill = list(Loading_fmt = "—", SE_fmt = "—", p_fmt = "—", is_sig = FALSE)
  )

  # Reorder columns and select for display
  display_cols <- c("Factor_Label", "Indicator_Label")
  if ("Loading_fmt_Female" %in% names(dt_wide)) {
    display_cols <- c(display_cols, "Loading_fmt_Female", "SE_fmt_Female", "p_fmt_Female")
  }
  if ("Loading_fmt_Male" %in% names(dt_wide)) {
    display_cols <- c(display_cols, "Loading_fmt_Male", "SE_fmt_Male", "p_fmt_Male")
  }

  dt_display <- dt_wide[, ..display_cols]

  # Build gt table
  gt_tbl <- dt_display |>
    gt(groupname_col = "Factor_Label") |>
    cols_label(
      Indicator_Label = "Indicator"
    )

  # Add column labels and spanners for sex groups
  if ("Loading_fmt_Female" %in% names(dt_display)) {
    gt_tbl <- gt_tbl |>
      cols_label(Loading_fmt_Female = "λ", SE_fmt_Female = "SE", p_fmt_Female = md("*p*")) |>
      tab_spanner(label = "Female", columns = c(Loading_fmt_Female, SE_fmt_Female, p_fmt_Female))
  }
  if ("Loading_fmt_Male" %in% names(dt_display)) {
    gt_tbl <- gt_tbl |>
      cols_label(Loading_fmt_Male = "λ", SE_fmt_Male = "SE", p_fmt_Male = md("*p*")) |>
      tab_spanner(label = "Male", columns = c(Loading_fmt_Male, SE_fmt_Male, p_fmt_Male))
  }

  # Bold significant rows (p < 0.05)
  if ("is_sig_Female" %in% names(dt_wide) && "Loading_fmt_Female" %in% names(dt_display)) {
    sig_rows_f <- which(dt_wide$is_sig_Female)
    if (length(sig_rows_f) > 0) {
      gt_tbl <- gt_tbl |>
        tab_style(style = cell_text(weight = "bold"),
                  locations = cells_body(columns = c(Loading_fmt_Female, SE_fmt_Female, p_fmt_Female),
                                         rows = sig_rows_f))
    }
  }
  if ("is_sig_Male" %in% names(dt_wide) && "Loading_fmt_Male" %in% names(dt_display)) {
    sig_rows_m <- which(dt_wide$is_sig_Male)
    if (length(sig_rows_m) > 0) {
      gt_tbl <- gt_tbl |>
        tab_style(style = cell_text(weight = "bold"),
                  locations = cells_body(columns = c(Loading_fmt_Male, SE_fmt_Male, p_fmt_Male),
                                         rows = sig_rows_m))
    }
  }

  gt_tbl |>
    tab_header(title = sprintf("Factor Loadings: %s", model_name)) |>
    fmt_markdown(columns = everything()) |>
    sub_missing(missing_text = "—") |>
    gt_pdf_style()
}

#' Format SEM residual variances table
#' @param residuals_dt data.table from extract_sem_residual_variances()
#' @param model_name Model name for title
#' @return gt table object
#' @export
format_sem_residuals_gt <- function(residuals_dt, model_name = "SEM") {
  if (is.null(residuals_dt) || nrow(residuals_dt) == 0) {
    return(gt(data.table(Message = "Residual variances not available")))
  }

  # Create display table
  display_dt <- residuals_dt[, .(
    Variable = VARIABLE,
    Sex = fifelse(is.na(GROUP) | GROUP == 0, "Pooled",
                  fifelse(GROUP == 1, "Female", "Male")),
    Variance = sprintf("%.3f", VARIANCE),
    SE = sprintf("%.3f", SE),
    CI_95 = sprintf("[%.3f, %.3f]", CI_LOWER, CI_UPPER)
  )]

  display_dt |>
    gt() |>
    cols_label(
      Variable = "Variable",
      Sex = "Sex",
      Variance = "Residual Var.",
      SE = "SE",
      CI_95 = "95% CI"
    ) |>
    tab_header(
      title = sprintf("Residual Variances: %s", model_name),
      subtitle = "Standardized residual variances (1 - R²)"
    ) |>
    tab_source_note("Lower values indicate more variance explained by latent factors") |>
    gt_pdf_style()
}

#' Format SEM covariate effects table
#' @param covar_dt data.table from extract_sem_covariate_effects()
#' @param model_name Model name for title
#' @param exclude_outcomes Outcomes to exclude (default: "MEM_s" due to factor collapse)
#' @return gt table object
#' @export
format_sem_covariates_gt <- function(covar_dt, model_name = "SEM",
                                      exclude_outcomes = "MEM_s") {
  if (is.null(covar_dt) || nrow(covar_dt) == 0) {
    return(gt(data.table(Message = "Covariate effects not available")))
  }

  # Filter out excluded outcomes (e.g., MEM_s due to factor collapse)
  if (length(exclude_outcomes) > 0) {
    covar_dt <- covar_dt[!OUTCOME %in% exclude_outcomes]
  }

  # Outcome labels (publication-ready)
  outcome_labels <- c(
    "g" = "General Cognition (<em>g</em>)",
    "MEM_s" = "Memory-Specific",
    "PRSP_s" = "Processing Speed-Specific",
    "HVR" = "Hippocampal-Ventricle Ratio",
    "HC" = "Hippocampal Volume",
    "HC_RES" = "Hippocampus (Residualized)",
    "LV" = "Lateral Ventricular Volume"
  )

  # Covariate labels (publication-ready)
  covariate_labels <- c(
    "AGE" = "Age",
    "AGE_sq" = "Age²",
    "ICV" = "Intracranial Volume",
    "IMDP" = "Deprivation Index",
    "EDUC" = "Education",
    "SEX" = "Sex",
    "SITE" = "Assessment Site"
  )

  # Create display table with proper labels
  dt <- copy(covar_dt)
  dt[, `:=`(
    Outcome_Label = fifelse(OUTCOME %in% names(outcome_labels),
                            outcome_labels[OUTCOME], OUTCOME),
    Covariate_Label = fifelse(COVARIATE %in% names(covariate_labels),
                               covariate_labels[COVARIATE], COVARIATE),
    Sex = fifelse(is.na(GROUP) | GROUP == 0, "Pooled",
                  fifelse(GROUP == 1, "Female", "Male")),
    # Separate beta and SE columns
    Beta_fmt = fifelse(is.na(BETA), "—", sprintf("%.3f", BETA)),
    SE_fmt = fifelse(is.na(SE), "—", sprintf("%.3f", SE)),
    # P-value formatting: scientific notation for p < 0.001
    p_fmt = fifelse(is.na(PVALUE), "—",
                    fifelse(PVALUE == 0, "< 2.2e-16",
                    fifelse(PVALUE < 0.001, sprintf("%.2e", PVALUE), sprintf("%.3f", PVALUE)))),
    # Flag significant for bolding
    is_sig = !is.na(PVALUE) & PVALUE < 0.05
  )]

  # Set outcome factor order: Brain measures first (HVR, HC-Res, HC, LV), then cognitive
  outcome_order <- c("Hippocampal-Ventricle Ratio", "Hippocampus (Residualized)",
                     "Hippocampal Volume", "Lateral Ventricular Volume",
                     "General Cognition (<em>g</em>)", "Processing Speed-Specific", "Memory-Specific")
  dt[, Outcome_Label := factor(Outcome_Label, levels = outcome_order)]

  # Reshape to wide format (separate columns for Female/Male)
  dt_wide <- dcast(
    dt,
    Outcome_Label + Covariate_Label ~ Sex,
    value.var = c("Beta_fmt", "SE_fmt", "p_fmt", "is_sig"),
    fill = list(Beta_fmt = "—", SE_fmt = "—", p_fmt = "—", is_sig = FALSE)
  )

  # Reorder columns and select for display
  display_cols <- c("Outcome_Label", "Covariate_Label")
  if ("Beta_fmt_Female" %in% names(dt_wide)) {
    display_cols <- c(display_cols, "Beta_fmt_Female", "SE_fmt_Female", "p_fmt_Female")
  }
  if ("Beta_fmt_Male" %in% names(dt_wide)) {
    display_cols <- c(display_cols, "Beta_fmt_Male", "SE_fmt_Male", "p_fmt_Male")
  }

  dt_display <- dt_wide[, ..display_cols]

  # Build gt table
  gt_tbl <- dt_display |>
    gt(groupname_col = "Outcome_Label") |>
    cols_label(
      Covariate_Label = "Covariate"
    )

  # Add column labels and spanners for sex groups
  if ("Beta_fmt_Female" %in% names(dt_display)) {
    gt_tbl <- gt_tbl |>
      cols_label(Beta_fmt_Female = "β", SE_fmt_Female = "SE", p_fmt_Female = md("*p*")) |>
      tab_spanner(label = "Female", columns = c(Beta_fmt_Female, SE_fmt_Female, p_fmt_Female))
  }
  if ("Beta_fmt_Male" %in% names(dt_display)) {
    gt_tbl <- gt_tbl |>
      cols_label(Beta_fmt_Male = "β", SE_fmt_Male = "SE", p_fmt_Male = md("*p*")) |>
      tab_spanner(label = "Male", columns = c(Beta_fmt_Male, SE_fmt_Male, p_fmt_Male))
  }

  # Bold significant coefficients
  if ("is_sig_Female" %in% names(dt_wide)) {
    sig_rows_f <- which(dt_wide$is_sig_Female)
    if (length(sig_rows_f) > 0) {
      gt_tbl <- gt_tbl |>
        tab_style(
          style = cell_text(weight = "bold"),
          locations = cells_body(columns = c(Beta_fmt_Female, SE_fmt_Female, p_fmt_Female), rows = sig_rows_f)
        )
    }
  }
  if ("is_sig_Male" %in% names(dt_wide)) {
    sig_rows_m <- which(dt_wide$is_sig_Male)
    if (length(sig_rows_m) > 0) {
      gt_tbl <- gt_tbl |>
        tab_style(
          style = cell_text(weight = "bold"),
          locations = cells_body(columns = c(Beta_fmt_Male, SE_fmt_Male, p_fmt_Male), rows = sig_rows_m)
        )
    }
  }

  gt_tbl |>
    tab_header(title = sprintf("Covariate Effects: %s", model_name)) |>
    fmt_markdown(columns = everything()) |>
    sub_missing(missing_text = "—") |>
    gt_pdf_style()
}

# =============================================================================
# GAMLSS Model Coefficient Tables
# =============================================================================

#' Format GAMLSS model coefficients table
#' @param coef_dt data.table from get_gamlss_model_summary()
#' @param title Optional title
#' @param include_site Whether to include site terms (default FALSE for transfer models)
#' @return gt table object
#' @export
format_gamlss_coef_gt <- function(coef_dt, title = "GAMLSS Model Coefficients",
                                   include_site = FALSE) {
  if (is.null(coef_dt) || nrow(coef_dt) == 0) {
    return(gt(data.table(Message = "GAMLSS coefficients not available")))
  }

  # Parameter descriptions (row groups)
  param_desc <- c(
    "mu" = "Location (μ)",
    "sigma" = "Scale (σ)",
    "nu" = "Shape (ν)",
    "tau" = "Shape (τ)"
  )

  # Term labels (publication-ready)
  term_labels <- c(
    "(Intercept)" = "Intercept",
    "Intercept" = "Intercept",
    "AGE" = "Age",
    "I(AGE^2)" = "Age²",
    "AGE_sq" = "Age²",
    "ICV" = "Intracranial Volume",
    "ICC" = "Intracranial Volume",
    "EDUC_num" = "Education (years)",
    "SITE" = "Assessment Site",
    "random(SITE)" = "Site (random)",
    "pb(AGE)" = "Age (smoothed)",
    "cs(AGE)" = "Age (cubic spline)"
  )

  dt <- copy(coef_dt)

  # Remove site rows if not including site (for transfer models)
  if (!include_site) {
    dt <- dt[!TERM %like% "SITE"]
  }

  # If no data remains after filtering
  if (nrow(dt) == 0) {
    return(gt(data.table(Message = "No coefficients to display")))
  }

  # Check if Distribution is uniform - if so, we'll mention it in subtitle instead
  unique_families <- unique(dt$FAMILY)
  single_family <- length(unique_families) == 1

  # Create labels
  dt[, `:=`(
    Parameter_Label = fifelse(PARAMETER %in% names(param_desc),
                               param_desc[PARAMETER], PARAMETER),
    Term_Label = fifelse(TERM %in% names(term_labels),
                          term_labels[TERM], TERM)
  )]

  # Separate Estimate and SE columns
  dt[, Estimate_fmt := fifelse(
    is.na(ESTIMATE) | is.nan(ESTIMATE) | is.infinite(ESTIMATE),
    "—",
    sprintf("%.4f", ESTIMATE)
  )]
  dt[, SE_fmt := fifelse(
    is.na(SE) | is.nan(SE) | is.infinite(SE),
    "—",
    sprintf("%.4f", SE)
  )]

  # Format p-values - handle NA and use scientific notation for small values
  dt[, p_fmt := fifelse(
    is.na(PVALUE) | is.nan(PVALUE) | is.infinite(PVALUE),
    "—",
    fifelse(PVALUE == 0, "< 2.2e-16",
    fifelse(PVALUE < 0.001, sprintf("%.2e", PVALUE), sprintf("%.3f", PVALUE)))
  )]
  dt[, is_sig := !is.na(PVALUE) & PVALUE < 0.05]

  # Reshape to wide format by sex with spanners
  dt_wide <- dcast(
    dt,
    Parameter_Label + Term_Label ~ SEX,
    value.var = c("Estimate_fmt", "SE_fmt", "p_fmt", "is_sig"),
    fill = list(Estimate_fmt = "—", SE_fmt = "—", p_fmt = "—", is_sig = FALSE)
  )

  # Rename columns for display
  setnames(dt_wide, c("Parameter_Label", "Term_Label"), c("Parameter", "Term"))

  # Build gt table with parameter as row group
  display_cols <- c("Parameter", "Term")
  if ("Estimate_fmt_Female" %in% names(dt_wide)) {
    display_cols <- c(display_cols, "Estimate_fmt_Female", "SE_fmt_Female", "p_fmt_Female")
  }
  if ("Estimate_fmt_Male" %in% names(dt_wide)) {
    display_cols <- c(display_cols, "Estimate_fmt_Male", "SE_fmt_Male", "p_fmt_Male")
  }

  gt_tbl <- dt_wide[, ..display_cols] |>
    gt(groupname_col = "Parameter") |>
    cols_label(Term = "Term")

  # Add spanners for sex (with unique IDs to avoid conflict with column names)
  if ("Estimate_fmt_Female" %in% names(dt_wide)) {
    gt_tbl <- gt_tbl |>
      cols_label(Estimate_fmt_Female = "Estimate", SE_fmt_Female = "SE", p_fmt_Female = md("*p*")) |>
      tab_spanner(label = "Female", columns = c(Estimate_fmt_Female, SE_fmt_Female, p_fmt_Female), id = "spanner_female")
  }
  if ("Estimate_fmt_Male" %in% names(dt_wide)) {
    gt_tbl <- gt_tbl |>
      cols_label(Estimate_fmt_Male = "Estimate", SE_fmt_Male = "SE", p_fmt_Male = md("*p*")) |>
      tab_spanner(label = "Male", columns = c(Estimate_fmt_Male, SE_fmt_Male, p_fmt_Male), id = "spanner_male")
  }

  # Bold significant rows (p < 0.05)
  if ("is_sig_Female" %in% names(dt_wide) && "Estimate_fmt_Female" %in% names(dt_wide)) {
    sig_rows_f <- which(dt_wide$is_sig_Female)
    if (length(sig_rows_f) > 0) {
      gt_tbl <- gt_tbl |>
        tab_style(style = cell_text(weight = "bold"),
                  locations = cells_body(columns = c(Estimate_fmt_Female, SE_fmt_Female, p_fmt_Female),
                                         rows = sig_rows_f))
    }
  }
  if ("is_sig_Male" %in% names(dt_wide) && "Estimate_fmt_Male" %in% names(dt_wide)) {
    sig_rows_m <- which(dt_wide$is_sig_Male)
    if (length(sig_rows_m) > 0) {
      gt_tbl <- gt_tbl |>
        tab_style(style = cell_text(weight = "bold"),
                  locations = cells_body(columns = c(Estimate_fmt_Male, SE_fmt_Male, p_fmt_Male),
                                         rows = sig_rows_m))
    }
  }

  # Build subtitle
  subtitle <- "GAMLSS distributional regression parameters"
  if (single_family) {
    subtitle <- sprintf("%s (%s distribution)", subtitle, unique_families)
  }

  gt_tbl |>
    tab_header(
      title = title,
      subtitle = subtitle
    ) |>
    sub_missing(missing_text = "—") |>
    gt_pdf_style()
}

#' Format GAMLSS model summary by ROI
#' @param roi ROI name (e.g., "HVR", "HC")
#' @param adj Adjustment method
#' @return gt table object
#' @export
format_gamlss_roi_summary_gt <- function(roi, adj = "NON") {
  coef_dt <- get_gamlss_model_summary(roi = roi, adj = adj)

  if (is.null(coef_dt)) {
    return(gt(data.table(Message = sprintf("GAMLSS coefficients for %s not available", roi))))
  }

  # ROI labels
  roi_labels <- c("HVR" = "Hippocampal-to-Ventricle Ratio",
                  "HC" = "Hippocampal Volume",
                  "LV" = "Lateral Ventricular Volume")

  title <- sprintf("GAMLSS Model: %s",
                   fifelse(roi %in% names(roi_labels), roi_labels[roi], roi))

  format_gamlss_coef_gt(coef_dt, title = title) |>
    gt_pdf_style()
}

# =============================================================================
# SEM Path Coefficients by Sex
# =============================================================================

#' Format SEM brain-cognition paths by sex as merged gt table
#'
#' Creates a table showing paths to both *g* and Speed_s with row groups by brain measure
#' and Sex spanners (Female/Male) with β, CI, p columns.
#'
#' @param sem_params SEM parameters list from load_analysis_data()
#' @param title Optional table title
#' @return gt table object
#' @export
format_sem_paths_by_sex_gt <- function(sem_params, title = "Brain-Cognition Path Coefficients by Sex") {
  if (is.null(sem_params)) {
    return(gt(data.table(Message = "SEM parameters not available")))
  }

  # Helper to extract sex-stratified paths (group 1 = Female, group 2 = Male)
  get_sex_path <- function(model, lhs_var, rhs_var, group_num) {
    params <- sem_params[[model]]
    if (is.null(params)) return(NULL)
    if (!is.data.table(params)) setDT(params)
    row <- params[lhs == lhs_var & rhs == rhs_var & op == "~" & group == group_num]
    if (nrow(row) == 0) return(NULL)
    row[1]
  }

  # Helper to safely extract values
  safe_val <- function(path, field) {
    if (is.null(path)) return(NA_real_)
    path[[field]]
  }

  # Helper to format CI
  fmt_ci <- function(lower, upper) {
    if (is.na(lower) || is.na(upper)) return("—")
    sprintf("[%.3f, %.3f]", lower, upper)
  }

  # Define brain measures in correct order
  brain_measures <- list(
    list(model = "HVR_COG", rhs = "HVR", label = "Hippocampal-to-Ventricle Ratio"),
    list(model = "HC_RES_COG", rhs = "HC_RES", label = "Hippocampus (Residualized)"),
    list(model = "HC_COG", rhs = "HC", label = "Hippocampus (Unadjusted)")
  )

  # Define cognitive outcomes
  cog_outcomes <- list(
    list(lhs = "g", path_label = "→ <em>g</em>"),
    list(lhs = "PRSP_s", path_label = "→ Speed<sub>s</sub>")
  )

  # Build data
  rows <- list()
  for (brain in brain_measures) {
    for (cog in cog_outcomes) {
      # Get Female and Male paths
      f_path <- get_sex_path(brain$model, cog$lhs, brain$rhs, 1)
      m_path <- get_sex_path(brain$model, cog$lhs, brain$rhs, 2)

      rows[[length(rows) + 1]] <- data.table(
        Brain_Measure = brain$label,
        Path = cog$path_label,
        F_beta = safe_val(f_path, "est.std"),
        F_ci = fmt_ci(safe_val(f_path, "ci.lower"), safe_val(f_path, "ci.upper")),
        F_p = safe_val(f_path, "pvalue"),
        M_beta = safe_val(m_path, "est.std"),
        M_ci = fmt_ci(safe_val(m_path, "ci.lower"), safe_val(m_path, "ci.upper")),
        M_p = safe_val(m_path, "pvalue")
      )
    }
  }

  dt <- rbindlist(rows)

  # Set brain measure as ordered factor
  dt[, Brain_Measure := factor(Brain_Measure, levels = c(
    "Hippocampal-to-Ventricle Ratio",
    "Hippocampus (Residualized)",
    "Hippocampus (Unadjusted)"
  ))]

  # Sort by factor order BEFORE computing row indices
  setorder(dt, Brain_Measure)
  dt[, row_id := .I]

  # Format p-values: scientific notation only for p < 0.001
  format_p_smart <- function(p) {
    if (is.na(p)) return("—")
    if (p == 0) return("< 2.2e-16")
    if (p < 0.001) return(sprintf("%.2e", p))
    sprintf("%.3f", p)
  }
  dt[, F_p_fmt := sapply(F_p, format_p_smart)]
  dt[, M_p_fmt := sapply(M_p, format_p_smart)]

  # Identify significant rows for bold styling (after ordering)
  female_sig_rows <- dt[F_p <= 0.05, row_id]
  male_sig_rows <- dt[M_p <= 0.05, row_id]

  # Build gt table
  gt_tbl <- dt[, .(Brain_Measure, Path, F_beta, F_ci, F_p_fmt, M_beta, M_ci, M_p_fmt)] |>
    gt(groupname_col = "Brain_Measure") |>
    cols_label(
      Path = "Outcome",
      F_beta = "β",
      F_ci = "95% CI",
      F_p_fmt = md("*p*"),
      M_beta = "β",
      M_ci = "95% CI",
      M_p_fmt = md("*p*")
    ) |>
    fmt_number(columns = c(F_beta, M_beta), decimals = 3) |>
    tab_spanner(label = "Female", columns = c(F_beta, F_ci, F_p_fmt), id = "spanner_female") |>
    tab_spanner(label = "Male", columns = c(M_beta, M_ci, M_p_fmt), id = "spanner_male") |>
    fmt_markdown(columns = c(Path)) |>
    tab_header(
      title = title,
      subtitle = "Multi-group SEM with sex-stratified parameter estimation"
    )

  # Apply bold styling to significant results
  if (length(female_sig_rows) > 0) {
    gt_tbl <- gt_tbl |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = c(F_beta, F_ci, F_p_fmt), rows = female_sig_rows)
      )
  }
  if (length(male_sig_rows) > 0) {
    gt_tbl <- gt_tbl |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = c(M_beta, M_ci, M_p_fmt), rows = male_sig_rows)
      )
  }

  gt_tbl |>
    gt_pdf_style()
}

#' Format HVR-ICV validation table
#'
#' Creates a gt table showing correlations between brain measures and ICV,
#' validating the self-normalizing property of HVR.
#'
#' @param hvr_icv_dt data.table with HVR-ICV validation data
#'   Must contain: VARIABLE, CORRELATION, CI_LOWER, CI_UPPER
#'   Optional: SAMPLE column for Primary/Sensitivity distinction
#' @param title Title for the table
#' @param subtitle Subtitle for the table
#' @return gt table object
#' @export
format_hvr_validation_gt <- function(hvr_icv_dt,
                                     title = "HVR Self-Normalizing Validation",
                                     subtitle = "Pearson correlation with intracranial volume (ICV)") {

  if (is.null(hvr_icv_dt) || nrow(hvr_icv_dt) == 0) {
    return(gt(data.table(Message = "Validation data not available")))
  }

  dt <- copy(hvr_icv_dt)

  # Filter to relevant variables
  dt <- dt[VARIABLE %in% c("HC", "HC_PRP", "HC_STX", "HC_RES", "HVR")]

  # Check if SAMPLE column exists
  has_sample_col <- "SAMPLE" %in% names(dt)
  if (!has_sample_col) {
    dt[, SAMPLE := "Primary"]
  }

  # Add readable labels
  dt[, Measure := fcase(
    VARIABLE == "HC", "Unadjusted",
    VARIABLE == "HC_PRP", "Proportions",
    VARIABLE == "HC_STX", "Stereotaxic",
    VARIABLE == "HC_RES", "Residualized",
    VARIABLE == "HVR", "Self-normalizing"
  )]

  # Add ROI grouping with full names
  dt[, ROI := fifelse(VARIABLE == "HVR", "Hippocampal-to-Ventricle Ratio", "Hippocampus")]

  # Format CI
  dt[, `95% CI` := sprintf("[%.3f, %.3f]", CI_LOWER, CI_UPPER)]

  # Order for display
  dt[, Measure := factor(Measure, levels = c(
    "Unadjusted", "Proportions", "Stereotaxic", "Residualized", "Self-normalizing"
  ))]
  dt[, SAMPLE := factor(SAMPLE, levels = c("Primary", "Sensitivity"))]
  setorder(dt, SAMPLE, Measure)

  # Build table based on whether multiple samples present
  if (dt[, uniqueN(SAMPLE)] > 1) {
    gt_tbl <- dt[, .(SAMPLE, ROI, Measure, `r(ICV)` = round(CORRELATION, 3), `95% CI`)] |>
      gt(groupname_col = "SAMPLE") |>
      tab_header(title = title, subtitle = subtitle)
  } else {
    gt_tbl <- dt[, .(ROI, Measure, `r(ICV)` = round(CORRELATION, 3), `95% CI`)] |>
      gt(groupname_col = "ROI") |>
      tab_header(title = title, subtitle = subtitle)
  }

  gt_tbl |> gt_pdf_style()
}

#' Format SEM covariate fields table
#'
#' Creates a gt table showing the UK Biobank field definitions for SEM covariates.
#'
#' @param fields_dt data.table with field definitions
#'   Must contain: Variable, Field, Description
#'   Optional: Notes column (will be converted to row footnotes)
#' @param title Title for the table
#' @return gt table object
#' @export
format_sem_covariates_fields_gt <- function(fields_dt,
                                            title = "UK Biobank Field Definitions for SEM Covariates") {

  if (is.null(fields_dt) || nrow(fields_dt) == 0) {
    return(gt(data.table(Message = "Field definitions not available")))
  }

  dt <- copy(fields_dt)

  # Extract notes before removing column
  has_notes <- "Notes" %in% names(dt)
  notes_list <- list()
  if (has_notes) {
    for (i in seq_len(nrow(dt))) {
      if (!is.na(dt$Notes[i]) && nchar(dt$Notes[i]) > 0) {
        notes_list[[dt$Variable[i]]] <- dt$Notes[i]
      }
    }
    dt[, Notes := NULL]
  }

  # Build gt table
  gt_tbl <- dt |>
    gt(rowname_col = "Variable") |>
    tab_header(title = title) |>
    tab_stubhead(label = "Variable")

  # Add row-specific footnotes for variables with notes
  for (var_name in names(notes_list)) {
    gt_tbl <- gt_tbl |>
      tab_footnote(
        footnote = notes_list[[var_name]],
        locations = cells_stub(rows = var_name)
      )
  }

  gt_tbl |> gt_pdf_style()
}
