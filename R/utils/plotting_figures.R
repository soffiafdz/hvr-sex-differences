# =============================================================================
# Manuscript Figure Functions
# =============================================================================
# Publication-ready figure generation for Quarto documents.
# These fig_* functions generate ggplot objects for specific analyses.
# Part of the R/utils/ module system.
# =============================================================================

#' Compute confidence interval for the difference between two estimates
#' Uses SE derived from individual CIs: SE = (CI_UPPER - CI_LOWER) / 3.92
#' Assumes independence between estimates for SE_diff = sqrt(SE_A^2 + SE_B^2)
#' @param est_a Numeric vector of estimates for condition A
#' @param ci_lo_a Lower CI bounds for condition A
#' @param ci_hi_a Upper CI bounds for condition A
#' @param est_b Numeric vector of estimates for condition B
#' @param ci_lo_b Lower CI bounds for condition B
#' @param ci_hi_b Upper CI bounds for condition B
#' @return data.table with DIFF, CI_LOWER_DIFF, CI_UPPER_DIFF columns
compute_difference_ci <- function(est_a, ci_lo_a, ci_hi_a, est_b, ci_lo_b, ci_hi_b) {
  se_a <- (ci_hi_a - ci_lo_a) / 3.92
  se_b <- (ci_hi_b - ci_lo_b) / 3.92
  diff <- est_a - est_b
  se_diff <- sqrt(se_a^2 + se_b^2)
  data.table(
    DIFF = diff,
    CI_LOWER_DIFF = diff - 1.96 * se_diff,
    CI_UPPER_DIFF = diff + 1.96 * se_diff
  )
}

#' Create a forest-plot panel for differences between two conditions
#' @param diff_data data.table with Y_LABEL, ADJ_LABEL, DIFF, CI_LOWER_DIFF, CI_UPPER_DIFF
#' @param y_order Factor levels for Y_LABEL ordering
#' @param adj_colors Named vector of colors for ADJ_LABEL
#' @param x_label X-axis label (default Delta*d)
#' @return ggplot object
plot_difference_panel <- function(diff_data, y_order, adj_colors,
                                  x_label = bquote(Delta * "d"),
                                  base_size = 10, show_legend = FALSE,
                                  legend_title = NULL) {
  diff_data[, Y_LABEL := factor(Y_LABEL, levels = y_order)]
  diff_data[, SIGNIFICANT := (CI_LOWER_DIFF > 0) | (CI_UPPER_DIFF < 0)]
  diff_data[, LINE_ALPHA := fifelse(SIGNIFICANT, 1.0, 0.4)]

  p <- ggplot(diff_data, aes(y = Y_LABEL, x = DIFF, color = ADJ_LABEL)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_LOWER_DIFF, xmax = CI_UPPER_DIFF, alpha = LINE_ALPHA),
                   height = 0, linewidth = 0.7) +
    scale_color_manual(values = adj_colors, name = legend_title) +
    scale_alpha_identity() +
    labs(y = NULL, x = x_label) +
    theme_publication(base_size = base_size) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  if (show_legend) {
    p <- p + theme(legend.position = "bottom")
  } else {
    p <- p + theme(legend.position = "none")
  }

  p
}

#' Plot effect sizes for HC, LV, and HVR (forest plot style)
#' Y-axis = ROI (Adjustment), shape = sample type
#' @param sex_diff Sex differences data from load_analysis_data()
#' @param include_matched Logical, whether to include matched sample comparison (default TRUE)
#' @return ggplot object
#' @export
fig_effect_sizes_by_adjustment <- function(sex_diff, include_matched = TRUE, base_size = 10) {
  if (is.null(sex_diff) || !"OVERALL" %in% names(sex_diff)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  }

  overall <- sex_diff$OVERALL
  if (is.null(overall) || nrow(overall) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No effect size data"))
  }

  # Get HC, LV, and HVR data
  hc_data <- overall[SIDE == "LR" & ROI == "HC"]
  lv_data <- overall[SIDE == "LR" & ROI == "LV"]
  hvr_data <- overall[SIDE == "LR" & ROI == "HVR" & ADJ == "NON"]

  # Full adjustment labels for Y-axis
  adj_full <- c("NON" = "Unadjusted", "PRP" = "Proportions", "STX" = "Stereotaxic", "RES" = "Residualized")

  # Y-axis labels: Full ROI name \n (Adjustment)
  hc_data[, Y_LABEL := paste0("Hippocampus\n(", adj_full[ADJ], ")")]
  hc_data[, ADJ_LABEL := adj_full[ADJ]]
  lv_data[, Y_LABEL := paste0("Lateral Ventricles\n(", adj_full[ADJ], ")")]
  lv_data[, ADJ_LABEL := adj_full[ADJ]]
  hvr_data <- copy(hvr_data)
  hvr_data[, Y_LABEL := "Hippocampal-to-Ventricle Ratio\n(Self-normalized)"]
  hvr_data[, ADJ_LABEL := "Self-normalized"]

  plot_data <- rbind(hc_data, lv_data, hvr_data, fill = TRUE)
  plot_data[, SAMPLE := "Full Sample"]

  # Add matched sample data if requested
  if (include_matched && "COMPARISON" %in% names(sex_diff)) {
    comp <- sex_diff$COMPARISON
    # HC matched
    matched_hc <- comp[ROI == "HC" & !is.na(ESTIMATE_MTCH),
                       .(ROI, ADJ, ESTIMATE = ESTIMATE_MTCH, CI_LOWER = CI_LOWER_MTCH, CI_UPPER = CI_UPPER_MTCH)]
    matched_hc[, Y_LABEL := paste0("Hippocampus\n(", adj_full[ADJ], ")")]
    matched_hc[, ADJ_LABEL := adj_full[ADJ]]
    # LV matched
    matched_lv <- comp[ROI == "LV" & !is.na(ESTIMATE_MTCH),
                       .(ROI, ADJ, ESTIMATE = ESTIMATE_MTCH, CI_LOWER = CI_LOWER_MTCH, CI_UPPER = CI_UPPER_MTCH)]
    matched_lv[, Y_LABEL := paste0("Lateral Ventricles\n(", adj_full[ADJ], ")")]
    matched_lv[, ADJ_LABEL := adj_full[ADJ]]
    # HVR matched
    matched_hvr <- comp[ROI == "HVR" & !is.na(ESTIMATE_MTCH),
                        .(ROI, ADJ, ESTIMATE = ESTIMATE_MTCH, CI_LOWER = CI_LOWER_MTCH, CI_UPPER = CI_UPPER_MTCH)]
    matched_hvr[, Y_LABEL := "Hippocampal-to-Ventricle Ratio\n(Self-normalized)"]
    matched_hvr[, ADJ_LABEL := "Self-normalized"]

    matched_data <- rbind(matched_hc, matched_lv, matched_hvr, fill = TRUE)
    matched_data[, SAMPLE := "Head-Size Matched"]
    plot_data <- rbind(plot_data, matched_data, fill = TRUE)
  }

  # Y-axis order: HVR at top, then HC, then LV (last level appears at top)
  y_order <- c("Lateral Ventricles\n(Residualized)", "Lateral Ventricles\n(Stereotaxic)",
               "Lateral Ventricles\n(Proportions)", "Lateral Ventricles\n(Unadjusted)",
               "Hippocampus\n(Residualized)", "Hippocampus\n(Stereotaxic)",
               "Hippocampus\n(Proportions)", "Hippocampus\n(Unadjusted)",
               "Hippocampal-to-Ventricle Ratio\n(Self-normalized)")
  plot_data[, Y_LABEL := factor(Y_LABEL, levels = y_order)]

  # Sample factor - Full Sample first (will appear above in dodge)
  plot_data[, SAMPLE := factor(SAMPLE, levels = c("Head-Size Matched", "Full Sample"))]

  # Adjustment colors (same for both samples)
  adj_colors_full <- c("Unadjusted" = "#E64B35", "Proportions" = "#F39B7F",
                       "Stereotaxic" = "#00A087", "Residualized" = "#3C5488", "Self-normalized" = "#9467BD")

  # Determine significance (CI doesn't cross 0)
  plot_data[, SIGNIFICANT := (CI_LOWER > 0) | (CI_UPPER < 0)]

  # Alpha based on significance (dim non-significant results)
  plot_data[, LINE_ALPHA := fifelse(SIGNIFICANT, 1.0, 0.5)]

  # Convert Y_LABEL to numeric for manual dodge
  plot_data[, Y_NUM := as.numeric(Y_LABEL)]

  # Manual vertical offset to prevent overlap (0.15 units apart)
  full_data <- copy(plot_data[SAMPLE == "Full Sample"])
  matched_data <- copy(plot_data[SAMPLE == "Head-Size Matched"])
  full_data[, Y_DODGE := Y_NUM + 0.15]
  matched_data[, Y_DODGE := Y_NUM - 0.15]

  # Assign colors (same palette for both)
  full_data[, COLOR := adj_colors_full[ADJ_LABEL]]
  matched_data[, COLOR := adj_colors_full[ADJ_LABEL]]

  # Combine for plotting with shape legend
  combined_data <- rbind(full_data, matched_data)

  p <- ggplot(combined_data, aes(y = Y_DODGE, x = ESTIMATE)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_LOWER, xmax = CI_UPPER, color = COLOR, alpha = LINE_ALPHA),
                   height = 0.2, linewidth = 0.6) +
    geom_point(aes(color = COLOR, alpha = LINE_ALPHA, shape = SAMPLE), size = 1.5) +
    scale_color_identity() +
    scale_alpha_identity() +
    scale_shape_manual(values = c("Full Sample" = 16, "Head-Size Matched" = 17),
                       name = "Sample") +
    scale_y_continuous(breaks = seq_along(y_order), labels = y_order) +
    labs(y = NULL, x = "Sex Difference (Cohen's d)\nPositive = Male < Female") +
    theme_publication(base_size = base_size) +
    theme(axis.text.y = element_text(size = base_size - 1),
          legend.position = "bottom") +
    guides(shape = guide_legend(override.aes = list(size = 2)))

  # Add difference panel if comparison data available
  if (include_matched && "COMPARISON" %in% names(sex_diff)) {
    comp <- sex_diff$COMPARISON
    diff_ci <- compute_difference_ci(
      comp$ESTIMATE_ALL, comp$CI_LOWER_ALL, comp$CI_UPPER_ALL,
      comp$ESTIMATE_MTCH, comp$CI_LOWER_MTCH, comp$CI_UPPER_MTCH
    )
    diff_data <- cbind(comp[, .(ROI, ADJ)], diff_ci)

    diff_data[ROI == "HC", Y_LABEL := paste0("Hippocampus\n(", adj_full[ADJ], ")")]
    diff_data[ROI == "HC", ADJ_LABEL := adj_full[ADJ]]
    diff_data[ROI == "LV", Y_LABEL := paste0("Lateral Ventricles\n(", adj_full[ADJ], ")")]
    diff_data[ROI == "LV", ADJ_LABEL := adj_full[ADJ]]
    diff_data[ROI == "HVR", Y_LABEL := "Hippocampal-to-Ventricle Ratio\n(Self-normalized)"]
    diff_data[ROI == "HVR", ADJ_LABEL := "Self-normalized"]

    diff_panel <- plot_difference_panel(diff_data, y_order, adj_colors_full,
                                         x_label = bquote(Delta * "d (Full - Matched)"),
                                         base_size = base_size)
    p <- p + diff_panel + plot_layout(widths = c(3, 1.5))
  }

  p
}

#' Plot effect sizes for lateral ventricles (forest/lollipop style)
#' @param sex_diff Sex differences data from load_analysis_data()
#' @return ggplot object
#' @export
fig_effect_sizes_lv <- function(sex_diff) {
  if (is.null(sex_diff) || !"OVERALL" %in% names(sex_diff)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  }
  lv_data <- sex_diff$OVERALL[SIDE == "LR" & ROI == "LV"]
  if (nrow(lv_data) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No LV data found"))

  adj_labels <- c("NON" = "Unadjusted", "PRP" = "Proportions", "STX" = "Stereotaxic", "RES" = "Residualized")
  lv_data[, ADJ_LABEL := adj_labels[ADJ]]
  # Use green palette (avoid blue/red used for sex) - consistent with HC figure
  adj_colors <- c("Unadjusted" = "#2E7D32", "Proportions" = "#66BB6A", "Stereotaxic" = "#A5D6A7", "Residualized" = "#1B5E20")
  adj_order <- c("RES", "STX", "PRP", "NON")  # Reversed for horizontal orientation
  label_order <- c("Residualized", "Stereotaxic", "Proportions", "Unadjusted")
  lv_data[, ADJ := factor(ADJ, levels = adj_order)]
  lv_data[, ADJ_LABEL := factor(ADJ_LABEL, levels = label_order)]

  # Forest/lollipop style with horizontal orientation (effect size on X, method on Y)
  ggplot(lv_data, aes(y = ADJ_LABEL, x = ESTIMATE, color = ADJ_LABEL)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_segment(aes(x = 0, xend = ESTIMATE, yend = ADJ_LABEL), linewidth = 1.2) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = CI_LOWER, xmax = CI_UPPER), height = 0.25, linewidth = 0.8) +
    scale_color_manual(values = adj_colors, guide = "none") +
    labs(y = NULL, x = "Sex Difference (Cohen's d)\nPositive = Females > Males") +
    theme_publication(base_size = 11) +
    theme(axis.text.y = element_text(size = 10, face = "bold"))
}

#' Plot age trajectory of sex differences (HC, LV, and HVR combined)
#' @param sex_diff Sex differences data from load_analysis_data()
#' @return ggplot object with ROI facets
#' @export
fig_age_trajectory <- function(sex_diff, base_size = 10) {
  if (is.null(sex_diff) || !"AGE_STRATIFIED" %in% names(sex_diff)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Age-stratified data not available"))
  }
  age_strat <- sex_diff$AGE_STRATIFIED
  if (is.null(age_strat) || nrow(age_strat) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Age-stratified data is empty"))
  }

  # Get HC, LV, and HVR data
  hc_data <- age_strat[SIDE == "LR" & ROI == "HC"]
  lv_data <- age_strat[SIDE == "LR" & ROI == "LV"]
  hvr_data <- age_strat[SIDE == "LR" & ROI == "HVR" & ADJ == "NON"]

  # HVR goes in HC facet
  hvr_data <- copy(hvr_data)
  hvr_data[, ADJ := "HVR"]
  hvr_data[, ROI_LABEL := "Hippocampus (HC) & HVR"]

  # Set ROI labels
  hc_data[, ROI_LABEL := "Hippocampus (HC) & HVR"]
  lv_data[, ROI_LABEL := "Lateral Ventricles (LV)"]

  age_data <- rbind(hc_data, lv_data, hvr_data, fill = TRUE)
  if (nrow(age_data) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No bilateral data found"))

  adj_labels <- c("NON" = "Unadjusted", "PRP" = "Proportions", "STX" = "Stereotaxic", "RES" = "Residualized", "HVR" = "HVR")
  age_data[, ADJ_LABEL := adj_labels[ADJ]]
  age_data[, AGE_MID := sapply(AGE_BIN, function(x) mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(x)), ",")))))]
  age_data[, SIGNIFICANT := (CI_LOWER > 0) | (CI_UPPER < 0)]

  adj_colors <- c("Unadjusted" = "#E64B35", "Proportions" = "#F39B7F", "Stereotaxic" = "#00A087", "Residualized" = "#3C5488", "HVR" = "#9467BD")
  age_data[, ADJ := factor(ADJ, levels = c("NON", "PRP", "STX", "RES", "HVR"))]
  age_data[, ADJ_LABEL := factor(ADJ_LABEL, levels = c("Unadjusted", "Proportions", "Stereotaxic", "Residualized", "HVR"))]
  age_data[, ROI_LABEL := factor(ROI_LABEL, levels = c("Hippocampus (HC) & HVR", "Lateral Ventricles (LV)"))]

  # Manual offset for each adjustment method to avoid overlap
  # Spread methods across ~5 years total (±2.5 years from center)
  method_offsets <- c("Unadjusted" = -2.0, "Proportions" = -1.0, "Stereotaxic" = 0, "Residualized" = 1.0, "HVR" = 2.0)
  age_data[, AGE_OFFSET := method_offsets[as.character(ADJ_LABEL)]]
  age_data[, AGE_PLOT := AGE_MID + AGE_OFFSET]

  # Linetype: HVR = solid, all others = dashed
  age_data[, LINE_TYPE := fifelse(ADJ == "HVR", "solid", "dashed")]

  # Fill color: white for non-significant (solid, no alpha)
  age_data[, FILL_COLOR := fifelse(SIGNIFICANT, as.character(ADJ_LABEL), "white")]

  # Colors including white for fill scale
  fill_colors <- c(adj_colors, "white" = "white")

  ggplot(age_data, aes(x = AGE_PLOT, y = ESTIMATE, color = ADJ_LABEL, group = ADJ_LABEL)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    # Ribbon for CI (faded) - use original AGE_MID for smooth ribbons
    geom_ribbon(aes(x = AGE_MID, ymin = CI_LOWER, ymax = CI_UPPER, fill = ADJ_LABEL), alpha = 0.08, color = NA) +
    # Lines: HVR solid, others dashed
    geom_line(aes(linetype = LINE_TYPE), linewidth = 0.8) +
    # Points with white fill for non-significant - NO alpha
    geom_point(aes(fill = FILL_COLOR), size = 2.5, shape = 21, stroke = 0.8) +
    facet_wrap(~ ROI_LABEL, scales = "fixed", ncol = 2) +
    scale_color_manual(values = adj_colors, drop = FALSE, name = "Method") +
    scale_fill_manual(values = fill_colors, drop = FALSE, guide = "none") +
    scale_linetype_identity() +
    labs(x = "Age (years)", y = "Sex Difference (Cohen's d)\nPositive = Male < Female") +
    theme_publication(base_size = base_size) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = base_size + 1)) +
    guides(color = guide_legend(title = "Method", nrow = 1))
}

#' Plot HVR centile curves
#' @param norm_tables Normative tables from load_norm_tables()
#' @param combined Logical, whether to show both sexes on same plot (default TRUE)
#' @return ggplot object
#' @export
fig_hvr_centiles <- function(norm_tables, combined = TRUE) {
  if (is.null(norm_tables)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  hvr_f <- copy(norm_tables$Female$HVR$NON$LR)
  hvr_m <- copy(norm_tables$Male$HVR$NON$LR)
  if (is.null(hvr_f) || is.null(hvr_m)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "HVR centile data not available"))

  hvr_f[, SEX := "Female"]; hvr_m[, SEX := "Male"]
  cent_data <- rbind(hvr_f, hvr_m)
  centiles <- c("p5", "p10", "p25", "p50", "p75", "p90", "p95")
  centile_labels <- c("p5" = "5th", "p10" = "10th", "p25" = "25th", "p50" = "Median", "p75" = "75th", "p90" = "90th", "p95" = "95th")
  available_centiles <- centiles[centiles %in% names(cent_data)]
  if (length(available_centiles) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No centile columns found"))
  cent_long <- melt(cent_data, id.vars = c("AGE", "SEX"), measure.vars = available_centiles, variable.name = "Centile", value.name = "Value")
  cent_long[, Centile := factor(Centile, levels = available_centiles, labels = centile_labels[available_centiles])]
  sex_colors <- c("Female" = "darkred", "Male" = "midnightblue")

  if (combined) {
    p <- ggplot(cent_long, aes(x = AGE, y = Value, color = SEX, linetype = Centile)) +
      geom_line(linewidth = 0.6, alpha = 0.8) +
      scale_color_manual(values = sex_colors) +
      scale_linetype_manual(values = c("5th" = "dotted", "10th" = "dashed", "25th" = "dotdash", "Median" = "solid", "75th" = "dotdash", "90th" = "dashed", "95th" = "dotted"), drop = FALSE) +
      labs(x = "Age (years)", y = "Hippocampal-Ventricle Ratio (HVR)", color = "Sex", linetype = "Centile") +
      theme_publication(base_size = 11) +
      theme(legend.position = "bottom", legend.box = "vertical") +
      guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2, ncol = 4))
  } else {
    p <- ggplot(cent_long, aes(x = AGE, y = Value, linetype = Centile)) +
      geom_line(aes(color = SEX), linewidth = 0.6) +
      facet_wrap(~ SEX, ncol = 2) +
      scale_color_manual(values = sex_colors, guide = "none") +
      scale_linetype_manual(values = c("5th" = "dotted", "10th" = "dashed", "25th" = "dotdash", "Median" = "solid", "75th" = "dotdash", "90th" = "dashed", "95th" = "dotted"), drop = FALSE) +
      labs(x = "Age (years)", y = "Hippocampal-Ventricle Ratio (HVR)", linetype = "Centile") +
      theme_publication(base_size = 11) +
      theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 12)) +
      guides(linetype = guide_legend(ncol = 4))
  }
  p
}

#' Combined normative centile panel: 3 ROIs x 3 columns (Female, Both, Male)
#'
#' Produces a single figure with rows for HVR, HC, and LV. The left and right
#' columns show sex-specific centile curves; the centre column overlays both
#' sexes for direct comparison.
#'
#' @param norm_tables Normative tables from load_norm_tables()
#' @return A patchwork composite ggplot
#' @export
fig_normative_centiles_panel <- function(norm_tables) {
  if (is.null(norm_tables)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  }

  require(patchwork)

  sex_colors <- c("Female" = "darkred", "Male" = "midnightblue")
  centiles <- c("p5", "p10", "p25", "p50", "p75", "p90", "p95")
  centile_labels <- c("p5" = "5th", "p10" = "10th", "p25" = "25th",
                      "p50" = "Median", "p75" = "75th", "p90" = "90th",
                      "p95" = "95th")
  lt_vals <- c("5th" = "dotted", "10th" = "dashed", "25th" = "dotdash",
               "Median" = "solid", "75th" = "dotdash", "90th" = "dashed",
               "95th" = "dotted")

  # ROI definitions: key, y-label
  roi_defs <- list(
    list(roi = "HVR", adj = "NON", ylab = "HVR"),
    list(roi = "HC",  adj = "NON", ylab = "Hippocampal Volume (mm\u00B3)"),
    list(roi = "LV",  adj = "NON", ylab = "Lateral Ventricle Volume (cc)")
  )

  # Helper: prepare long data for one ROI
  prep_data <- function(def) {
    dt_f <- copy(norm_tables$Female[[def$roi]][[def$adj]]$LR)
    dt_m <- copy(norm_tables$Male[[def$roi]][[def$adj]]$LR)
    if (is.null(dt_f) || is.null(dt_m)) return(NULL)
    dt_f[, SEX := "Female"]; dt_m[, SEX := "Male"]
    cent_data <- rbind(dt_f, dt_m)
    avail <- centiles[centiles %in% names(cent_data)]
    long <- melt(cent_data, id.vars = c("AGE", "SEX"), measure.vars = avail,
                 variable.name = "Centile", value.name = "Value")
    long[, Centile := factor(Centile, levels = avail,
                             labels = centile_labels[avail])]
    long
  }

  # Base theme for all panels
  base_theme <- theme_publication(base_size = 9) +
    theme(plot.margin = margin(2, 4, 2, 4))

  # Line widths: median thickest, inner centiles medium, outer centiles thin
  lw_vals <- c("5th" = 0.4, "10th" = 0.5, "25th" = 0.6,
               "Median" = 1.1, "75th" = 0.6, "90th" = 0.5, "95th" = 0.4)

  # Helper: single-sex panel (Female or Male)
  panel_single <- function(long, sex, ylab, show_y = TRUE, show_x = FALSE) {
    d <- long[SEX == sex]
    p <- ggplot(d, aes(x = AGE, y = Value, linetype = Centile, linewidth = Centile)) +
      geom_line(color = sex_colors[sex]) +
      scale_linetype_manual(values = lt_vals, drop = FALSE) +
      scale_linewidth_manual(values = lw_vals, guide = "none") +
      scale_x_continuous(breaks = seq(50, 80, 10)) +
      labs(x = if (show_x) "Age (years)" else NULL,
           y = if (show_y) ylab else NULL) +
      base_theme +
      theme(legend.position = "none")
    if (!show_y) p <- p + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank())
    p
  }

  # Helper: overlaid both-sex panel (carries the legend aesthetics)
  panel_both <- function(long, ylab, show_x = FALSE) {
    ggplot(long, aes(x = AGE, y = Value, linetype = Centile, color = SEX,
                     linewidth = Centile)) +
      geom_line(alpha = 0.8) +
      scale_color_manual(values = sex_colors) +
      scale_linetype_manual(values = lt_vals, drop = FALSE) +
      scale_linewidth_manual(values = lw_vals, guide = "none") +
      scale_x_continuous(breaks = seq(50, 80, 10)) +
      labs(x = if (show_x) "Age (years)" else NULL, y = NULL) +
      base_theme +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.margin = margin(0, 0, 0, 0)) +
      guides(linetype = guide_legend(title = "Centile", ncol = 7,
                                     override.aes = list(color = "grey30")),
             color = guide_legend(title = "Sex", order = 1))
  }

  # Build all 9 panels (3 ROIs x 3 columns)
  panels <- list()
  for (i in seq_along(roi_defs)) {
    def <- roi_defs[[i]]
    long <- prep_data(def)
    if (is.null(long)) next
    is_bottom <- (i == length(roi_defs))
    panels <- c(panels, list(
      panel_single(long, "Female", def$ylab, show_y = TRUE, show_x = is_bottom),
      panel_both(long, def$ylab, show_x = is_bottom),
      panel_single(long, "Male", def$ylab, show_y = FALSE, show_x = is_bottom)
    ))
  }

  if (length(panels) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No normative data available"))
  }

  # Column header labels
  lab <- function(txt) {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = txt,
                        fontface = "bold", size = 3.5) + theme_void()
  }

  # Assemble: header row + 3 ROI rows + guide area
  # Layout: 4 rows (header + 3 ROI) x 3 columns, plus guide_area spanning bottom
  design <- "ABC
             DEF
             GHI
             JKL
             mmm"

  combined <- lab("Female") + lab("All") + lab("Male") +
    panels[[1]] + panels[[2]] + panels[[3]] +
    panels[[4]] + panels[[5]] + panels[[6]] +
    panels[[7]] + panels[[8]] + panels[[9]] +
    guide_area() +
    plot_layout(design = design,
                heights = c(0.8, 5, 5, 5, 1.5),
                guides = "collect")

  combined
}

#' Plot SEM forest plot for structural paths
#' @param sem_params SEM parameters from load_analysis_data()
#' @param models_to_include Character vector of model name patterns to include
#' @param cog_outcomes Character vector of cognitive outcomes to include
#' @return ggplot object
#' @export
fig_sem_forest <- function(sem_params, models_to_include = c("HC_COG", "HC_RES_COG", "HVR_COG"), cog_outcomes = c("g", "PRSP_s")) {
  if (is.null(sem_params)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  brain_labels <- c("HC" = "HC", "HVR" = "HVR", "LV" = "LV", "HC_RES" = "HC (Res)")
  cog_labels <- c("g" = "g", "MEM_s" = "Memory", "PRSP_s" = "Speed")

  paths_list <- list()
  for (mod_name in names(sem_params)) {
    if (!any(sapply(models_to_include, function(pat) grepl(pat, mod_name)))) next
    params <- sem_params[[mod_name]]; if (is.null(params)) next
    setDT(params)
    struct_paths <- params[op == "~" & lhs %in% cog_outcomes & rhs %in% c("HC", "HVR", "LV", "HC_RES")]
    if (nrow(struct_paths) > 0) { struct_paths[, MODEL := gsub("_POOLED$", "", mod_name)]; paths_list[[mod_name]] <- struct_paths }
  }
  if (length(paths_list) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No structural paths found"))

  all_paths <- rbindlist(paths_list, fill = TRUE)
  all_paths[, Sex := fifelse(is.na(group), "Combined", fifelse(group == 1, "Female", "Male"))]
  all_paths[, Brain_Label := brain_labels[rhs]][, Cog_Label := cog_labels[lhs]]
  all_paths[is.na(Brain_Label), Brain_Label := rhs][is.na(Cog_Label), Cog_Label := lhs]

  # Create Y-axis label with arrow symbol
  all_paths[, Y_LABEL := paste0(Brain_Label, " \u2192 ", Cog_Label)]
  all_paths[is.na(pvalue), pvalue := 1][, SIGN := pvalue < 0.05]
  all_paths <- all_paths[!is.na(est.std)]
  if (nrow(all_paths) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No valid path estimates"))

  # Define Y-axis order (paths grouped logically) with arrow symbol
  y_order <- c("HVR \u2192 g", "HVR \u2192 Speed",
               "HC (Res) \u2192 g", "HC (Res) \u2192 Speed",
               "HC \u2192 g", "HC \u2192 Speed")
  all_paths[, Y_LABEL := factor(Y_LABEL, levels = rev(y_order))]
  all_paths[, Sex := factor(Sex, levels = c("Combined", "Female", "Male"))]
  sex_colors <- c("Female" = "darkred", "Male" = "midnightblue", "Combined" = "#7F7F7F")

  # Significance-based styling: white fill for n.s., NO alpha
  # Use actual color values for fill
  all_paths[, FILL_COLOR := fifelse(SIGN, sex_colors[as.character(Sex)], "white")]

  # Sort data by Y_LABEL and Sex to ensure consistent dodge order
  setorder(all_paths, Y_LABEL, Sex)

  pd <- position_dodge(width = 0.7)

  ggplot(all_paths, aes(y = Y_LABEL, x = est.std, color = Sex, group = Sex)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # Error bars FIRST (so points with white fill are on top)
    geom_errorbarh(aes(xmin = ci.lower, xmax = ci.upper),
                   height = 0.25, linewidth = 0.6, position = pd) +
    # Points on top with explicit fill colors - NO alpha
    geom_point(aes(fill = FILL_COLOR),
               size = 1.8, shape = 21, stroke = 0.6, position = pd) +
    scale_color_manual(values = sex_colors, name = "Sex") +
    scale_fill_identity() +
    labs(y = NULL, x = bquote("Standardized Path Coefficient (" * beta * ")")) +
    theme_publication(base_size = 12) +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 11, face = "bold")) +
    guides(color = guide_legend(title = "Sex",
                                override.aes = list(fill = c("#7F7F7F", "darkred", "midnightblue"),
                                                    shape = 21, size = 3)))
}

#' Plot hemisphere comparison (forest plot style)
#' @param sex_diff Sex differences data from load_analysis_data()
#' @return ggplot object
#' @export
fig_hemisphere_comparison <- function(sex_diff, base_size = 10) {
  if (is.null(sex_diff) || !"OVERALL" %in% names(sex_diff)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  # Get HC and LV with RES adjustment, HVR with NON (self-normalizing)
  hc_lv_data <- sex_diff$OVERALL[SIDE %in% c("L", "R") & ADJ == "RES" & ROI %in% c("HC", "LV")]
  hvr_data <- sex_diff$OVERALL[SIDE %in% c("L", "R") & ADJ == "NON" & ROI == "HVR"]
  hemi_data <- rbind(hc_lv_data, hvr_data, fill = TRUE)
  if (nrow(hemi_data) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Hemisphere data not available"))

  # Use orange/teal palette (avoid blue/red used for sex)
  side_colors <- c("Left" = "#FF8C00", "Right" = "#008080", "white" = "white")
  hemi_data[, SIDE_LABEL := fifelse(SIDE == "L", "Left", "Right")]
  hemi_data[, SIDE_LABEL := factor(SIDE_LABEL, levels = c("Left", "Right"))]

  # Y-axis = ROI full names (HVR at top, HC middle, LV bottom)
  roi_labels <- c("HC" = "Hippocampus", "LV" = "Lateral Ventricles", "HVR" = "Hippocampal-to-Ventricle Ratio")
  hemi_data[, Y_LABEL := roi_labels[ROI]]
  hemi_data[is.na(Y_LABEL), Y_LABEL := ROI]
  # Order: LV bottom, HC middle, HVR top (last level appears at top in ggplot Y-axis)
  roi_order <- c("Lateral Ventricles", "Hippocampus", "Hippocampal-to-Ventricle Ratio")
  hemi_data[, Y_LABEL := factor(Y_LABEL, levels = roi_order)]

  # Significance-based styling
  hemi_data[, SIGNIFICANT := (CI_LOWER > 0) | (CI_UPPER < 0)]
  # Alpha based on significance
  hemi_data[, LINE_ALPHA := fifelse(SIGNIFICANT, 1.0, 0.4)]

  # Sort data to ensure consistent dodge order
  setorder(hemi_data, Y_LABEL, SIDE_LABEL)

  pd <- position_dodge(width = 0.6)

  # Forest plot style - errorbars with color distinction only
  p <- ggplot(hemi_data, aes(y = Y_LABEL, x = ESTIMATE, color = SIDE_LABEL, group = SIDE_LABEL)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_LOWER, xmax = CI_UPPER, alpha = LINE_ALPHA),
                   height = 0.25, linewidth = 0.7, position = pd) +
    scale_color_manual(values = side_colors[1:2], name = "Hemisphere") +
    scale_alpha_identity() +
    labs(y = NULL, x = "Effect Size (Cohen's d)\nPositive = Male < Female") +
    theme_publication(base_size = base_size) +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = base_size, face = "bold"))

  # Add difference panel (Right - Left) colored by sex advantage direction
  left_data <- hemi_data[SIDE == "L"]
  right_data <- hemi_data[SIDE == "R"]
  if (nrow(left_data) > 0 && nrow(right_data) > 0) {
    setkey(left_data, ROI)
    setkey(right_data, ROI)
    diff_ci <- compute_difference_ci(
      right_data$ESTIMATE, right_data$CI_LOWER, right_data$CI_UPPER,
      left_data$ESTIMATE, left_data$CI_LOWER, left_data$CI_UPPER
    )
    diff_data <- cbind(right_data[, .(ROI, Y_LABEL)], diff_ci)

    # Determine sex advantage from bilateral overall effects
    # Positive bilateral d = female advantage; negative = male advantage
    bilateral <- sex_diff$OVERALL[SIDE == "LR"]
    sex_adv <- bilateral[ROI %in% c("HC", "LV", "HVR"), .(ROI, BILATERAL_D = ESTIMATE)]
    # HC/LV use RES, HVR uses NON
    sex_adv <- rbind(
      bilateral[SIDE == "LR" & ROI %in% c("HC", "LV") & ADJ == "RES", .(ROI, BILATERAL_D = ESTIMATE)],
      bilateral[SIDE == "LR" & ROI == "HVR" & ADJ == "NON", .(ROI, BILATERAL_D = ESTIMATE)]
    )
    diff_data <- merge(diff_data, sex_adv, by = "ROI", all.x = TRUE)
    sex_colors <- c("Female" = "darkred", "Male" = "midnightblue")
    diff_data[, ADJ_LABEL := fifelse(BILATERAL_D > 0, "Female", "Male")]

    diff_panel <- plot_difference_panel(diff_data, roi_order, sex_colors,
                                         x_label = bquote(Delta * "d (Right - Left)"),
                                         base_size = base_size,
                                         show_legend = TRUE,
                                         legend_title = "Sex advantage")
    # Use guides="collect" to move legends to bottom, widen diff panel for readability
    p <- p + diff_panel + plot_layout(widths = c(3, 2), guides = "collect") &
      theme(legend.position = "bottom")
  }

  p
}

#' Plot sensitivity comparison (forest plot matching Fig 2 style)
#' Y-axis = ROI (Adjustment), shape = sample type, colors = adjustment
#' @param sex_diff Sex differences data from load_analysis_data()
#' @return ggplot object
#' @export
fig_sensitivity_comparison <- function(sex_diff, base_size = 10) {
  if (is.null(sex_diff) || !"SENS_COMPARISON" %in% names(sex_diff)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  sens_data <- copy(sex_diff$SENS_COMPARISON)
  if ("SIDE" %in% names(sens_data)) sens_data <- sens_data[SIDE == "LR"]
  if (nrow(sens_data) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Sensitivity data not available"))

  # Full adjustment labels for Y-axis (matching Fig 2)
  adj_full <- c("NON" = "Unadjusted", "PRP" = "Proportions", "STX" = "Stereotaxic", "RES" = "Residualized")

  # Create Y-axis labels and ADJ_LABEL for colors
  sens_data[ROI == "HC", Y_LABEL := paste0("Hippocampus\n(", adj_full[ADJ], ")")]
  sens_data[ROI == "HC", ADJ_LABEL := adj_full[ADJ]]
  sens_data[ROI == "LV", Y_LABEL := paste0("Lateral Ventricles\n(", adj_full[ADJ], ")")]
  sens_data[ROI == "LV", ADJ_LABEL := adj_full[ADJ]]
  sens_data[ROI == "HVR", Y_LABEL := "Hippocampal-to-Ventricle Ratio\n(Self-normalized)"]
  sens_data[ROI == "HVR", ADJ_LABEL := "Self-normalized"]

  # Melt data for estimates
  id_vars <- c("ROI", "ADJ", "Y_LABEL", "ADJ_LABEL")
  plot_data <- melt(sens_data, id.vars = id_vars,
                    measure.vars = c("ESTIMATE_PRIMARY", "ESTIMATE_SENS"),
                    variable.name = "Sample", value.name = "ESTIMATE")
  plot_data[, Sample := fifelse(Sample == "ESTIMATE_PRIMARY", "Primary", "Sensitivity")]

  # Add CI columns if available
  if (all(c("CI_LOWER_PRIMARY", "CI_UPPER_PRIMARY", "CI_LOWER_SENS", "CI_UPPER_SENS") %in% names(sens_data))) {
    ci_lower <- melt(sens_data, id.vars = id_vars,
                     measure.vars = c("CI_LOWER_PRIMARY", "CI_LOWER_SENS"),
                     variable.name = "Sample", value.name = "CI_LOWER")
    ci_lower[, Sample := fifelse(Sample == "CI_LOWER_PRIMARY", "Primary", "Sensitivity")]

    ci_upper <- melt(sens_data, id.vars = id_vars,
                     measure.vars = c("CI_UPPER_PRIMARY", "CI_UPPER_SENS"),
                     variable.name = "Sample", value.name = "CI_UPPER")
    ci_upper[, Sample := fifelse(Sample == "CI_UPPER_PRIMARY", "Primary", "Sensitivity")]

    plot_data <- merge(plot_data, ci_lower[, .(ROI, ADJ, Sample, CI_LOWER)], by = c("ROI", "ADJ", "Sample"))
    plot_data <- merge(plot_data, ci_upper[, .(ROI, ADJ, Sample, CI_UPPER)], by = c("ROI", "ADJ", "Sample"))
  }

  # Y-axis order: HVR at top, then HC, then LV (last level appears at top)
  y_order <- c("Lateral Ventricles\n(Residualized)", "Lateral Ventricles\n(Stereotaxic)",
               "Lateral Ventricles\n(Proportions)", "Lateral Ventricles\n(Unadjusted)",
               "Hippocampus\n(Residualized)", "Hippocampus\n(Stereotaxic)",
               "Hippocampus\n(Proportions)", "Hippocampus\n(Unadjusted)",
               "Hippocampal-to-Ventricle Ratio\n(Self-normalized)")
  plot_data[, Y_LABEL := factor(Y_LABEL, levels = y_order)]

  # Sample factor - Primary first (will appear above in dodge)
  plot_data[, Sample := factor(Sample, levels = c("Sensitivity", "Primary"))]

  # Adjustment colors (same for both samples)
  adj_colors_primary <- c("Unadjusted" = "#E64B35", "Proportions" = "#F39B7F",
                          "Stereotaxic" = "#00A087", "Residualized" = "#3C5488", "Self-normalized" = "#9467BD")

  # Significance-based styling
  if ("CI_LOWER" %in% names(plot_data) && "CI_UPPER" %in% names(plot_data)) {
    plot_data[, SIGNIFICANT := (CI_LOWER > 0) | (CI_UPPER < 0)]
  } else {
    plot_data[, SIGNIFICANT := TRUE]
  }

  # Alpha based on significance (dim non-significant results)
  plot_data[, LINE_ALPHA := fifelse(SIGNIFICANT, 1.0, 0.5)]

  # Convert Y_LABEL to numeric for manual dodge
  plot_data[, Y_NUM := as.numeric(Y_LABEL)]

  # Manual vertical offset to prevent overlap (0.15 units apart)
  primary_data <- copy(plot_data[Sample == "Primary"])
  sens_data_plot <- copy(plot_data[Sample == "Sensitivity"])
  primary_data[, Y_DODGE := Y_NUM + 0.15]
  sens_data_plot[, Y_DODGE := Y_NUM - 0.15]

  # Assign colors (same palette for both)
  primary_data[, COLOR := adj_colors_primary[ADJ_LABEL]]
  sens_data_plot[, COLOR := adj_colors_primary[ADJ_LABEL]]

  # Combine for plotting with shape legend
  combined_data <- rbind(primary_data, sens_data_plot)

  # Forest plot matching Fig 2 style - tiny points with different shapes
  p <- ggplot(combined_data, aes(y = Y_DODGE, x = ESTIMATE)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_LOWER, xmax = CI_UPPER, color = COLOR, alpha = LINE_ALPHA),
                   height = 0.2, linewidth = 0.6) +
    geom_point(aes(color = COLOR, alpha = LINE_ALPHA, shape = Sample), size = 1.5) +
    scale_color_identity() +
    scale_alpha_identity() +
    scale_shape_manual(values = c("Primary" = 16, "Sensitivity" = 17),
                       name = "Sample") +
    scale_y_continuous(breaks = seq_along(y_order), labels = y_order) +
    labs(y = NULL, x = "Effect Size (Cohen's d)\nPositive = Male < Female") +
    theme_publication(base_size = base_size) +
    theme(axis.text.y = element_text(size = base_size - 1),
          legend.position = "bottom") +
    guides(shape = guide_legend(override.aes = list(size = 2)))

  # Add difference panel
  sens_data_orig <- copy(sex_diff$SENS_COMPARISON)
  if ("SIDE" %in% names(sens_data_orig)) sens_data_orig <- sens_data_orig[SIDE == "LR"]

  if (all(c("CI_LOWER_PRIMARY", "CI_UPPER_PRIMARY", "CI_LOWER_SENS", "CI_UPPER_SENS") %in% names(sens_data_orig))) {
    diff_ci <- compute_difference_ci(
      sens_data_orig$ESTIMATE_PRIMARY, sens_data_orig$CI_LOWER_PRIMARY, sens_data_orig$CI_UPPER_PRIMARY,
      sens_data_orig$ESTIMATE_SENS, sens_data_orig$CI_LOWER_SENS, sens_data_orig$CI_UPPER_SENS
    )
    diff_data <- cbind(sens_data_orig[, .(ROI, ADJ)], diff_ci)

    diff_data[ROI == "HC", Y_LABEL := paste0("Hippocampus\n(", adj_full[ADJ], ")")]
    diff_data[ROI == "HC", ADJ_LABEL := adj_full[ADJ]]
    diff_data[ROI == "LV", Y_LABEL := paste0("Lateral Ventricles\n(", adj_full[ADJ], ")")]
    diff_data[ROI == "LV", ADJ_LABEL := adj_full[ADJ]]
    diff_data[ROI == "HVR", Y_LABEL := "Hippocampal-to-Ventricle Ratio\n(Self-normalized)"]
    diff_data[ROI == "HVR", ADJ_LABEL := "Self-normalized"]

    diff_panel <- plot_difference_panel(diff_data, y_order, adj_colors_primary,
                                         x_label = bquote(Delta * "d (Primary - Sensitivity)"),
                                         base_size = base_size)
    p <- p + diff_panel + plot_layout(widths = c(3, 1))
  }

  p
}

#' Plot combined GAMLSS validation (calibration) with facets
#' @param gamlss_test GAMLSS TEST data from load_gamlss_calibration()
#' @return ggplot object with HVR and HC (Residualized) faceted
#' @export
fig_gamlss_validation_combined <- function(gamlss_test) {
  if (is.null(gamlss_test)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  roi_config <- list(list(roi = "HVR", adj = "NON", label = "HVR"), list(roi = "HC", adj = "RES", label = "HC (Residualized)"))
  cal_list <- list()
  for (cfg in roi_config) {
    for (sex in c("Female", "Male")) {
      cal_data <- gamlss_test[[sex]][[cfg$roi]][[cfg$adj]]$LR$CENT_VALID
      if (!is.null(cal_data)) { cal_data <- copy(cal_data); cal_data[, SEX := sex]; cal_data[, MEASURE := cfg$label]; cal_list[[paste(cfg$roi, sex)]] <- cal_data }
    }
  }
  if (length(cal_list) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Calibration data not available"))
  cal_combined <- rbindlist(cal_list)
  cal_combined[, MEASURE := factor(MEASURE, levels = c("HVR", "HC (Residualized)"))]
  sex_colors <- c("Female" = "darkred", "Male" = "midnightblue")

  ggplot(cal_combined, aes(x = CENT_exp, y = CENT_obs, color = SEX)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray70", linewidth = 0.3, alpha = 0.4) +
    geom_line(linewidth = 0.4, alpha = 0.5) +
    geom_point(size = 1.5, shape = 1, alpha = 0.7) +
    facet_wrap(~ MEASURE, ncol = 2) +
    scale_color_manual(values = sex_colors, name = "Sex") +
    coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    labs(x = "Expected Centile (%)", y = "Observed Centile (%)") +
    theme_publication(base_size = 11) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 12))
}

#' Plot overall GAMLSS validation (calibration) for a ROI
#' @param gamlss_test GAMLSS TEST data from load_gamlss_calibration()
#' @param roi ROI code ("HC", "HVR", "LV")
#' @return ggplot object
#' @export
fig_gamlss_validation_overall <- function(gamlss_test, roi) {
  if (is.null(gamlss_test)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  adj <- if (roi == "HVR") "NON" else "RES"
  cal_list <- list()
  for (sex in c("Female", "Male")) {
    cal_data <- gamlss_test[[sex]][[roi]][[adj]]$LR$CENT_VALID
    if (!is.null(cal_data)) { cal_data <- copy(cal_data); cal_data[, SEX := sex]; cal_list[[sex]] <- cal_data }
  }
  if (length(cal_list) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Calibration data not available"))
  cal_combined <- rbindlist(cal_list)
  sex_colors <- c("Female" = "darkred", "Male" = "midnightblue")
  roi_labels <- c("HC" = "Hippocampus", "HVR" = "HVR", "LV" = "Lateral Ventricles")

  ggplot(cal_combined, aes(x = CENT_exp, y = CENT_obs, color = SEX)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray70", linewidth = 0.3, alpha = 0.4) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = sex_colors, name = "Sex") +
    coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    labs(title = paste0(roi_labels[roi], " Model Calibration"), x = "Expected Centile (%)", y = "Observed Centile (%)") +
    theme_publication(base_size = 11) +
    theme(legend.position = "bottom")
}

#' Plot sex-stratified GAMLSS validation (faceted by adjustment)
#' @param gamlss_test GAMLSS TEST data from load_gamlss_calibration()
#' @param roi ROI code ("HC", "HVR", "LV")
#' @param sex Sex ("Female" or "Male")
#' @return ggplot object
#' @export
fig_gamlss_validation_faceted <- function(gamlss_test, roi, sex) {
  if (is.null(gamlss_test)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  sex_data <- gamlss_test[[sex]][[roi]]
  if (is.null(sex_data)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))

  adj_labels <- c("NON" = "Unadjusted", "PRP" = "Proportions", "STX" = "Stereotaxic", "RES" = "Residualized")
  cal_list <- list()
  for (adj in names(sex_data)) {
    cal_data <- sex_data[[adj]]$LR$CENT_VALID
    if (!is.null(cal_data)) { cal_data <- copy(cal_data); cal_data[, ADJ := adj_labels[adj]]; cal_list[[adj]] <- cal_data }
  }
  if (length(cal_list) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Calibration data not available"))
  cal_combined <- rbindlist(cal_list)
  roi_labels <- c("HC" = "Hippocampus", "HVR" = "HVR", "LV" = "Lateral Ventricles")
  sex_color <- if (sex == "Female") "darkred" else "midnightblue"

  ggplot(cal_combined, aes(x = CENT_exp, y = CENT_obs)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_point(size = 2.5, alpha = 0.8, color = sex_color) +
    geom_line(color = sex_color, linewidth = 0.8) +
    facet_wrap(~ ADJ, ncol = 2) +
    coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    labs(title = paste0(roi_labels[roi], " Model Calibration (", sex, ")"), subtitle = "Faceted by adjustment method",
         x = "Expected Centile (%)", y = "Observed Centile (%)") +
    theme_publication(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
}

#' Plot age trajectory for lateral ventricles (supplementary)
#' @param sex_diff Sex differences data from load_analysis_data()
#' @return ggplot object
#' @export
fig_age_trajectory_lv <- function(sex_diff) {
  if (is.null(sex_diff) || !"AGE_STRATIFIED" %in% names(sex_diff)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  }

  age_strat <- sex_diff$AGE_STRATIFIED
  lv_data <- age_strat[SIDE == "LR" & ROI == "LV"]

  if (nrow(lv_data) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No LV data"))
  }

  adj_labels <- c("NON" = "Unadjusted", "PRP" = "Proportions",
                  "STX" = "Stereotaxic", "RES" = "Residualized")
  lv_data[, ADJ_LABEL := adj_labels[ADJ]]

  lv_data[, AGE_MID := sapply(AGE_BIN, function(x) {
    mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(x)), ","))))
  })]
  lv_data[, SIGNIFICANT := (CI_LOWER > 0) | (CI_UPPER < 0)]

  adj_colors <- c(
    "Unadjusted" = "#E64B35",
    "Proportions" = "#F39B7F",
    "Stereotaxic" = "#00A087",
    "Residualized" = "#3C5488"
  )

  lv_data[, ADJ_LABEL := factor(ADJ_LABEL, levels = names(adj_colors))]
  pd <- position_dodge(width = 2.5)  # Increased dodge for less overlap

  # Fill color: white for non-significant, method color for significant
  lv_data[, FILL_COLOR := fifelse(SIGNIFICANT, as.character(ADJ_LABEL), "white")]

  ggplot(lv_data, aes(x = AGE_MID, y = ESTIMATE, color = ADJ_LABEL, group = ADJ_LABEL)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = CI_LOWER, ymax = CI_UPPER, fill = ADJ_LABEL),
                alpha = 0.08, color = NA, position = pd) +
    # Dashed lines for LV (consistent with HC adjustment methods)
    geom_line(linetype = "dashed", linewidth = 0.8, position = pd) +
    # Points with white fill for non-significant - NO alpha
    geom_point(aes(fill = FILL_COLOR), size = 3.5, shape = 21, stroke = 1, position = pd) +
    scale_color_manual(values = adj_colors, name = "Adjustment") +
    scale_fill_manual(values = c(adj_colors, "white" = "white"), guide = "none") +
    labs(x = "Age (years)", y = "Sex Difference (Cohen's d)\nPositive = Male < Female",
         caption = "White-filled points: 95% CI crosses zero (n.s.).") +
    theme_publication(base_size = 10) +
    theme(legend.position = "none",
          plot.caption = element_text(size = 8, hjust = 0, face = "italic"))
}

#' Plot distribution overlap for HC adjustment methods and HVR
#' @param brain_data Brain volume data from load_brain_volumes()
#' @param matched_brain_data Optional matched sample brain data for comparison
#' @param show_matched Whether to show matched sample distributions (default FALSE)
#' @return ggplot object
#' @export
fig_distribution_overlap <- function(brain_data, matched_brain_data = NULL, show_matched = FALSE) {
  if (is.null(brain_data)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  }

  # Helper function to prepare data
  prepare_dist_data <- function(data, sample_label) {
    plot_data <- data[SIDE == "LR" & SUBFIELD == "total"]
    hc_data <- plot_data[, .(SEX, ADJ, VAL = HC)]
    hc_data[, MEASURE := paste0("HC_", ADJ)]
    hvr_data <- plot_data[ADJ == "NON", .(SEX, ADJ = "HVR", VAL = HVR)]
    hvr_data[, MEASURE := "HVR"]
    combined <- rbind(hc_data, hvr_data[, .(SEX, ADJ, VAL, MEASURE)])
    combined[, SAMPLE := sample_label]
    combined
  }

  measure_order <- c("HC_NON", "HC_PRP", "HC_STX", "HC_RES", "HVR")
  measure_labels <- c(
    "HC_NON" = "HC (Unadjusted)",
    "HC_PRP" = "HC (Proportions)",
    "HC_STX" = "HC (Stereotaxic)",
    "HC_RES" = "HC (Residualized)",
    "HVR" = "HVR"
  )

  full_data <- prepare_dist_data(brain_data, "Full Sample")
  sex_colors <- c("Female" = "darkred", "Male" = "midnightblue")

  if (show_matched && !is.null(matched_brain_data) && nrow(matched_brain_data) > 0) {
    matched_data <- prepare_dist_data(matched_brain_data, "Head-Size Matched")
    combined <- rbind(full_data, matched_data)
    combined <- combined[MEASURE %in% measure_order]
    combined[, MEASURE := factor(MEASURE, levels = measure_order, labels = measure_labels[measure_order])]
    combined[, SAMPLE := factor(SAMPLE, levels = c("Full Sample", "Head-Size Matched"))]
    combined[, VAL_z := scale(VAL), by = MEASURE]

    ggplot(combined, aes(x = VAL_z, fill = SEX, color = SEX)) +
      geom_density(alpha = 0.4, linewidth = 0.5) +
      facet_grid(MEASURE ~ SAMPLE, scales = "free_y") +
      scale_fill_manual(values = sex_colors, name = "Sex") +
      scale_color_manual(values = sex_colors, guide = "none") +
      labs(x = "Standardized Value (z-score)", y = "Density") +
      theme_publication(base_size = 10) +
      theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 8))
  } else {
    combined <- full_data
    combined <- combined[MEASURE %in% measure_order]
    combined[, MEASURE := factor(MEASURE, levels = measure_order, labels = measure_labels[measure_order])]
    combined[, VAL_z := scale(VAL), by = MEASURE]

    ggplot(combined, aes(x = VAL_z, fill = SEX, color = SEX)) +
      geom_density(alpha = 0.4, linewidth = 0.5) +
      facet_wrap(~ MEASURE, ncol = 5, scales = "free_y") +
      scale_fill_manual(values = sex_colors, name = "Sex") +
      scale_color_manual(values = sex_colors, guide = "none") +
      labs(x = "Standardized Value (z-score)", y = "Density") +
      theme_publication(base_size = 10) +
      theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 8))
  }
}

#' Plot generic GAMLSS validation/calibration
#' @param valid_data Validation data with Expected and Observed columns
#' @param title Plot title
#' @return ggplot object
#' @export
fig_gamlss_validation <- function(valid_data, title = "Model Calibration") {
  if (is.null(valid_data) || nrow(valid_data) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Validation data not available"))
  }

  ggplot(valid_data, aes(x = Expected, y = Observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "#E64B35") +
    labs(
      title = title,
      x = "Expected Proportion",
      y = "Observed Proportion"
    ) +
    coord_equal() +
    theme_publication(base_size = 11)
}

#' Plot z-score trajectories
#' @param norm_tables Normative tables
#' @param sex Sex to plot ("Female" or "Male")
#' @param roi ROI to plot ("HC", "HVR", "LV")
#' @return ggplot object
#' @export
fig_zscore_trajectories <- function(norm_tables, sex = "Female", roi = "HC") {
  if (is.null(norm_tables)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  }

  # Get centile data for residualized
  cent_data <- norm_tables[[sex]][[roi]]$RES$LR

  if (is.null(cent_data)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Centile data not available"))
  }

  # Calculate z-scores from centiles (approximate using median and spread)
  # p50 is median (z=0), p2.5 is ~z=-2, p97.5 is ~z=+2
  cent_data[, z_neg2 := (p2.5 - p50) / (p50 - p2.5) * 2]  # normalized

  centiles <- c("p2.5", "p10", "p25", "p50", "p75", "p90", "p97.5")
  cent_long <- melt(cent_data, id.vars = "AGE", measure.vars = centiles,
                    variable.name = "Centile", value.name = "Value")

  ggplot(cent_long, aes(x = AGE, y = Value, linetype = Centile)) +
    geom_line(color = ifelse(sex == "Female", "darkred", "midnightblue"), linewidth = 0.6) +
    scale_linetype_manual(values = c(
      "p2.5" = "dotted", "p10" = "dashed", "p25" = "dotdash",
      "p50" = "solid",
      "p75" = "dotdash", "p90" = "dashed", "p97.5" = "dotted"
    )) +
    labs(
      title = paste0(roi, " (", sex, ")"),
      x = "Age (years)",
      y = "Volume (residualized)",
      linetype = "Centile"
    ) +
    theme_publication(base_size = 11) +
    theme(legend.position = "right")
}

#' Plot HC centile curves by sex (faceted)
#' @param norm_tables Normative tables from load_norm_tables()
#' @return ggplot object with sex facets
#' @export
fig_hc_centiles_sexcomp <- function(norm_tables) {
  if (is.null(norm_tables)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  }

  # Get HC centiles (unadjusted with ICV covariate in GAMLSS, bilateral)
  hc_f <- norm_tables$Female$HC$NON$LR
  hc_m <- norm_tables$Male$HC$NON$LR

  if (is.null(hc_f) || is.null(hc_m)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "HC data not available"))
  }

  hc_f[, SEX := "Female"]
  hc_m[, SEX := "Male"]

  cent_data <- rbind(hc_f, hc_m)

  # Consistent centile levels with HVR (F3)
  centiles <- c("p5", "p10", "p25", "p50", "p75", "p90", "p95")
  centile_labels <- c("p5" = "5th", "p10" = "10th", "p25" = "25th",
                      "p50" = "Median", "p75" = "75th", "p90" = "90th", "p95" = "95th")
  centiles <- centiles[centiles %in% names(cent_data)]
  cent_long <- melt(cent_data, id.vars = c("AGE", "SEX"), measure.vars = centiles,
                    variable.name = "Centile", value.name = "Value")
  cent_long[, Centile := factor(Centile, levels = centiles, labels = centile_labels[centiles])]

  sex_colors <- c("Female" = "darkred", "Male" = "midnightblue")

  ggplot(cent_long, aes(x = AGE, y = Value, linetype = Centile, color = SEX)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ SEX, ncol = 2) +
    scale_color_manual(values = sex_colors, guide = "none") +
    scale_linetype_manual(
      values = c(
        "5th" = "dotted", "10th" = "dashed", "25th" = "dotdash",
        "Median" = "solid",
        "75th" = "dotdash", "90th" = "dashed", "95th" = "dotted"
      ), drop = FALSE
    ) +
    labs(
      x = "Age (years)",
      y = "Hippocampal Volume (mm³)",
      linetype = "Centile"
    ) +
    theme_publication(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 11)
    )
}

#' Plot LV centile curves by sex (faceted)
#' @param norm_tables Normative tables from load_norm_tables()
#' @return ggplot object with sex facets
#' @export
fig_lv_centiles_sexcomp <- function(norm_tables) {
  if (is.null(norm_tables)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available"))
  }

  # Get LV centiles (unadjusted, bilateral)
  lv_f <- norm_tables$Female$LV$NON$LR
  lv_m <- norm_tables$Male$LV$NON$LR

  if (is.null(lv_f) || is.null(lv_m)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "LV data not available"))
  }

  lv_f[, SEX := "Female"]
  lv_m[, SEX := "Male"]

  cent_data <- rbind(lv_f, lv_m)

  # Consistent centile levels with HVR (F3) and HC (FS3)
  centiles <- c("p5", "p10", "p25", "p50", "p75", "p90", "p95")
  centile_labels <- c("p5" = "5th", "p10" = "10th", "p25" = "25th",
                      "p50" = "Median", "p75" = "75th", "p90" = "90th", "p95" = "95th")
  centiles <- centiles[centiles %in% names(cent_data)]
  cent_long <- melt(cent_data, id.vars = c("AGE", "SEX"), measure.vars = centiles,
                    variable.name = "Centile", value.name = "Value")
  cent_long[, Centile := factor(Centile, levels = centiles, labels = centile_labels[centiles])]

  sex_colors <- c("Female" = "darkred", "Male" = "midnightblue")

  ggplot(cent_long, aes(x = AGE, y = Value, linetype = Centile, color = SEX)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ SEX, ncol = 2) +
    scale_color_manual(values = sex_colors, guide = "none") +
    scale_linetype_manual(
      values = c(
        "5th" = "dotted", "10th" = "dashed", "25th" = "dotdash",
        "Median" = "solid",
        "75th" = "dotdash", "90th" = "dashed", "95th" = "dotted"
      ), drop = FALSE
    ) +
    labs(
      x = "Age (years)",
      y = "Lateral Ventricle Volume (cc)",
      linetype = "Centile"
    ) +
    theme_publication(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 11)
    )
}
