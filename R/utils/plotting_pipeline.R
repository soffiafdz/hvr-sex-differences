# =============================================================================
# Pipeline Plotting Functions
# =============================================================================
# Exploratory and pipeline-specific plotting functions (plot_* functions).
# Part of the R/utils/ module system.
# =============================================================================

#' Create GAMLSS centile plot
#' @param obs_data Observed data
#' @param cent_data Centile predictions
#' @param age_var Age variable name
#' @param value_var Value variable name
#' @param sex_var Sex variable name
#' @param centiles Vector of centile columns
#' @param linetypes Named vector of linetypes
#' @return ggplot object
plot_gamlss_centiles <- function(
    obs_data, cent_data,
    age_var = "AGE", value_var = "VAL_scl",
    sex_var = "Sex", centiles = c("p10", "p25", "p50", "p75", "p90"),
    linetypes = NULL) {
  if (is.null(linetypes)) {
    linetypes <- c(
      "p10" = "dashed",
      "p25" = "dotdash",
      "p50" = "solid",
      "p75" = "dotdash",
      "p90" = "dashed"
    )
  }

  sex_colors <- get_palette("sex")

  # Melt centile data
  cent_long <- data.table::melt(
    cent_data,
    measure.vars = centiles,
    variable.name = "Centiles",
    value.name = "VAL_scl"
  )

  ggplot2::ggplot() +
    theme_publication(use_markdown = TRUE) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[age_var]],
        y = .data[[value_var]],
        colour = .data[[sex_var]]
      ),
      data = obs_data,
      alpha = 0.05,
      size = 0.3,
      shape = 21
    ) +
    ggplot2::scale_colour_manual(values = sex_colors) +
    ggplot2::geom_smooth(
      ggplot2::aes(
        x = .data[[age_var]],
        y = VAL_scl,
        linetype = Centiles
      ),
      data = cent_long[cent_long[[sex_var]] == "Male", ],
      se = FALSE,
      colour = sex_colors["Male"],
      method = "loess",
      span = 0.35,
      linewidth = 0.6
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(
        x = .data[[age_var]],
        y = VAL_scl,
        linetype = Centiles
      ),
      data = cent_long[cent_long[[sex_var]] == "Female", ],
      se = FALSE,
      colour = sex_colors["Female"],
      method = "loess",
      span = 0.35,
      linewidth = 0.6
    ) +
    ggplot2::scale_linetype_manual(values = linetypes)
}

#' Create GAMLSS centile curve plot
#' @param norm_table Data table with normative centiles
#' @param obs_data Optional observed data to overlay
#' @param centiles Vector of centiles to plot
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param x_label X-axis label
#' @param y_label Y-axis label
#' @return ggplot object
plot_gamlss_curves <- function(
    norm_table,
    obs_data = NULL,
    centiles = c("p2.5", "p10", "p25", "p50", "p75", "p90", "p97.5"),
    title = "Normative Centiles",
    subtitle = NULL,
    x_label = "Age (years)",
    y_label = "Volume (cc)") {
  # Filter available centiles
  centiles <- centiles[centiles %in% names(norm_table)]

  # Melt normative data
  cent_long.dt <- data.table::melt(
    norm_table,
    id.vars = "AGE",
    measure.vars = centiles,
    variable.name = "Centile",
    value.name = "Value"
  )

  # Define line types
  linetypes <- c(
    "p2.5" = "dotted", "p5" = "dotted",
    "p10" = "dashed", "p25" = "dotdash",
    "p50" = "solid",
    "p75" = "dotdash", "p90" = "dashed",
    "p95" = "dotted", "p97.5" = "dotted"
  )

  # Get sex from data if available
  sex_col <- if ("SEX" %in% names(norm_table)) {
    unique(norm_table$SEX)[1]
  } else {
    "Male"
  }
  line_color <- get_palette("sex")[sex_col]

  # Base plot
  p <- ggplot2::ggplot()

  # Add observed data if provided
  if (!is.null(obs_data) && nrow(obs_data) > 0) {
    p <- p +
      ggplot2::geom_point(
        data = obs_data,
        ggplot2::aes(x = AGE, y = VAL),
        alpha = 0.1, size = 0.5,
        color = line_color
      )
  }

  # Add centile curves
  p <- p +
    ggplot2::geom_line(
      data = cent_long.dt,
      ggplot2::aes(
        x = AGE, y = Value,
        group = Centile, linetype = Centile
      ),
      color = line_color,
      linewidth = 0.8
    ) +
    ggplot2::scale_linetype_manual(
      values = linetypes,
      labels = function(x) gsub("p", "", x)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label,
      linetype = "Centile"
    ) +
    theme_publication(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      legend.position = "right",
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )

  p
}

#' Create GAMLSS centile curve plot with prediction intervals
#'
#' Enhanced version of plot_gamlss_curves that adds 95% prediction intervals
#' to show individual-level uncertainty, not just population centiles.
#'
#' @param norm_table Data table with normative centiles (must include p2.5 and p97.5)
#' @param obs_data Optional observed data to overlay as scatter
#' @param centiles Vector of centiles to plot as lines
#' @param show_prediction_interval Whether to shade the 95% prediction interval
#' @param show_scatter Whether to show individual data points
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param x_label X-axis label
#' @param y_label Y-axis label
#' @return ggplot object
#' @export
plot_gamlss_curves_with_pi <- function(
    norm_table,
    obs_data = NULL,
    centiles = c("p10", "p25", "p50", "p75", "p90"),
    show_prediction_interval = TRUE,
    show_scatter = TRUE,
    title = "Normative Centiles with 95% Prediction Interval",
    subtitle = NULL,
    x_label = "Age (years)",
    y_label = "Volume (cc)") {

  # Check for prediction interval bounds
  has_pi <- all(c("p2.5", "p97.5") %in% names(norm_table))

  # Filter available centiles
  centiles <- centiles[centiles %in% names(norm_table)]

  # Melt normative data for centile lines
  cent_long.dt <- data.table::melt(
    norm_table,
    id.vars = "AGE",
    measure.vars = centiles,
    variable.name = "Centile",
    value.name = "Value"
  )

  # Define line types
  linetypes <- c(
    "p5" = "dotted", "p10" = "dashed", "p25" = "dotdash",
    "p50" = "solid",
    "p75" = "dotdash", "p90" = "dashed", "p95" = "dotted"
  )

  # Get sex from data if available
  sex_col <- if ("SEX" %in% names(norm_table)) {
    unique(norm_table$SEX)[1]
  } else {
    "Male"
  }
  line_color <- get_palette("sex")[sex_col]
  fill_color <- adjustcolor(line_color, alpha.f = 0.15)

  # Base plot
  p <- ggplot2::ggplot()

  # Add 95% prediction interval ribbon (p2.5 to p97.5)
  if (show_prediction_interval && has_pi) {
    p <- p +
      ggplot2::geom_ribbon(
        data = norm_table,
        ggplot2::aes(x = AGE, ymin = p2.5, ymax = p97.5),
        fill = line_color,
        alpha = 0.15
      )
  }

  # Add observed data if provided and requested
  if (show_scatter && !is.null(obs_data) && nrow(obs_data) > 0) {
    p <- p +
      ggplot2::geom_point(
        data = obs_data,
        ggplot2::aes(x = AGE, y = VAL),
        alpha = 0.08, size = 0.4,
        color = line_color
      )
  }

  # Add centile curves
  p <- p +
    ggplot2::geom_line(
      data = cent_long.dt,
      ggplot2::aes(
        x = AGE, y = Value,
        group = Centile, linetype = Centile
      ),
      color = line_color,
      linewidth = 0.8
    ) +
    ggplot2::scale_linetype_manual(
      values = linetypes,
      labels = function(x) gsub("p", "", x)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = if (is.null(subtitle) && show_prediction_interval && has_pi) {
        "Shaded region = 95% prediction interval; lines = population centiles"
      } else {
        subtitle
      },
      x = x_label,
      y = y_label,
      linetype = "Centile"
    ) +
    theme_publication(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      legend.position = "right",
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )

  p
}

#' Create GAMLSS validation plot (centile calibration)
#' @param validation_data Data table with expected and observed centiles
#' @param title Plot title
#' @return ggplot object
plot_gamlss_validation <- function(
    validation_data,
    title = "Centile Calibration") {
  ggplot2::ggplot(
    validation_data,
    ggplot2::aes(x = CENT_exp, y = CENT_obs)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", color = "grey50", linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 3, alpha = 0.7, color = "#0072B2") +
    ggplot2::geom_line(color = "#0072B2", linewidth = 1) +
    ggplot2::labs(
      title = title,
      subtitle = "Perfect calibration shown by dashed line",
      x = "Expected Centile (%)",
      y = "Observed Centile (%)"
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    theme_publication(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )
}

#' Create longitudinal stability plot
#' @param centile_data Data table with T1 and T2 centiles
#' @param title Plot title
#' @return ggplot object
plot_longitudinal_stability <- function(
    centile_data,
    title = "Longitudinal Centile Stability") {
  ggplot2::ggplot(
    centile_data,
    ggplot2::aes(x = CENT_t1, y = CENT_t2)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", color = "grey50", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      alpha = 0.4, size = 2, color = "#0072B2"
    ) +
    ggplot2::geom_smooth(
      method = "lm", se = TRUE,
      color = "#D55E00", fill = "#D55E00", alpha = 0.2
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Stable centiles should fall along the diagonal",
      x = "Centile at Time 1 (%)",
      y = "Centile at Time 2 (%)",
      caption = sprintf(
        "r = %.3f | Mean change = %.2f",
        cor(centile_data$CENT_t1, centile_data$CENT_t2, use = "complete.obs"),
        mean(centile_data$CENT_t2 - centile_data$CENT_t1, na.rm = TRUE)
      )
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    theme_publication(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      plot.caption = ggplot2::element_text(size = 10, hjust = 0),
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )
}

#' Create faceted centile plot by adjustment method
#' @param norm_tables_list Named list of norm tables by adjustment
#' @param obs_data_list Named list of obs data by adjustment
#' @param roi_label ROI label for title
#' @param sex Sex label for title
#' @param sex_color Color for sex-specific plotting
#' @param y_label Y-axis label
#' @return ggplot object
#' @export
plot_gamlss_faceted_by_adjustment <- function(
    norm_tables_list,
    obs_data_list,
    roi_label,
    sex,
    sex_color,
    y_label = "Volume (cc)") {
  # Combine data across adjustments
  norm_combined.dt <- rbindlist(norm_tables_list, fill = TRUE)
  obs_combined.dt <- rbindlist(obs_data_list, fill = TRUE)

  # Get centile columns
  cent_cols <- grep("^p[0-9.]+$", names(norm_combined.dt), value = TRUE)
  cent_cols <- cent_cols[cent_cols %in% c(
    "p2.5", "p10", "p25", "p50", "p75", "p90", "p97.5"
  )]

  # Melt for plotting
  cent_long.dt <- data.table::melt(
    norm_combined.dt,
    id.vars = c("AGE", "ADJ"),
    measure.vars = cent_cols,
    variable.name = "Centile",
    value.name = "Value"
  )

  # Define aesthetics
  linetypes <- c(
    "p2.5" = "dotted", "p10" = "dashed", "p25" = "dotdash",
    "p50" = "solid", "p75" = "dotdash", "p90" = "dashed",
    "p97.5" = "dotted"
  )

  # Create faceted plot
  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = obs_combined.dt,
      ggplot2::aes(x = AGE, y = VAL),
      alpha = 0.05, size = 0.2, color = sex_color
    ) +
    ggplot2::geom_line(
      data = cent_long.dt,
      ggplot2::aes(
        x = AGE, y = Value,
        linetype = Centile, group = Centile
      ),
      color = sex_color, linewidth = 0.6
    ) +
    ggplot2::facet_wrap(~ ADJ, ncol = 2, scales = "free_y") +
    ggplot2::scale_linetype_manual(
      values = linetypes,
      labels = function(x) gsub("p", "", x)
    ) +
    ggplot2::labs(
      title = sprintf("Normative Centiles: %s (%s)", roi_label, sex),
      subtitle = "Faceted by adjustment method",
      x = "Age (years)",
      y = y_label,
      linetype = "Centile"
    ) +
    theme_publication(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "grey95", color = "grey70"
      ),
      legend.position = "bottom",
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )
}

#' Create faceted sex comparison plot
#' @param norm_tables_list Named list (by adj_sex) of norm tables
#' @param obs_data_list Named list (by adj_sex) of obs data
#' @param roi_label ROI label for title
#' @param sex_colors Named vector of sex colors
#' @param y_label Y-axis label
#' @return ggplot object
#' @export
plot_gamlss_sexcomp_faceted <- function(
    norm_tables_list,
    obs_data_list,
    roi_label,
    sex_colors,
    y_label = "Volume (cc)") {
  # Combine all data
  norm_combined.dt <- rbindlist(norm_tables_list, fill = TRUE)
  obs_combined.dt <- rbindlist(obs_data_list, fill = TRUE)

  # Get centile columns (use fewer for clarity)
  cent_cols <- grep("^p[0-9.]+$", names(norm_combined.dt), value = TRUE)
  cent_cols <- cent_cols[cent_cols %in% c("p25", "p50", "p75")]

  # Melt for plotting
  cent_long.dt <- data.table::melt(
    norm_combined.dt,
    id.vars = c("AGE", "ADJ", "SEX"),
    measure.vars = cent_cols,
    variable.name = "Centile",
    value.name = "Value"
  )

  # Define aesthetics
  linetypes <- c("p25" = "dashed", "p50" = "solid", "p75" = "dashed")

  # Create faceted comparison plot
  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = obs_combined.dt,
      ggplot2::aes(x = AGE, y = VAL, color = SEX),
      alpha = 0.03, size = 0.2
    ) +
    ggplot2::geom_line(
      data = cent_long.dt,
      ggplot2::aes(
        x = AGE, y = Value,
        color = SEX, linetype = Centile,
        group = interaction(SEX, Centile)
      ),
      linewidth = 0.7
    ) +
    ggplot2::facet_grid(SEX ~ ADJ, scales = "free_y") +
    ggplot2::scale_color_manual(values = sex_colors, name = "Sex") +
    ggplot2::scale_linetype_manual(
      values = linetypes,
      labels = function(x) gsub("p", "", x),
      name = "Centile"
    ) +
    ggplot2::labs(
      title = sprintf("%s Normative Centiles", roi_label),
      subtitle = "Rows: Sex | Columns: Adjustment method",
      x = "Age (years)",
      y = y_label
    ) +
    theme_publication(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "grey95", color = "grey70"
      ),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 2)
    )
}

#' Create faceted validation plot for both sexes
#' @param valid_data_list Named list of validation data (by adj_sex)
#' @param roi_label ROI label for title
#' @param sex_colors Named vector of sex colors
#' @return ggplot object
#' @export
plot_validation_sexcomp_faceted <- function(
    valid_data_list, roi_label, sex_colors) {
  valid_combined.dt <- rbindlist(valid_data_list)

  ggplot2::ggplot(
    valid_combined.dt,
    ggplot2::aes(x = CENT_exp, y = CENT_obs, color = SEX)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", color = "grey50", linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 2, alpha = 0.7) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_grid(SEX ~ ADJ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    ggplot2::scale_color_manual(values = sex_colors, name = "Sex") +
    ggplot2::labs(
      title = sprintf("Centile Calibration: %s", roi_label),
      subtitle = "Rows: Sex | Columns: Adjustment method",
      x = "Expected Centile (%)",
      y = "Observed Centile (%)"
    ) +
    theme_publication(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "grey95", color = "grey70"
      ),
      legend.position = "bottom",
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )
}

#' Create faceted validation plot (single sex)
#' @param valid_data_list Named list of validation data by adjustment
#' @param roi_label ROI label for title
#' @param sex Sex label for title
#' @return ggplot object
#' @export
plot_validation_faceted <- function(valid_data_list, roi_label, sex) {
  valid_combined.dt <- rbindlist(valid_data_list)

  ggplot2::ggplot(
    valid_combined.dt,
    ggplot2::aes(x = CENT_exp, y = CENT_obs)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", color = "grey50", linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 2, alpha = 0.7, color = "#0072B2") +
    ggplot2::geom_line(color = "#0072B2", linewidth = 0.8) +
    ggplot2::facet_wrap(~ ADJ, ncol = 2) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100)) +
    ggplot2::labs(
      title = sprintf("Centile Calibration: %s (%s)", roi_label, sex),
      subtitle = "Faceted by adjustment method",
      x = "Expected Centile (%)",
      y = "Observed Centile (%)"
    ) +
    theme_publication(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      panel.grid.major = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )
}

#' Create faceted stability plot with summary metrics
#' @param stability_data_list Named list of stability summary data (by adj_sex)
#' @param roi_label ROI label for title
#' @param sex_colors Named vector of sex colors
#' @return ggplot object
#' @export
plot_stability_faceted <- function(
    stability_data_list, roi_label, sex_colors) {
  stab_combined.dt <- rbindlist(stability_data_list)

  # Reshape to long format for faceting by metric
  stab_long.dt <- data.table::melt(
    stab_combined.dt,
    id.vars = c("SEX", "ADJ"),
    measure.vars = c("CORR", "MEAN_CHG", "SD_CHG", "CROSS"),
    variable.name = "METRIC",
    value.name = "VALUE"
  )

  # Create labels for metrics
  metric_labels <- c(
    CORR = "Correlation (r)",
    MEAN_CHG = "Mean Change",
    SD_CHG = "SD Change",
    CROSS = "Threshold Crossing"
  )

  stab_long.dt[, METRIC_LABEL := metric_labels[METRIC]]

  ggplot2::ggplot(
    stab_long.dt,
    ggplot2::aes(x = ADJ, y = VALUE, fill = SEX)
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.7
    ) +
    ggplot2::facet_wrap(~ METRIC_LABEL, scales = "free_y", ncol = 2) +
    ggplot2::scale_fill_manual(values = sex_colors, name = "Sex") +
    ggplot2::labs(
      title = sprintf("Longitudinal Stability: %s", roi_label),
      subtitle = "Test-retest metrics across adjustment methods",
      x = "Adjustment Method",
      y = "Value"
    ) +
    theme_publication(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "grey95", color = "grey70"
      ),
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      ),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

#' Create z-score standardized centile plot for HC across adjustment methods
#'
#' Standardizes each adjustment method to z-scores (mean=0, sd=1) so they can
#' be plotted on the same scale for direct visual comparison.
#'
#' @param obs_data_list Named list of observed data by adjustment method
#' @param roi_label ROI label for title
#' @param sex_colors Named vector of sex colors
#' @param show_points Whether to show individual data points
#' @return ggplot object
#' @export
plot_zscore_trajectories <- function(
    obs_data_list,
    roi_label = "Hippocampal Volume",
    sex_colors = c(Female = "darkred", Male = "midnightblue"),
    show_points = FALSE) {

  # Combine and standardize within each adjustment method
  obs_combined.dt <- rbindlist(obs_data_list, fill = TRUE)

  # Standardize within each ADJ method

  obs_combined.dt[, ZSCORE := scale(VAL), by = ADJ]

  # Calculate smoothed trajectories by sex and adjustment
  smooth_data.dt <- obs_combined.dt[, {
    # Fit loess for smoothing
    fit <- loess(ZSCORE ~ AGE, span = 0.3)
    age_seq <- seq(min(AGE), max(AGE), length.out = 100)
    pred <- predict(fit, newdata = data.frame(AGE = age_seq))
    list(AGE = age_seq, ZSCORE_SMOOTH = pred)
  }, by = .(ADJ, SEX)]

  # Create plot
  p <- ggplot2::ggplot()

  if (show_points) {
    p <- p + ggplot2::geom_point(
      data = obs_combined.dt,
      ggplot2::aes(x = AGE, y = ZSCORE, color = SEX),
      alpha = 0.02, size = 0.3
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_line(
      data = smooth_data.dt,
      ggplot2::aes(x = AGE, y = ZSCORE_SMOOTH, color = SEX),
      linewidth = 1
    ) +
    ggplot2::facet_wrap(~ ADJ, ncol = 2) +
    ggplot2::scale_color_manual(values = sex_colors, name = "Sex") +
    ggplot2::scale_y_continuous(
      breaks = c(-2, -1, 0, 1, 2),
      labels = c("-2 SD", "-1 SD", "Mean", "+1 SD", "+2 SD")
    ) +
    ggplot2::labs(
      title = sprintf("%s: Z-Score Trajectories by Adjustment Method", roi_label),
      subtitle = "Standardized within each method for direct comparison",
      x = "Age (years)",
      y = "Z-Score"
    ) +
    theme_publication(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      legend.position = "bottom"
    )

  return(p)
}

#' Create combined HC vs HVR sex difference trajectory plot
#'
#' Shows age trajectories of sex difference effect sizes for HC (multiple
#' adjustment methods) and HVR (self-normalizing) on the same scale.
#'
#' @param effect_data Data.table with AGE_BIN, ROI, ADJ, ESTIMATE, CI_LOWER, CI_UPPER
#' @param highlight_hvr Whether to emphasize HVR line
#' @return ggplot object
#' @export
plot_hc_hvr_effect_comparison <- function(effect_data, highlight_hvr = TRUE) {

  # Create combined label for HC adjustment methods
  plot_data <- copy(effect_data)
  plot_data[, METHOD := fifelse(
    ROI == "HVR", "HVR (self-normalizing)",
    sprintf("HC (%s)", ADJ_LABEL)
  )]

  # Set line properties
  plot_data[, LINE_WIDTH := fifelse(ROI == "HVR", 1.2, 0.7)]
  plot_data[, LINE_ALPHA := fifelse(ROI == "HVR", 1.0, 0.7)]

  # Color palette: HVR distinct, HC methods in gradient
  method_colors <- c(
    "HVR (self-normalizing)" = "#009E73",
    "HC (Unadjusted)" = "#E69F00",
    "HC (Proportions)" = "#56B4E9",
    "HC (Stereotaxic)" = "#CC79A7",
    "HC (Residualized)" = "#0072B2"
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = AGE_BIN, y = ESTIMATE, color = METHOD, group = METHOD)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = CI_LOWER, ymax = CI_UPPER, fill = METHOD),
      alpha = 0.15, color = NA
    ) +
    ggplot2::geom_line(ggplot2::aes(linewidth = LINE_WIDTH, alpha = LINE_ALPHA)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = method_colors, name = "Measure") +
    ggplot2::scale_fill_manual(values = method_colors, guide = "none") +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_alpha_identity() +
    ggplot2::labs(
      title = "Sex Difference Trajectories: HC vs HVR",
      subtitle = "Effect size (Cohen's d) across age; positive = females higher",
      x = "Age Group",
      y = "Cohen's d"
    ) +
    theme_publication(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

#' Create simple HVR centile plot (self-normalizing, both sexes)
#'
#' Clean plot showing HVR normative centiles for both sexes side-by-side.
#' HVR is inherently self-normalizing (ratio measure), so no head-size adjustment is needed.
#' Uses raw HVR values (0-1 scale) for clinical interpretability.
#'
#' @param norm_tables_list List with Female and Male normative tables
#' @param obs_data_list List with Female and Male observed data
#' @param sex_colors Named vector of sex colors
#' @return ggplot object
#' @export
plot_hvr_simple <- function(
    norm_tables_list,
    obs_data_list,
    sex_colors = c(Female = "darkred", Male = "midnightblue")) {

  # Combine normative tables
  norm_combined.dt <- rbindlist(
    lapply(names(norm_tables_list), function(sex) {
      dt <- copy(norm_tables_list[[sex]])
      dt[, SEX := sex]
      dt
    }),
    fill = TRUE
  )

  # Combine observed data
  obs_combined.dt <- rbindlist(obs_data_list, fill = TRUE)

  # Get centile columns
  cent_cols <- grep("^p[0-9.]+$", names(norm_combined.dt), value = TRUE)
  cent_cols <- cent_cols[cent_cols %in% c("p2.5", "p10", "p25", "p50", "p75", "p90", "p97.5")]

  # Melt for plotting
  cent_long.dt <- data.table::melt(
    norm_combined.dt,
    id.vars = c("AGE", "SEX"),
    measure.vars = cent_cols,
    variable.name = "Centile",
    value.name = "Value"
  )

  # Define line types
  linetypes <- c(
    "p2.5" = "dotted", "p10" = "dashed", "p25" = "dotdash",
    "p50" = "solid", "p75" = "dotdash", "p90" = "dashed",
    "p97.5" = "dotted"
  )

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = obs_combined.dt,
      ggplot2::aes(x = AGE, y = VAL, color = SEX),
      alpha = 0.03, size = 0.3
    ) +
    ggplot2::geom_line(
      data = cent_long.dt,
      ggplot2::aes(
        x = AGE, y = Value, color = SEX,
        linetype = Centile, group = interaction(SEX, Centile)
      ),
      linewidth = 0.7
    ) +
    ggplot2::facet_wrap(~ SEX, ncol = 2) +
    ggplot2::scale_color_manual(values = sex_colors, guide = "none") +
    ggplot2::scale_linetype_manual(
      values = linetypes,
      labels = function(x) gsub("p", "", x),
      name = "Centile"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.2, 1.0),
      breaks = seq(0.2, 1.0, by = 0.1)
    ) +
    ggplot2::labs(
      title = "Hippocampal-Ventricle Ratio: Normative Centiles by Sex",
      subtitle = "HVR = HC / (HC + LV); higher values indicate better structural integrity",
      x = "Age (years)",
      y = "HVR"
    ) +
    theme_publication(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 11, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      legend.position = "bottom"
    )
}

#' Create HVR vs HC_RES scatter plot with marginal distributions
#'
#' Shows the relationship between HVR and residualized hippocampal volume,
#' demonstrating that HVR captures unique information beyond HC_RES.
#'
#' @param data Data.table with HC_RES (or HC_z), HVR (or HVR_z), and SEX columns
#' @param hc_col Name of the HC column (default "HC_z")
#' @param hvr_col Name of the HVR column (default "HVR_z")
#' @param sex_colors Named vector of sex colors
#' @param show_marginals Whether to include marginal density plots
#' @param annotate_r Whether to add correlation annotation
#' @return ggplot object (or patchwork object if show_marginals = TRUE)
#' @export
plot_hvr_vs_hcres <- function(
    data,
    hc_col = "HC_z",
    hvr_col = "HVR_z",
    sex_colors = c(Female = "darkred", Male = "midnightblue"),
    show_marginals = TRUE,
    annotate_r = TRUE) {

  # Calculate correlations
  r_overall <- cor(data[[hc_col]], data[[hvr_col]], use = "complete.obs")
  r_female <- cor(data[SEX == "Female"][[hc_col]], data[SEX == "Female"][[hvr_col]], use = "complete.obs")
  r_male <- cor(data[SEX == "Male"][[hc_col]], data[SEX == "Male"][[hvr_col]], use = "complete.obs")

  # Main scatter plot
  p_main <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[hc_col]], y = .data[[hvr_col]], color = SEX)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.5) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.2) +
    ggplot2::scale_color_manual(values = sex_colors, name = "Sex") +
    ggplot2::labs(
      x = "Residualized Hippocampal Volume (z-score)",
      y = "Hippocampal-Ventricle Ratio (z-score)"
    ) +
    theme_publication(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  # Add correlation annotation
  if (annotate_r) {
    p_main <- p_main +
      ggplot2::annotate(
        "text", x = Inf, y = Inf,
        label = sprintf("r = %.2f (overall)\nFemale: %.2f | Male: %.2f",
                       r_overall, r_female, r_male),
        hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "bold"
      )
  }

  if (!show_marginals) {
    return(p_main +
      ggplot2::labs(
        title = "HVR vs Residualized HC: Unique Information",
        subtitle = "Moderate correlation indicates HVR captures distinct variance"
      ))
  }

  # Marginal density for HC (top)
  p_top <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[hc_col]], fill = SEX, color = SEX)) +
    ggplot2::geom_density(alpha = 0.3, linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = sex_colors) +
    ggplot2::scale_color_manual(values = sex_colors) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  # Marginal density for HVR (right)
  p_right <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[hvr_col]], fill = SEX, color = SEX)) +
    ggplot2::geom_density(alpha = 0.3, linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = sex_colors) +
    ggplot2::scale_color_manual(values = sex_colors) +
    ggplot2::coord_flip() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  # Combine using patchwork
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- (p_top + patchwork::plot_spacer() + p_main + p_right) +
      patchwork::plot_layout(
        ncol = 2, nrow = 2,
        widths = c(4, 1),
        heights = c(1, 4)
      ) +
      patchwork::plot_annotation(
        title = "HVR vs Residualized HC: Unique Information",
        subtitle = "Moderate correlation indicates HVR captures distinct variance beyond HC_RES",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 10)
        )
      )
    return(combined)
  } else {
    # Fallback without patchwork
    return(p_main +
      ggplot2::labs(
        title = "HVR vs Residualized HC: Unique Information",
        subtitle = "Moderate correlation indicates HVR captures distinct variance"
      ))
  }
}

#' Create incremental R² visualization
#'
#' Bar plot showing R² for HC_RES only vs HC_RES + HVR models,
#' highlighting the incremental variance explained by HVR.
#'
#' @param r2_base R² from base model (HC_RES + covariates)
#' @param r2_full R² from full model (HC_RES + HVR + covariates)
#' @param title Plot title
#' @return ggplot object
#' @export
plot_incremental_r2 <- function(
    r2_base,
    r2_full,
    title = "HVR Adds Unique Variance Beyond HC_RES") {

  delta_r2 <- r2_full - r2_base

  plot_data <- data.frame(
    Model = factor(c("HC_RES only", "HC_RES + HVR"),
                  levels = c("HC_RES only", "HC_RES + HVR")),
    R2 = c(r2_base, r2_full),
    Component = c("HC_RES", "Total")
  )

  # Create stacked bar showing increment
  increment_data <- data.frame(
    Model = factor("HC_RES + HVR", levels = c("HC_RES only", "HC_RES + HVR")),
    ymin = r2_base,
    ymax = r2_full
  )

  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = plot_data,
      ggplot2::aes(x = Model, y = R2),
      fill = c("#0072B2", "#009E73"),
      width = 0.6
    ) +
    ggplot2::geom_segment(
      data = increment_data,
      ggplot2::aes(x = 1.4, xend = 1.6, y = ymin, yend = ymin),
      linetype = "dashed", color = "grey50"
    ) +
    ggplot2::geom_segment(
      data = increment_data,
      ggplot2::aes(x = 1.5, xend = 1.5, y = ymin, yend = ymax),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"), type = "closed"),
      color = "#D55E00", linewidth = 1
    ) +
    ggplot2::annotate(
      "text", x = 1.7, y = (r2_base + r2_full) / 2,
      label = sprintf("ΔR² = %.4f\n(+%.1f%%)", delta_r2, 100 * delta_r2 / r2_base),
      hjust = 0, size = 3.5, fontface = "bold", color = "#D55E00"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.01),
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Base R² = %.4f | Full R² = %.4f", r2_base, r2_full),
      x = "Model",
      y = "Variance Explained (R²)"
    ) +
    theme_publication(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      panel.grid.major.x = ggplot2::element_blank()
    )
}

#' Create distribution overlap plot showing sex differences
#'
#' Visualizes male and female distributions with overlap region highlighted.
#' Annotates with Cohen's d effect size and percent overlap.
#'
#' @param data Data.table with VAL (value) and SEX columns
#' @param effect_size Cohen's d effect size (positive = males larger)
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param x_label X-axis label
#' @param sex_colors Named vector of sex colors
#' @param show_overlap Whether to shade the overlap region
#' @return ggplot object
#' @export
plot_distribution_overlap <- function(
    data,
    effect_size = NULL,
    title = "Distribution by Sex",
    subtitle = NULL,
    x_label = "Value",
    sex_colors = c(Female = "darkred", Male = "midnightblue"),
    show_overlap = TRUE) {

  # Calculate effect size if not provided
  if (is.null(effect_size)) {
    effect_size <- effsize::cohen.d(data$VAL ~ data$SEX, na.rm = TRUE)$estimate
  }

  # Calculate percent overlap (Cohen's U3 complement)
  # For normal distributions: overlap = 2 * pnorm(-abs(d)/2)
  overlap_pct <- 2 * pnorm(-abs(effect_size) / 2) * 100

  # Get density data for each sex
  dens_f <- density(data[SEX == "Female", VAL], na.rm = TRUE)
  dens_m <- density(data[SEX == "Male", VAL], na.rm = TRUE)

  # Create data frames for plotting
  dens_f_df <- data.frame(x = dens_f$x, y = dens_f$y, SEX = "Female")
  dens_m_df <- data.frame(x = dens_m$x, y = dens_m$y, SEX = "Male")
  dens_combined <- rbind(dens_f_df, dens_m_df)

  # Create base plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_density(
      data = data,
      ggplot2::aes(x = VAL, fill = SEX, color = SEX),
      alpha = 0.3,
      linewidth = 0.8
    ) +
    ggplot2::scale_fill_manual(values = sex_colors, name = "Sex") +
    ggplot2::scale_color_manual(values = sex_colors, name = "Sex")

  # Add effect size annotation
  if (is.null(subtitle)) {
    subtitle <- sprintf(
      "Cohen's d = %.2f | Distribution overlap = %.1f%%",
      effect_size, overlap_pct
    )
  }

  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = "Density"
    ) +
    theme_publication(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      legend.position = "bottom"
    )

  p
}

#' Create faceted distribution overlap plot by adjustment method
#'
#' Shows how sex differences in distributions change across head-size
#' adjustment methods. Demonstrates that HC reverses while HVR is stable.
#'
#' @param data Data.table with VAL, SEX, ROI, and ADJ columns
#' @param effect_sizes Data.table with ROI, ADJ, ESTIMATE columns for Cohen's d
#' @param title Plot title
#' @param sex_colors Named vector of sex colors
#' @param roi_order Order of ROIs for faceting
#' @param adj_order Order of adjustments for faceting
#' @return ggplot object
#' @export
plot_distribution_overlap_faceted <- function(
    data,
    effect_sizes,
    title = "Sex Differences Across Adjustment Methods",
    sex_colors = c(Female = "darkred", Male = "midnightblue"),
    roi_order = c("HC", "LV", "HVR"),
    adj_order = c("NON", "PRP", "STX", "RES")) {

  # Set factor levels
  data <- copy(data)
  data[, ROI := factor(ROI, levels = roi_order)]
  data[, ADJ := factor(ADJ, levels = adj_order)]

  # Merge effect sizes for annotation
  effect_sizes <- copy(effect_sizes)
  effect_sizes[, d_label := sprintf("d = %.2f", ESTIMATE)]

  # Create base plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = VAL, fill = SEX, color = SEX)) +
    ggplot2::geom_density(alpha = 0.3, linewidth = 0.6) +
    ggplot2::facet_grid(ROI ~ ADJ, scales = "free") +
    ggplot2::scale_fill_manual(values = sex_colors, name = "Sex") +
    ggplot2::scale_color_manual(values = sex_colors, name = "Sex")

  # Add effect size labels
  p <- p +
    ggplot2::geom_text(
      data = effect_sizes,
      ggplot2::aes(x = Inf, y = Inf, label = d_label),
      hjust = 1.1, vjust = 1.5,
      inherit.aes = FALSE,
      size = 3, fontface = "bold"
    )

  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = "Positive d = males larger | Rows: Brain measure | Columns: Adjustment method",
      x = "Value",
      y = "Density"
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      legend.position = "bottom",
      axis.text = ggplot2::element_text(size = 8)
    )

  p
}

#' Create HVR-focused distribution comparison plot
#'
#' Highlights that HVR shows consistent female advantage regardless of
#' whether HC shows male or female advantage depending on adjustment.
#'
#' @param data Data.table with VAL, SEX, ROI, ADJ columns
#' @param effect_sizes Data.table with effect sizes for annotation
#' @param sex_colors Named vector of sex colors
#' @return ggplot object
#' @export
plot_hvr_vs_hc_distributions <- function(
    data,
    effect_sizes,
    sex_colors = c(Female = "darkred", Male = "midnightblue")) {

  # Create comparison data: HC with different adjustments + HVR (NON only)
  hc_data <- data[ROI == "HC"]
  hvr_data <- data[ROI == "HVR" & ADJ == "NON"]
  hvr_data[, ADJ := "HVR"]  # Rename for clarity

  # Combine
  plot_data <- rbind(
    hc_data[, .(VAL, SEX, METHOD = paste0("HC (", ADJ, ")"))],
    hvr_data[, .(VAL, SEX, METHOD = "HVR")]
  )

  # Set method order
  method_order <- c("HC (NON)", "HC (PRP)", "HC (STX)", "HC (RES)", "HVR")
  plot_data[, METHOD := factor(METHOD, levels = method_order)]

  # Merge effect sizes
  hc_effects <- effect_sizes[ROI == "HC", .(METHOD = paste0("HC (", ADJ, ")"), ESTIMATE)]
  hvr_effects <- effect_sizes[ROI == "HVR" & ADJ == "NON", .(METHOD = "HVR", ESTIMATE)]
  all_effects <- rbind(hc_effects, hvr_effects)
  all_effects[, d_label := sprintf("d = %.2f", ESTIMATE)]
  all_effects[, METHOD := factor(METHOD, levels = method_order)]

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = VAL, fill = SEX, color = SEX)) +
    ggplot2::geom_density(alpha = 0.35, linewidth = 0.7) +
    ggplot2::facet_wrap(~ METHOD, scales = "free_x", ncol = 5) +
    ggplot2::scale_fill_manual(values = sex_colors, name = "Sex") +
    ggplot2::scale_color_manual(values = sex_colors, name = "Sex") +
    ggplot2::geom_text(
      data = all_effects,
      ggplot2::aes(x = Inf, y = Inf, label = d_label),
      hjust = 1.1, vjust = 1.5,
      inherit.aes = FALSE,
      size = 3.5, fontface = "bold"
    ) +
    ggplot2::labs(
      title = "Head-Size Adjustment Reverses HC Sex Differences; HVR is Stable",
      subtitle = "HC shows male advantage (d > 0) when unadjusted but female advantage (d < 0) when residualized; HVR consistently shows female advantage",
      x = "Volume / Ratio",
      y = "Density"
    ) +
    theme_publication(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      legend.position = "bottom"
    )

  p
}
