# =============================================================================
# Manuscript Export Utilities
# =============================================================================
# Functions for exporting publication-ready figures and tables from pipeline
# scripts using config-defined naming. These functions bridge the pipeline
# analysis scripts with the Quarto document figure/table generation functions
# in reports-src/_common.R.
# =============================================================================

library(here)

#' Check if _common.R is loaded
#'
#' @return TRUE if _common.R functions are available
.common_loaded <- function() {
  exists("fig_distribution_overlap", mode = "function")
}

#' Load _common.R if not already loaded
#'
#' @return TRUE if successful
.ensure_common_loaded <- function() {
  if (!.common_loaded()) {
    common_path <- here("reports-src/_common.R")
    if (!file.exists(common_path)) {
      warning("_common.R not found at: ", common_path, call. = FALSE)
      return(FALSE)
    }
    source(common_path, local = FALSE)
  }
  TRUE
}

#' Export manuscript figures for a script
#'
#' Generates all figures defined in config for the specified script.
#' Uses the fig_* functions from _common.R.
#'
#' @param script_name Name of the script (e.g., "sex_differences", "gamlss")
#' @param data_env Environment or list containing required data objects
#' @param dpi Resolution for PNG output (default from config)
#' @return Invisible list of generated figure paths
#' @export
export_manuscript_figures <- function(script_name, data_env = NULL, dpi = NULL) {
  if (!.ensure_common_loaded()) {
    warning("Cannot export figures: _common.R not available", call. = FALSE)
    return(invisible(list()))
  }

  # Get figures for this script
  fig_names <- get_script_figures(script_name)
  if (length(fig_names) == 0) {
    message("No figures defined for script: ", script_name)
    return(invisible(list()))
  }

  # Get output settings
  if (is.null(dpi)) {
    dpi <- get_config("output", "plots_dpi", default = 300)
  }
  fig_dir <- get_output_path("figures")

  # Figure function mapping - maps config names to _common.R function names
  fig_functions <- list(
    fig1_distribution_overlap = "fig_distribution_overlap",
    fig2_effect_sizes = "fig_effect_sizes_by_adjustment",
    fig3_age_trajectory = "fig_age_trajectory",
    fig4_hvr_centiles = "fig_hvr_centiles",
    fig5_sem_forest = "fig_sem_forest",
    sfig_hc_centiles_female = function(env) fig_zscore_trajectories(env$norm_tables, sex = "Female", roi = "HC"),
    sfig_hc_centiles_male = function(env) fig_zscore_trajectories(env$norm_tables, sex = "Male", roi = "HC"),
    sfig_lv_centiles = "fig_lv_centiles_sexcomp",
    sfig_calibration_hvr = function(env) fig_gamlss_validation_overall(env$gamlss_calib, "HVR"),
    sfig_calibration_hc = function(env) fig_gamlss_validation_overall(env$gamlss_calib, "HC"),
    sfig_calibration_lv = function(env) fig_gamlss_validation_overall(env$gamlss_calib, "LV"),
    sfig_hemisphere_comparison = "fig_hemisphere_comparison",
    sfig_sensitivity_comparison = "fig_sensitivity_comparison"
  )

  # Data argument mapping (which data object each named function needs)
  fig_data_args <- list(
    fig1_distribution_overlap = "brain_data",
    fig2_effect_sizes = "sex_diff",
    fig3_age_trajectory = "sex_diff",
    fig4_hvr_centiles = c("norm_tables", "brain_data"),
    fig5_sem_forest = "sem_params",
    sfig_lv_centiles = "norm_tables",
    sfig_hemisphere_comparison = "sex_diff",
    sfig_sensitivity_comparison = "sex_diff"
  )

  generated <- list()

  for (fig_name in fig_names) {
    fig_config <- get_figure_config(fig_name)
    output_path <- file.path(fig_dir, paste0(fig_name, ".png"))

    tryCatch({
      # Get the function
      func_ref <- fig_functions[[fig_name]]

      if (is.null(func_ref)) {
        warning("No function mapping for figure: ", fig_name, call. = FALSE)
        next
      }

      # Generate the plot
      if (is.function(func_ref)) {
        # Custom function with data_env parameter
        p <- func_ref(data_env)
      } else {
        # Named function - get data arguments
        func <- get(func_ref, mode = "function")
        data_args <- fig_data_args[[fig_name]]

        if (is.null(data_args) || is.null(data_env)) {
          p <- func()
        } else if (length(data_args) == 1) {
          p <- func(data_env[[data_args]])
        } else {
          # Multiple arguments
          args <- lapply(data_args, function(x) data_env[[x]])
          p <- do.call(func, args)
        }
      }

      # Save
      ggplot2::ggsave(
        output_path, p,
        width = fig_config$width,
        height = fig_config$height,
        dpi = dpi
      )

      generated[[fig_name]] <- output_path
      message("  Generated: ", fig_name)

    }, error = function(e) {
      warning("Failed to generate ", fig_name, ": ", e$message, call. = FALSE)
    })
  }

  invisible(generated)
}

#' Export manuscript tables for a script
#'
#' Saves gt tables defined in config for the specified script.
#' Tables are saved in both HTML and TEX formats.
#'
#' @param script_name Name of the script (e.g., "sex_differences", "sem")
#' @param tables_list Named list of gt table objects to export
#' @return Invisible list of generated table paths
#' @export
export_manuscript_tables <- function(script_name, tables_list) {
  # Get tables for this script
  tbl_names <- get_script_tables(script_name)
  if (length(tbl_names) == 0) {
    message("No tables defined for script: ", script_name)
    return(invisible(list()))
  }

  tbl_dir <- get_output_path("tables")
  generated <- list()

  for (tbl_name in tbl_names) {
    if (!tbl_name %in% names(tables_list)) {
      warning("Table not provided: ", tbl_name, call. = FALSE)
      next
    }

    gt_table <- tables_list[[tbl_name]]

    tryCatch({
      html_path <- file.path(tbl_dir, paste0(tbl_name, ".html"))
      tex_path <- file.path(tbl_dir, paste0(tbl_name, ".tex"))

      gt::gtsave(gt_table, html_path)
      gt::gtsave(gt_table, tex_path)

      generated[[tbl_name]] <- list(html = html_path, tex = tex_path)
      message("  Generated: ", tbl_name)

    }, error = function(e) {
      warning("Failed to generate ", tbl_name, ": ", e$message, call. = FALSE)
    })
  }

  invisible(generated)
}

#' Export all manuscript assets for a script
#'
#' Convenience function that exports both figures and tables for a script.
#'
#' @param script_name Name of the script
#' @param data_env Environment or list containing data for figures
#' @param tables_list Named list of gt table objects
#' @return Invisible list with figures and tables paths
#' @export
export_manuscript_assets <- function(script_name, data_env = NULL, tables_list = NULL) {
  message("Exporting manuscript assets for: ", script_name)

  result <- list(figures = list(), tables = list())

  if (!is.null(data_env)) {
    result$figures <- export_manuscript_figures(script_name, data_env)
  }

  if (!is.null(tables_list)) {
    result$tables <- export_manuscript_tables(script_name, tables_list)
  }

  invisible(result)
}
