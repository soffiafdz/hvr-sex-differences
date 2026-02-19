# =============================================================================
# Configuration Management
# =============================================================================
# Functions for loading and accessing pipeline configuration
# =============================================================================

library(yaml)
library(here)

.config <- NULL

#' Load pipeline configuration
#'
#' @param config_file Path to YAML configuration file
#' @return List containing configuration
#' @export
load_config <- function(config_file = here("config/pipeline_config.yaml")) {
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file, call. = FALSE)
  }

  .config <<- read_yaml(config_file)
  invisible(.config)
}

#' Get configuration value
#'
#' @param ... Path to configuration value (e.g., "general", "random_seed")
#' @param default Default value if not found
#' @return Configuration value
#' @export
get_config <- function(..., default = NULL) {
  if (is.null(.config)) {
    load_config()
  }

  path <- list(...)
  value <- .config

  for (key in path) {
    if (is.null(value[[key]])) {
      if (!is.null(default)) {
        return(default)
      }
      stop("Configuration key not found: ", paste(path, collapse = " -> "),
           call. = FALSE)
    }
    value <- value[[key]]
  }

  return(value)
}

#' Get data path from configuration
#'
#' @param ... Path keys (e.g., "processed", "covars_fst")
#' @return Full path to data file
#' @export
get_data_path <- function(...) {
  rel_path <- get_config("data", ...)
  here(rel_path)
}

#' Get script setting
#'
#' @param script_name Name of script (e.g., "gamlss")
#' @param setting Setting name (e.g., "redo_plots")
#' @param default Default value
#' @export
get_script_setting <- function(script_name, setting, default = NULL) {
  get_config("scripts", script_name, setting, default = default)
}

#' Get parameter value
#'
#' @param ... Path to parameter
#' @param default Default value
#' @export
get_parameter <- function(..., default = NULL) {
  get_config("parameters", ..., default = default)
}

#' Get number of cores for parallel processing
#'
#' @export
get_n_cores <- function() {
  n_cores <- get_config("general", "n_cores", default = 1)
  max_cores <- parallel::detectCores() - 1
  min(n_cores, max_cores)
}

#' Get random seed
#'
#' @export
get_seed <- function() {
  get_config("general", "random_seed", default = 1618)
}

#' Set random seed from configuration
#'
#' @export
set_seed <- function() {

  seed <- get_seed()
  set.seed(seed)
  invisible(seed)
}

#' Validate configuration has required fields
#'
#' Checks that all required configuration sections and keys exist.
#' Call after load_config() to ensure config is complete.
#'
#' @param config Configuration list (default: loaded config)
#' @return TRUE if valid, throws error otherwise
#' @export
validate_config <- function(config = NULL) {
  if (is.null(config)) {
    if (is.null(.config)) {
      load_config()
    }
    config <- .config
  }

  # Required top-level sections

required_sections <- c("general", "data", "parameters", "scripts")
  missing_sections <- setdiff(required_sections, names(config))
  if (length(missing_sections) > 0) {
    stop("Missing required config sections: ",
         paste(missing_sections, collapse = ", "), call. = FALSE)
  }

  # Required general settings
  required_general <- c("random_seed", "n_cores")
  missing_general <- setdiff(required_general, names(config$general))
  if (length(missing_general) > 0) {
    stop("Missing required general settings: ",
         paste(missing_general, collapse = ", "), call. = FALSE)
  }

  # Required data paths
  required_data <- c("raw", "processed", "models")
  missing_data <- setdiff(required_data, names(config$data))
  if (length(missing_data) > 0) {
    stop("Missing required data sections: ",
         paste(missing_data, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)
}

#' Get education years mapping
#'
#' Returns the education level to years mapping from config.
#' Centralizes access to avoid duplication across scripts.
#'
#' @return Named numeric vector mapping education levels to years
#' @export
get_education_years_map <- function() {
  edu_map <- get_parameter("education")
  if (is.null(edu_map)) {
    stop("Education mapping not found in config parameters.education", call. = FALSE)
  }
  sapply(edu_map, as.numeric)
}

#' Get output path
#'
#' Returns the path for output files (figures, tables, documents).
#'
#' @param type Type of output: "figures", "tables", or "documents"
#' @param filename Optional filename to append
#' @return Full path to output directory or file
#' @export
get_output_path <- function(type = c("figures", "tables", "documents"), filename = NULL) {
  type <- match.arg(type)
  base_path <- get_config("output", type, default = file.path("outputs", type))
  full_path <- here(base_path)

  # Create directory if it doesn't exist
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE, showWarnings = FALSE)
  }

  if (!is.null(filename)) {
    file.path(full_path, filename)
  } else {
    full_path
  }
}

#' Get manuscript figure configuration
#'
#' Returns the configuration for a manuscript figure (dimensions, script source).
#'
#' @param fig_name Name of the figure (e.g., "fig1_distribution_overlap")
#' @return List with width, height, script, or NULL if not found
#' @export
get_figure_config <- function(fig_name) {
  fig_config <- get_config("manuscript", "figures", fig_name, default = NULL)
  if (is.null(fig_config)) {
    warning("Figure not defined in config: ", fig_name, call. = FALSE)
    return(list(width = 8, height = 6))  # defaults

  }
  fig_config
}

#' Get manuscript table configuration
#'
#' Returns the configuration for a manuscript table.
#'
#' @param tbl_name Name of the table (e.g., "table1_demographics")
#' @return List with script source, or NULL if not found
#' @export
get_table_config <- function(tbl_name) {
  get_config("manuscript", "tables", tbl_name, default = NULL)
}

#' Get all manuscript figures for a script
#'
#' Returns the names of all figures that should be generated by a given script.
#'
#' @param script_name Name of the script (e.g., "sex_differences", "gamlss")
#' @return Character vector of figure names
#' @export
get_script_figures <- function(script_name) {
  all_figs <- get_config("manuscript", "figures", default = list())
  fig_names <- names(all_figs)

  # Filter to figures for this script
  script_figs <- sapply(fig_names, function(fn) {
    all_figs[[fn]]$script == script_name
  })

  fig_names[script_figs]
}

#' Get all manuscript tables for a script
#'
#' Returns the names of all tables that should be generated by a given script.
#'
#' @param script_name Name of the script (e.g., "sex_differences", "sem")
#' @return Character vector of table names
#' @export
get_script_tables <- function(script_name) {
  all_tbls <- get_config("manuscript", "tables", default = list())
  tbl_names <- names(all_tbls)

  # Filter to tables for this script
  script_tbls <- sapply(tbl_names, function(tn) {
    all_tbls[[tn]]$script == script_name
  })

  tbl_names[script_tbls]
}

#' Get standard factor levels for categorical variables
#'
#' Returns the expected factor levels for common categorical variables
#' to ensure consistency across the pipeline.
#'
#' @param variable_name Name of the variable (e.g., "SEX", "ADJ_METHOD")
#' @return Character vector of expected levels in order
#' @export
get_factor_levels <- function(variable_name) {
  # Standard factor levels used throughout the pipeline
  factor_levels <- list(
    SEX = c("Female", "Male"),
    ADJ_METHOD = c("NON", "PRP", "STX", "RES"),
    ADJ_LABEL = c("Unadjusted", "Proportions", "Stereotaxic", "Residuals"),
    SIDE = c("L", "R", "LR"),
    SEGMENTATION = c("CRS", "LNG"),
    COHORT = c("ALL", "MTCH", "SENS")
  )

  if (!variable_name %in% names(factor_levels)) {
    warning("Unknown variable name: ", variable_name,
            ". No standard factor levels defined.", call. = FALSE)
    return(NULL)
  }

  factor_levels[[variable_name]]
}
