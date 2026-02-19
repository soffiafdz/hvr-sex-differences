# =============================================================================
# Data I/O Utilities
# =============================================================================
# Functions for reading and writing data with consistent error handling
# =============================================================================

library(here)
library(data.table)

source(here("R/utils/logging.R"))

#' Check if required files exist
#'
#' @param files Character vector of file paths
#' @param stop_on_missing Throw error if files missing
#' @return Logical indicating all files exist
#' @export
check_files_exist <- function(files, stop_on_missing = TRUE) {
  missing <- files[!file.exists(files)]

  if (length(missing) > 0) {
    msg <- sprintf(
      "Missing required file(s):\n%s",
      paste("  -", missing, collapse = "\n")
    )

    if (stop_on_missing) {
      log_error(msg)
      stop(msg, call. = FALSE)
    } else {
      log_warn(msg)
      return(FALSE)
    }
  }

  log_debug("All required files exist (%d files)", length(files))
  return(TRUE)
}

#' Read RDS file with error handling
#'
#' @param file_path Path to RDS file
#' @param description Optional description for logging
#' @return Data object from RDS file
#' @export
read_rds_safe <- function(file_path, description = NULL) {
  if (!file.exists(file_path)) {
    log_error("RDS file not found: %s", file_path)
    stop("File not found: ", file_path, call. = FALSE)
  }

  desc <- if (!is.null(description)) description else basename(file_path)
  log_debug("Reading RDS: %s", desc)

  tryCatch(
    readRDS(file_path),
    error = function(e) {
      log_error("Failed to read RDS file %s: %s", file_path, e$message)
      stop(e)
    }
  )
}

#' Write RDS file with error handling
#'
#' @param object Object to save
#' @param file_path Path to save to
#' @param description Optional description for logging
#' @export
write_rds_safe <- function(object, file_path, description = NULL) {
  # Create directory if needed
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    log_debug("Creating directory: %s", dir_path)
    dir.create(dir_path, recursive = TRUE)
  }

  desc <- if (!is.null(description)) description else basename(file_path)
  log_debug("Writing RDS: %s", desc)

  tryCatch(
    saveRDS(object, file_path),
    error = function(e) {
      log_error("Failed to write RDS file %s: %s", file_path, e$message)
      stop(e)
    }
  )

  log_info("Saved: %s", file_path)
  invisible(file_path)
}

#' Read FST file with error handling
#'
#' @param file_path Path to FST file
#' @param as_data_table Return as data.table
#' @param description Optional description for logging
#' @return Data from FST file
#' @export
read_fst_safe <- function(file_path, as_data_table = TRUE, description = NULL) {
  if (!requireNamespace("fst", quietly = TRUE)) {
    stop("Package 'fst' is required", call. = FALSE)
  }

  if (!file.exists(file_path)) {
    log_error("FST file not found: %s", file_path)
    stop("File not found: ", file_path, call. = FALSE)
  }

  desc <- if (!is.null(description)) description else basename(file_path)
  log_debug("Reading FST: %s", desc)

  tryCatch(
    fst::read_fst(file_path, as.data.table = as_data_table),
    error = function(e) {
      log_error("Failed to read FST file %s: %s", file_path, e$message)
      stop(e)
    }
  )
}

#' Write FST file with error handling
#'
#' @param x Data to save
#' @param file_path Path to save to
#' @param compress Compression level (0-100)
#' @param description Optional description for logging
#' @export
write_fst_safe <- function(x, file_path, compress = 85, description = NULL) {
  if (!requireNamespace("fst", quietly = TRUE)) {
    stop("Package 'fst' is required", call. = FALSE)
  }

  # Create directory if needed
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    log_debug("Creating directory: %s", dir_path)
    dir.create(dir_path, recursive = TRUE)
  }

  desc <- if (!is.null(description)) description else basename(file_path)
  log_debug("Writing FST: %s", desc)

  tryCatch(
    fst::write_fst(x, file_path, compress = compress),
    error = function(e) {
      log_error("Failed to write FST file %s: %s", file_path, e$message)
      stop(e)
    }
  )

  log_info("Saved: %s", file_path)
  invisible(file_path)
}

#' Read CSV file with error handling
#'
#' @param file_path Path to CSV file
#' @param ... Additional arguments for fread
#' @param description Optional description for logging
#' @export
read_csv_safe <- function(file_path, ..., description = NULL) {
  if (!file.exists(file_path)) {
    log_error("CSV file not found: %s", file_path)
    stop("File not found: ", file_path, call. = FALSE)
  }

  desc <- if (!is.null(description)) description else basename(file_path)
  log_debug("Reading CSV: %s", desc)

  tryCatch(
    fread(file_path, ...),
    error = function(e) {
      log_error("Failed to read CSV file %s: %s", file_path, e$message)
      stop(e)
    }
  )
}

#' Write CSV file with error handling
#'
#' @param x Data to save
#' @param file_path Path to save to
#' @param ... Additional arguments for fwrite
#' @param description Optional description for logging
#' @export
write_csv_safe <- function(x, file_path, ..., description = NULL) {
  # Create directory if needed
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    log_debug("Creating directory: %s", dir_path)
    dir.create(dir_path, recursive = TRUE)
  }

  desc <- if (!is.null(description)) description else basename(file_path)
  log_debug("Writing CSV: %s", desc)

  tryCatch(
    fwrite(x, file_path, ...),
    error = function(e) {
      log_error("Failed to write CSV file %s: %s", file_path, e$message)
      stop(e)
    }
  )

  log_info("Saved: %s", file_path)
  invisible(file_path)
}

#' Check if output file needs to be regenerated
#'
#' @param output_path Path to output file
#' @param input_paths Vector of input file paths
#' @param force_regenerate Force regeneration regardless of modification times
#' @return Logical indicating if regeneration is needed
#' @export
needs_regeneration <- function(output_path,
                               input_paths = NULL,
                               force_regenerate = FALSE) {
  if (force_regenerate) {
    log_debug("Force regeneration: %s", basename(output_path))
    return(TRUE)
  }

  if (!file.exists(output_path)) {
    log_debug("Output does not exist: %s", basename(output_path))
    return(TRUE)
  }

  if (is.null(input_paths)) {
    log_debug("Output exists, no input tracking: %s", basename(output_path))
    return(FALSE)
  }

  # Check modification times
  output_mtime <- file.mtime(output_path)
  input_mtimes <- sapply(input_paths[file.exists(input_paths)], file.mtime)

  if (any(input_mtimes > output_mtime)) {
    log_debug("Inputs newer than output: %s", basename(output_path))
    return(TRUE)
  }

  log_debug("Output up-to-date: %s", basename(output_path))
  return(FALSE)
}

#' Create directory if it doesn't exist
#'
#' @param dir_path Directory path
#' @export
ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    log_debug("Creating directory: %s", dir_path)
    dir.create(dir_path, recursive = TRUE)
  }
  invisible(dir_path)
}

# =============================================================================
# Analysis Data Loaders
# =============================================================================
# High-level functions for loading analysis outputs with validation

#' Load primary analysis data with validation
#' Returns a list with all commonly needed datasets
#' @return List with demog_data, sex_diff, sem_fit, sem_params, cfa_minv, gamlss
#' @export
load_analysis_data <- function() {
  # Define required files
  required_files <- list(
    demog = get_data_path("demographics", "data"),
    sex_diff = get_data_path("processed", "sex_differences"),
    sem_fit = get_data_path("models", "results", "sem_fit_measures"),
    sem_params = get_data_path("models", "results", "sem_params"),
    cfa_minv = get_data_path("models", "results", "cfa-minv_cog")
  )

  # Optional files (GAMLSS)
  optional_files <- list(
    gamlss = get_data_path("models", "fit", "gamlss")
  )

  # Check required files exist
  missing <- sapply(required_files, function(f) !file.exists(f))
  if (any(missing)) {
    stop("Required data files not found:\n  ",
         paste(names(required_files)[missing], collapse = "\n  "),
         call. = FALSE)
  }

  # Load required data
  data <- list(
    demog_data = read_rds_safe(required_files$demog),
    sex_diff = read_rds_safe(required_files$sex_diff),
    sem_fit = read_rds_safe(required_files$sem_fit),
    sem_params = read_rds_safe(required_files$sem_params),
    cfa_minv = read_rds_safe(required_files$cfa_minv)
  )

  # Load optional data (GAMLSS summaries for normative modeling tables)
  if (file.exists(optional_files$gamlss)) {
    data$gamlss <- read_rds_safe(optional_files$gamlss)
  }

  # Validate critical nested structures
  if (!"CRS" %in% names(data$demog_data)) {
    stop("demog_data missing expected 'CRS' element", call. = FALSE)
  }
  if (!"ALL" %in% names(data$demog_data$CRS)) {
    stop("demog_data$CRS missing expected 'ALL' element", call. = FALSE)
  }
  if (!"OVERALL" %in% names(data$sex_diff)) {
    stop("sex_diff missing expected 'OVERALL' element", call. = FALSE)
  }

  data
}

#' Extract sample sizes from demographics data
#' @param demog_data Demographics data list (from load_analysis_data)
#' @return List with n_total, n_female, n_male, n_matched_pairs, etc.
#' @export
extract_sample_sizes <- function(demog_data) {
  all_demog <- demog_data$CRS$ALL

  sizes <- list(
    n_total = nrow(all_demog),
    n_female = sum(all_demog$SEX == "Female"),
    n_male = sum(all_demog$SEX == "Male")
  )

  # Add matched sample if available
  if ("MTCH" %in% names(demog_data$CRS)) {
    matched_demog <- demog_data$CRS$MTCH
    sizes$n_matched_total <- nrow(matched_demog)
    sizes$n_matched_pairs <- nrow(matched_demog) / 2
  }

  # Add sensitivity sample if available
  if ("SENS" %in% names(demog_data$CRS)) {
    sens_demog <- demog_data$CRS$SENS
    sizes$n_sens <- nrow(sens_demog)
    sizes$n_fcodes <- sizes$n_sens - sizes$n_total
  }

  sizes
}

#' Load HVR vs HC comparison analysis results
#' @return List with correlation_comparison, age_sex_interactions, age_stratified,
#'         standardized_effects, memory_comparison
#' @export
load_hvr_comparison <- function() {
  hvr_path <- get_data_path("models", "results", "hvr_hc_comparison")

  if (!file.exists(hvr_path)) {
    warning("HVR comparison data not found. Run R/run_hvr_comparison.R first.",
            call. = FALSE)
    return(NULL)
  }

  data <- read_rds_safe(hvr_path)

  # Validate structure (including incremental R² and sex-stratified results)
  expected_elements <- c("correlation_comparison", "age_sex_interactions",
                         "age_stratified", "standardized_effects", "memory_comparison",
                         "incremental_r2", "incremental_r2_by_sex")
  missing <- setdiff(expected_elements, names(data))
  if (length(missing) > 0) {
    # Only warn for critical elements, not optional ones
    critical_missing <- intersect(missing, c("correlation_comparison", "age_sex_interactions"))
    if (length(critical_missing) > 0) {
      warning("HVR comparison data missing critical elements: ", paste(critical_missing, collapse = ", "),
              call. = FALSE)
    }
  }

  data
}

#' Load normative centile tables from GAMLSS analysis
#' @return data.table with centile values by age, sex, ROI, and adjustment method
#' @export
load_normative_tables <- function() {
  norm_path <- get_data_path("processed", "norm_tables")

  if (!file.exists(norm_path)) {
    warning("Normative tables not found at: ", norm_path,
            "\nRun R/08_normative_tables.R first.", call. = FALSE)
    return(NULL)
  }

  read_rds_safe(norm_path)
}

#' Alias for load_normative_tables
#' @return list with normative tables by sex/ROI/ADJ
#' @export
load_norm_tables <- function() {
  load_normative_tables()
}

#' Get centile table for specific ROI, adjustment method, and sex
#' @param norm_tables Normative tables from load_normative_tables()
#' @param roi ROI code (e.g., "HVR", "HC", "LV")
#' @param adj Adjustment method (e.g., "RES", "PRP", "STX", "NON")
#' @param side Laterality (e.g., "LR", "L", "R")
#' @param sex Sex ("Female" or "Male")
#' @return data.table with age and centile columns, or NULL if not found
#' @export
get_centile_table <- function(norm_tables, roi, adj, side, sex) {
  if (is.null(norm_tables)) return(NULL)

  # Normalize sex to title case
  sex_key <- tools::toTitleCase(tolower(sex))

  # Navigate nested structure: norm_tables[[sex]][[roi]][[adj]][[side]]
  tryCatch({
    result <- norm_tables[[sex_key]][[roi]][[adj]][[side]]
    if (is.null(result)) {
      warning("Centile table not found for: ", sex_key, "/", roi, "/", adj, "/", side,
              call. = FALSE)
    }
    result
  }, error = function(e) {
    warning("Error accessing centile table: ", e$message, call. = FALSE)
    NULL
  })
}

#' Load SEM fit indices for all models
#' @return data.table with fit indices for each model
#' @export
load_sem_fit <- function() {
  fit_path <- get_data_path("models", "results", "sem_fit_measures")

  if (!file.exists(fit_path)) {
    warning("SEM fit measures not found. Run R/11_sem_analysis.R first.",
            call. = FALSE)
    return(NULL)
  }

  fit_list <- read_rds_safe(fit_path)

  # Combine into single data.table
  fit_dt <- rbindlist(fit_list, idcol = "MODEL")

  # Select key columns
  key_cols <- c("MODEL", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr")
  available_cols <- intersect(key_cols, names(fit_dt))

  fit_dt[, ..available_cols]
}

#' Load CFA model fits and extract fit measures
#' @return data.table with fit indices for CFA models
#' @export
load_cfa_fit <- function() {
  cfa_path <- get_data_path("models", "fit", "cfa_cog")

  if (!file.exists(cfa_path)) {
    warning("CFA fits not found. Run R/06_cognitive_factors.R first.",
            call. = FALSE)
    return(NULL)
  }

  cfa_fits <- read_rds_safe(cfa_path)

  # Extract fit measures from each model
  fit_list <- lapply(names(cfa_fits), function(name) {
    fit <- cfa_fits[[name]]
    if (inherits(fit, "lavaan")) {
      fm <- lavaan::fitMeasures(fit, c("cfi.robust", "tli.robust",
                                        "rmsea.robust", "rmsea.ci.lower.robust",
                                        "rmsea.ci.upper.robust", "srmr"))
      data.table(
        MODEL = name,
        CFI = fm["cfi.robust"],
        TLI = fm["tli.robust"],
        RMSEA = fm["rmsea.robust"],
        RMSEA_CI_LOWER = fm["rmsea.ci.lower.robust"],
        RMSEA_CI_UPPER = fm["rmsea.ci.upper.robust"],
        SRMR = fm["srmr"]
      )
    } else {
      NULL
    }
  })

  rbindlist(fit_list[!sapply(fit_list, is.null)])
}

#' Load CFA measurement invariance results
#' @return data.table with invariance test results
#' @export
load_measurement_invariance <- function() {
  minv_path <- get_data_path("models", "results", "cfa-minv_cog")

  if (!file.exists(minv_path)) {
    warning("Measurement invariance results not found.", call. = FALSE)
    return(NULL)
  }

  read_rds_safe(minv_path)
}

#' Load cognitive factor test-retest reliability data
#' @return data.table with ICC and correlation for each factor
#' @export
load_cog_reliability <- function() {
  rel_path <- get_data_path("processed", "cog_reliability")

  if (!file.exists(rel_path)) {
    warning("Cognitive reliability data not found. Run R/06_cognitive_factors.R first.",
            call. = FALSE)
    return(NULL)
  }

  read_rds_safe(rel_path)
}

#' Load brain volume data for figures
#' @return data.table with brain volumes
#' @export
load_brain_volumes <- function() {
  hc_hvr_path <- get_data_path("processed", "hc_hvr_adjusted")
  if (!file.exists(hc_hvr_path)) {
    warning("Brain volume data not found")
    return(NULL)
  }
  hc_hvr <- read_rds_safe(hc_hvr_path)
  # Return cross-sectional ALL data from LPP_CNN network
  hc_hvr$LPP_CNN$CRS$ALL
}

#' Load matched sample brain volume data
#' @param demog_data Demographics data containing matched sample EIDs
#' @return data.table with brain volumes for matched sample only
#' @export
load_matched_brain_volumes <- function(demog_data) {
  brain_data <- load_brain_volumes()
  if (is.null(brain_data) || is.null(demog_data)) return(NULL)

  # Get matched sample EIDs
  mtch_demog <- demog_data$CRS$MTCH
  if (is.null(mtch_demog)) return(NULL)

  # Filter brain data to matched sample
  brain_data[EID %in% mtch_demog$EID]
}

#' Load GAMLSS model diagnostics
#' @return list with validation data
#' @export
load_gamlss_diagnostics <- function() {
  diag_path <- get_data_path("models", "diagnostics", "gamlss")
  if (!file.exists(diag_path)) {
    warning("GAMLSS diagnostics not found")
    return(NULL)
  }
  read_rds_safe(diag_path)
}

#' Load GAMLSS calibration data
#' @return List with calibration data by ROI and sex
#' @export
load_gamlss_calibration <- function() {
  # Try site model first
  gamlss_path <- get_data_path("models", "fit", "gamlss_site")
  if (file.exists(gamlss_path)) {
    gamlss <- read_rds_safe(gamlss_path)
    if ("TEST" %in% names(gamlss) && !is.null(gamlss$TEST)) {
      return(gamlss$TEST)
    }
  }

  # Fall back to transfer model if site model has no calibration data
  gamlss_path <- get_data_path("models", "fit", "gamlss")
  if (!file.exists(gamlss_path)) {
    warning("GAMLSS fits not found")
    return(NULL)
  }
  gamlss <- read_rds_safe(gamlss_path)
  if (!"TEST" %in% names(gamlss) || is.null(gamlss$TEST)) return(NULL)
  gamlss$TEST
}

# =============================================================================
# SEM Validation and Extraction Utilities
# =============================================================================

# SEM Group Encoding
# lavaan assigns group numbers alphabetically, so for SEX:
# group == 1 -> "Female" (F comes before M)
# group == 2 -> "Male"
SEM_GROUP_FEMALE <- 1L
SEM_GROUP_MALE <- 2L

#' Validate SEM group encoding from model parameters
#' Checks that group labels in the data match expected alphabetical encoding.
#' @param sem_params SEM parameters data.table with 'group' and optionally 'group.label' columns
#' @return TRUE if validation passes, otherwise stops with informative error
#' @export
validate_sem_groups <- function(sem_params) {
  if (is.null(sem_params) || nrow(sem_params) == 0) {
    warning("Cannot validate SEM groups: sem_params is NULL or empty", call. = FALSE)
    return(FALSE)
  }

  # Check if group.label column exists (lavaan adds this)
  if ("group.label" %in% names(sem_params)) {
    group_labels <- unique(sem_params[, .(group, group.label)])
    setorder(group_labels, group)

    expected_labels <- data.table(
      group = c(1L, 2L),
      expected = c("Female", "Male")
    )

    # Validate encoding
    if (nrow(group_labels) >= 2) {
      if (group_labels[group == 1, group.label] != "Female") {
        stop(sprintf(
          "SEM group encoding mismatch! Expected group 1 = 'Female', got '%s'. Factor level ordering may have changed.",
          group_labels[group == 1, group.label]
        ), call. = FALSE)
      }
      if (group_labels[group == 2, group.label] != "Male") {
        stop(sprintf(
          "SEM group encoding mismatch! Expected group 2 = 'Male', got '%s'. Factor level ordering may have changed.",
          group_labels[group == 2, group.label]
        ), call. = FALSE)
      }
    }
  }

  TRUE
}

#' Safe extraction of SEM parameters with group validation
#' @param params_dt Data.table of SEM parameters
#' @param filter_expr Quoted filter expression (use quote())
#' @param validate_groups Whether to validate group encoding (default TRUE on first call)
#' @return Filtered data.table row(s), or empty data.table if no match
#' @export
safe_extract <- function(params_dt, filter_expr, validate_groups = TRUE) {
  if (is.null(params_dt) || nrow(params_dt) == 0) {
    warning("safe_extract: params_dt is NULL or empty", call. = FALSE)
    return(data.table())
  }

  # Validate groups on first extraction
  if (validate_groups && !exists(".sem_groups_validated", envir = .GlobalEnv)) {
    validate_sem_groups(params_dt)
    assign(".sem_groups_validated", TRUE, envir = .GlobalEnv)
  }

  result <- tryCatch(
    params_dt[eval(filter_expr)],
    error = function(e) {
      warning(sprintf("safe_extract failed: %s", e$message), call. = FALSE)
      data.table()
    }
  )

  if (nrow(result) == 0) {
    warning(sprintf("safe_extract: no rows matched filter"), call. = FALSE)
  }

  result
}

#' Validate sex_diff data structure has expected columns and values
#' @param sex_diff Sex differences data from pipeline
#' @return TRUE if valid, otherwise stops with error
#' @export
validate_sex_diff_structure <- function(sex_diff) {
  required_cols <- c("SIDE", "ROI", "ADJ", "ESTIMATE", "ROI_LABEL", "ADJ_LABEL")
  expected_adj_labels <- get_factor_levels("ADJ_LABEL")

  if (!"OVERALL" %in% names(sex_diff)) {
    stop("sex_diff missing 'OVERALL' element", call. = FALSE)
  }

  validate_columns(sex_diff$OVERALL, required_cols, "sex_diff$OVERALL")

  actual_labels <- unique(sex_diff$OVERALL$ADJ_LABEL)
  missing_labels <- setdiff(expected_adj_labels, actual_labels)
  if (length(missing_labels) > 0) {
    warning("sex_diff$OVERALL$ADJ_LABEL missing expected values: ",
            paste(missing_labels, collapse = ", "))
  }

  TRUE
}

# =============================================================================
# SEM Full Parameter Extraction
# =============================================================================

#' Load SEM full parameters (factor loadings, residuals, covariate effects)
#' @param model_name Optional: specific model to extract (e.g., "HC_COG", "HVR_COG")
#' @return data.table with full parameter estimates
#' @export
load_sem_full_params <- function(model_name = NULL) {
  params_path <- get_data_path("models", "results", "sem_params")
  if (!file.exists(params_path)) {
    warning("SEM parameters file not found: ", params_path)
    return(NULL)
  }

  params.lst <- read_rds_safe(params_path)

  if (!is.null(model_name)) {
    if (!model_name %in% names(params.lst)) {
      warning("Model not found: ", model_name)
      return(NULL)
    }
    return(params.lst[[model_name]])
  }

  # Combine all models
  rbindlist(params.lst, idcol = "MODEL", fill = TRUE)
}

#' Extract SEM measurement model parameters (factor loadings)
#' @param params_dt Data.table from load_sem_full_params
#' @param model_filter Optional model name filter
#' @return data.table with factor loadings
#' @export
extract_sem_factor_loadings <- function(params_dt, model_filter = NULL) {
  if (is.null(params_dt)) return(NULL)

  # Filter to factor loadings (op == "=~")
  loadings <- params_dt[op == "=~"]

  if (!is.null(model_filter)) {
    loadings <- loadings[MODEL == model_filter]
  }

  # Format output
  loadings[, .(
    MODEL = MODEL,
    FACTOR = lhs,
    INDICATOR = rhs,
    GROUP = if ("group" %in% names(loadings)) group else NA_integer_,
    LOADING = est.std,
    SE = se,
    CI_LOWER = ci.lower,
    CI_UPPER = ci.upper,
    PVALUE = pvalue
  )]
}

#' Extract SEM residual variances
#' @param params_dt Data.table from load_sem_full_params
#' @param model_filter Optional model name filter
#' @return data.table with residual variances
#' @export
extract_sem_residual_variances <- function(params_dt, model_filter = NULL) {
  if (is.null(params_dt)) return(NULL)

  # Filter to residual variances (op == "~~" and lhs == rhs)
  residuals <- params_dt[op == "~~" & lhs == rhs]

  if (!is.null(model_filter)) {
    residuals <- residuals[MODEL == model_filter]
  }

  # Format output
  residuals[, .(
    MODEL = MODEL,
    VARIABLE = lhs,
    GROUP = if ("group" %in% names(residuals)) group else NA_integer_,
    VARIANCE = est.std,
    SE = se,
    CI_LOWER = ci.lower,
    CI_UPPER = ci.upper,
    PVALUE = pvalue
  )]
}

#' Extract SEM covariate effects
#' @param params_dt Data.table from load_sem_full_params
#' @param model_filter Optional model name filter
#' @return data.table with covariate effects
#' @export
extract_sem_covariate_effects <- function(params_dt, model_filter = NULL) {
  if (is.null(params_dt)) return(NULL)

  # Define covariates (not brain measures or cognitive factors)
  brain_vars <- c("HC", "LV", "HVR", "HC_RES")
  cog_vars <- c("g", "MEM_s", "PRSP_s", "COG", "MEM", "PRSP")
  covariate_list <- c("AGE", "AGE_sq", "ICV", "IMDP", "EDUC", "SEX", "SITE")

  # Filter to regressions where RHS is a covariate
  covar_effects <- params_dt[
    op == "~" &
    rhs %in% covariate_list
  ]

  if (!is.null(model_filter)) {
    covar_effects <- covar_effects[MODEL == model_filter]
  }

  # Format output
  covar_effects[, .(
    MODEL = MODEL,
    OUTCOME = lhs,
    COVARIATE = rhs,
    GROUP = if ("group" %in% names(covar_effects)) group else NA_integer_,
    BETA = est.std,
    SE = se,
    CI_LOWER = ci.lower,
    CI_UPPER = ci.upper,
    PVALUE = pvalue
  )]
}

#' Get complete SEM model summary for supplementary table
#' @param model_name Model name (e.g., "HVR_COG")
#' @return List with loadings, residuals, covariate_effects, and structural_paths
#' @export
get_sem_model_summary <- function(model_name) {
  params <- load_sem_full_params(model_name)
  if (is.null(params)) return(NULL)

  # Add MODEL column if not present (single model case)
  if (!"MODEL" %in% names(params)) {
    params[, MODEL := model_name]
  }

  # Define brain and cognitive variables for structural paths
  brain_vars <- c("HC", "LV", "HVR", "HC_RES")
  cog_vars <- c("g", "MEM_s", "PRSP_s")

  list(
    loadings = extract_sem_factor_loadings(params),
    residuals = extract_sem_residual_variances(params),
    covariate_effects = extract_sem_covariate_effects(params),
    structural_paths = params[
      op == "~" &
      lhs %in% cog_vars &
      rhs %in% brain_vars,
      .(
        MODEL = MODEL,
        OUTCOME = lhs,
        PREDICTOR = rhs,
        GROUP = if ("group" %in% names(params)) group else NA_integer_,
        BETA = est.std,
        SE = se,
        CI_LOWER = ci.lower,
        CI_UPPER = ci.upper,
        PVALUE = pvalue
      )
    ]
  )
}

# =============================================================================
# SEM Path Extraction Helpers (for manuscript inline R)
# =============================================================================

#' Extract a single SEM structural path (pooled model)
#' @param sem_params Named list of SEM parameter data.tables
#' @param model Model name key in sem_params (e.g., "HC_COG_POOLED")
#' @param lhs_var Left-hand side variable (outcome, e.g., "g")
#' @param rhs_var Right-hand side variable (predictor, e.g., "HVR")
#' @return Single-row data.table or NULL
#' @export
get_sem_path <- function(sem_params, model, lhs_var, rhs_var) {
  params <- sem_params[[model]]
  if (is.null(params)) return(NULL)
  setDT(params)
  row <- params[lhs == lhs_var & rhs == rhs_var & op == "~"]
  if (nrow(row) == 0) return(NULL)
  row[1]
}

#' Extract a single SEM structural path for a specific sex group
#' @param sem_params Named list of SEM parameter data.tables
#' @param model Model name key (e.g., "HC_COG")
#' @param lhs_var Left-hand side variable
#' @param rhs_var Right-hand side variable
#' @param group_num Group number (1 = Female, 2 = Male)
#' @return Single-row data.table or NULL
#' @export
get_sem_sex_path <- function(sem_params, model, lhs_var, rhs_var, group_num) {
  params <- sem_params[[model]]
  if (is.null(params)) return(NULL)
  setDT(params)
  row <- params[lhs == lhs_var & rhs == rhs_var & op == "~" & group == group_num]
  if (nrow(row) == 0) return(NULL)
  row[1]
}

#' Extract age-squared coefficient from SEM model
#' @param sem_params Named list of SEM parameter data.tables
#' @param model_name Model name key (e.g., "HVR_COG_POOLED")
#' @param outcome Outcome variable (e.g., "HVR", "g")
#' @return Single-row data.table or list with NAs
#' @export
get_sem_age_sq_coef <- function(sem_params, model_name, outcome) {
  params <- sem_params[[model_name]]
  if (is.null(params)) return(list(est.std = NA, pvalue = NA))
  setDT(params)
  row <- params[lhs == outcome & rhs == "AGE_sq" & op == "~"]
  if (nrow(row) == 0) return(list(est.std = NA, pvalue = NA))
  row[1]
}

# =============================================================================
# GAMLSS Model Coefficient Extraction
# =============================================================================

#' Load GAMLSS model fits
#' @param use_site_model Whether to use site-controlled models (default TRUE)
#' @return List containing GAMLSS fits
#' @export
load_gamlss_fits <- function(use_site_model = TRUE) {
  if (use_site_model) {
    gamlss_path <- get_data_path("models", "fit", "gamlss_site")
    if (!file.exists(gamlss_path)) {
      gamlss_path <- get_data_path("models", "fit", "gamlss")
    }
  } else {
    gamlss_path <- get_data_path("models", "fit", "gamlss")
  }

  if (!file.exists(gamlss_path)) {
    warning("GAMLSS fits not found")
    return(NULL)
  }

  read_rds_safe(gamlss_path)
}

#' Extract GAMLSS model coefficients
#' @param gamlss_fit A fitted gamlss model object
#' @param model_data Optional data.table used to fit the model. Required for SE
#'   extraction when models are loaded from RDS (the model's internal data
#'   reference is broken after deserialization).
#' @return data.table with coefficients for mu, sigma, nu, tau
#' @export
extract_gamlss_coefficients <- function(gamlss_fit, model_data = NULL) {
  if (is.null(gamlss_fit) || !inherits(gamlss_fit, "gamlss")) {
    return(NULL)
  }

  # Restore the model's data reference for vcov/summary to work.
  # GAMLSS stores the data variable name in fit$call$data; after RDS
  # deserialization that variable no longer exists in the calling environment.
  if (!is.null(model_data)) {
    data_name <- tryCatch(deparse(gamlss_fit$call$data), error = function(e) NULL)
    if (!is.null(data_name) && nchar(data_name) > 0) {
      assign(data_name, model_data, envir = environment())
    }
  }

  # Parse SE from summary() text output — the most reliable method for
  # deserialized GAMLSS objects with smooth (cs/pb) terms.
  parse_se_from_summary <- function(fit) {
    summ_text <- tryCatch(
      capture.output(suppressWarnings(summary(fit, type = "qr"))),
      error = function(e) NULL
    )
    if (is.null(summ_text)) return(list())

    result <- list()
    current_param <- NULL
    in_coef_block <- FALSE

    for (line in summ_text) {
      if (grepl("^Mu Coefficients:", line)) {
        current_param <- "mu"; in_coef_block <- TRUE; next
      } else if (grepl("^Sigma Coefficients:", line)) {
        current_param <- "sigma"; in_coef_block <- TRUE; next
      } else if (grepl("^Nu Coefficients:", line)) {
        current_param <- "nu"; in_coef_block <- TRUE; next
      } else if (grepl("^Tau Coefficients:", line)) {
        current_param <- "tau"; in_coef_block <- TRUE; next
      }

      if (in_coef_block && (grepl("^---", line) || grepl("^Signif", line))) {
        in_coef_block <- FALSE; next
      }
      if (in_coef_block && grepl("Estimate\\s+Std\\. Error", line)) next

      if (in_coef_block && !is.null(current_param)) {
        parts <- trimws(line)
        if (nchar(parts) > 0) {
          tokens <- strsplit(parts, "\\s+")[[1]]
          # Remove significance stars (e.g., "*", "**", "***", ".")
          tokens <- tokens[!grepl("^[*\\.]+$", tokens)]
          if (length(tokens) >= 4) {
            term_name <- tokens[1]
            se_val <- suppressWarnings(as.numeric(tokens[3]))
            # P-value parsing: handle multiple formats
            # Standard format: Term  Estimate  Std.Error  t-value  Pr(>|t|)
            #                  [1]   [2]       [3]        [4]      [5]
            # P-value can be: "2.07e-08", "<2e-16", or "< 2e-16" (two tokens)
            pval <- NA_real_
            if (length(tokens) >= 5) {
              if (tokens[5] == "<" && length(tokens) >= 6) {
                # "< 2e-16" format (two tokens): treat as machine epsilon
                pval <- 2.2e-16
              } else if (grepl("^<", tokens[5])) {
                # "<2e-16" format (one token): extract number after <
                pval <- suppressWarnings(as.numeric(sub("^<", "", tokens[5])))
                if (is.na(pval)) pval <- 2.2e-16
              } else {
                pval <- suppressWarnings(as.numeric(tokens[5]))
              }
            }
            if (!is.na(se_val)) {
              if (is.null(result[[paste0(current_param, "_se")]])) {
                result[[paste0(current_param, "_se")]] <- c()
                result[[paste0(current_param, "_pval")]] <- c()
              }
              result[[paste0(current_param, "_se")]] <- c(result[[paste0(current_param, "_se")]],
                                           setNames(se_val, term_name))
              result[[paste0(current_param, "_pval")]] <- c(result[[paste0(current_param, "_pval")]],
                                           setNames(pval, term_name))
            }
          }
        }
      }
    }
    result
  }

  # Get all SE values by parsing summary output
  se_lookup <- parse_se_from_summary(gamlss_fit)

  # Helper to create coefficient data.table with SE and p-value from parsed summary
  make_coef_dt <- function(fit, what) {
    coef_vals <- tryCatch(coef(fit, what = what), error = function(e) NULL)
    if (is.null(coef_vals) || length(coef_vals) == 0) return(NULL)

    term_names <- names(coef_vals)
    se_key <- paste0(what, "_se")
    pval_key <- paste0(what, "_pval")

    # Match SE values from parsed summary by term name (handles random/smooth
    # terms that appear in coef() but not in summary())
    se_vals <- if (!is.null(se_lookup[[se_key]])) {
      vapply(term_names, function(tn) {
        if (tn %in% names(se_lookup[[se_key]])) se_lookup[[se_key]][[tn]]
        else NA_real_
      }, numeric(1))
    } else {
      rep(NA_real_, length(term_names))
    }

    # Match p-values from parsed summary
    pval_vals <- if (!is.null(se_lookup[[pval_key]])) {
      vapply(term_names, function(tn) {
        if (tn %in% names(se_lookup[[pval_key]])) se_lookup[[pval_key]][[tn]]
        else NA_real_
      }, numeric(1))
    } else {
      rep(NA_real_, length(term_names))
    }

    data.table(
      PARAMETER = what,
      TERM = term_names,
      ESTIMATE = as.numeric(coef_vals),
      SE = as.numeric(se_vals),
      PVALUE = as.numeric(pval_vals)
    )
  }

  # Extract coefficients for each parameter
  coef_list <- list()

  mu_dt <- make_coef_dt(gamlss_fit, "mu")
  if (!is.null(mu_dt)) coef_list$mu <- mu_dt

  sigma_dt <- make_coef_dt(gamlss_fit, "sigma")
  if (!is.null(sigma_dt)) coef_list$sigma <- sigma_dt

  nu_dt <- make_coef_dt(gamlss_fit, "nu")
  if (!is.null(nu_dt)) coef_list$nu <- nu_dt

  tau_dt <- make_coef_dt(gamlss_fit, "tau")
  if (!is.null(tau_dt)) coef_list$tau <- tau_dt

  if (length(coef_list) == 0) return(NULL)

  rbindlist(coef_list, fill = TRUE)
}

#' Get GAMLSS model summary for all ROI/sex combinations
#' @param roi ROI name (e.g., "HVR", "HC")
#' @param adj Adjustment method (e.g., "NON", "RES")
#' @param side Hemisphere ("LR" for bilateral)
#' @param use_site_model Whether to use site-controlled models (default FALSE
#'   for transfer models without random site effect)
#' @return data.table with model coefficients and metadata
#' @export
get_gamlss_model_summary <- function(roi = NULL, adj = "NON", side = "LR",
                                     use_site_model = FALSE) {
  # Suppress all GAMLSS output (models print summary on access)
  invisible(capture.output({
    gamlss.lst <- load_gamlss_fits(use_site_model = use_site_model)
  }, type = "output"))

  if (is.null(gamlss.lst) || !"FINAL" %in% names(gamlss.lst)) {
    warning("GAMLSS FINAL models not found")
    return(NULL)
  }

  results <- list()

  for (sex in c("Female", "Male")) {
    # Navigate the nested structure (suppress any print output)
    sex_models <- suppressWarnings(suppressMessages({
      gamlss.lst$FINAL[[sex]]
    }))
    if (is.null(sex_models)) next

    rois_to_check <- if (!is.null(roi)) roi else names(sex_models)

    for (r in rois_to_check) {
      if (is.null(sex_models[[r]])) next

      adj_models <- sex_models[[r]][[adj]]
      if (is.null(adj_models)) next

      model_info <- adj_models[[side]]
      if (is.null(model_info) || is.null(model_info$FIT)) next

      # Extract coefficients (suppress any print output from GAMLSS models)
      coefs <- suppressWarnings(suppressMessages({
        invisible(capture.output({
          result <- extract_gamlss_coefficients(model_info$FIT, model_data = model_info$DATA)
        }, type = "output"))
        result
      }))
      if (is.null(coefs)) next

      # Add metadata
      coefs[, `:=`(
        SEX = sex,
        ROI = r,
        ADJ = adj,
        SIDE = side,
        FAMILY = model_info$FAMILY
      )]

      results[[paste(sex, r, adj, side, sep = "_")]] <- coefs
    }
  }

  if (length(results) == 0) return(NULL)

  rbindlist(results, fill = TRUE)
}
