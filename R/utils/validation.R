# =============================================================================
# Validation Utilities
# =============================================================================
# Functions for validating data and parameters
# =============================================================================

library(data.table)
source(here::here("R/utils/logging.R"))

#' Validate data.table has required columns
#'
#' @param data_table data.table to validate
#' @param required_cols Character vector of required column names
#' @param data_name Name of data for error message
#' @export
validate_columns <- function(data_table, required_cols, data_name = "data") {
  missing_cols <- setdiff(required_cols, names(data_table))

  if (length(missing_cols) > 0) {
    msg <- sprintf(
      "Missing required columns in %s: %s",
      data_name,
      paste(missing_cols, collapse = ", ")
    )
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  log_debug("All required columns present in %s", data_name)
  invisible(TRUE)
}

#' Validate data structure is not empty
#'
#' @param x Data structure to validate (data.table, data.frame, or list)
#' @param data_name Name of data for error message
#' @export
validate_not_empty <- function(x, data_name = "data") {
  if (is.data.frame(x)) {
    if (nrow(x) == 0) {
      msg <- sprintf("Dataset %s is empty", data_name)
      log_error(msg)
      stop(msg, call. = FALSE)
    }
    log_debug("Dataset %s has %d rows", data_name, nrow(x))
  } else if (is.list(x)) {
    if (length(x) == 0) {
      msg <- sprintf("List %s is empty", data_name)
      log_error(msg)
      stop(msg, call. = FALSE)
    }
    log_debug("List %s has %d elements", data_name, length(x))
  } else {
    if (length(x) == 0) {
      msg <- sprintf("Object %s is empty", data_name)
      log_error(msg)
      stop(msg, call. = FALSE)
    }
    log_debug("Object %s has length %d", data_name, length(x))
  }
  invisible(TRUE)
}

#' Validate numeric columns are in expected range
#'
#' @param data_table data.table
#' @param col_name Column name
#' @param min_val Minimum expected value
#' @param max_val Maximum expected value
#' @param na_allowed Are NAs allowed
#' @export
validate_num_range <- function(
    data_table, col_name, min_val = -Inf, max_val = Inf, na_allowed = TRUE) {
  if (!col_name %in% names(data_table)) {
    msg <- sprintf("Column %s not found", col_name)
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  values <- data_table[[col_name]]

  # Check NAs
  n_na <- sum(is.na(values))
  if (n_na > 0 && !na_allowed) {
    msg <- sprintf("Column %s contains %d NAs (not allowed)", col_name, n_na)
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  # Check range
  values_no_na <- values[!is.na(values)]
  if (length(values_no_na) > 0) {
    if (any(values_no_na < min_val) || any(values_no_na > max_val)) {
      actual_min <- min(values_no_na)
      actual_max <- max(values_no_na)
      msg <- sprintf(
        "Column %s out of range [%.2f, %.2f], actual range [%.2f, %.2f]",
        col_name, min_val, max_val, actual_min, actual_max
      )
      log_error(msg)
      stop(msg, call. = FALSE)
    }
  }

  log_debug(
    "Column %s validated (range: [%.2f, %.2f])", col_name, min_val, max_val
  )
  invisible(TRUE)
}

#' Validate categorical column has expected levels
#'
#' @param data_table data.table
#' @param col_name Column name
#' @param expected_levels Expected factor levels or unique values
#' @param allow_extra Allow extra levels not in expected
#' @export
validate_categorical <- function(
    data_table, col_name, expected_levels, allow_extra = FALSE) {
  if (!col_name %in% names(data_table)) {
    msg <- sprintf("Column %s not found", col_name)
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  actual_levels <- unique(data_table[[col_name]])
  actual_levels <- actual_levels[!is.na(actual_levels)]

  missing_levels <- setdiff(expected_levels, actual_levels)
  extra_levels <- setdiff(actual_levels, expected_levels)

  if (length(missing_levels) > 0) {
    log_warn(
      "Column %s missing expected levels: %s",
      col_name, paste(missing_levels, collapse = ", ")
    )
  }

  if (length(extra_levels) > 0 && !allow_extra) {
    msg <- sprintf(
      "Column %s has unexpected levels: %s",
      col_name, paste(extra_levels, collapse = ", ")
    )
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  log_debug("Column %s validated (%d levels)", col_name, length(actual_levels))
  invisible(TRUE)
}

#' Validate age values are reasonable
#'
#' @param data_table data.table
#' @param age_col Name of age column
#' @param min_age Minimum reasonable age
#' @param max_age Maximum reasonable age
#' @export
validate_age <- function(
    data_table, age_col = "AGE", min_age = 18, max_age = 120, ...) {
  validate_num_range(
    data_table, age_col,
    min_val = min_age, max_val = max_age, ...
  )
}

#' Validate IDs are unique
#'
#' @param data_table data.table
#' @param id_cols Character vector of ID column names
#' @export
validate_unique_ids <- function(data_table, id_cols) {
  if (!all(id_cols %in% names(data_table))) {
    missing <- setdiff(id_cols, names(data_table))
    msg <- sprintf("ID columns not found: %s", paste(missing, collapse = ", "))
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  n_before <- nrow(data_table)
  n_unique <- nrow(unique(data_table[, ..id_cols]))

  if (n_unique < n_before) {
    n_duplicates <- n_before - n_unique
    msg <- sprintf(
      "Found %d duplicate IDs in columns: %s",
      n_duplicates, paste(id_cols, collapse = ", ")
    )
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  log_debug("IDs validated (%d unique)", n_unique)
  invisible(TRUE)
}

#' Check package is installed
#'
#' @param package Package name
#' @export
validate_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    msg <- sprintf("Required package not installed: %s", package)
    log_error(msg)
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate required packages are installed
#'
#' @param packages Character vector of package names
#' @export
validate_packages <- function(packages) {
  missing <- character(0)

  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }

  if (length(missing) > 0) {
    msg <- sprintf(
      "Missing required packages: %s\nInstall with: install.packages(c(%s))",
      paste(missing, collapse = ", "),
      paste0("'", missing, "'", collapse = ", ")
    )
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  log_debug("All required packages installed (%d packages)", length(packages))
  invisible(TRUE)
}
