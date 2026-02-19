# =============================================================================
# Logging Utilities
# =============================================================================
# Functions for consistent logging throughout the pipeline
# Author: Pipeline Refactoring
# Date: 2024
# =============================================================================

library(here)

# Initialize logger
.log_file <- NULL
.log_level <- "INFO"

#' Initialize logging system
#'
#' @param log_dir Directory for log files
#' @param log_level Minimum level to log (DEBUG, INFO, WARN, ERROR)
#' @param append Append to existing log file
#' @export
init_logger <- function(log_dir = here("logs"),
                       log_level = "INFO",
                       append = FALSE) {
  # Create log directory with error handling
  if (!dir.exists(log_dir)) {
    if (!dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)) {
      warning("Failed to create log directory: ", log_dir, call. = FALSE)
    }
  }

  # Set log file with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  .log_file <<- here(log_dir, paste0("pipeline_", timestamp, ".log"))
  .log_level <<- log_level

  # Write header
  if (!append) {
    writeLines(
      c(
        paste0(rep("=", 80), collapse = ""),
        paste("Pipeline Execution Log"),
        paste("Started:", Sys.time()),
        paste("R Version:", R.version.string),
        paste0(rep("=", 80), collapse = ""),
        ""
      ),
      .log_file
    )
  }

  log_info("Logger initialized: %s", .log_file)
  invisible(.log_file)
}

#' Get current log file path
#' @export
get_log_file <- function() {
  .log_file
}

#' Write log message
#'
#' @param level Log level (DEBUG, INFO, WARN, ERROR)
#' @param message Message to log
#' @param ... Arguments for sprintf formatting
#' @export
log_message <- function(level, message, ...) {
  if (is.null(.log_file)) {
    init_logger()
  }

  # Check if level should be logged
  levels <- c("DEBUG" = 1, "INFO" = 2, "WARN" = 3, "ERROR" = 4)
  if (levels[level] < levels[.log_level]) {
    return(invisible(NULL))
  }

  # Format message
  if (length(list(...)) > 0) {
    message <- sprintf(message, ...)
  }

  # Replace glue-style placeholders
  message <- gsub("\\{([^}]+)\\}", "\\1", message)

  # Create log entry
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- sprintf("[%s] %s: %s", timestamp, level, message)

  # Write to file
  cat(entry, "\n", file = .log_file, append = TRUE)

  # Also print to console for INFO and above
  if (levels[level] >= levels["INFO"]) {
    message(entry)
  }

  invisible(NULL)
}

#' Log debug message
#' @export
log_debug <- function(message, ...) {
  log_message("DEBUG", message, ...)
}

#' Log info message
#' @export
log_info <- function(message, ...) {
  log_message("INFO", message, ...)
}

#' Log warning message
#' @export
log_warn <- function(message, ...) {
  log_message("WARN", message, ...)
}

#' Log error message
#' @export
log_error <- function(message, ...) {
  log_message("ERROR", message, ...)
}

#' Log script start
#' @export
log_script_start <- function(script_name) {
  log_info(paste0(rep("-", 80), collapse = ""))
  log_info("Starting script: %s", script_name)
  log_info(paste0(rep("-", 80), collapse = ""))
}

#' Log script end
#' @export
log_script_end <- function(script_name, success = TRUE) {
  status <- ifelse(success, "COMPLETED", "FAILED")
  log_info("Script %s: %s", status, script_name)
  log_info(paste0(rep("-", 80), collapse = ""))
}

#' Log section
#' @export
log_section <- function(section_name) {
  log_info("")
  log_info(">>> %s", section_name)
}

#' Time a code block and log duration
#' @export
log_time <- function(description, expr) {
  log_info("Starting: %s", description)
  start_time <- Sys.time()

  result <- tryCatch(
    expr,
    error = function(e) {
      log_error("Error in %s: %s", description, e$message)
      stop(e)
    }
  )

  duration <- difftime(Sys.time(), start_time, units = "secs")
  log_info("Completed: %s (%.2f seconds)", description, as.numeric(duration))

  invisible(result)
}

#' Rotate log files
#'
#' Archives old log files and optionally deletes logs older than max_days.
#' Call periodically to prevent log directory from growing unbounded.
#'
#' @param log_dir Directory containing log files
#' @param max_files Maximum number of log files to keep
#' @param max_days Maximum age of log files in days (NULL = no age limit)
#' @param archive_dir Directory for archived logs (NULL = delete old logs)
#' @return Number of files rotated/deleted
#' @export
rotate_logs <- function(log_dir = here("logs"),
                        max_files = 50,
                        max_days = 30,
                        archive_dir = NULL) {
  if (!dir.exists(log_dir)) {
    return(0)
  }

  # Get all log files

log_files <- list.files(log_dir, pattern = "\\.log$", full.names = TRUE)
  if (length(log_files) == 0) {
    return(0)
  }

  # Get file info
  file_info <- file.info(log_files)
  file_info$path <- log_files
  file_info <- file_info[order(file_info$mtime, decreasing = TRUE), ]

  files_to_remove <- character(0)

  # Remove files exceeding max count
  if (nrow(file_info) > max_files) {
    files_to_remove <- c(files_to_remove, file_info$path[(max_files + 1):nrow(file_info)])
  }

  # Remove files older than max_days
  if (!is.null(max_days)) {
    cutoff_time <- Sys.time() - as.difftime(max_days, units = "days")
    old_files <- file_info$path[file_info$mtime < cutoff_time]
    files_to_remove <- c(files_to_remove, old_files)
  }

  files_to_remove <- unique(files_to_remove)

  if (length(files_to_remove) == 0) {
    return(0)
  }

  # Archive or delete with error handling
  if (!is.null(archive_dir)) {
    if (!dir.exists(archive_dir)) {
      if (!dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)) {
        warning("Failed to create archive directory: ", archive_dir, call. = FALSE)
        return(invisible(0L))
      }
    }
    archived <- 0L
    for (f in files_to_remove) {
      archive_path <- file.path(archive_dir, basename(f))
      if (file.rename(f, archive_path)) {
        archived <- archived + 1L
      } else {
        warning("Failed to archive log file: ", f, call. = FALSE)
      }
    }
    message(sprintf("Archived %d log files to %s", archived, archive_dir))
  } else {
    removed <- suppressWarnings(file.remove(files_to_remove))
    if (any(!removed)) {
      warning("Failed to remove some log files", call. = FALSE)
    }
    message(sprintf("Removed %d old log files", sum(removed)))
  }

  invisible(length(files_to_remove))
}

#' Get log file summary
#'
#' Returns summary statistics about log files in the log directory.
#'
#' @param log_dir Directory containing log files
#' @return data.frame with log file information
#' @export
get_log_summary <- function(log_dir = here("logs")) {
  if (!dir.exists(log_dir)) {
    return(data.frame())
  }

  log_files <- list.files(log_dir, pattern = "\\.log$", full.names = TRUE)
  if (length(log_files) == 0) {
    return(data.frame())
  }

  file_info <- file.info(log_files)
  data.frame(
    file = basename(log_files),
    size_kb = round(file_info$size / 1024, 1),
    modified = file_info$mtime,
    age_days = round(as.numeric(difftime(Sys.time(), file_info$mtime, units = "days")), 1)
  )
}

#' Set logging verbosity level
#'
#' Allows runtime adjustment of logging verbosity.
#'
#' @param level Log level: "DEBUG", "INFO", "WARN", "ERROR"
#' @export
set_log_level <- function(level) {
  valid_levels <- c("DEBUG", "INFO", "WARN", "ERROR")
  level <- toupper(level)
  if (!level %in% valid_levels) {
    stop("Invalid log level. Must be one of: ", paste(valid_levels, collapse = ", "),
         call. = FALSE)
  }
  .log_level <<- level
  message(sprintf("Log level set to: %s", level))
  invisible(level)
}
