#!/usr/bin/env Rscript

# =============================================================================
# Master Pipeline Orchestration
# =============================================================================
# Executes the complete head-size sex effects analysis pipeline
#
# Usage:
#   Rscript R/00_run_pipeline.R [options]
#
# Options:
#   --scripts=01,02,11  Run specific scripts (comma-separated)
#   --from=21           Run from script 21 onwards
#   --to=41             Run up to script 41
#   --force             Force regeneration of all outputs
#   --parallel          Run independent scripts in parallel (where possible)
#   --help              Show this help message
#
# Examples:
#   Rscript R/00_run_pipeline.R                    # Run full pipeline
#   Rscript R/00_run_pipeline.R --from=21 --to=41  # Run scripts 21-41
#   Rscript R/00_run_pipeline.R --scripts=51,61    # Run only scripts 51,61
#   Rscript R/00_run_pipeline.R --force            # Force regenerate all
# =============================================================================

# --- Setup ---
suppressPackageStartupMessages({
  library(here)
  library(yaml)
})

# Source utilities
# Note: plotting.R and tables.R were split into focused modules during refactoring
util_modules <- c(
  "logging.R", "config.R", "data_io.R", "validation.R", "statistics.R",
  # Plotting modules (split from plotting.R)
  "plotting_core.R", "plotting_figures.R", "plotting_pipeline.R",
  # Table modules (split from tables.R)
  "tables_core.R", "tables_data.R", "tables_gt.R", "tables_normative.R",
  # Other utilities
  "formatting.R", "gamlss.R", "export.R"
)

for (util in util_modules) {
  util_path <- here("R/utils", util)
  if (!file.exists(util_path)) {
    stop(sprintf("Required utility module not found: %s", util_path), call. = FALSE)
  }
  source(util_path)
}

rm(util, util_path, util_modules)

# --- Parse Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

# Default options
options <- list(
  scripts = NULL,
  from = NULL,
  to = NULL,
  force = FALSE,
  parallel = FALSE,
  help = FALSE
)

# Parse arguments
if (length(args) > 0) {
  for (arg in args) {
    if (arg == "--help" || arg == "-h") {
      options$help <- TRUE
    } else if (arg == "--force" || arg == "-f") {
      options$force <- TRUE
    } else if (arg == "--parallel" || arg == "-p") {
      options$parallel <- TRUE
    } else if (grepl("^--scripts=", arg)) {
      scripts_str <- sub("^--scripts=", "", arg)
      options$scripts <- as.integer(strsplit(scripts_str, ",")[[1]])
    } else if (grepl("^--from=", arg)) {
      options$from <- as.integer(sub("^--from=", "", arg))
    } else if (grepl("^--to=", arg)) {
      options$to <- as.integer(sub("^--to=", "", arg))
    }
  }
}

# Show help
if (options$help) {
  cat(readLines(here("R/00_run_pipeline.R"))[1:20], sep = "\n")
  quit(status = 0)
}

# --- Initialize ---
init_logger(log_level = get_config("general", "log_level", default = "INFO"))
log_info(paste0(rep("=", 80), collapse = ""))
log_info("PIPELINE EXECUTION START")
log_info(paste0(rep("=", 80), collapse = ""))
log_info("")

# Load configuration
config <- load_config()
log_info("Configuration loaded: %s", here("config/pipeline_config.yaml"))

# Validate packages
log_section("Validating Required Packages")
required_packages <- get_config("packages", "required")
validate_packages(required_packages)

# Set random seed
set_seed()
log_info("Random seed set: %d", get_seed())

# --- Define Pipeline Scripts ---
pipeline_scripts <- list(
  list(
    number = 1,
    name = "01_parse_volumes",
    description = "Parse preprocessed volumes",
    depends_on = NULL,
    output_key = "icc_scale"
  ),
  list(
    number = 2,
    name = "02_quality_control",
    description = "Quality control filtering",
    depends_on = c(1),
    output_key = "qc_combined"
  ),
  list(
    number = 3,
    name = "03_parse_covariates",
    description = "Parse UK Biobank covariates",
    depends_on = NULL,
    output_key = "covars_fst"
  ),
  list(
    number = 4,
    name = "04_clean_ses_cognitive",
    description = "Clean SES and cognitive data",
    depends_on = c(3),
    output_key = "cog_tests"
  ),
  list(
    number = 5,
    name = "05_adjust_headsize",
    description = "Head-size adjustment and matching",
    depends_on = c(1, 2, 3),
    output_key = "hc_hvr_adjusted"
  ),
  list(
    number = 6,
    name = "06_cognitive_factors",
    description = "Cognitive factor analysis",
    depends_on = c(4, 5),
    output_key = "lat-cog_values"
  ),
  list(
    number = 7,
    name = "07_demographics",
    description = "Generate demographic tables (uses table utilities)",
    depends_on = c(5, 6),
    output_key = NULL,  # Produces tables, not data
    uses_utils = c("tables")
  ),
  list(
    number = 8,
    name = "08_normative_tables",
    description = "GAMLSS normative modeling (uses plotting utilities)",
    depends_on = c(5),
    output_key = "gamlss_fits",
    uses_utils = c("plotting", "tables")
  ),
  list(
    number = 9,
    name = "09_sex_differences",
    description = "Sex differences analysis (uses plotting utilities)",
    depends_on = c(5),
    output_key = "sex_differences",
    uses_utils = c("plotting")
  ),
  list(
    number = 10,
    name = "10_hvr_comparison",
    description = "HVR vs HC comparison analysis",
    depends_on = c(5, 6),
    output_key = "hvr_hc_comparison",
    uses_utils = c("plotting", "statistics")
  ),
  list(
    number = 11,
    name = "11_sem_analysis",
    description = "Structural equation modeling (brain-cognition associations)",
    depends_on = c(4, 5, 6),
    output_key = "sem_analysis",
    uses_utils = c("plotting", "tables")
  ),
  list(
    number = 12,
    name = "12_manuscript_objects",
    description = "Pre-compute manuscript environment for Quarto rendering",
    depends_on = c(7, 8, 9, 10, 11),
    output_key = "manuscript_env",
    uses_utils = c("tables")
  )
)

# --- Filter Scripts to Run ---
scripts_to_run <- pipeline_scripts

if (!is.null(options$scripts)) {
  scripts_to_run <- Filter(function(s) s$number %in% options$scripts, pipeline_scripts)
  log_info("Running selected scripts: %s",
           paste(sapply(scripts_to_run, function(s) s$number), collapse = ", "))
} else if (!is.null(options$from) || !is.null(options$to)) {
  from_num <- if (!is.null(options$from)) options$from else 1
  to_num <- if (!is.null(options$to)) options$to else 99
  scripts_to_run <- Filter(function(s) s$number >= from_num && s$number <= to_num,
                           pipeline_scripts)
  log_info("Running scripts %d to %d", from_num, to_num)
}

# --- Check Dependencies ---
check_dependencies <- function(script, completed_scripts) {
  if (is.null(script$depends_on)) {
    return(TRUE)
  }

  missing_deps <- setdiff(script$depends_on, completed_scripts)
  if (length(missing_deps) > 0) {
    log_warn("Script %02d depends on scripts: %s (not yet run)",
             script$number, paste(missing_deps, collapse = ", "))
    return(FALSE)
  }

  return(TRUE)
}

# --- Run Individual Script ---
run_script <- function(script, force = FALSE) {
  script_num <- sprintf("%02d", script$number)
  script_path <- here("R", paste0(script$name, ".R"))

  log_info("")
  log_info(paste0(rep("=", 80), collapse = ""))
  log_info("SCRIPT %s: %s", script_num, script$description)
  log_info(paste0(rep("=", 80), collapse = ""))

  # Check if output exists and is up-to-date
  # Config structure: scripts.force_regenerate.script_name (not scripts.script_name.force_regenerate)
  if (!force && !is.null(script$output_key)) {
    output_path <- get_data_path("processed", script$output_key)
    script_short_name <- gsub("^\\d+_", "", script$name)
    force_regen <- get_script_setting("force_regenerate", script_short_name, default = FALSE)
    if (file.exists(output_path) && !force_regen) {
      log_info("Output exists, skipping: %s", basename(output_path))
      return(list(success = TRUE, skipped = TRUE))
    }
  }

  # Check script exists
  if (!file.exists(script_path)) {
    log_error("Script not found: %s", script_path)
    return(list(success = FALSE, error = "Script not found"))
  }

  # Run script
  start_time <- Sys.time()
  result <- tryCatch({
    source(script_path, local = new.env())
    list(success = TRUE, skipped = FALSE)
  }, error = function(e) {
    log_error("Script failed: %s", e$message)
    list(success = FALSE, error = e$message)
  })

  duration <- difftime(Sys.time(), start_time, units = "secs")
  if (result$success) {
    status <- if (result$skipped) "SKIPPED" else "COMPLETED"
    log_info("Script %s: %s (%.1f seconds)", status, script_num, as.numeric(duration))
  } else {
    log_error("Script FAILED: %s", script_num)
  }

  return(result)
}

# --- Main Execution Loop ---
log_section("Pipeline Execution Plan")
log_info("Total scripts to run: %d", length(scripts_to_run))
for (script in scripts_to_run) {
  log_info("  %02d: %s", script$number, script$description)
}
log_info("")

completed_scripts <- integer(0)
failed_scripts <- list()

for (script in scripts_to_run) {
  # Check dependencies
  if (!check_dependencies(script, completed_scripts)) {
    log_error("Skipping script %02d due to missing dependencies", script$number)
    next
  }

  # Run script
  result <- run_script(script, force = options$force)

  if (result$success) {
    completed_scripts <- c(completed_scripts, script$number)
  } else {
    failed_scripts[[length(failed_scripts) + 1]] <- list(
      script = script,
      error = result$error
    )

    # Stop pipeline on error
    log_error("Pipeline halted due to script failure")
    break
  }
}

# --- Summary ---
log_info("")
log_info(paste0(rep("=", 80), collapse = ""))
log_info("PIPELINE EXECUTION SUMMARY")
log_info(paste0(rep("=", 80), collapse = ""))
log_info("Completed scripts: %d", length(completed_scripts))
if (length(failed_scripts) > 0) {
  log_error("Failed scripts: %d", length(failed_scripts))
  for (failed in failed_scripts) {
    log_error("  %02d: %s - %s",
              failed$script$number,
              failed$script$name,
              failed$error)
  }
}
log_info("")
log_info("Log file: %s", get_log_file())
log_info(paste0(rep("=", 80), collapse = ""))

# Exit with appropriate status
if (length(failed_scripts) > 0) {
  quit(status = 1)
} else {
  quit(status = 0)
}
