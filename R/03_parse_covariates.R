#!/usr/bin/env Rscript

# =============================================================================
# Parse UK Biobank Covariates
# =============================================================================
# Extracts and parses covariates from UK Biobank parquet file
#
# Inputs:
#   - data/raw/fields.yaml: Selected field IDs
#   - data/raw/fields.tsv: Full UK Biobank field metadata
#   - data/raw/ukb676574: UK Biobank parquet data
#
# Outputs:
#   - data/derivatives/metadata.rds: Field metadata
#   - data/fst/ukb_covars.fst: Parsed covariates (FST format)
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
# library(stringr)
library(duckdb)
library(yaml)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("03_parse_covariates.R")

# Load configuration
config <- load_config()

# ----- Input Files -----
log_section("Checking requirements")
ukb_dataset_file <- get_data_path("raw", "ukb_dataset")
fields_yaml_file <- get_data_path("raw", "fields_yaml")
fields_tsv_file <- get_data_path("raw", "fields_tsv")

check_files_exist(c(fields_yaml_file, fields_tsv_file, ukb_dataset_file))

# ----- UK Biobank Field Metadata -----
log_section("Loading UK Biobank field metadata")

# Selected fields from YAML
log_info("Reading selected fields from YAML")
sel_fields.dt <- fields_yaml_file |>
  read_yaml() |>
  unlist() |>
  as.data.table(keep.rownames = TRUE) |>
  setnames(c("COLNAME", "FIELD"))
log_info("Loaded %d selected fields", nrow(sel_fields.dt))

# Parse full UKB field metadata
cols.lst <- list(
  ORIG = c(
    "field_id", "title", "num_participants", "value_type", "base_type",
    "instanced", "arrayed", "sexed", "units", "main_category",
    "instance_min", "instance_max", "array_min", "array_max", "notes"
  ),
  NEW = c(
    "FIELD", "TITLE", "N", "TYPE_value", "TYPE_base",
    "INSTANCED", "ARRAYED", "SEXED", "UNITS", "CATEGORY",
    "INSTANCE_min", "INSTANCE_max", "ARRAY_min", "ARRAY_max", "NOTES"
  )
)

ukb_fields.dt <- read_csv_safe(
  fields_tsv_file,
  select = cols.lst$ORIG,
  col.names = cols.lst$NEW,
  description = "UK Biobank field definitions"
)

# Inner join
metadata.dt <- ukb_fields.dt[
  sel_fields.dt,
  on = "FIELD"
] |>
  setcolorder("COLNAME", after = "FIELD")

# Separate notes
notes.dt <- metadata.dt[, .(FIELD, NOTES)]
metadata.dt[, NOTES := NULL]
rm(cols.lst, sel_fields.dt, ukb_fields.dt)

# ----- Database Querying -----
log_section("Querying UK Biobank parquet database")
con <- dbConnect(duckdb())

# Filtering
all_cols.v <- dbGetQuery(
  con,
  sprintf("DESCRIBE SELECT * FROM read_parquet('%s')", ukb_dataset_file)
)$column_name

# Regex pattern for columns of interest
cols.re <- metadata.dt$FIELD |>
  sprintf(fmt = "f\\.%i\\.") |>
  paste(collapse = "|") |>
  paste0("^f\\.eid$|", ... = _)

# Filter columns on REGEX
sel_cols.v <- sprintf('"%s"', grep(cols.re, all_cols.v, value = TRUE))
log_info("Selected %d columns from parquet file", length(sel_cols.v))

# Query
log_time("Querying parquet database", {
  covars.subdt <- dbGetQuery(
    con,
    sprintf(
      "SELECT %s FROM read_parquet('%s')",
      paste(sel_cols.v, collapse = ", "),
      ukb_dataset_file
    )
  )
  setDT(covars.subdt)
  covars.subdt
})

log_info(
  "Retrieved %d rows, %d columns", nrow(covars.subdt), ncol(covars.subdt)
)
dbDisconnect(con, shutdown = TRUE)
rm(
  all_cols.v,
  cols.re,
  con,
  sel_cols.v,
  fields_yaml_file,
  fields_tsv_file,
  ukb_dataset_file
)

# ----- Data Cleaning -----
log_section("Cleaning and renaming columns")
orig.dt <- copy(covars.subdt)
setnames(covars.subdt, "f.eid", "EID", skip_absent = TRUE)
cols.v <- names(covars.subdt)
last_col <- "EID"

n_fields <- metadata.dt[, .N]
log_debug("Processing %d fields from metadata", n_fields)
for (i in metadata.dt[, seq_len(.N)]) {
  field <- metadata.dt[i, FIELD]
  name_col <- metadata.dt[i, COLNAME]
  instanced <- metadata.dt[i, INSTANCED]
  arrayed <- metadata.dt[i, ARRAYED]

  # Relevant columns
  cols_orig.v <- sprintf("f\\.%i\\.", field) |> grep(cols.v, value = TRUE)
  if (length(cols_orig.v) == 0) next

  # Rename
  cols_new.v <- sprintf("f\\.%i", field) |> sub(name_col, cols_orig.v)

  # Keep sessions?
  cols_new.v <- if (instanced) {
    sub("\\.(\\d)\\.", "_ses\\1.", cols_new.v)
  } else {
    sub("\\.\\d(?=\\.)", "", cols_new.v, perl = TRUE)
  }

  # Delete unneeded arrays
  cols_new.v <- if (arrayed) {
    sub("\\.(\\d*)$", "_\\1", cols_new.v)
  } else {
    sub("\\.\\d*$", "", cols_new.v)
  }

  # Rename columns
  setnames(covars.subdt, cols_orig.v, cols_new.v)
  setcolorder(covars.subdt, cols_new.v, after = last_col)
  last_col <- cols_new.v[length(cols_new.v)]

  # Aggregate diagnoses: Family History | ICD[9,10]
  if (name_col %like% "ILL|ICD") {
    if (instanced) {
      for (j in metadata.dt[i, seq(INSTANCE_min, INSTANCE_max)]) {
        out_col <- sprintf("%s_ses%i", name_col, j)
        src_cols.v <- sprintf("ses%i", j) |> grep(cols_new.v, value = TRUE)
        covars.subdt[
          , (out_col) := do.call(paste, c(.SD, sep = ",")),
          .SDcols = src_cols.v
        ] |> suppressWarnings()
        # Clean up NA strings from paste: remove ",NA" (trailing), "NA," (leading), and pure "NA"
        covars.subdt[, (out_col) := gsub("(^NA,|,NA$|,NA(?=,))", "", get(out_col), perl = TRUE)]
        covars.subdt[get(out_col) %in% c("NA", ""), (out_col) := NA_character_]
        setcolorder(covars.subdt, out_col, before = src_cols.v[1])
        covars.subdt[, (src_cols.v) := NULL] |> suppressWarnings()
      }
      rm(j)
    } else {
      out_col <- name_col
      src_cols.v <- cols_new.v
      covars.subdt[
        , (out_col) := do.call(paste, c(.SD, sep = ",")),
        .SDcols = src_cols.v
      ] |> suppressWarnings()
      # Clean up NA strings from paste: remove ",NA" (trailing), "NA," (leading), and pure "NA"
      covars.subdt[, (out_col) := gsub("(^NA,|,NA$|,NA(?=,))", "", get(out_col), perl = TRUE)]
      covars.subdt[get(out_col) %in% c("NA", ""), (out_col) := NA_character_]
      setcolorder(covars.subdt, out_col, before = src_cols.v[1])
      covars.subdt[, (src_cols.v) := NULL] |> suppressWarnings()
    }
    last_col <- out_col
    rm(out_col, src_cols.v)
  }

  rm(i, field, name_col, instanced, arrayed, cols_orig.v, cols_new.v)
}

covars.dt <- copy(covars.subdt)
rm(covars.subdt, cols.v, last_col)

log_info(
  "Final covariate data: %d rows, %d columns",
  nrow(covars.dt),
  ncol(covars.dt)
)

# ----- Validate SEX column if present (CL4) -----
if ("SEX" %in% names(covars.dt)) {
  log_section("Validating SEX column")
  expected_sex_levels <- get_factor_levels("SEX")

  # Check for unexpected values
  actual_sex_values <- unique(covars.dt$SEX)
  actual_sex_values <- actual_sex_values[!is.na(actual_sex_values)]

  # Handle numeric encoding (UK Biobank uses 0=Female, 1=Male)
  if (all(actual_sex_values %in% c(0, 1))) {
    log_info("Converting SEX from numeric (0/1) to character (Female/Male)")
    covars.dt[, SEX := fifelse(SEX == 0, "Female", "Male")]
    actual_sex_values <- unique(covars.dt$SEX)
    actual_sex_values <- actual_sex_values[!is.na(actual_sex_values)]
  }

  unexpected_sex <- setdiff(actual_sex_values, expected_sex_levels)
  if (length(unexpected_sex) > 0) {
    log_warn("Unexpected SEX values found: %s", paste(unexpected_sex, collapse = ", "))
    log_warn("Expected: %s", paste(expected_sex_levels, collapse = ", "))
  }

  # Set as factor with standard levels
  covars.dt[, SEX := factor(SEX, levels = expected_sex_levels)]

  n_female <- covars.dt[SEX == "Female", .N]
  n_male <- covars.dt[SEX == "Male", .N]
  n_na <- covars.dt[is.na(SEX), .N]
  log_info("SEX distribution: Female=%d, Male=%d, NA=%d", n_female, n_male, n_na)
}

# ----- Outputs -----
log_section("Saving outputs")

# Validate outputs before saving
validate_not_empty(metadata.dt, "field metadata")
validate_columns(metadata.dt, c("FIELD", "COLNAME", "TITLE"), "metadata")
log_info("Field metadata: %d fields", nrow(metadata.dt))

validate_not_empty(notes.dt, "field notes")
log_info("Field notes: %d entries", nrow(notes.dt))

validate_not_empty(covars.dt, "covariates")
validate_columns(covars.dt, "EID", "covariates")
log_info("Covariates data: %d subjects", covars.dt[!duplicated(EID), .N])

write_rds_safe(
  metadata.dt,
  get_data_path("processed", "metadata"),
  description = "Field metadata"
)

write_rds_safe(
  notes.dt,
  get_data_path("processed", "notes"),
  description = "Field notes"
)

write_fst_safe(
  covars.dt,
  get_data_path("processed", "covars_fst"),
  compress = 85,
  description = "UK Biobank covariates"
)

log_script_end("03_parse_covariates.R", success = TRUE)
