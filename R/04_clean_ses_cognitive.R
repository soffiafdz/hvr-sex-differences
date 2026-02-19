#!/usr/bin/env Rscript

# =============================================================================
# Clean SES and Cognitive Data
# =============================================================================
# Processes socioeconomic status (SES) and cognitive test data
#
# Inputs:
#   - data/fst/ukb_covars.fst: UK Biobank covariates
#   - data/derivatives/metadata.rds: Field metadata
#   - data/derivatives/notes.rds: Field notes
#
# Outputs:
#   - data/derivatives/ses_data.rds: SES scores by country
#   - data/derivatives/cog-tests_metadata.rds: Cognitive test metadata
#   - data/derivatives/cog-tests.rds: Cleaned cognitive test data
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
library(stringr)
library(progress)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("04_clean_ses_cognitive.R")

# Load configuration
config <- load_config()

# ----- Input Files -----
log_section("Loading input data")
metadata_file <- get_data_path("processed", "metadata")
notes_file <- get_data_path("processed", "notes")
covars_file <- get_data_path("processed", "covars_fst")

if (!check_files_exist(
  c(covars_file, metadata_file, notes_file),
  stop_on_missing = FALSE
)) {
  log_warn("Running prerequisite script")
  source(here("R/03_parse_covariates.R"))
}
metadata.dt <- read_rds_safe(metadata_file, description = "Field metadata")
notes.dt <- read_rds_safe(notes_file, description = "Field notes")
covars.dt <- read_fst_safe(
  covars_file,
  as_data_table = TRUE, description = "UK Biobank covariates"
)
rm(covars_file, metadata_file, notes_file)

# ----- Socioeconomic Status (SES) -----
log_section("Processing SES data")
cols <- covars.dt |>
  names() |>
  grep(pattern = "EID|*_score_*")

ses.dt <- covars.dt[, ..cols] |>
  melt(id = 1, value = "SCORE") |>
  suppressWarnings() |>
  na.omit() |>
  (\(data_table) {
    data_table[
      , c("ITEM", "COUNTRY") := tstrsplit(variable, "_score_", fixed = TRUE)
    ][
      , let(variable = NULL, COUNTRY = toupper(COUNTRY))
    ]
  })() |>
  dcast(EID + COUNTRY ~ ITEM, value.var = "SCORE")

# Fix incorrect coding for Scotland & Wales
wrong_educ_eid <- ses.dt[is.na(EMPL), .(EID)]

really_sct <- ses.dt[
  wrong_educ_eid,
  on = "EID"
][
  !is.na(EDUC)
][
  "WLS",
  on = "COUNTRY"
][
  , .(EID, COUNTRY = "SCT", EDUC)
]

really_wls <- ses.dt[
  wrong_educ_eid,
  on = "EID"
][
  !is.na(EDUC)
][
  "SCT",
  on = "COUNTRY"
][
  , .(EID, COUNTRY = "WLS", EDUC)
]

ses.subdt <- ses.dt[wrong_educ_eid, on = "EID"][is.na(EDUC)][, EDUC := NULL]
sct.dt <- really_sct[ses.subdt["SCT", on = "COUNTRY"], on = .(EID, COUNTRY)]
wls.dt <- really_wls[ses.subdt["WLS", on = "COUNTRY"], on = .(EID, COUNTRY)]

ses.dt <- ses.dt |>
  na.omit() |>
  rbind(sct.dt, wls.dt, use.names = TRUE)

rm(cols, ses.subdt, sct.dt, wls.dt, really_sct, really_wls, wrong_educ_eid)
log_info("SES data processed: %d rows", nrow(ses.dt))

# ----- Cognitive Tests -----
log_section("Processing cognitive test data")
# Functions for parsing cognitive data
ext_cols.fn <- function(
    data_table,
    pattern,
    include_eid = TRUE,
    ignore_case = FALSE) {
  if (include_eid) pattern <- sprintf("EID|%s", pattern)
  data_table[, .SD, .SDcols = patterns(pattern, ignore.case = ignore_case)]
}

ext_sess.fn <- function(
    data_table,
    id_vars = "EID",
    ses.re = "ses\\d",
    ses_name = "SESSION") {
  data_table <- melt(data_table, id = id_vars) |> suppressWarnings()
  data_table[, (ses_name) := str_extract(variable, ses.re)]
  data_table[, variable := str_remove(variable, paste0("_", ses.re))]
  data_table <- dcast(data_table, ... ~ variable, value.var = "value")
  return(data_table)
}

ext_arrays.fn <- function(
    data_table,
    id_vars = "EID",
    include_sess = TRUE,
    sess_name = "SESSION",
    array_name = "ARRAY",
    array.re = "\\d+$") {
  if (include_sess) id_vars <- c(id_vars, sess_name)
  data_table <- melt(data_table, id = id_vars)
  data_table[
    , (array_name) := variable |> str_extract(array.re) |> as.integer()
  ]
  data_table[, variable := str_remove(variable, paste0("_", array.re))]
  data_table <- dcast(data_table, ... ~ variable, value.var = "value")
  return(data_table)
}

sep_online.fn <- function(
    data_table,
    id_vars = "EID",
    online.re = "onl",
    online_name = "ONLINE",
    include_sess = TRUE,
    sess_name = "SESSION",
    include_array = TRUE,
    array_name = "ARRAY") {
  if (include_sess) id_vars <- c(id_vars, sess_name)
  if (include_array) id_vars <- c(id_vars, array_name)
  data_table <- melt(data_table, id = id_vars)
  data_table[, (online_name) := fifelse(variable %like% online.re, TRUE, FALSE)]
  data_table[, variable := str_remove(variable, paste0("_", online.re))]
  data_table <- dcast(data_table, ... ~ variable, value.var = "value")
  return(data_table)
}

# Extend metadata for cognitive scores
# Filter by known cognitive test patterns instead of hardcoded row index
# This makes the code robust to changes in metadata ordering
cog_test_patterns <- c(
  "PRS", "NUM", "PRMEM", "FLINT", "MATS", "TRLS", "TOWER",
  "REACT", "SYM", "VOCAB", "DSS", "TMT"  # Include common cognitive test abbreviations
)
cog_pattern <- paste(cog_test_patterns, collapse = "|")
cog_meta.dt <- metadata.dt[
  COLNAME %like% cog_pattern,
  .(TITLE, UNITS, COLNAME, INSTANCED, ARRAYED)
]

if (nrow(cog_meta.dt) == 0) {
  log_warn("No cognitive test metadata found matching patterns: %s", cog_pattern)
  log_warn("Falling back to row-based extraction (rows 44:N)")
  cog_meta.dt <- metadata.dt[44:.N, .(TITLE, UNITS, COLNAME, INSTANCED, ARRAYED)]
}
log_debug("Extracted %d cognitive test metadata rows", nrow(cog_meta.dt))

# Fill Units
invisible({
  cog_meta.dt[TITLE %like% "incorrect", UNITS := "error-count"]
  cog_meta.dt[TITLE %like% "digits|letters", UNITS := "correct-answers"]
  cog_meta.dt[TITLE %like% "viewed", UNITS := "attempts"]
  cog_meta.dt[COLNAME %in% c("VOCAB_lv", "PRMEM_res"), UNITS := "score"]
  cog_meta.dt[TITLE %like% "Uncertainty", UNITS := "uncertainty"]
  cog_meta.dt[TITLE %like% "completion", UNITS := "completion"]
})

# Disassociate online from COLNAME
invisible({
  cog_meta.dt[, ONLINE := fifelse(COLNAME %like% "onl", "UNIT_online", "UNIT")]
  cog_meta.dt[, COLNAME := str_remove(COLNAME, "_onl")]
  cog_meta.dt <- dcast(cog_meta.dt, ... ~ ONLINE, value.var = "UNITS")
})

# Standardizing time units
invisible({
  cog_meta.dt[, CONV_time := !(is.na(UNIT_online) | UNIT == UNIT_online)]
  cog_meta.dt[is.na(CONV_time), CONV_time := FALSE]
})

t_conv.dt <-
  cog_meta.dt[UNIT != UNIT_online, .SD, .SDcols = patterns("COL|UNI")] |>
  melt(id = "COLNAME") |>
  (\(DT) {
    DT[
      , let(
        ONLINE = fifelse(variable %like% "onl", TRUE, FALSE),
        CONV_time = fcase(
          value %like% "dec", .1,
          value %like% "mil", .001,
          default = 1
        )
      )
    ][
      , c("variable", "value") := NULL
    ]
  })()

# Pattern from column name
cog_meta.dt[, PATTERN := str_extract(COLNAME, "^[A-Z]*")] |> invisible()

# Test names and COL pattern
cog_meta.dt <- data.table(
  NAME = c(
    "Pair_matching",
    "Fluid_intelligence",
    "Trail_making",
    "Symbol_digit_substitution",
    "Numeric_memory",
    "Matrix_pattern_recognition",
    "Reaction_time",
    "Tower_rearranging",
    "Prospective_memory"
  ),
  PATTERN = c(
    "PRS",
    "FLINT",
    "TRLS",
    "SYM",
    "NUM",
    "MATS",
    "REACT",
    "TOWER",
    "PRMEM"
  )
)[cog_meta.dt, on = "PATTERN"]

# Create list with data.tables for cognitive scores
pb <- progress_bar$new(
  format = "Processing tests | :what [:bar] :current/:total",
  total  = cog_meta.dt[!duplicated(NAME), .N],
  clear  = FALSE,
  width  = 75
)

cog_tests.lst <- list()
for (cog_test in unique(cog_meta.dt$NAME)) {
  test_name <- str_replace_all(cog_test, "_", " ")
  pb$tick(tokens = list(what = test_name))
  rm(test_name)

  # Get specific metadata
  submetadata <- cog_meta.dt[
    cog_test,
    on = "NAME",
    .(
      PATTERN = unique(PATTERN),
      INSTANCED = any(INSTANCED),
      ARRAYED = any(ARRAYED),
      ONLINE = any(!is.na(UNIT_online)),
      CONV_T = any(CONV_time)
    )
  ]

  # Get column names
  columnames <- cog_meta.dt[cog_test, on = "NAME", COLNAME]

  # Extract preliminary data.table
  prelim.subdt <- ext_cols.fn(covars.dt, submetadata$PATTERN)

  # Disaggregate columns by type
  cols_num <- prelim.subdt[, names(.SD), .SDcols = is.numeric]
  cols_chr <- c("EID", prelim.subdt[, names(.SD), .SDcols = is.character])
  cols_fct <- c("EID", prelim.subdt[, names(.SD), .SDcols = is.factor])
  sublist <- list(
    prelim.subdt[, ..cols_num],
    prelim.subdt[, ..cols_chr],
    prelim.subdt[, ..cols_fct]
  )

  # Remove empty data.tables
  sublist <- sublist[sapply(sublist, \(DT) length(DT) > 1)]

  # Extract sessions
  if (submetadata$INSTANCED) sublist <- lapply(sublist, ext_sess.fn)

  # Extract arrays
  if (submetadata$ARRAYED) {
    sublist <- lapply(
      sublist,
      ext_arrays.fn,
      include_sess = submetadata$INSTANCED
    )
  }

  # Separate online tests
  if (submetadata$ONLINE) {
    sublist <- lapply(
      sublist,
      sep_online.fn,
      include_sess = submetadata$INSTANCED,
      include_array = submetadata$ARRAYED
    )
  } else {
    lapply(
      sublist,
      \(data_table) {
        data_table[, ONLINE := FALSE]
        setcolorder(data_table, "ONLINE", after = "SESSION")
        setkey(data_table, "EID", "SESSION", "ONLINE")
      }
    )
  }

  # Standardize time
  if (submetadata$CONV_T) {
    for (columname in columnames) {
      if (columname %in% t_conv.dt$COLNAME) {
        t_conv.subdt <- t_conv.dt[columname, on = "COLNAME"]
        sublist[[1]][, (columname) := as.numeric(get(columname))]
        # ONLINE
        if (sublist[[1]][ONLINE == TRUE, .N > 0]) {
          sublist[[1]][
            ONLINE == TRUE,
            (columname) := .SD * t_conv.subdt[ONLINE == TRUE, CONV_time],
            .SDcols = columname
          ]
        }
        # OFFLINE
        if (sublist[[1]][ONLINE == FALSE, .N > 0]) {
          sublist[[1]][
            ONLINE == FALSE,
            (columname) := .SD * t_conv.subdt[ONLINE == FALSE, CONV_time],
            .SDcols = columname
          ]
        }
      }
    }
  }

  # Merge back data.tables
  prelim.subdt <- Reduce(
    \(data_table1, data_table2) merge(data_table1, data_table2, all = TRUE),
    sublist
  )

  # Remove rows where all test columns are NA
  cog_tests.lst[[cog_test]] <- prelim.subdt[
    !(rowSums(is.na(prelim.subdt[, ..columnames])) == length(columnames))
  ]
}

rm(
  prelim.subdt, sublist, cog_test, cols_chr, cols_fct, cols_num,
  columname, columnames, submetadata, t_conv.subdt, t_conv.dt
)

# Parse notes
notes.dt <-
  metadata.dt[
    , FIELD,
    keyby = "COLNAME"
  ][
    cog_meta.dt[, NAME, keyby = "COLNAME"]
  ][
    notes.dt,
    on = "FIELD", .(NAME, COLNAME, NOTES), nomatch = NULL
  ]

# Clinical history
# Filter out people with ICD_10 specific Dx
dx_pattern <- get_config("parameters", "exclusions", "icd10_patterns")
icd_10.dt <- covars.dt[ICD_10 %like% dx_pattern, .(EID)] |> setkey(EID)
cog_tests.lst <- lapply(cog_tests.lst, function(DT) DT[!icd_10.dt])
rm(icd_10.dt, dx_pattern)

# ----- Specific Cleaning -----
# Matrices (Abstract reasoning/Problem solving)
cog_tests.lst$Matrix_pattern_recognition <- merge(
  x = cog_tests.lst$Matrix[
    !is.na(ARRAY),
    .(MATS_time_mean = -mean(MATS_time)),
    keyby = .(EID, SESSION, ONLINE)
  ],
  y = cog_tests.lst$Matrix[
    is.na(ARRAY),
    .(MATS_corr_try = MATS_corr / MATS_try, MATS_corr, MATS_try),
    keyby = .(EID, SESSION, ONLINE)
  ]
)

# Trail making (Processing speed/executive function)
cols <- cog_tests.lst$Trail[, names(.SD), .SDcols = is.numeric][-1]
for (column in cols) {
  cog_tests.lst[["Trail_making"]][is.na(get(column)), (column) := 0]
}
rm(column, cols)

cog_tests.lst$Trail[
  is.na(TRLS_complete),
  TRLS_complete := fcase(
    TRLS_alnum_time == 0 | TRLS_num_time == 0, "Abandoned",
    default = "Completed"
  )
]

cog_tests.lst$Trail[
  , let(
    TRLS_complete_ord = factor(
      TRLS_complete,
      levels = c("Abandoned", "Timed-out due to inactivity", "Completed"),
      ordered = TRUE
    ),
    TRLS_complete_bin = fifelse(TRLS_complete == "Completed", 1, 0)
  )
][
  , TRLS_complete := NULL
]

# Reverse time
cog_tests.lst$Trail[
  TRLS_num_time > 0,
  TRLS_num_time := -TRLS_num_time,
  .(EID, SESSION, ONLINE)
]

cog_tests.lst$Trail[
  TRLS_alnum_time > 0,
  TRLS_alnum_time := -TRLS_alnum_time,
  .(EID, SESSION, ONLINE)
]

# Reaction time (Processing speed)
cog_tests.lst$Reaction[, REACT := -REACT]

# Pair matching (Memory/Visual-spatial memory)
cog_tests.lst$Pair_matching[
  PRS_time > 0,
  PRS_inc_time := -mean(PRS_inc / PRS_time),
  .(EID, SESSION, ONLINE)
]

cog_tests.lst$Pair_matching[
  ,
  PRS_wt_inc_time := -sum(PRS_inc * PRS_time) / sum(PRS_time),
  .(EID, SESSION, ONLINE)
]

cog_tests.lst$Pair_matching[
  ,
  PRS_wt_time_inc := -sum(PRS_inc * PRS_time) / sum(PRS_inc),
  .(EID, SESSION, ONLINE)
]

cog_tests.lst$Pair_matching[
  ,
  let(PRS_mean_inc = -mean(PRS_inc), PRS_mean_time = -mean(PRS_time)),
  .(EID, SESSION, ONLINE)
]

cog_tests.lst$Pair_matching[, c("ARRAY", "PRS_inc", "PRS_time") := NULL]
cog_tests.lst$Pair_matching <- unique(cog_tests.lst$Pair_matching)
setkey(cog_tests.lst$Pair_matching, "EID", "SESSION", "ONLINE")

# Tower rearranging (Planning/Executive function)
cog_tests.lst$Tower[, TOWER_corr_try := TOWER_corr / TOWER_try]

# Symbol-digit substitution (Processing speed/Attention)
cog_tests.lst$Symbol[, SYM_corr_try := SYM_corr / SYM_try]

# Prospective memory (Memory/Prospective memory)
cog_tests.lst$Prospective[
  ,
  PRMEM_res_n := fcase(
    PRMEM_res %like% "incorrect", 0,
    PRMEM_res %like% "second", 1,
    PRMEM_res %like% "first", 2
  )
]

# ----- Outputs -----
log_section("Saving outputs")

# Validate outputs before saving
validate_not_empty(ses.dt, "SES data")
validate_columns(ses.dt, c("EID", "COUNTRY"), "SES data")
log_info("SES data: %d rows across %d countries",
  nrow(ses.dt), ses.dt[, uniqueN(COUNTRY)]
)

validate_not_empty(cog_meta.dt, "cognitive metadata")
log_info("Cognitive metadata: %d test variables", nrow(cog_meta.dt))

n_cog_tests <- length(cog_tests.lst)
if (n_cog_tests == 0) {
  log_error("No cognitive tests processed")
  stop("Failed to process cognitive tests", call. = FALSE)
}
log_info("Cognitive tests: %d tests processed", n_cog_tests)

write_rds_safe(
  ses.dt,
  get_data_path("processed", "ses_data"),
  description = "SES data"
)

write_rds_safe(
  cog_meta.dt,
  get_data_path("processed", "cog_metadata"),
  description = "Cognitive test metadata"
)

write_rds_safe(
  cog_tests.lst,
  get_data_path("processed", "cog_tests"),
  description = "Cognitive test data"
)

log_script_end("04_clean_ses_cognitive.R", success = TRUE)
