#!/usr/bin/env Rscript

# =============================================================================
# Demographics Table data
# =============================================================================
# Generates demographic comparison tables for
# Cross-sectional and Longitudinal samples
#
# Inputs:
#   - data/fst/ukb_covars.fst: UK Biobank covariates
#   - data/derivatives/hc-hvr_adj.rds: Head-size adjusted data
#   - data/derivatives/lat-cog_vals.rds: Cognitive latent factor scores
#
# Outputs:
#   - data/derivatives/demog_data.rds
#   - data/derivatives/demog_t-tests.rds
#   - outputs/tables/demographics/*.{html,tex}: Publication tables
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
library(gt)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))
source(here("R/utils/tables_core.R"))
source(here("R/utils/export.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("07_demographics.R")

validate_packages(c("data.table", "gt"))

# Load configuration
config <- load_config()
set_seed()

# ----- Constants -----
NETWORK <- get_script_setting("demographics", "network")
log_info("Network: %s", NETWORK)

BILATERAL <- get_script_setting("demographics", "bilateral")
log_info(fifelse(
  BILATERAL,
  "Bilateral: tables will contain data for L/R hemispheres",
  "Bilateral: tables will show total volumes (sum of both hemispheres)"
))

# ----- Input Files -----
log_section("Loading input data")

# Define input paths
covars.path <- get_data_path("processed", "covars_fst")
hc_hvr.path <- get_data_path("processed", "hc_hvr_adjusted")
cog_scores.path <- get_data_path("processed", "lat-cog_values")

# Check and load covariates
if (!check_files_exist(covars.path, stop_on_missing = FALSE)) {
  log_warn("Covariates file not found, running 03_parse_covariates.R")
  source(here("R", "03_parse_covariates.R"))
}
covars.dt <- read_fst_safe(
  covars.path,
  as_data_table = TRUE, description = "Covariates"
)
validate_not_empty(covars.dt, "covariates")
validate_columns(covars.dt, c("EID", "EDUC_lv_ses2_0"), "covariates")
rm(covars.path)

# Check and load head-size adjusted data
if (!check_files_exist(hc_hvr.path, stop_on_missing = FALSE)) {
  log_warn("Head-size adjusted data not found, running 05_adjust_headsize.R")
  source(here("R", "05_adjust_headsize.R"))
}
hc_hvr.lst <- read_rds_safe(
  hc_hvr.path,
  description = "Head-size adjusted data"
)
validate_not_empty(hc_hvr.lst, "head-size adjusted data")
if (!NETWORK %in% names(hc_hvr.lst)) {
  log_error("Network '%s' not found. Available: %s", NETWORK, paste(names(hc_hvr.lst), collapse = ", "))
  stop(sprintf("Invalid network specification: %s", NETWORK), call. = FALSE)
}
rm(hc_hvr.path)

# Check and load cognitive scores
if (!check_files_exist(cog_scores.path, stop_on_missing = FALSE)) {
  log_warn("Cognitive scores not found, running 06_cognitive_factors.R")
  source(here("R", "06_cognitive_factors.R"))
}
lat_cog.lst <- read_rds_safe(
  cog_scores.path,
  description = "Cognitive factor scores"
)
validate_not_empty(lat_cog.lst, "cognitive scores")
rm(cog_scores.path)

# ----- Data Cleaning -----
log_section("Preparing demographics data")

# Get education mapping from config (centralized in config.R)
edu_years_map <- get_education_years_map()
educ.dt <- covars.dt[, .(EDUC_char = as.character(EDUC_lv_ses2_0)), EID]
educ.dt[, EDUC_num := edu_years_map[EDUC_char]] |> invisible()
setkey(educ.dt, EID)

# Total and bilateral HC/LV (unadjusted) & merge with cognitive factors
demog.lst <- Map(
  \(crs_lvl.lst, cog.dt) {
    lapply(crs_lvl.lst, \(mtch_lvl) {
      group_cols.v <- c("EID", "MATCH", "SEX", "INST", "AGE", "ICC")
      if (!"MATCH" %in% names(mtch_lvl)) group_cols.v <- group_cols.v[-2]

      # Whether to use bilateral or total volumes
      rois.v <- c("HC", "LV")
      if (BILATERAL) {
        hc_hvr.subdt <- mtch_lvl[!"LR", on = "SIDE", -"HVR"] |>
          melt(measure = rois.v) |>
          dcast(... ~ variable + SIDE)
        rois.v <- hc_hvr.subdt[, names(.SD), .SDcols = patterns("(L|R)$")]
      } else {
        hc_hvr.subdt <- mtch_lvl["LR", on = "SIDE", -c("SIDE", "HVR")]
      }

      # If subfield data available, show only Total volumes
      if ("SUBFIELD" %in% names(hc_hvr.subdt)) {
        hc_hvr.subdt <- hc_hvr.subdt["total", on = "SUBFIELD", -"SUBFIELD"]
      }

      setkey(hc_hvr.subdt, ADJ)

      hc_hvr.subdt["NON", .SD, by = group_cols.v, .SDcols = rois.v] |>
        setkey(EID, INST) |>
        merge(cog.dt[, -"SEX"], all.x = TRUE) |>
        merge(educ.dt[, -"EDUC_char"], all.x = TRUE) |>
        setcolorder("EDUC_num", after = "AGE") |>
        setcolorder("COG", before = "MEM")
    })
  }, hc_hvr.lst[[NETWORK]], lat_cog.lst
)
rm(covars.dt, hc_hvr.lst, lat_cog.lst)

demog.lst$LNG <- demog.lst$LNG$ALL

# ----- England Subsample -----
# Load SES data to get COUNTRY for England filtering
# (needed for SEM analyses which use IMDP, only comparable within England)
ses_path <- get_data_path("processed", "ses_data")
if (check_files_exist(ses_path, stop_on_missing = FALSE)) {
  ses.dt <- read_rds_safe(ses_path, "SES data for country filtering")
  eng_eids <- ses.dt[COUNTRY == "ENG", EID]
  demog.lst$CRS$ENG <- demog.lst$CRS$ALL[EID %in% eng_eids]
  log_info("England subsample: %d subjects (%.1f%% of full sample)",
           nrow(demog.lst$CRS$ENG),
           100 * nrow(demog.lst$CRS$ENG) / nrow(demog.lst$CRS$ALL))
  rm(ses.dt, eng_eids, ses_path)
} else {
  log_warn("SES data not found - England subsample not created")
}

# Group comparisons
ttests.lst <- list()
for (s in c("ALL", "MTCH", "ENG")) {
  if (!s %in% names(demog.lst$CRS)) next
  ttests.lst$CRS[[s]] <- data.table()
  for (v in names(demog.lst$CRS[[s]][, .SD, .SDcols = AGE:PRSP])) {
    tt <- demog.lst$CRS[[s]][, t.test(get(v) ~ SEX, na.rm = TRUE)]
    ttests.lst$CRS[[s]] <- rbind(
      ttests.lst$CRS[[s]],
      data.table(
        X = v,
        Tstat = tt$statistic,
        DF = tt$parameter,
        Pval = tt$p.value
      )
    ) |> invisible()
  }
}
rm(s, v, tt)

# Longitudinal sample
demog.lst$LNG[
  , TIME_diff := floor((AGE[INST == "ses-3"] - AGE[INST == "ses-2"]) * 12), EID
] |> invisible()

setcolorder(demog.lst$LNG, "TIME_diff", after = "AGE")

ttests.lst[["LNG"]] <- data.table()
fu_cols <- names(demog.lst$LNG[, .SD, .SDcols = ICC:PRSP])
for (v in names(demog.lst$LNG[, .SD, .SDcols = AGE:TIME_diff])) {
  tt <- demog.lst$LNG[INST %like% "2", t.test(get(v) ~ SEX, na.rm = TRUE)]
  ttests.lst$LNG <- rbind(
    ttests.lst$LNG,
    data.table(
      X = v,
      INST = "BL",
      Tstat = tt$statistic,
      DF = tt$parameter,
      Pval = tt$p.value
    )
  )
  if (v %in% fu_cols) {
    tt <- demog.lst$LNG[INST %like% "3", t.test(get(v) ~ SEX, na.rm = TRUE)]
    ttests.lst$LNG <- rbind(
      ttests.lst$LNG,
      data.table(
        X = v,
        INST = "FU",
        Tstat = tt$statistic,
        DF = tt$parameter,
        Pval = tt$p.value
      )
    )
  }
}
rm(v, tt)

ttests.lst$LNG["TIME_diff", on = "X", INST := "FU"] |> invisible()

# ----- Save Results for Manuscript -----
log_section("Saving outputs")

# Validate outputs before saving
validate_not_empty(demog.lst$CRS$ALL, "cross-sectional demographics (ALL)")
validate_not_empty(demog.lst$CRS$MTCH, "cross-sectional demographics (MTCH)")
if ("ENG" %in% names(demog.lst$CRS)) {
  validate_not_empty(demog.lst$CRS$ENG, "cross-sectional demographics (ENG)")
}
validate_not_empty(demog.lst$LNG, "longitudinal demographics")
validate_not_empty(ttests.lst$CRS$ALL, "cross-sectional t-tests (ALL)")
validate_not_empty(ttests.lst$CRS$MTCH, "cross-sectional t-tests (MTCH)")
if ("ENG" %in% names(ttests.lst$CRS)) {
  validate_not_empty(ttests.lst$CRS$ENG, "cross-sectional t-tests (ENG)")
}
validate_not_empty(ttests.lst$LNG, "longitudinal t-tests")

log_info("Cross-sectional ALL: %d subjects", demog.lst$CRS$ALL[!duplicated(EID), .N])
if ("ENG" %in% names(demog.lst$CRS)) {
  log_info("Cross-sectional ENG: %d subjects", demog.lst$CRS$ENG[!duplicated(EID), .N])
}
log_info("Cross-sectional MTCH: %d subjects", demog.lst$CRS$MTCH[!duplicated(EID), .N])
log_info("Longitudinal: %d subjects", demog.lst$LNG[!duplicated(EID), .N])

# Save demographics data and t-test results for use in manuscript.qmd
write_rds_safe(
  demog.lst,
  get_data_path("demographics", "data"),
  description = "Demographics data"
)

write_rds_safe(
  ttests.lst,
  get_data_path("demographics", "t-tests"),
  description = "Demographics t-test results"
)

log_info("Demographics data saved")

# ===========================================================================
# Create Publication-Ready Tables
# ===========================================================================
log_section("Creating publication tables")

# Get settings
create_tables <- get_script_setting(
  "demographics", "create_tables", default = TRUE
)

if (create_tables) {
  tables_dir <- get_output_path("tables")
  ensure_directory(tables_dir)

  # ----- Cross-Sectional Demographics (ALL) -----
  log_info("Creating cross-sectional demographics table (ALL)")

  crs_all.dt <- demog.lst$CRS$ALL
  crs_all_ttest.dt <- ttests.lst$CRS$ALL

  # Calculate summary statistics by sex
  # Note: Current bifactor model uses COG (g-factor), MEM (memory), PRSP (processing speed)
  summ_cols <- c("AGE", "EDUC_num", "ICC", "HC_L", "HC_R", "LV_L", "LV_R",
    "COG", "MEM", "PRSP")

  # Filter to available columns
  summ_cols <- summ_cols[summ_cols %in% names(crs_all.dt)]

  crs_summ.dt <- crs_all.dt[,
    {
      res <- list(SEX = "Total", N = .N)
      for (col in summ_cols) {
        if (col %in% names(.SD)) {
          res[[col]] <- format_mean_sd(
            get(col),
            is_icc = col == "ICC",
            use_cog_format = col %in% c("COG", "MEM", "PRSP")
          )
        }
      }
      res
    }
  ]

  for (sex in c("Female", "Male")) {
    sex_summ <- crs_all.dt[SEX == sex,
      {
        res <- list(SEX = sex, N = .N)
        for (col in summ_cols) {
          if (col %in% names(.SD)) {
            res[[col]] <- format_mean_sd(
              get(col),
              is_icc = col == "ICC",
              use_cog_format = col %in% c("COG", "MEM", "PRSP")
            )
          }
        }
        res
      }
    ]
    crs_summ.dt <- rbind(crs_summ.dt, sex_summ)
  }

  # Add p-values
  crs_summ.dt[, P := ""]
  for (var in summ_cols) {
    if (var %in% crs_all_ttest.dt$X) {
      pval <- crs_all_ttest.dt[X == var, Pval]
      crs_summ.dt[SEX == "Male", (var) := paste0(get(var), " ",
        format_pval(pval, stars = TRUE))]
    }
  }

  # Create gt table
  crs_all_gt_base <- crs_summ.dt |>
    gt(rowname_col = "SEX") |>
    tab_header(
      title = "Demographics: Cross-Sectional Sample (All Subjects)",
      subtitle = sprintf("N = %s", format(.N <- crs_all.dt[, uniqueN(EID)],
        big.mark = ","))
    ) |>
    cols_label(
      N = "N",
      AGE = "Age (years)",
      EDUC_num = "Education (years)",
      ICC = "TIV (cc)",
      HC_L = "HC Left (cc)",
      HC_R = "HC Right (cc)",
      LV_L = "LV Left (cc)",
      LV_R = "LV Right (cc)",
      COG = "g-Factor",
      MEM = "Memory",
      PRSP = "Processing Speed"
    ) |>
    tab_source_note("Values are Mean (SD); TIV = Total Intracranial Volume; HC = Hippocampus; LV = Lateral Ventricles; *** p < 0.001, ** p < 0.01, * p < 0.05")

  # HTML version
  crs_all_gt_html <- crs_all_gt_base |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) |>
    tab_options(
      table.font.size = px(11),
      heading.background.color = "#4A90E2",
      heading.title.font.size = px(14),
      column_labels.font.weight = "bold",
      column_labels.background.color = "#F0F0F0",
      data_row.padding = px(3),
      table.border.top.style = "solid",
      table.border.bottom.style = "solid"
    )

  gt::gtsave(crs_all_gt_html,
    file.path(tables_dir, "demographics_crosssectional_all.html"))

  # LaTeX version
  crs_all_gt_tex <- crs_all_gt_base |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) |>
    tab_options(
      latex.use_longtable = TRUE,
      table.font.size = px(10),
      table.border.top.style = "solid",
      table.border.bottom.style = "solid",
      column_labels.border.top.style = "solid",
      column_labels.border.bottom.style = "solid"
    )

  tex_path <- file.path(tables_dir, "demographics_crosssectional_all.tex")
  gt::gtsave(crs_all_gt_tex, tex_path)
  fix_latex_source_notes(tex_path, table_type = "summary")
  wrap_latex_table(tex_path)

  log_info("  - Cross-sectional ALL table created")

  # ----- Cross-Sectional Demographics (MTCH) -----
  log_info("Creating cross-sectional demographics table (MTCH)")

  crs_mtch.dt <- demog.lst$CRS$MTCH
  crs_mtch_ttest.dt <- ttests.lst$CRS$MTCH

  # Calculate summary statistics
  crs_mtch_summ.dt <- crs_mtch.dt[,
    {
      res <- list(SEX = "Total", N = .N)
      for (col in summ_cols) {
        if (col %in% names(.SD)) {
          res[[col]] <- format_mean_sd(
            get(col),
            is_icc = col == "ICC",
            use_cog_format = col %in% c("COG", "MEM", "PRSP")
          )
        }
      }
      res
    }
  ]

  for (sex in c("Female", "Male")) {
    sex_summ <- crs_mtch.dt[SEX == sex,
      {
        res <- list(SEX = sex, N = .N)
        for (col in summ_cols) {
          if (col %in% names(.SD)) {
            res[[col]] <- format_mean_sd(
              get(col),
              is_icc = col == "ICC",
              use_cog_format = col %in% c("COG", "MEM", "PRSP")
            )
          }
        }
        res
      }
    ]
    crs_mtch_summ.dt <- rbind(crs_mtch_summ.dt, sex_summ)
  }

  # Add p-values
  for (var in summ_cols) {
    if (var %in% crs_mtch_ttest.dt$X) {
      pval <- crs_mtch_ttest.dt[X == var, Pval]
      crs_mtch_summ.dt[SEX == "Male", (var) := paste0(get(var), " ",
        format_pval(pval, stars = TRUE))]
    }
  }

  # Create gt table
  crs_mtch_gt_base <- crs_mtch_summ.dt |>
    gt(rowname_col = "SEX") |>
    tab_header(
      title = "Demographics: Cross-Sectional Sample (Age-Matched)",
      subtitle = sprintf("N = %s", format(.N <- crs_mtch.dt[, uniqueN(EID)],
        big.mark = ","))
    ) |>
    cols_label(
      N = "N",
      AGE = "Age (years)",
      EDUC_num = "Education (years)",
      ICC = "TIV (cc)",
      HC_L = "HC Left (cc)",
      HC_R = "HC Right (cc)",
      LV_L = "LV Left (cc)",
      LV_R = "LV Right (cc)",
      COG = "g-Factor",
      MEM = "Memory",
      PRSP = "Processing Speed"
    ) |>
    tab_source_note("Values are Mean (SD); TIV = Total Intracranial Volume; HC = Hippocampus; LV = Lateral Ventricles; *** p < 0.001, ** p < 0.01, * p < 0.05")

  # HTML version
  crs_mtch_gt_html <- crs_mtch_gt_base |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) |>
    tab_options(
      table.font.size = px(11),
      heading.background.color = "#4A90E2",
      heading.title.font.size = px(14),
      column_labels.font.weight = "bold",
      column_labels.background.color = "#F0F0F0",
      data_row.padding = px(3),
      table.border.top.style = "solid",
      table.border.bottom.style = "solid"
    )

  gt::gtsave(crs_mtch_gt_html,
    file.path(tables_dir, "demographics_crosssectional_matched.html"))

  # LaTeX version
  crs_mtch_gt_tex <- crs_mtch_gt_base |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) |>
    tab_options(
      latex.use_longtable = TRUE,
      table.font.size = px(10),
      table.border.top.style = "solid",
      table.border.bottom.style = "solid",
      column_labels.border.top.style = "solid",
      column_labels.border.bottom.style = "solid"
    )

  tex_path <- file.path(tables_dir, "demographics_crosssectional_matched.tex")
  gt::gtsave(crs_mtch_gt_tex, tex_path)
  fix_latex_source_notes(tex_path, table_type = "summary")
  wrap_latex_table(tex_path)

  log_info("  - Cross-sectional MTCH table created")

  # ----- Longitudinal Demographics -----
  log_info("Creating longitudinal demographics table")

  lng.dt <- demog.lst$LNG
  lng_ttest.dt <- ttests.lst$LNG

  # Baseline and Follow-up summaries
  # Note: Current bifactor model uses COG (g-factor), MEM (memory), PRSP (processing speed)
  lng_summ_cols <- c("AGE", "TIME_diff", "EDUC_num", "ICC", "HC_L", "HC_R",
    "LV_L", "LV_R", "COG", "MEM", "PRSP")
  lng_summ_cols <- lng_summ_cols[lng_summ_cols %in% names(lng.dt)]

  # Baseline
  lng_bl_summ.dt <- lng.dt[INST == "ses-2",
    {
      res <- list(TIMEPOINT = "Baseline", SEX = "Total", N = uniqueN(EID))
      for (col in lng_summ_cols) {
        if (col %in% names(.SD) && col != "TIME_diff") {
          res[[col]] <- format_mean_sd(
            get(col),
            is_icc = col == "ICC",
            use_cog_format = col %in% c("COG", "MEM", "PRSP")
          )
        }
      }
      res
    }
  ]

  for (sex in c("Female", "Male")) {
    sex_summ <- lng.dt[INST == "ses-2" & SEX == sex,
      {
        res <- list(TIMEPOINT = "Baseline", SEX = sex, N = uniqueN(EID))
        for (col in lng_summ_cols) {
          if (col %in% names(.SD) && col != "TIME_diff") {
            res[[col]] <- format_mean_sd(
              get(col),
              is_icc = col == "ICC",
              use_cog_format = col %in% c("COG", "MEM", "PRSP")
            )
          }
        }
        res
      }
    ]
    lng_bl_summ.dt <- rbind(lng_bl_summ.dt, sex_summ)
  }

  # Add BL p-values
  for (var in lng_summ_cols) {
    if (var %in% lng_ttest.dt[INST == "BL", X]) {
      pval <- lng_ttest.dt[INST == "BL" & X == var, Pval]
      lng_bl_summ.dt[SEX == "Male", (var) := paste0(get(var), " ",
        format_pval(pval, stars = TRUE))]
    }
  }

  # Follow-up
  lng_fu_summ.dt <- lng.dt[INST == "ses-3",
    {
      res <- list(TIMEPOINT = "Follow-up", SEX = "Total", N = uniqueN(EID))
      for (col in lng_summ_cols) {
        if (col %in% names(.SD)) {
          res[[col]] <- format_mean_sd(
            get(col),
            is_icc = col == "ICC",
            use_cog_format = col %in% c("COG", "MEM", "PRSP")
          )
        }
      }
      res
    }
  ]

  for (sex in c("Female", "Male")) {
    sex_summ <- lng.dt[INST == "ses-3" & SEX == sex,
      {
        res <- list(TIMEPOINT = "Follow-up", SEX = sex, N = uniqueN(EID))
        for (col in lng_summ_cols) {
          if (col %in% names(.SD)) {
            res[[col]] <- format_mean_sd(
              get(col),
              is_icc = col == "ICC",
              use_cog_format = col %in% c("COG", "MEM", "PRSP")
            )
          }
        }
        res
      }
    ]
    lng_fu_summ.dt <- rbind(lng_fu_summ.dt, sex_summ)
  }

  # Add FU p-values
  for (var in lng_summ_cols) {
    if (var %in% lng_ttest.dt[INST == "FU", X]) {
      pval <- lng_ttest.dt[INST == "FU" & X == var, Pval]
      lng_fu_summ.dt[SEX == "Male", (var) := paste0(get(var), " ",
        format_pval(pval, stars = TRUE))]
    }
  }

  # Combine
  lng_summ.dt <- rbind(lng_bl_summ.dt, lng_fu_summ.dt, fill = TRUE)

  # Create gt table
  lng_gt_base <- lng_summ.dt |>
    gt(groupname_col = "TIMEPOINT", rowname_col = "SEX") |>
    tab_header(
      title = "Demographics: Longitudinal Sample",
      subtitle = sprintf("N = %s subjects",
        format(.N <- lng.dt[, uniqueN(EID)], big.mark = ","))
    ) |>
    cols_label(
      N = "N",
      AGE = "Age (years)",
      TIME_diff = "Time (months)",
      EDUC_num = "Education (years)",
      ICC = "TIV (cc)",
      HC_L = "HC Left (cc)",
      HC_R = "HC Right (cc)",
      LV_L = "LV Left (cc)",
      LV_R = "LV Right (cc)",
      COG = "g-Factor",
      MEM = "Memory",
      PRSP = "Processing Speed"
    ) |>
    tab_source_note("Values are Mean (SD); TIV = Total Intracranial Volume; HC = Hippocampus; LV = Lateral Ventricles; *** p < 0.001, ** p < 0.01, * p < 0.05")

  # HTML version
  lng_gt_html <- lng_gt_base |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) |>
    tab_options(
      table.font.size = px(11),
      heading.background.color = "#4A90E2",
      heading.title.font.size = px(14),
      column_labels.font.weight = "bold",
      column_labels.background.color = "#F0F0F0",
      data_row.padding = px(3),
      table.border.top.style = "solid",
      table.border.bottom.style = "solid"
    )

  gt::gtsave(lng_gt_html,
    file.path(tables_dir, "demographics_longitudinal.html"))

  # LaTeX version
  lng_gt_tex <- lng_gt_base |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) |>
    tab_options(
      latex.use_longtable = TRUE,
      table.font.size = px(10),
      table.border.top.style = "solid",
      table.border.bottom.style = "solid",
      column_labels.border.top.style = "solid",
      column_labels.border.bottom.style = "solid"
    )

  tex_path <- file.path(tables_dir, "demographics_longitudinal.tex")
  gt::gtsave(lng_gt_tex, tex_path)
  fix_latex_source_notes(tex_path, table_type = "summary")
  wrap_latex_table(tex_path)

  log_info("  - Longitudinal table created")
  log_info("Demographic tables saved to: %s", tables_dir)
}

log_script_end("07_demographics.R", success = TRUE)
