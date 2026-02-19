#!/usr/bin/env Rscript

# =============================================================================
# Cognitive Factor Analysis
# =============================================================================
# Performs confirmatory factor analysis on cognitive test data
#
# Inputs:
#   - data/derivatives/ses_data.rds: SES scores
#   - data/derivatives/cog-tests.rds: Cognitive test data
#   - data/derivatives/cog-tests_metadata.rds: Test metadata
#   - data/derivatives/hc-hvr_adj.rds: Head-size adjusted data
#
# Outputs:
#   - data/derivatives/cog_cfa_fits.rds: CFA model fits
#   - data/derivatives/cog_cfa_vals.rds: Latent cognitive factor scores
#   - data/derivatives/cog_cfa_msr-inv.rds: Measurement invariance fits
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
library(lavaan)
library(psych)  # For ICC reliability calculation

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("06_cognitive_factors.R")

# Load configuration
config <- load_config()
set_seed()

# ----- Constants -----
SES_SUBSET <- get_script_setting("cfa", "ses_subset")
REFIT_CFA <- get_script_setting("cfa", "refit_models")
log_info("SES subset (England only): %s", SES_SUBSET)
log_info("Refit CFA models: %s", REFIT_CFA)

# ----- Input Files -----
log_section("Loading input data")

ses.path <- get_data_path("processed", "ses_data")
cog_meta.path <- get_data_path("processed", "cog_metadata")
cog_tests.path <- get_data_path("processed", "cog_tests")
ses_cog_paths.v <- c(ses.path, cog_meta.path, cog_tests.path)
if (!check_files_exist(ses_cog_paths.v, stop_on_missing = FALSE)) {
  ses_cog_paths.v |>
    sapply(file.exists) |>
    (\(bool) names(bool)[!bool])() |>
    paste(collapse = " & ") |>
    log_warn(message = "File(s) not found: %s - running prerequisite script")
  source(here("R", "04_clean_ses_cognitive.R"))
}

if (!exists("ses.dt")) ses.dt <- read_rds_safe(ses.path, "SES data")
if (!exists("cog_tests.lst")) {
  cog_tests.lst <- read_rds_safe(cog_tests.path, "Cognitive tests")
}
if (!exists("cog_meta.dt")) {
  cog_meta.dt <- read_rds_safe(cog_meta.path, "Cognitive metadata")
}
rm(ses.path, cog_meta.path, cog_tests.path, ses_cog_paths.v)

hc_hvr.path <- get_data_path("processed", "hc_hvr_adjusted")
if (!check_files_exist(hc_hvr.path, stop_on_missing = FALSE)) {
  log_warn("File not found: %s - running prerequisite script", hc_hvr.path)
  source(here("R", "05_adjust_headsize.R"))
}
if (!exists("hc_hvr.lst")) {
  hc_hvr.lst <- read_rds_safe(hc_hvr.path, "Head-size adjusted data")
}
validate_not_empty(hc_hvr.lst, "head-size adjusted data")
rm(hc_hvr.path)

# ----- Settings -----
# Initialize CFA settings that are needed across conditional blocks
estim <- get_parameter("cfa", "estimator")
if (get_parameter("cfa", "scaled_scores")) {
  scaled <- "SCL"
} else {
  scaled <- "RAW"
}

# ----- Functions -----
ext_tests.fn <- function(
    cog_list.lst,
    cog_tests.i,
    id_vars.v = c("EID", "SESSION", "ONLINE"),
    scaling = FALSE) {
  cog_tests.dt <- cog_list.lst[cog_tests.i] |>
    lapply(\(data_table) {
      data_table[, EID := as.character(EID)]
      data_table[, .SD, .SDcols = is.numeric, keyby = id_vars.v] |>
        melt(id = id_vars.v) |>
        suppressWarnings()
    }) |>
    rbindlist() |>
    (\(data_table) data_table[, EID := as.integer(EID)])()

  if ("ONLINE" %in% id_vars.v) {
    cog_tests.dt[(ONLINE), variable := paste0(variable, "_online")]
    cog_tests.dt[, ONLINE := NULL]
  }

  cog_tests.dt <- cog_tests.dt[!is.na(value)] |>
    unique() |>
    dcast(... ~ variable)

  if (scaling) {
    cols <- names(cog_tests.dt)[!names(cog_tests.dt) %in% id_vars.v]
    if ("SESSION" %in% id_vars.v) {
      cog_tests.dt[, (cols) := lapply(.SD, scale), by = SESSION, .SDcols = cols]
    } else {
      cog_tests.dt[, (cols) := lapply(.SD, scale), .SDcols = cols]
    }
  }

  return(cog_tests.dt)
}

# ----- Data Preparation -----
log_section("Preparing cognitive data")
cog_clean.path <- here("data/derivatives/cog_clean-data.rds")

if (check_files_exist(cog_clean.path, stop_on_missing = FALSE)) {
  if (!exists("cog.lst")) {
    cog.lst <- read_rds_safe(cog_clean.path, "Cleaned cognitive data")
  }
} else {
  log_info("Extracting and cleaning cognitive tests")
  cog.lst <- list()
  cog.lst$Memory <- ext_tests.fn(cog_tests.lst, c(6, 4, 9))
  cog.lst$Proc_speed <- ext_tests.fn(cog_tests.lst, c(5, 2, 8))
  cog.lst$Reas_Exec <- ext_tests.fn(cog_tests.lst, c(1:3, 7))
  write_rds_safe(cog.lst, cog_clean.path, "Cleaned cognitive data")
}

data.lst <- tests.lst <- list()
# It's irrelevant which dataset to use, since both have the same subjects
data.lst$RAW$CRS <- hc_hvr.lst[[1]]$CRS$ALL[!duplicated(EID), .(EID, INST, SEX)]
data.lst$RAW$LNG <- hc_hvr.lst[[1]]$LNG$ALL[, .(EID, INST, SEX)] |> unique()

if (SES_SUBSET) {
  data.lst$RAW$CRS <- ses.dt[
    "ENG",
    on = "COUNTRY",
    .(EID)
  ][
    data.lst$RAW$CRS,
    on = "EID",
    nomatch = NULL
  ]
}

# Cognitive tests included in factors
tests.lst$MEM <- get_parameter("cognitive", "memory")
tests.lst$PRSP <- get_parameter("cognitive", "proc_speed")
tests.lst$EXEC <- get_parameter("cognitive", "exec_func")

# Extract data.tables for each cognitive domain
ids.v <- c("EID", "SESSION")
subcog.lst <- Map(
  \(items, data.subdt) {
    cols <- c(ids.v, items)
    data.subdt[
      SESSION %like% "2|3", ..cols
    ][
      , SESSION := sub("ses", "ses-", SESSION)
    ] |>
      setnames("SESSION", "INST") |>
      setkey(EID, INST)
  },
  tests.lst,
  cog.lst
)

data.lst$RAW <- lapply(
  data.lst$RAW,
  \(data.subdt) {
    Reduce(
      \(x, y) merge(x, y, all = TRUE, by = c("EID", "INST")),
      subcog.lst
    ) |> merge(data.subdt, all.y = TRUE)
  }
)

# Separate baseline & follow-up
data.lst$RAW$LNG <- data.lst$RAW$LNG |>
  setcolorder("SEX", after = "INST") |>
  melt(1:3) |>
  (\(data.subdt) {
    data.subdt[, .(
      EID,
      SEX,
      VAR = paste0(variable, "_t", fifelse(INST %like% "2", "1", "2")),
      value
    )]
  })() |>
  dcast(... ~ VAR, value.var = "value")

# Scaling
tests_cols <- unlist(tests.lst)

data.lst$SCL <- copy(data.lst$RAW)
data.lst$SCL$CRS[
  ,
  (tests_cols) := lapply(.SD, scale),
  .SDcols = tests_cols
] |> invisible()

tests_cols_lng <- tests_cols |>
  rep(each = 2) |>
  paste(1:2, sep = "_t")

data.lst$SCL$LNG[
  ,
  (tests_cols_lng) := lapply(.SD, scale),
  .SDcols = tests_cols_lng
] |> invisible()

# ----- CFA Modeling -----
# Measurement invariance by sex
log_section("Measurement invariance")
minv_fits.path <- get_data_path("models", "fit", "cfa-minv_cog")
minv_ran <- FALSE
if (any(!check_files_exist(minv_fits.path, stop_on_missing = FALSE), REFIT_CFA)) {
  log_info("Testing measurement invariance by sex")
  models.lst <- list()
  models.lst$CRS <- get_data_path("models", "definition", "cfa_cog")
  models.lst$LNG <- get_data_path("models", "definition", "cfa-lng_cog")
  models.lst |>
    unlist() |>
    check_files_exist(stop_on_missing = TRUE) |>
    invisible()

  ## Notes:
  # 3Factor model showed Nonpositive definite (NPD) covariance matrix
  # Memory and Executive Function were highly correlated r > .85
  # Bifactor model with Executive function items loading only on g
  mod <- models.lst$CRS |>
    readLines() |>
    paste(collapse = "\n")

  log_info("Estimator for CFA: %s", estim)
  log_info("Using %s cognitive scores", scaled)

  data.dt <- data.lst[[scaled]]$CRS

  invar.lst <- list(
    METRIC = "loadings",
    SCALAR = c("loadings", "intercepts"),
    STRICT = c("loadings", "intercepts", "residuals")
  )

  minv_fits.lst <- list()
  log_info("Estimating: Configural invariance by Sex")
  minv_fits.lst$CONFIG <- cfa(
    mod, data.dt,
    estimator = estim, missing = "fiml", group = "SEX"
  )
  if (!inspect(minv_fits.lst$CONFIG, "converged")) {
    msg <- "CFA estimation failed for Configural invariance between Sex"
    log_error(msg)
    stop(msg, call. = FALSE)
  }
  log_info("Estimating: Metric invariance by Sex")
  minv_fits.lst$METRIC <- cfa(
    mod, data.dt,
    estimator = estim, missing = "fiml",
    group = "SEX", group.equal = invar.lst$METRIC
  )
  if (!inspect(minv_fits.lst$METRIC, "converged")) {
    msg <- "CFA estimation failed for Metric invariance between Sex"
    log_error(msg)
    stop(msg, call. = FALSE)
  }
  log_info("Estimating: Scalar invariance by Sex")
  minv_fits.lst$SCALAR <- cfa(
    mod, data.dt,
    estimator = estim, missing = "fiml",
    group = "SEX", group.equal = invar.lst$SCALAR
  )
  if (!inspect(minv_fits.lst$SCALAR, "converged")) {
    msg <- "CFA estimation failed for Scalar invariance between Sex"
    log_error(msg)
    stop(msg, call. = FALSE)
  }
  log_info("Estimating: Strict invariance by Sex")
  minv_fits.lst$STRICT <- cfa(
    mod, data.dt,
    estimator = estim, missing = "fiml",
    group = "SEX", group.equal = invar.lst$STRICT
  )
  if (!inspect(minv_fits.lst$STRICT, "converged")) {
    msg <- "CFA estimation failed for Strict invariance between Sex"
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  log_info("Saving measurement invariance fits")
  write_rds_safe(minv_fits.lst, minv_fits.path, "Measurement invariance fits")
  minv_ran <- TRUE
  rm(mod)
} else {
  log_info("Loading existing measurement invariance fits")
  minv_fits.lst <- read_rds_safe(minv_fits.path, "Measurement invariance fits")
}

# Evaluating Measurement invariance results
minv_res.path <- get_data_path("models", "results", "cfa-minv_cog")
if (any(!check_files_exist(minv_res.path, stop_on_missing = FALSE), minv_ran)) {
  log_info("Evaluating measurement invariance test results")
  # Extract robust fit indices
  minv_test.dt <- rbindlist(
    lapply(names(minv_fits.lst), \(name) {
      m <- minv_fits.lst[[name]]
      data.table(
        LEVEL         = name,
        CFI_robust    = fitMeasures(m, "cfi.robust"),
        TLI_robust    = fitMeasures(m, "tli.robust"),
        RMSEA_robust  = fitMeasures(m, "rmsea.robust"),
        SRMR          = fitMeasures(m, "srmr")
      )
    })
  )

  invisible({
    # Compute differences (Δ) from previous level
    minv_test.dt[, DeltaCFI := c(NA, diff(CFI_robust))]
    minv_test.dt[, DeltaTLI := c(NA, diff(TLI_robust))]
    minv_test.dt[, DeltaRMSEA := c(NA, diff(RMSEA_robust))]
    minv_test.dt[, DeltaSRMR := c(NA, diff(SRMR))]

    # Get measurement invariance thresholds from config
    # References: Chen (2007); Cheung & Rensvold (2002)
    thresh_cfi <- get_parameter("cfa", "invariance_thresholds", "delta_cfi", default = 0.01)
    thresh_rmsea <- get_parameter("cfa", "invariance_thresholds", "delta_rmsea", default = 0.015)
    thresh_srmr <- get_parameter("cfa", "invariance_thresholds", "delta_srmr", default = 0.03)

    log_debug("Invariance thresholds: CFI=%s, RMSEA=%s, SRMR=%s",
              thresh_cfi, thresh_rmsea, thresh_srmr)

    # Flag invariance according to conventional rules
    # ΔCFI & ΔTLI ≤ thresh_cfi
    # ΔRMSEA ≤ thresh_rmsea
    # ΔSRMR ≤ thresh_srmr (typically 0.03 for Metric, 0.01 for Scalar)
    minv_test.dt[, HoldsCFI := c(TRUE, abs(DeltaCFI[-1]) <= thresh_cfi)]
    minv_test.dt[, HoldsRMSEA := c(TRUE, abs(DeltaRMSEA[-1]) <= thresh_rmsea)]
    minv_test.dt[, HoldsSRMR := c(TRUE, abs(DeltaSRMR[-1]) <= thresh_srmr)]

    # Optionally, create a summary column combining all rules
    minv_test.dt[, InvarianceHolds := HoldsCFI & HoldsRMSEA & HoldsSRMR]
    setkey(minv_test.dt, LEVEL)
  })

  # Format and save results with column names matching _common.R expectations
  minv_results.dt <- minv_test.dt[, .(
    MODEL = LEVEL,
    CFI = CFI_robust,
    TLI = TLI_robust,
    RMSEA = RMSEA_robust,
    SRMR = SRMR,
    DELTA_CFI = DeltaCFI,
    DELTA_TLI = DeltaTLI,
    DELTA_RMSEA = DeltaRMSEA,
    DELTA_SRMR = DeltaSRMR,
    InvarianceHolds = InvarianceHolds
  )]
  setkey(minv_results.dt, MODEL)

  log_info("Saving measurement invariance results")
  write_rds_safe(minv_results.dt, minv_res.path, "Measurement invariance results")
} else {
  log_info("Loading measurement invariance test results")
  minv_test.dt <- read_rds_safe(minv_res.path, "Measurement invariance results")
  # Handle backward compatibility: old files use MODEL, new use LEVEL
  if ("MODEL" %in% names(minv_test.dt) && !"LEVEL" %in% names(minv_test.dt)) {
    setnames(minv_test.dt, "MODEL", "LEVEL")
    setkey(minv_test.dt, LEVEL)
  }
}

# Use row matching instead of keyed indexing for robustness
minv.res <- fcase(
  minv_test.dt[LEVEL == "STRICT", InvarianceHolds], "Strict",
  minv_test.dt[LEVEL == "SCALAR", InvarianceHolds], "Scalar",
  minv_test.dt[LEVEL == "METRIC", InvarianceHolds], "Metric",
  default = NA
)

if (is.na(minv.res)) {
  log_info("Measurement invariance across groups (Sex) failed at all levels")
} else {
  log_info(
    "Measurement invariance across groups (Sex) held up to: %s level",
    minv.res
  )
}

# Final definite CFA specification
log_section("Final CFA modeling")
cfa_fits.path <- get_data_path("models", "fit", "cfa_cog")
minv_param <- get_parameter("cfa", "group_invariance")
if (minv.res != minv_param) {
  log_warn(
    paste(
      "Group invariance parameter (%s)",
      "is inconsistent with measurement invariance test result (%s)"
    ),
    minv_param, minv.res
  )
}

if (any(!check_files_exist(cfa_fits.path, stop_on_missing = FALSE), REFIT_CFA)) {
  log_info("Fitting CFA models")

  if (!exists("models.lst")) {
    # Get model file paths from config
    models.lst <- list()
    models.lst$CRS <- get_data_path("models", "definition", "cfa_cog")
    models.lst$LNG <- get_data_path("models", "definition", "cfa-lng_cog")
    models.lst |>
      unlist() |>
      check_files_exist(stop = TRUE) |>
      invisible()
  }

  # estim and scaled are now initialized at script start

  ref_fits.lst <- Map(
    \(model, data.subdt, data_name) {
      modfile <- model |>
        readLines() |>
        paste(collapse = "\n")
      cohort <- c(
        CRS = "Cross-sectional cohort",
        LNG = "Longitudinal cohort"
      )[data_name]
      log_info("Estimating CFA for %s", cohort)
      cfa(
        model = modfile,
        data = data.subdt,
        estimator = estim,
        missing = "fiml",
        group = "SEX",
        group.equal = invar.lst[toupper(minv_param)]
      )
    },
    models.lst,
    data.lst[[scaled]],
    names(data.lst[[scaled]])
  )

  if (ref_fits.lst |> sapply(inspect, "converged") |> (\(bool) any(!bool))()) {
    msg <- ref_fits.lst |>
      sapply(inspect, "converged") |>
      (\(bool) names(bool)[!bool])() |>
      paste(collapse = " & ") |>
      sprintf(fmt = "Model(s): %s did not converge")
    log_error(msg)
    stop(msg, call. = FALSE)
  }

  write_rds_safe(ref_fits.lst, cfa_fits.path, "CFA model fits")
} else {
  log_info("Loading existing CFA fits")
  ref_fits.lst <- read_rds_safe(cfa_fits.path, "CFA model fits")
}

# Extract latent factors
lat_cog.lst <- Map(
  \(data.subdt, cfa.fit) {
    if ("INST" %in% names(data.subdt)) {
      lat_scores <- lavPredict(cfa.fit) |> lapply(as.data.table)
      lat_cog.subdt <- rbind(
        cbind(data.subdt[SEX %like% "F", .(EID, INST, SEX)], lat_scores$Female),
        cbind(data.subdt[SEX %like% "M", .(EID, INST, SEX)], lat_scores$Male)
      ) |> na.omit()
      setnames(lat_cog.subdt, "g", "COG")
      setkey(lat_cog.subdt, EID, INST)
      lat_cog.subdt
    } else {
      lat_scores <- lavPredict(cfa.fit) |> lapply(as.data.table)
      lat_cog.subdt <- rbind(
        cbind(data.subdt[SEX %like% "F", .(EID, SEX)], lat_scores$Female),
        cbind(data.subdt[SEX %like% "M", .(EID, SEX)], lat_scores$Male)
      ) |> melt(id = 1:2)
      lat_cog.subdt[, c("variable", "INST") := tstrsplit(variable, "_t")]
      # Convert timepoint (1 or 2) to session label (ses-2 or ses-3)
      lat_cog.subdt[, INST := fifelse(as.integer(INST) == 1L, "ses-2", "ses-3")]
      lat_cog.subdt <- dcast(lat_cog.subdt, EID + INST + SEX ~ variable)
      setnames(lat_cog.subdt, "g", "COG")
      setkey(lat_cog.subdt, EID, INST)
      lat_cog.subdt
    }
  },
  data.lst[[scaled]],
  ref_fits.lst
)

# ----- Output -----
log_section("Saving outputs")

# Validate outputs before saving
validate_not_empty(lat_cog.lst$CRS, "cross-sectional cognitive scores")
validate_not_empty(lat_cog.lst$LNG, "longitudinal cognitive scores")
validate_columns(lat_cog.lst$CRS, c("EID", "SEX", "COG", "MEM", "PRSP"), "CRS cognitive scores")
validate_columns(lat_cog.lst$LNG, c("EID", "INST", "SEX", "COG", "MEM", "PRSP"), "LNG cognitive scores")

log_info("Cross-sectional: %d subjects with cognitive scores", lat_cog.lst$CRS[!duplicated(EID), .N])
log_info("Longitudinal: %d subjects with cognitive scores", lat_cog.lst$LNG[!duplicated(EID), .N])

# ----- Test-Retest Reliability -----
log_section("Cognitive Factor Test-Retest Reliability")

# Reshape longitudinal data to wide format for reliability computation
lng_wide.dt <- dcast(
  lat_cog.lst$LNG,
  EID + SEX ~ INST,
  value.var = c("COG", "MEM", "PRSP")
)

# Compute test-retest correlations and ICCs for each factor
cog_reliability.lst <- list()

for (factor in c("COG", "MEM", "PRSP")) {
  t1_col <- paste0(factor, "_ses-2")
  t2_col <- paste0(factor, "_ses-3")

  # Get complete cases for this factor
  complete_idx <- complete.cases(lng_wide.dt[, .SD, .SDcols = c(t1_col, t2_col)])
  n_complete <- sum(complete_idx)

  if (n_complete < 10) {
    log_warn("Insufficient data for %s reliability (n=%d)", factor, n_complete)
    next
  }

  t1_scores <- lng_wide.dt[complete_idx, get(t1_col)]
  t2_scores <- lng_wide.dt[complete_idx, get(t2_col)]

  # Pearson correlation
  r_test <- cor.test(t1_scores, t2_scores)

  # ICC(2,1) - two-way random, single measures, absolute agreement
  # Using psych::ICC for correct calculation and confidence intervals
  icc_mat <- cbind(t1_scores, t2_scores)
  icc_result <- psych::ICC(icc_mat, missing = FALSE, alpha = 0.05, lmer = FALSE)

  # Extract ICC(2,1) - "ICC2" row (two-way random, single measures)
  icc_value <- icc_result$results$ICC[2]
  icc_lower <- icc_result$results$`lower bound`[2]
  icc_upper <- icc_result$results$`upper bound`[2]

  # Mean difference and SD of differences (for Bland-Altman style reporting)
  diff_scores <- t2_scores - t1_scores
  mean_diff <- mean(diff_scores)
  sd_diff <- sd(diff_scores)

  # Time between assessments (if available in data)
  # For now, just report the number of subjects

  cog_reliability.lst[[factor]] <- data.table(
    FACTOR = factor,
    N = n_complete,
    R = r_test$estimate,
    R_CI_LOWER = r_test$conf.int[1],
    R_CI_UPPER = r_test$conf.int[2],
    R_PVALUE = r_test$p.value,
    ICC = icc_value,
    ICC_CI_LOWER = icc_lower,
    ICC_CI_UPPER = icc_upper,
    MEAN_DIFF = mean_diff,
    SD_DIFF = sd_diff
  )

  log_info(
    "%s: r = %.3f [%.3f, %.3f], ICC = %.3f [%.3f, %.3f], n = %d",
    factor, r_test$estimate, r_test$conf.int[1], r_test$conf.int[2],
    icc_value, icc_lower, icc_upper, n_complete
  )
}

cog_reliability.dt <- rbindlist(cog_reliability.lst)

# Interpret ICC values (Koo & Li, 2016 guidelines)
cog_reliability.dt[, INTERPRETATION := fcase(
  ICC < 0.50, "Poor",
  ICC < 0.75, "Moderate",
  ICC < 0.90, "Good",
  ICC >= 0.90, "Excellent"
)]

log_info("Test-retest reliability summary:")
log_info("  COG (g-factor): ICC = %.3f (%s)",
         cog_reliability.dt[FACTOR == "COG", ICC],
         cog_reliability.dt[FACTOR == "COG", INTERPRETATION])
log_info("  MEM (memory-specific): ICC = %.3f (%s)",
         cog_reliability.dt[FACTOR == "MEM", ICC],
         cog_reliability.dt[FACTOR == "MEM", INTERPRETATION])
log_info("  PRSP (processing speed-specific): ICC = %.3f (%s)",
         cog_reliability.dt[FACTOR == "PRSP", ICC],
         cog_reliability.dt[FACTOR == "PRSP", INTERPRETATION])

write_rds_safe(
  lat_cog.lst,
  get_data_path("processed", "lat-cog_values"),
  description = "Latent cognitive factor scores"
)

write_rds_safe(
  cog_reliability.dt,
  get_data_path("processed", "cog_reliability"),
  description = "Cognitive factor test-retest reliability"
)

log_script_end("06_cognitive_factors.R", success = TRUE)
