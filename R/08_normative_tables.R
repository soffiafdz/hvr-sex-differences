#!/usr/bin/env Rscript
# =============================================================================
# GAMLSS Normative Tables
# =============================================================================
# Fits GAMLSS models and generates normative centile tables for brain volumes
# with validation of model fit quality and longitudinal stability
#
# Inputs:
#   - data/derivatives/hc-hvr_adj.rds: Head-size adjusted volumes
#   - data/fst/ukb_covars.fst: Covariates (age, education, etc.)
#   - config/pipeline_config.yaml: Configuration parameters
#
# Outputs:
#   Models:
#   - models/fits/gamlss.rds: GAMLSS model fits & validation
#   - models/diagnostics/gamlss.rds: Diagnostic summaries
#
#   Data:
#   - data/derivatives/norm_tables.rds: Normative centile tables
#   - data/derivatives/z_scores.rds: Z-scores from GAMLSS models
#
#   Tables (HTML & LaTeX):
#   - outputs/tables/gamlss/normative/*.{html,tex}: Full age-specific centiles
#   - outputs/tables/gamlss/summary/summary_*.{html,tex}: Binned summary tables
#   - outputs/tables/gamlss/summary/*_standalone.tex: Compilable LaTeX documents
#
#   Plots:
#   - outputs/figures/gamlss/centiles/*.png: Centile curves
#   - outputs/figures/gamlss/validation/*.png: Calibration plots
#   - outputs/figures/gamlss/stability/*.png: Stability plots
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
library(progress)
library(gamlss)
library(gamlss.add)
library(rsample)
library(gt)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))
source(here("R/utils/tables_core.R"))
source(here("R/utils/tables_normative.R"))
source(here("R/utils/plotting_core.R"))
source(here("R/utils/plotting_pipeline.R"))
source(here("R/utils/export.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("08_normative_tables.R")

# Validate required packages
validate_packages(c("gamlss", "gamlss.add", "rsample", "gt", "data.table"))

# Load configuration
config <- load_config()
set_seed()

# ----- Constants -----
REDO_TABLES <- get_script_setting("gamlss", "redo_tables", default = TRUE)
REDO_PLOTS <- get_script_setting("gamlss", "redo_plots", default = TRUE)
log_info("Redo tables: %s, Redo plots: %s", REDO_TABLES, REDO_PLOTS)

NETWORK <- get_script_setting("gamlss", "network")
log_info("Network: %s", NETWORK)

ADJS <- get_parameter("adjustment_methods")
log_info("Head-size adjustment methods: %s", paste(ADJS, collapse = ", "))

ROIS <- get_parameter("rois")
log_info("Regions of Interest: %s", paste(ROIS, collapse = ", "))

# ----- Input Files -----
log_section("Loading input data")

# Check required input files exist
covars.path <- get_data_path("processed", "covars_fst")
hc_hvr.path <- get_data_path("processed", "hc_hvr_adjusted")

if (!check_files_exist(covars.path, stop_on_missing = FALSE)) {
  log_warn("Covariates file not found, running 03_parse_covariates.R")
  source(here("R", "03_parse_covariates.R"))
}
covars.dt <- read_fst_safe(
  covars.path,
  as_data_table = TRUE, description = "Covariates"
)
rm(covars.path)

# Validate covariates
validate_not_empty(covars.dt, "covariates")
validate_columns(
  covars.dt,
  c("EID", "EDUC_lv_ses2_0"),
  "covariates"
)

if (!check_files_exist(hc_hvr.path, stop_on_missing = FALSE)) {
  log_warn("Head-size adjusted data not found, running 05_adjust_headsize.R")
  source(here("R", "05_adjust_headsize.R"))
}
hc_hvr.lst <- read_rds_safe(
  hc_hvr.path,
  description = "Head-size adjusted data"
)

# Validate structure
validate_not_empty(hc_hvr.lst, "head-size adjusted data")
if (!is.list(hc_hvr.lst)) {
  log_error("Expected list structure from hc-hvr_adj.rds, got %s", class(hc_hvr.lst))
  stop("Invalid data structure: expected list", call. = FALSE)
}
if (!NETWORK %in% names(hc_hvr.lst)) {
  log_error("Network '%s' not found. Available: %s", NETWORK, paste(names(hc_hvr.lst), collapse = ", "))
  stop(sprintf("Invalid network specification: %s", NETWORK), call. = FALSE)
}

hc_hvr.lst <- hc_hvr.lst[[NETWORK]]

# Validate contents
invisible({
  hc_hvr.lst |> lapply(lapply, \(hc_hvr.subdt) {
    cols <- c("INST", "SEX", "AGE", "ICC", "SIDE", "ADJ", "HC", "LV", "HVR")
    validate_not_empty(hc_hvr.subdt, "imaging data")
    validate_columns(hc_hvr.subdt, cols, "imaging data")
    validate_categorical(hc_hvr.subdt, "INST", c("ses-2", "ses-3"))
    validate_categorical(hc_hvr.subdt, "SEX", c("Female", "Male"))
    validate_age(hc_hvr.subdt, min_age = 40, max_age = 85, na_allowed = FALSE)
    validate_num_range(hc_hvr.subdt, "ICC", min = 1000, max = 2600, na = FALSE)
    validate_categorical(hc_hvr.subdt, "SIDE", c("L", "R", "LR"))
    validate_categorical(hc_hvr.subdt, "ADJ", names(ADJS))
    validate_num_range(hc_hvr.subdt, "HC", na_allowed = FALSE)
    validate_num_range(hc_hvr.subdt, "LV", na_allowed = FALSE)
    validate_num_range(hc_hvr.subdt, "HVR")
  })
})

rm(hc_hvr.path)

# ----- Data Cleaning -----
log_section("Preparing data for GAMLSS")

# Get education mapping from config (centralized in config.R)
edu_years_map <- get_education_years_map()
educ.dt <- covars.dt[, .(EDUC_char = as.character(EDUC_lv_ses2_0)), EID]
educ.dt[, EDUC_num := edu_years_map[EDUC_char]] |> invisible()
setkey(educ.dt, EID)

# Extract scanner site (Assessment Centre) for random effect
site.dt <- covars.dt[, .(SITE = as.factor(SITE_ses2)), EID]
site.dt <- site.dt[!is.na(SITE)]
setkey(site.dt, EID)
log_info("Site data: %d subjects with valid scanner site", nrow(site.dt))

# Summed sides for HC/HVR & make long format
hc_hvr.lst <- lapply(
  hc_hvr.lst,
  \(crs_lvl) {
    # Ignore matched subcohorts
    data.subdt <- crs_lvl$ALL

    # Remove redundant HVR(ICC) == HVR(Unadj)
    data.subdt["STX", on = "ADJ", HVR := NA]

    # If subfield data; use only total HC/LV volumes
    if ("SUBFIELD" %in% names(data.subdt)) {
      data.subdt <- data.subdt["total", on = "SUBFIELD", -"SUBFIELD"]
    }

    # Merge Education and Site data, remove NAs, and convert to long format
    data.subdt <- data.subdt[educ.dt][site.dt, nomatch = NULL] |>
      melt(
        na.rm = TRUE,
        measure = names(ROIS),
        value = "VAL",
        variable = "ROI"
      ) |>
      na.omit()

    # Final cleaning and output
    data.subdt[, SEX := factor(as.character(SEX))]
    data.subdt[, SIDE := factor(SIDE, levels = c("L", "R", "LR"))]
    data.subdt[, ADJ := factor(ADJ)]
    data.subdt[, SITE := factor(SITE)]
    setcolorder(data.subdt, c("EDUC_char", "EDUC_num", "SITE"), after = "SEX")
    setkey(data.subdt, EID, INST)
    data.subdt
  }
)

# ----- GAMLSS Fitting -----
log_section("GAMLSS modeling (stratified by Sex)")
sexes.v <- levels(hc_hvr.lst$CRS$SEX)
group_cols.v <- c("SEX", "ROI", "ADJ", "SIDE")

# Create sorting table for iteration
sort.dt <- unique(hc_hvr.lst$CRS[
  , lapply(.SD, as.character),
  .SDcols = group_cols.v
])
setkey(sort.dt)

# Validate sort.dt structure
validate_not_empty(sort.dt, "model iteration table")
validate_columns(sort.dt, group_cols.v, "model iteration table")
log_info("Model combinations to fit: %d", nrow(sort.dt))

# GAMLSS family names (defined here so available in all code paths)
fam_names.v <- c(
  BE = "Beta",
  BCCG = "Box-Cox Cole & Green",
  L_NO = "Gaussian (logit-transformed)",
  NO = "Gaussian"
)

gamlss.path <- get_data_path("models", "fit", "gamlss")
if (check_files_exist(gamlss.path, stop_on_missing = FALSE)) {
  log_info("Loading existing GAMLSS data/fits")
  gamlss.lst <- read_rds_safe(gamlss.path, "GAMLSS data & model-fits")

  # Validate loaded data structure
  if (!is.list(gamlss.lst) || !all(c("DATA", "TRAIN", "TEST", "FINAL") %in% names(gamlss.lst))) {
    log_error("Invalid GAMLSS data structure: missing required components")
    stop("Invalid GAMLSS data file: regenerate by deleting and rerunning", call. = FALSE)
  }
  validate_not_empty(gamlss.lst$DATA$CRS, "loaded GAMLSS cross-sectional data")
  log_info("GAMLSS data validated: %d subjects", gamlss.lst$DATA$CRS[!duplicated(EID), .N])
} else {
  # ---- Data setup ----
  log_info("Setting up GAMLSS data (Cross-sectional)")
  gamlss.lst <- list()
  gamlss.lst$DATA <- copy(hc_hvr.lst)
  gamlss.dt <- gamlss.lst$DATA$CRS

  # Validate GAMLSS data
  validate_not_empty(gamlss.dt, "GAMLSS cross-sectional data")
  validate_columns(
    gamlss.dt,
    c("EID", "INST", "SEX", "AGE", "ICC", "SIDE", "ADJ", "ROI", "VAL", "EDUC_num", "SITE"),
    "GAMLSS data"
  )
  validate_categorical(gamlss.dt, "SEX", c("Female", "Male"))
  validate_categorical(gamlss.dt, "SIDE", c("L", "R", "LR"))
  validate_categorical(gamlss.dt, "ADJ", names(ADJS))
  validate_categorical(gamlss.dt, "ROI", names(ROIS))
  validate_age(gamlss.dt, min_age = 40, max_age = 85, na_allowed = FALSE)

  log_info(
    "Cross-sectional data: %d Females & %d Males",
    gamlss.dt[!duplicated(EID)]["Female", on = "SEX", .N],
    gamlss.dt[!duplicated(EID)]["Male", on = "SEX", .N]
  )

  # ---- Subject-level train/test split (by sex) ----
  log_info("Splitting Train/Test cohorts (80%/20%) by sex and age (5y bins)")
  sexes.v |>
    lapply(\(sex) {
      data.subdt <- copy(gamlss.dt[!duplicated(EID)])
      train.dt <- data.subdt[
        sex,
        on = "SEX",
        .(EID, AGE_bin = cut(
          AGE,
          breaks = seq(floor(min(AGE)), ceiling(max(AGE)) + 1, by = 5),
          right = FALSE
        ))
      ] |>
        initial_split(.8, "AGE_bin") |>
        training() |>
        setkey(EID)
    }) |>
    lapply(\(train.dt) {
      gamlss.dt[train.dt, SPLIT := "train"]
      gamlss.dt[is.na(SPLIT), SPLIT := "test"]
    }) |>
    invisible()

  split_n.dt <- gamlss.dt[!duplicated(EID), .N, keyby = .(SPLIT, SEX)]

  # Validate train/test split has data for both sexes
  validate_categorical(gamlss.dt, "SPLIT", c("train", "test"))
  if (split_n.dt[.("train", "Female"), N] == 0 || split_n.dt[.("train", "Male"), N] == 0) {
    log_error("Training set is empty for one or both sexes")
    stop("Invalid train/test split: empty training set", call. = FALSE)
  }
  if (split_n.dt[.("test", "Female"), N] == 0 || split_n.dt[.("test", "Male"), N] == 0) {
    log_error("Test set is empty for one or both sexes")
    stop("Invalid train/test split: empty test set", call. = FALSE)
  }

  log_info(
    "Training: %d Females & %d Males; Testing: %d Females & %d Males",
    split_n.dt[.("train", "Female"), N], split_n.dt[.("train", "Male"), N],
    split_n.dt[.("test", "Female"), N], split_n.dt[.("test", "Male"), N]
  )
  rm(split_n.dt)

  # ---- Helper functions ----
  save_fit.fn <- function(mod.lst, mod_name, mod.fit) {
    comp.subdt <- data.table(
      MOD   = mod_name,
      AIC   = mod.fit$aic,
      BIC   = BIC(mod.fit),
      DF    = mod.fit$df.fit,
      GDEV  = deviance(mod.fit)
    )
    if (is.null(mod.lst$COMP)) {
      mod.lst$COMP <- comp.subdt
    } else {
      mod.lst$COMP <- rbind(mod.lst$COMP, comp.subdt)
    }
    setkey(mod.lst$COMP, MOD)
    mod.lst[[mod_name]] <- mod.fit
    mod.lst
  }

  pick_best_mod.fn <- function(comp.dt, mod_colname = "MOD", bic_crit = 10) {
    # Heuristics:
    # 1) Use AIC.
    # 2) If model is drastically worse in BIC (+10); reconsider.
    # Burnham, Kenneth P. & Anderson, David R. (2002).
    # Model Selection and Multimodel Inference:
    # A Practical Information-Theoretic Approach (2nd ed.).
    best_mod <- comp.dt[which.min(AIC), MOD]
    bic_diff <- comp.dt[, .SD[best_mod, on = mod_colname, BIC] - min(BIC)]
    if (bic_diff > bic_crit) {
      return(list(MOD = comp.dt[which.min(BIC), MOD], WARN = TRUE))
    }
    list(MOD = best_mod, WARN = FALSE)
  }

  # ---- GAMLSS families for each ROI ----
  families.v <- get_parameter("gamlss", "families")

  # ---- Training ----
  log_info("Fitting GAMLSS models (training)")
  pb <- progress_bar$new(
    format = "GAMLSS (training) | :what [:bar] :current/:total\n",
    # 2 model fits BE & LOGIT-NORM | BCCG & NORM
    total = sort.dt[, 2 * .N], clear = FALSE, width = 100, show_after = 0
  )
  pb$tick(0)
  setkeyv(gamlss.dt, group_cols.v)
  for (i in seq_len(nrow(sort.dt))) {
    row.dt <- sort.dt[i]
    sex <- row.dt$SEX
    roi <- row.dt$ROI
    adj <- row.dt$ADJ
    side <- row.dt$SIDE

    # --- Data preparation ---
    gamlss.subdt <- gamlss.dt[row.dt][
      VAL > 0, .(EID, AGE, EDUC_num, ICC, SITE, VAL),
      keyby = SPLIT
    ] |> na.omit()

    gamlss_train.subdt <- gamlss.subdt["train", -"SPLIT"]

    if (!nrow(gamlss_train.subdt)) {
      pb$tick(2)
      next
    }

    # ---- Training ----
    gamlss.lst$TRAIN[[sex]][[roi]][[adj]][[side]] <- list()
    sublist.lst <- gamlss.lst$TRAIN[[sex]][[roi]][[adj]][[side]]
    # Formula includes random effect for scanner site (Assessment Centre)
    formula_base <- if (roi != "HVR" && adj == "NON") {
      VAL ~ cs(AGE) + ICC + EDUC_num + random(SITE)
    } else {
      VAL ~ cs(AGE) + EDUC_num + random(SITE)
    }
    for (fam in families.v[[roi]]) {
      if (!fam %in% names(fam_names.v)) {
        log_warn("%s: %s family not recognized; skipping", roi, fam)
        pb$tick(1)
        next
      }
      pb$tick(tokens = list(what = sprintf(
        "%s — %s (%s; %s): %s", sex, roi, side, adj, fam_names.v[[fam]]
      )))
      if (fam == "NO") {
        # -- Normal model --
        output.fit <- gamlss(
          formula = formula_base,
          sigma.formula = ~ cs(AGE),
          data = gamlss_train.subdt,
          family = NO()
        )
      } else if (fam == "L_NO") {
        # -- Logit-Normal model --
        eps <- 1e-6
        gamlss_train.subdt[
          , VAL_logit := VAL |> pmin(1 - eps) |> pmax(eps) |> qlogis()
        ]
        output.fit <- gamlss(
          formula = VAL_logit ~ cs(AGE) + EDUC_num + random(SITE),
          sigma.formula = ~ cs(AGE),
          data = gamlss_train.subdt,
          family = NO()
        )
      } else if (fam == "BCCG") {
        # -- Box-Cox Cole & Green  model --
        output.fit <- gamlss(
          formula = formula_base,
          sigma.formula = ~ cs(AGE),
          nu.formula = ~ cs(AGE),
          data = gamlss_train.subdt,
          family = BCCG()
        )
      } else if (fam == "BE") {
        # -- Beta Model --
        output.fit <- gamlss(
          formula = formula_base,
          sigma.formula = ~ cs(AGE),
          nu.formula = ~ cs(AGE),
          data = gamlss_train.subdt,
          family = BE()
        )
      }
      sublist.lst <- save_fit.fn(sublist.lst, fam, output.fit)
    }

    # ---- Model selection: AIC/BIC ----
    comp_result.lst <- pick_best_mod.fn(sublist.lst$COMP)
    if (comp_result.lst$WARN) {
      log_warn(
        "%s — %s (%s; %s): best AIC model has BIC > 10 above the minimum BIC",
        sex, roi, side, adj
      )
    }
    best_mod <- comp_result.lst$MOD
    best_mod.fit <- sublist.lst[[best_mod]]
    sublist.lst$BEST <- list(MOD = best_mod, FIT = best_mod.fit)
    gamlss.lst$TRAIN[[sex]][[roi]][[adj]][[side]] <- sublist.lst

    log_info(
      "%s — %s (%s; %s): Best model - %s",
      sex, roi, side, adj, fam_names.v[[best_mod]]
    )
  }
  rm(fam, formula_base)

  # ---- Constrain bilateral chosen model ----
  for (
    i in sort.dt |>
      unique() |>
      nrow() |>
      seq_len()
  ) {
    row.dt <- sort.dt[i]
    sex <- row.dt$SEX
    roi <- row.dt$ROI
    adj <- row.dt$ADJ

    sublist.lst <- gamlss.lst$TRAIN[[sex]][[roi]][[adj]]
    left_mod <- sublist.lst$L$BEST$MOD
    right_mod <- sublist.lst$R$BEST$MOD

    if (left_mod == right_mod) next

    log_warn(
      "%s — %s (%s): Selected models L (%s) != R (%s)",
      sex, roi, adj, left_mod, right_mod
    )

    fit_msrs.subdt <- sublist.lst[c("L", "R")] |>
      lapply("[[", "COMP") |>
      rbindlist(idcol = "SIDE") |>
      setkey(MOD, SIDE)

    # Sum of AIC for best L$MOD == R$MOD
    aic_same <- fit_msrs.subdt[, sum(AIC), MOD][which.min(V1), V1]

    # Sum of AIC for min(L$MOD) + min(R$MOD)
    aic_free <- sum(
      fit_msrs.subdt[.(left_mod, "L"), AIC],
      fit_msrs.subdt[.(right_mod, "R"), AIC]
    )

    if ((aic_same - aic_free) > 10) {
      log_info(
        paste(
          "%s — %s (%s):",
          "Different families for L/R significantly improves fit",
          "-- %.0f vs %.0f"
        ),
        sex, roi, adj, aic_free, aic_same
      )
    } else {
      log_info(
        paste(
          "%s — %s (%s): Different families for L/R fit (%.0f vs %.0f)",
          "is not significantly better (+10 AIC);",
          "Constraining to same family."
        ),
        sex, roi, adj, aic_free, aic_same
      )
      mod <- fit_msrs.subdt[, sum(AIC), MOD][which.min(V1), MOD]
      sublist.lst$L$BEST <- list(MOD = mod, FIT = sublist.lst$L[[mod]])
      sublist.lst$R$BEST <- list(MOD = mod, FIT = sublist.lst$R[[mod]])
      gamlss.lst$TRAIN[[sex]][[roi]][[adj]] <- sublist.lst
    }
  }

  # --- Extract model comparison summary ---
  log_section("Model comparison across families")

  # Extract all model comparison tables (all fitted models, not just best)
  model_comparison.dt <- sort.dt[, {
    train.lst <- gamlss.lst$TRAIN[[SEX]][[ROI]][[ADJ]][[SIDE]]
    if (is.null(train.lst) || is.null(train.lst$COMP)) NULL else train.lst$COMP
  }, by = .(SEX, ROI, ADJ, SIDE)]

  # Calculate delta AIC/BIC relative to best model within each group
  if (nrow(model_comparison.dt) > 0) {
    model_comparison.dt[
      , `:=`(
        DELTA_AIC = AIC - min(AIC),
        DELTA_BIC = BIC - min(BIC),
        BEST_AIC = MOD[which.min(AIC)],
        BEST_BIC = MOD[which.min(BIC)]
      ),
      by = .(SEX, ROI, ADJ, SIDE)
    ]

    # Log summary
    n_comparisons <- model_comparison.dt[
      , .N, by = .(SEX, ROI, ADJ, SIDE)
    ][
      , .N
    ]
    n_aic_bic_agree <- model_comparison.dt[
      , .(AGREE = BEST_AIC[1] == BEST_BIC[1]),
      by = .(SEX, ROI, ADJ, SIDE)
    ][AGREE == TRUE, .N]

    log_info("Total model comparisons: %d combinations", n_comparisons)
    log_info(
      "AIC/BIC agreement: %d/%d (%.1f%%)",
      n_aic_bic_agree, n_comparisons,
      100 * n_aic_bic_agree / n_comparisons
    )
    log_info(
      "Mean models fitted per combination: %.1f",
      model_comparison.dt[, .N, by = .(SEX, ROI, ADJ, SIDE)][, mean(N)]
    )
  }

  # ---- Testing: Hold-out centile validation ----
  log_info("Evaluating best selected GAMLSS models (testing)")
  cent_thresholds.v <- get_parameter("gamlss", "summary_centiles") * 100
  pb <- progress_bar$new(
    format = "GAMLSS (testing) | :what [:bar] :current/:total\n",
    total = sort.dt[, .N], clear = FALSE, width = 100, show_after = 0
  )
  pb$tick(0)
  for (i in seq_len(nrow(sort.dt))) {
    row.dt <- sort.dt[i]
    sex <- row.dt$SEX
    roi <- row.dt$ROI
    adj <- row.dt$ADJ
    side <- row.dt$SIDE

    # --- Data preparation ---
    gamlss.subdt <- gamlss.dt[row.dt][
      VAL > 0, .(EID, AGE, EDUC_num, ICC, SITE, VAL),
      keyby = SPLIT
    ] |> na.omit()

    # Training dataset needs to be in environment for GAMLSS
    gamlss_train.subdt <- gamlss.subdt["train", -"SPLIT"]
    gamlss_test.subdt <- gamlss.subdt["test", -"SPLIT"]

    if (!nrow(gamlss_train.subdt) || !nrow(gamlss_test.subdt)) {
      pb$tick()
      next
    }

    # Selected best model
    best_mod.lst <- gamlss.lst$TRAIN[[sex]][[roi]][[adj]][[side]]$BEST

    # Progress bar
    pb$tick(tokens = list(what = sprintf(
      "%s — %s (%s; %s): %s",
      sex, roi, side, adj, fam_names.v[[best_mod.lst$MOD]]
    )))

    # --- Diagnostic: Check for potential extrapolation ---
    # Compare training vs test covariate ranges
    train_age_range <- range(gamlss_train.subdt$AGE)
    test_age_range <- range(gamlss_test.subdt$AGE)
    train_icc_range <- range(gamlss_train.subdt$ICC)
    test_icc_range <- range(gamlss_test.subdt$ICC)

    # Flag extrapolation if test is outside training range
    age_extrap <- test_age_range[1] < train_age_range[1] ||
      test_age_range[2] > train_age_range[2]
    icc_extrap <- test_icc_range[1] < train_icc_range[1] ||
      test_icc_range[2] > train_icc_range[2]

    if ((age_extrap || icc_extrap) && best_mod.lst$MOD %in% c("BCCG", "BE")) {
      log_warn(
        paste(
          "%s %s (%s; %s): Test set may extrapolate",
          "(Age: [%.1f-%.1f] vs [%.1f-%.1f], ICC: [%.0f-%.0f] vs [%.0f-%.0f])"
        ),
        sex, roi, side, adj,
        test_age_range[1], test_age_range[2],
        train_age_range[1], train_age_range[2],
        test_icc_range[1], test_icc_range[2],
        train_icc_range[1], train_icc_range[2]
      )
    }

    # --- Predicted centiles ---
    if (best_mod.lst$MOD == "NO") {
      # -- Selected model: Gaussian --
      gamlss_test.subdt[
        , CENT_pred := pNO(
          VAL,
          mu = predict(best_mod.lst$FIT, new = .SD, type = "response"),
          sigma = predict(
            best_mod.lst$FIT, "sigma", new = .SD, type = "response"
          )
        ) * 100
      ] |> invisible()
    } else if (best_mod.lst$MOD == "L_NO") {
      # -- Selected model: Logit-Normal --
      eps <- 1e-6
      invisible({
        gamlss_test.subdt[
          , VAL_logit := VAL |> pmin(1 - eps) |> pmax(eps) |> qlogis()
        ]
        gamlss_test.subdt[
          , CENT_pred := pNO(
            VAL_logit,
            mu = predict(best_mod.lst$FIT, new = .SD, type = "response"),
            sigma = predict(
              best_mod.lst$FIT, "sigma", new = .SD, type = "response"
            )
          ) * 100
        ]
      })
    } else if (best_mod.lst$MOD == "BE") {
      # -- Selected model: Beta --
      # Capture warnings about re-fitting for safe predictions
      refit_warning <- FALSE
      withCallingHandlers(
        {
          gamlss_test.subdt[
            , CENT_pred := pBE(
              VAL,
              mu = predict(best_mod.lst$FIT, new = .SD, type = "response"),
              sigma = predict(
                best_mod.lst$FIT, "sigma", new = .SD, type = "response"
              )
            ) * 100
          ] |> invisible()
        },
        warning = function(w) {
          if (grepl("discrepancy.*re-fit", w$message)) {
            refit_warning <<- TRUE
            log_debug(
              "%s %s (%s; %s): GAMLSS refit warning during prediction",
              sex, roi, side, adj
            )
          }
          invokeRestart("muffleWarning")
        }
      )
    } else if (best_mod.lst$MOD == "BCCG") {
      # -- Selected model: Box-Cox Cole & Green --
      # Capture warnings about re-fitting for safe predictions
      refit_warning <- FALSE
      withCallingHandlers(
        {
          gamlss_test.subdt[
            , CENT_pred := pBCCG(
              VAL,
              mu = predict(best_mod.lst$FIT, new = .SD, type = "response"),
              sigma = predict(
                best_mod.lst$FIT, "sigma", new = .SD, type = "response"
              ),
              nu = predict(
                best_mod.lst$FIT, "nu", new = .SD, type = "response"
              )
            ) * 100
          ] |> invisible()
        },
        warning = function(w) {
          if (grepl("discrepancy.*re-fit", w$message)) {
            refit_warning <<- TRUE
            log_debug(
              "%s %s (%s; %s): GAMLSS refit warning during prediction",
              sex, roi, side, adj
            )
          }
          invokeRestart("muffleWarning")
        }
      )
    }

    # Initialize refit_warning for other distributions
    if (!exists("refit_warning")) refit_warning <- FALSE

    # --- NA check for predicted centiles (CL3) ---
    n_na_cent <- sum(is.na(gamlss_test.subdt$CENT_pred))
    if (n_na_cent > 0) {
      log_warn(
        "%s %s (%s; %s): %d NA values in predicted centiles (%.1f%%)",
        sex, roi, side, adj, n_na_cent,
        100 * n_na_cent / nrow(gamlss_test.subdt)
      )
    }

    cent_calibration.v <- sapply(
      cent_thresholds.v,
      \(t) mean(gamlss_test.subdt$CENT_pred < t) * 100
    )

    mae <- mean(abs(cent_calibration.v - cent_thresholds.v))

    calibration.dt <- data.table(
      CENT_exp = cent_thresholds.v,
      CENT_obs = round(cent_calibration.v, 2),
      DIFF = round(cent_calibration.v - cent_thresholds.v, 2)
    )
    gamlss.lst$TEST[[sex]][[roi]][[adj]][[side]] <- list(
      MOD = best_mod.lst$MOD,
      FIT = best_mod.lst$FIT,
      CENT_VALID = calibration.dt,
      MAE = mae,
      INTERPRETATION = fcase(
        mae < 2, "Excellent",
        mae < 5, "Good",
        default = "Poor"
      ),
      REFIT_WARNING = refit_warning,
      AGE_EXTRAP = age_extrap,
      ICC_EXTRAP = icc_extrap
    )
  }

  # --- Diagnostic summary ---
  log_section("Hold-out validation diagnostics")

  # Extract validation diagnostics iterating row-wise over the definition grid
  test_results.dt <- sort.dt[, {
    res <- gamlss.lst$TEST[[SEX]][[ROI]][[ADJ]][[SIDE]]

    # If the model failed/doesn't exist, return NULL to drop the row
    # Otherwise, extract the specific scalar fields you need
    if (is.null(res)) NULL else list(
      MOD            = res$MOD,
      MAE            = res$MAE,
      INTERPRETATION = res$INTERPRETATION,
      REFIT_WARNING  = res$REFIT_WARNING,
      AGE_EXTRAP     = res$AGE_EXTRAP,
      ICC_EXTRAP     = res$ICC_EXTRAP
    )
  }, by = .(SEX, ROI, ADJ, SIDE)] # 'by' is crucial: it forces scalar processing

  # Optional: Extract calibration data.tables for detailed diagnostics
  # This creates a long-format data.table with all calibration results
  calibration_details.dt <- sort.dt[, {
    res <- gamlss.lst$TEST[[SEX]][[ROI]][[ADJ]][[SIDE]]
    if (is.null(res) || is.null(res$CENT_VALID)) NULL else res$CENT_VALID
  }, by = .(SEX, ROI, ADJ, SIDE)]

  # Summary statistics
  n_models <- nrow(test_results.dt)
  n_refit <- sum(test_results.dt$REFIT_WARNING, na.rm = TRUE)
  n_age_extrap <- sum(test_results.dt$AGE_EXTRAP, na.rm = TRUE)
  n_icc_extrap <- sum(test_results.dt$ICC_EXTRAP, na.rm = TRUE)
  n_excellent <- sum(
    test_results.dt$INTERPRETATION == "Excellent", na.rm = TRUE
  )
  n_good <- sum(test_results.dt$INTERPRETATION == "Good", na.rm = TRUE)
  n_poor <- sum(test_results.dt$INTERPRETATION == "Poor", na.rm = TRUE)

  log_info("Total models validated: %d", n_models)
  log_info(
    "Model quality: %d Excellent (%.1f%%), %d Good (%.1f%%), %d Poor (%.1f%%)",
    n_excellent, 100 * n_excellent / n_models,
    n_good, 100 * n_good / n_models,
    n_poor, 100 * n_poor / n_models
  )
  log_info(
    "GAMLSS refit warnings: %d (%.1f%%)",
    n_refit, 100 * n_refit / n_models
  )
  log_info(
    "Potential extrapolation: Age=%d (%.1f%%), ICC=%d (%.1f%%)",
    n_age_extrap, 100 * n_age_extrap / n_models,
    n_icc_extrap, 100 * n_icc_extrap / n_models
  )

  # Check correlation between refit warnings and validation quality
  if (n_refit > 0) {
    mae_refit <- mean(
      test_results.dt[REFIT_WARNING == TRUE]$MAE, na.rm = TRUE
    )
    mae_no_refit <- mean(
      test_results.dt[REFIT_WARNING == FALSE]$MAE, na.rm = TRUE
    )
    log_info(
      "Mean MAE: with refit warning = %.2f, without = %.2f",
      mae_refit, mae_no_refit
    )
  }

  # ---- Final Models ----
  log_info("Fitting final GAMLSS on the full dataset")
  pb <- progress_bar$new(
    format = "GAMLSS (Final) | :what [:bar] :current/:total\n",
    total = sort.dt[, .N], clear = FALSE, width = 100, show_after = 0
  )
  pb$tick(0)
  for (i in seq_len(nrow(sort.dt))) {
    row.dt <- sort.dt[i]
    sex <- row.dt$SEX
    roi <- row.dt$ROI
    adj <- row.dt$ADJ
    side <- row.dt$SIDE

    # --- Data preparation ---
    gamlss_final.subdt <- gamlss.dt[row.dt][
      VAL > 0, .(EID, AGE, EDUC_num, ICC, SITE, VAL)
    ] |> na.omit()

    if (!nrow(gamlss_final.subdt)) {
      pb$tick()
      next
    }

    # Selected best model
    best_mod.lst <- gamlss.lst$TRAIN[[sex]][[roi]][[adj]][[side]]$BEST

    pb$tick(tokens = list(what = sprintf(
      "%s — %s (%s; %s): %s",
      sex, roi, side, adj, fam_names.v[[best_mod.lst$MOD]]
    )))

    # Formula includes random effect for scanner site (Assessment Centre)
    formula_base <- if (roi != "HVR" && adj == "NON") {
      VAL ~ cs(AGE) + ICC + EDUC_num + random(SITE)
    } else {
      VAL ~ cs(AGE) + EDUC_num + random(SITE)
    }

    if (best_mod.lst$MOD == "NO") {
      # -- Normal model --
      output.fit <- gamlss(
        formula = formula_base,
        sigma.formula = ~ cs(AGE),
        data = gamlss_final.subdt,
        family = NO()
      )
    } else if (best_mod.lst$MOD == "L_NO") {
      # -- Logit-Normal model --
      eps <- 1e-6
      gamlss_final.subdt[
        , VAL_logit := VAL |> pmin(1 - eps) |> pmax(eps) |> qlogis()
      ]
      output.fit <- gamlss(
        formula = VAL_logit ~ cs(AGE) + EDUC_num + random(SITE),
        sigma.formula = ~ cs(AGE),
        data = gamlss_final.subdt,
        family = NO()
      )
    } else if (best_mod.lst$MOD == "BCCG") {
      # -- Box-Cox Cole & Green  model --
      output.fit <- gamlss(
        formula = formula_base,
        sigma.formula = ~ cs(AGE),
        nu.formula = ~ cs(AGE),
        data = gamlss_final.subdt,
        family = BCCG()
      )
    } else if (best_mod.lst$MOD == "BE") {
      # -- Beta Model --
      output.fit <- gamlss(
        formula = formula_base,
        sigma.formula = ~ cs(AGE),
        nu.formula = ~ cs(AGE),
        data = gamlss_final.subdt,
        family = BE()
      )
    }

    log_info(
      "%s — %s (%s; %s): %s -- AIC: %.0f - BIC: %.0f",
      sex, roi, side, adj,
      fam_names.v[[best_mod.lst$MOD]], output.fit$aic, BIC(output.fit)
    )

    # Definitive list
    gamlss.lst$FINAL[[sex]][[roi]][[adj]][[side]] <- list(
      FAMILY  = best_mod.lst$MOD,
      DATA    = gamlss_final.subdt,
      FIT     = output.fit,
      AIC     = output.fit$aic,
      BIC     = BIC(output.fit)
    )
    rm(output.fit)
  }

  # ---- Temporal validation (Longitudinal stability) ----
  # ---- Data setup ----
  log_info("Setting up GAMLSS data (Longitudinal)")
  gamlss.dt <- gamlss.lst$DATA$LNG

  # Validate longitudinal data
  validate_not_empty(gamlss.dt, "longitudinal data")
  validate_columns(
    gamlss.dt,
    c("EID", "INST", "SEX", "AGE", "ICC", "SIDE", "ADJ", "ROI", "VAL"),
    "longitudinal data"
  )
  validate_categorical(gamlss.dt, "SEX", c("Female", "Male"))
  validate_categorical(gamlss.dt, "INST", c("ses-2", "ses-3"))

  n.dt <- gamlss.dt[!duplicated(EID), .N, keyby = SEX]
  log_info(
    "Longitudinal data -- N: %d Females & %d Males",
    n.dt["Female", N], n.dt["Male", N]
  )

  time_lag.dt <- gamlss.dt[
    , .(Time = unique(AGE)),
    by = .(EID, INST, SEX)
  ][
    , .(LAG = diff(Time) * 12),
    by = .(EID, SEX)
  ][
    , .(AVG = mean(LAG), SD = sd(LAG)),
    keyby = SEX
  ]
  log_info(
    paste(
      "Longitudinal data -- Follow-up (months):",
      "F - %.1f (%.1f) & M - %.1f (%.1f)"
    ),
    time_lag.dt["Female", AVG], time_lag.dt["Female", SD],
    time_lag.dt["Male", AVG], time_lag.dt["Male", SD]
  )
  rm(n.dt, time_lag.dt)

  log_info("Normative curves validation")
  pb <- progress_bar$new(
    format = "GAMLSS (longitudinal) | :what [:bar] :current/:total\n",
    total = sort.dt[, .N], clear = FALSE, width = 100, show_after = 0
  )
  pb$tick(0)
  setkeyv(gamlss.dt, group_cols.v)
  for (i in seq_len(nrow(sort.dt))) {
    row.dt <- sort.dt[i]
    sex <- row.dt$SEX
    roi <- row.dt$ROI
    adj <- row.dt$ADJ
    side <- row.dt$SIDE

    # --- Data preparation ---
    gamlss.subdt <- gamlss.dt[row.dt][
      VAL > 0,
      .(AGE, EDUC_num, ICC, SITE, VAL),
      .(EID, TP = fifelse(INST %like% 2, "T1", "T2"))
    ] |> na.omit()

    if (!nrow(gamlss.subdt)) {
      pb$tick()
      next
    }

    # Final model
    final_mod.lst <- gamlss.lst$FINAL[[sex]][[roi]][[adj]][[side]]
    gamlss_final.subdt <- final_mod.lst$DATA

    pb$tick(tokens = list(what = sprintf(
      "%s — %s (%s; %s): %s",
      sex, roi, side, adj, fam_names.v[[final_mod.lst$FAMILY]]
    )))

    # --- Diagnostic: Check for potential extrapolation ---
    # Compare training (final model data) vs longitudinal covariate ranges
    train_age_range <- range(gamlss_final.subdt$AGE)
    lng_age_range <- range(gamlss.subdt$AGE)
    train_icc_range <- range(gamlss_final.subdt$ICC)
    lng_icc_range <- range(gamlss.subdt$ICC)

    # Flag extrapolation if longitudinal data outside training range
    lng_age_extrap <- lng_age_range[1] < train_age_range[1] ||
      lng_age_range[2] > train_age_range[2]
    lng_icc_extrap <- lng_icc_range[1] < train_icc_range[1] ||
      lng_icc_range[2] > train_icc_range[2]

    if ((lng_age_extrap || lng_icc_extrap) &&
        final_mod.lst$FAMILY %in% c("BCCG", "BE")
    ) {
      log_warn(
        paste(
          "%s %s (%s; %s): Longitudinal data may extrapolate",
          "(Age: [%.1f-%.1f] vs [%.1f-%.1f], ICC: [%.0f-%.0f] vs [%.0f-%.0f])"
        ),
        sex, roi, side, adj,
        lng_age_range[1], lng_age_range[2],
        train_age_range[1], train_age_range[2],
        lng_icc_range[1], lng_icc_range[2],
        train_icc_range[1], train_icc_range[2]
      )
    }

    # --- Predicted centiles ---
    if (final_mod.lst$FAMILY == "NO") {
      # -- Best model: Gaussian --
      gamlss.subdt[
        ,
        CENT := pNO(
          VAL,
          mu = predict(final_mod.lst$FIT, new = .SD, type = "response"),
          sigma = predict(
            final_mod.lst$FIT, "sigma", new = .SD, type = "response"
          )
        ) * 100,
        .SDcols = -"TP"
      ] |> invisible()
    } else if (final_mod.lst$FAMILY == "L_NO") {
      # -- Best model: Logit-Normal --
      eps <- 1e-6
      invisible({
        gamlss.subdt[
          , VAL_logit := VAL |> pmin(1 - eps) |> pmax(eps) |> qlogis()
        ]
        gamlss.subdt[
          ,
          CENT := pNO(
            VAL_logit,
            mu = predict(final_mod.lst$FIT, new = .SD, type = "response"),
            sigma = predict(
              final_mod.lst$FIT, "sigma", new = .SD, type = "response"
            )
          ) * 100,
          .SDcols = -"TP"
        ]
      })
    } else if (final_mod.lst$FAMILY == "BE") {
      # -- Best model: Beta --
      # Capture warnings about re-fitting for safe predictions
      lng_refit_warning <- FALSE
      withCallingHandlers(
        {
          gamlss.subdt[
            ,
            CENT := pBE(
              VAL,
              mu = predict(final_mod.lst$FIT, new = .SD, type = "response"),
              sigma = predict(
                final_mod.lst$FIT, "sigma", new = .SD, type = "response"
              )
            ) * 100,
            .SDcols = -"TP"
          ] |> invisible()
        },
        warning = function(w) {
          if (grepl("discrepancy.*re-fit", w$message)) {
            lng_refit_warning <<- TRUE
            log_debug(
              "%s %s (%s; %s): GAMLSS refit warning (longitudinal)",
              sex, roi, side, adj
            )
          }
          invokeRestart("muffleWarning")
        }
      )
    } else if (final_mod.lst$FAMILY == "BCCG") {
      # -- Best model: Box-Cox Cole & Green --
      # Capture warnings about re-fitting for safe predictions
      lng_refit_warning <- FALSE
      withCallingHandlers(
        {
          gamlss.subdt[
            ,
            CENT := pBCCG(
              VAL,
              mu = predict(final_mod.lst$FIT, new = .SD, type = "response"),
              sigma = predict(
                final_mod.lst$FIT, "sigma", new = .SD, type = "response"
              ),
              nu = predict(
                final_mod.lst$FIT, "nu", new = .SD, type = "response"
              )
            ) * 100,
            .SDcols = -"TP"
          ] |> invisible()
        },
        warning = function(w) {
          if (grepl("discrepancy.*re-fit", w$message)) {
            lng_refit_warning <<- TRUE
            log_debug(
              "%s %s (%s; %s): GAMLSS refit warning (longitudinal)",
              sex, roi, side, adj
            )
          }
          invokeRestart("muffleWarning")
        }
      )
    }

    # Initialize refit warning for other distributions
    if (!exists("lng_refit_warning")) lng_refit_warning <- FALSE

    # --- NA check for predicted centiles (CL3) ---
    n_na_cent <- sum(is.na(gamlss.subdt$CENT))
    if (n_na_cent > 0) {
      log_warn(
        "%s %s (%s; %s): %d NA values in longitudinal centiles (%.1f%%)",
        sex, roi, side, adj, n_na_cent,
        100 * n_na_cent / nrow(gamlss.subdt)
      )
    }

    cent.dt <- dcast(gamlss.subdt, EID ~ TP, value.var = "CENT")
    cent.dt[, CHANGE := T2 - T1]
    cent.dt[, CROSS_05 := (T1 < 5 & T2 > 5) | (T1 > 5 & T2 < 5)]
    cent.dt[, CROSS_95 := (T1 < 95 & T2 > 95) | (T1 > 95 & T2 < 95)]
    setnames(cent.dt, c("T1", "T2"), paste0("CENT_t", 1:2))

    cent.lst <- list(
      CENTILES = cent.dt,
      CORR = cent.dt[, cor(CENT_t1, CENT_t2)],
      MEAN_CHANGE = cent.dt[, mean(CHANGE)],
      SD_CHANGE = cent.dt[, sd(CHANGE)],
      CROSS_PROP = cent.dt[, mean(CROSS_05 | CROSS_95)],
      REFIT_WARNING = lng_refit_warning,
      AGE_EXTRAP = lng_age_extrap,
      ICC_EXTRAP = lng_icc_extrap
    )
    cent.lst$INTERPRETATION <- fifelse(
      cent.lst$CORR > 0.8 & abs(cent.lst$MEAN_CHANGE) < 2,
      "Good longitudinal stability",
      "May not capture individual trajectories"
    )

    gamlss.lst$FINAL[[sex]][[roi]][[adj]][[side]]$TEMP_VAL <- cent.lst
  }

  # --- Diagnostic summary ---
  log_section("Longitudinal stability diagnostics")

  # Extract temporal validation diagnostics
  stability_results.dt <- sort.dt[, {
    temp.lst <- gamlss.lst$FINAL[[SEX]][[ROI]][[ADJ]][[SIDE]]$TEMP_VAL
    if (is.null(temp.lst)) NULL else list(
      CORR = temp.lst$CORR,
      MEAN_CHANGE = temp.lst$MEAN_CHANGE,
      INTERPRETATION = temp.lst$INTERPRETATION,
      REFIT_WARNING = temp.lst$REFIT_WARNING,
      AGE_EXTRAP = temp.lst$AGE_EXTRAP,
      ICC_EXTRAP = temp.lst$ICC_EXTRAP
    )
  }, by = .(SEX, ROI, ADJ, SIDE)]

  # Optional: Extract longitudinal centile data.tables for detailed diagnostics
  # This creates a long-format data.table with all individual trajectories
  longitudinal_centiles.dt <- sort.dt[, {
    temp.lst <- gamlss.lst$FINAL[[SEX]][[ROI]][[ADJ]][[SIDE]]$TEMP_VAL
    if (is.null(temp.lst) || is.null(temp.lst$CENTILES)) {
      NULL
    } else {
      temp.lst$CENTILES
    }
  }, by = .(SEX, ROI, ADJ, SIDE)]

  # Derive centile-specific stability metrics
  if (nrow(longitudinal_centiles.dt) > 0) {
    centile_stability.dt <- longitudinal_centiles.dt[, .(
      N_SUBJECTS = .N,
      MEAN_T1 = mean(CENT_t1, na.rm = TRUE),
      MEAN_T2 = mean(CENT_t2, na.rm = TRUE),
      MEAN_CHANGE = mean(CHANGE, na.rm = TRUE),
      SD_CHANGE = sd(CHANGE, na.rm = TRUE),
      CORR_T1_T2 = cor(CENT_t1, CENT_t2, use = "complete.obs"),
      N_CROSS_05 = sum(CROSS_05, na.rm = TRUE),
      PCT_CROSS_05 = 100 * mean(CROSS_05, na.rm = TRUE),
      N_CROSS_95 = sum(CROSS_95, na.rm = TRUE),
      PCT_CROSS_95 = 100 * mean(CROSS_95, na.rm = TRUE)
    ), by = .(SEX, ROI, ADJ, SIDE)]

    # Define centile ranges for stratified analysis
    longitudinal_centiles.dt[
      , CENTILE_RANGE := fcase(
        CENT_t1 < 10, "Bottom_decile",
        CENT_t1 < 25, "Lower_quartile",
        CENT_t1 < 75, "Middle_range",
        CENT_t1 < 90, "Upper_quartile",
        default = "Top_decile"
      )
    ]

    # Stratified metrics by baseline centile range
    centile_stratified.dt <- longitudinal_centiles.dt[, .(
      N_SUBJECTS = .N,
      MEAN_CHANGE = mean(CHANGE, na.rm = TRUE),
      SD_CHANGE = sd(CHANGE, na.rm = TRUE),
      CORR_T1_T2 = cor(CENT_t1, CENT_t2, use = "complete.obs"),
      PCT_CROSS_ANY = 100 * mean(CROSS_05 | CROSS_95, na.rm = TRUE)
    ), by = .(SEX, ROI, ADJ, SIDE, CENTILE_RANGE)]

    log_info(
      "Centile stability: mean crossing 5th = %.1f%%, 95th = %.1f%%",
      centile_stability.dt[, mean(PCT_CROSS_05)],
      centile_stability.dt[, mean(PCT_CROSS_95)]
    )
  }

  # Summary statistics
  n_lng_models <- nrow(stability_results.dt)
  n_lng_refit <- sum(stability_results.dt$REFIT_WARNING, na.rm = TRUE)
  n_lng_age_extrap <- sum(stability_results.dt$AGE_EXTRAP, na.rm = TRUE)
  n_lng_icc_extrap <- sum(stability_results.dt$ICC_EXTRAP, na.rm = TRUE)
  n_good_stability <- sum(
    stability_results.dt$INTERPRETATION == "Good longitudinal stability",
    na.rm = TRUE
  )
  n_poor_stability <- n_lng_models - n_good_stability

  log_info("Total models validated (longitudinal): %d", n_lng_models)
  log_info(
    "Stability quality: %d Good (%.1f%%), %d Poor (%.1f%%)",
    n_good_stability, 100 * n_good_stability / n_lng_models,
    n_poor_stability, 100 * n_poor_stability / n_lng_models
  )
  log_info(
    "GAMLSS refit warnings: %d (%.1f%%)",
    n_lng_refit, 100 * n_lng_refit / n_lng_models
  )
  log_info(
    "Potential extrapolation: Age=%d (%.1f%%), ICC=%d (%.1f%%)",
    n_lng_age_extrap, 100 * n_lng_age_extrap / n_lng_models,
    n_lng_icc_extrap, 100 * n_lng_icc_extrap / n_lng_models
  )

  # Check correlation between refit warnings and stability quality
  if (n_lng_refit > 0) {
    corr_refit <- mean(
      stability_results.dt[REFIT_WARNING == TRUE]$CORR, na.rm = TRUE
    )
    corr_no_refit <- mean(
      stability_results.dt[REFIT_WARNING == FALSE]$CORR, na.rm = TRUE
    )
    log_info(
      "Mean correlation: with refit warning = %.3f, without = %.3f",
      corr_refit, corr_no_refit
    )
  }

  # Store summary data.tables in gamlss.lst for easy access
  gamlss.lst$SUMMARIES <- list(
    # Model comparison (all fitted models with AIC/BIC)
    MODEL_COMPARISON = if (exists("model_comparison.dt")) {
      model_comparison.dt
    } else {
      NULL
    },

    # Hold-out validation summaries
    VALIDATION_RESULTS = if (exists("test_results.dt")) {
      test_results.dt
    } else {
      NULL
    },
    CALIBRATION_DETAILS = if (exists("calibration_details.dt")) {
      calibration_details.dt
    } else {
      NULL
    },

    # Longitudinal stability summaries
    STABILITY_OVERALL = stability_results.dt,
    CENTILE_STABILITY = if (exists("centile_stability.dt")) {
      centile_stability.dt
    } else {
      NULL
    },
    CENTILE_STRATIFIED = if (exists("centile_stratified.dt")) {
      centile_stratified.dt
    } else {
      NULL
    },
    LONGITUDINAL_CENTILES = if (nrow(longitudinal_centiles.dt) > 0) {
      longitudinal_centiles.dt
    } else {
      NULL
    }
  )

  log_info("Saving GAMLSS data")
  n_summaries <- sum(!sapply(gamlss.lst$SUMMARIES, is.null))
  log_info("  - %d summary data.tables included", n_summaries)
  write_rds_safe(gamlss.lst, gamlss.path, "GAMLSS data & fits")
  rm(roi, adj, side, sex, gamlss.subdt)
}

# ----- Generate Normative Tables -----
log_info("Normative tables")
age_range <- get_parameter("gamlss", "age_range")
age_seq <- seq(age_range[1], age_range[2], by = 1)
# centiles <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
centiles <- get_parameter("gamlss", "centiles")
norm_tables.lst <- list()
pb <- progress_bar$new(
  format = "GAMLSS (longitudinal) | :what [:bar] :current/:total\n",
  total = sort.dt[, .N], clear = FALSE, width = 100, show_after = 0
)
for (i in seq_len(nrow(sort.dt))) {
  row.dt <- sort.dt[i]
  sex <- row.dt$SEX
  roi <- row.dt$ROI
  adj <- row.dt$ADJ
  side <- row.dt$SIDE

  # Final model
  final_mod.lst <- gamlss.lst$FINAL[[sex]][[roi]][[adj]][[side]]
  if (is.null(final_mod.lst)) {
    pb$tick()
    next
  }

  mod.fit <- final_mod.lst$FIT
  gamlss_final.subdt <- final_mod.lst$DATA

  pb$tick(tokens = list(what = sprintf(
    "%s — %s (%s; %s): %s",
    sex, roi, side, adj, fam_names.v[[final_mod.lst$FAMILY]]
  )))

  # --- Diagnostic: Check for potential extrapolation ---
  # Compare prediction grid vs training data ranges
  train_age_range <- range(gamlss_final.subdt$AGE)
  pred_age_range <- range(age_seq)
  train_icc <- mean(gamlss_final.subdt$ICC)
  train_educ <- mean(gamlss_final.subdt$EDUC_num)

  # Flag extrapolation if prediction grid extends beyond training range
  norm_age_extrap <- pred_age_range[1] < train_age_range[1] ||
    pred_age_range[2] > train_age_range[2]

  if (norm_age_extrap && final_mod.lst$FAMILY %in% c("BCCG", "BE")) {
    log_warn(
      paste(
        "%s %s (%s; %s): Prediction grid may extrapolate",
        "(Age: [%.1f-%.1f] vs training [%.1f-%.1f])"
      ),
      sex, roi, side, adj,
      pred_age_range[1], pred_age_range[2],
      train_age_range[1], train_age_range[2]
    )
  }

  # Generate centile predictions (Extract subsample N)
  # First create table with ages that exist in data AND are in age_seq
  cent_table.dt <- gamlss_final.subdt[
    order(AGE), .N, .(AGE = round(AGE))
  ][AGE %in% age_seq]

  # Create prediction data matching cent_table.dt rows exactly
  # Use most frequent site for prediction (random effect will be ~0 for population average)
  # Transfer models (no SITE) don't have SITE column
  has_site <- "SITE" %in% names(gamlss_final.subdt)
  if (has_site) {
    mode_site <- gamlss_final.subdt[, .N, SITE][order(-N), SITE[1]]
    new.subdt <- data.table(
      AGE = cent_table.dt$AGE,
      EDUC_num = mean(gamlss_final.subdt$EDUC_num),
      ICC = mean(gamlss_final.subdt$ICC),
      SITE = mode_site
    )
  } else {
    new.subdt <- data.table(
      AGE = cent_table.dt$AGE,
      EDUC_num = mean(gamlss_final.subdt$EDUC_num),
      ICC = mean(gamlss_final.subdt$ICC)
    )
  }
  if (final_mod.lst$FAMILY == "NO") {
    # Gaussian distribution
    cent_table.dt[
      , paste0("p", centiles * 100) := lapply(centiles, \(p) {
        qNO(
          p = p,
          mu = predict(final_mod.lst$FIT, new = new.subdt, type = "response"),
          sigma = predict(
            final_mod.lst$FIT, "sigma", new = new.subdt, type = "response"
          )
        )
      })
    ]
  } else if (final_mod.lst$FAMILY == "L_NO") {
    # Logit-Normal distribution
    cent_table.dt[
      , paste0("p", centiles * 100) := lapply(centiles, \(p) {
        qNO(
          p = p,
          mu = predict(final_mod.lst$FIT, new = new.subdt, type = "response"),
          sigma = predict(
            final_mod.lst$FIT, "sigma", new = new.subdt, type = "response"
          )
        ) |> plogis()
      })
    ]
  } else if (final_mod.lst$FAMILY == "BE") {
    # Beta distribution
    # Capture warnings about re-fitting for safe predictions
    norm_refit_warning <- FALSE
    withCallingHandlers(
      {
        cent_table.dt[
          , paste0("p", centiles * 100) := lapply(centiles, \(p) {
            qBE(
              p = p,
              mu = predict(
                final_mod.lst$FIT, new = new.subdt, type = "response"
              ),
              sigma = predict(
                final_mod.lst$FIT, "sigma", new = new.subdt, type = "response"
              )
            )
          })
        ]
      },
      warning = function(w) {
        if (grepl("discrepancy.*re-fit", w$message)) {
          norm_refit_warning <<- TRUE
          log_debug(
            "%s %s (%s; %s): GAMLSS refit warning (norm table)",
            sex, roi, side, adj
          )
        }
        invokeRestart("muffleWarning")
      }
    )
  } else if (final_mod.lst$FAMILY == "BCCG") {
    # Box-Cox Cole & Green distribution
    # Capture warnings about re-fitting for safe predictions
    norm_refit_warning <- FALSE
    withCallingHandlers(
      {
        cent_table.dt[
          , paste0("p", centiles * 100) := lapply(centiles, \(p) {
            qBCCG(
              p = p,
              mu = predict(
                final_mod.lst$FIT, new = new.subdt, type = "response"
              ),
              sigma = predict(
                final_mod.lst$FIT, "sigma", new = new.subdt, type = "response"
              ),
              nu = predict(
                final_mod.lst$FIT, "nu", new = new.subdt, type = "response"
              )
            )
          })
        ]
      },
      warning = function(w) {
        if (grepl("discrepancy.*re-fit", w$message)) {
          norm_refit_warning <<- TRUE
          log_debug(
            "%s %s (%s; %s): GAMLSS refit warning (norm table)",
            sex, roi, side, adj
          )
        }
        invokeRestart("muffleWarning")
      }
    )
  }

  # Initialize refit warning for other distributions
  if (!exists("norm_refit_warning")) norm_refit_warning <- FALSE

  # Store diagnostics in model structure
  gamlss.lst$FINAL[[sex]][[roi]][[adj]][[side]]$NORM_DIAG <- list(
    REFIT_WARNING = norm_refit_warning,
    AGE_EXTRAP = norm_age_extrap,
    PRED_AGE_RANGE = pred_age_range,
    TRAIN_AGE_RANGE = train_age_range
  )

  # Add metadata for table generation
  cent_table.dt[, FAMILY := fam_names.v[final_mod.lst$FAMILY]]

  # Store in list
  norm_tables.lst[[sex]][[roi]][[adj]][[side]] <- cent_table.dt
}

# --- Diagnostic summary ---
log_section("Normative table generation diagnostics")

# Extract diagnostics
norm_diag.dt <- sort.dt[, {
  diag.lst <- gamlss.lst$FINAL[[SEX]][[ROI]][[ADJ]][[SIDE]]$NORM_DIAG
  if (is.null(diag.lst)) NULL else list(
    REFIT_WARNING = diag.lst$REFIT_WARNING,
    AGE_EXTRAP = diag.lst$AGE_EXTRAP
  )
}, by = .(SEX, ROI, ADJ, SIDE)]

# Summary statistics
n_norm_models <- nrow(norm_diag.dt)
n_norm_refit <- sum(norm_diag.dt$REFIT_WARNING, na.rm = TRUE)
n_norm_age_extrap <- sum(norm_diag.dt$AGE_EXTRAP, na.rm = TRUE)

log_info("Total normative tables generated: %d", n_norm_models)
log_info(
  "GAMLSS refit warnings: %d (%.1f%%)",
  n_norm_refit, 100 * n_norm_refit / n_norm_models
)
log_info(
  "Potential age extrapolation: %d (%.1f%%)",
  n_norm_age_extrap, 100 * n_norm_age_extrap / n_norm_models
)

# Extract as data.table
norm_tables.dt <- sort.dt[
  , {
    temp.subdt <- norm_tables.lst[[SEX]][[ROI]][[ADJ]][[SIDE]]
    if (is.null(temp.subdt)) NULL else temp.subdt
  },
  by = .(SEX, ROI, ADJ, SIDE)
]

# Create and save publication-ready normative tables
log_section("Creating publication-ready normative tables")

# Get normative combinations from config (clinical use only)
# Other adjustments (PRP, STX, RES) are used for sex differences analysis but
# not provided as clinical norms to align with manuscript recommendation of HVR
norm_combos <- get_parameter("gamlss", "normative_combinations")
if (is.null(norm_combos)) {
  # Fallback: all combinations if not specified
  norm_combos <- lapply(seq_len(nrow(sort.dt)), function(i) {
    list(roi = sort.dt[i, ROI], adj = sort.dt[i, ADJ])
  })
  log_warn("No normative_combinations in config; generating all tables")
} else {

  log_info("Generating normative tables for clinical use:")
  for (combo in norm_combos) {
    log_info("  - %s (%s adjustment)", combo$roi, combo$adj)
  }
}

# Filter sort.dt to only normative combinations
norm_sort.dt <- rbindlist(lapply(norm_combos, function(combo) {
  sort.dt[ROI == combo$roi & ADJ == combo$adj]
}))

# Get table directory from config (html and tex in same directory)
table_dir <- file.path(get_output_path("tables"), "gamlss_normative")
needs_tables <- REDO_TABLES || !dir.exists(table_dir)
if (needs_tables) {
  ensure_directory(table_dir)

  pb <- progress_bar$new(
    format = "Creating tables | :what [:bar] :current/:total\n",
    total = norm_sort.dt[, .N], clear = FALSE, width = 100, show_after = 0
  )
  pb$tick(0)

  for (i in seq_len(nrow(norm_sort.dt))) {
    row.dt <- norm_sort.dt[i]
    sex <- row.dt$SEX
    roi <- row.dt$ROI
    adj <- row.dt$ADJ
    side <- row.dt$SIDE

    if (is.null(norm_tables.lst[[sex]][[roi]][[adj]][[side]])) {
      pb$tick()
      next
    }

    norm_table.dt <- norm_tables.lst[[sex]][[roi]][[adj]][[side]]

    pb$tick(tokens = list(what = sprintf(
      "%s — %s (%s; %s)", sex, roi, side, adj
    )))

    # Separate title and subtitle
    side_label <- fcase(
      side == "L", "Left",
      side == "R", "Right",
      default = "Bilateral"
    )
    table_title <- ROIS[roi]

    # For HVR && Proportions, omit adjustment label (HVR is inherently a proportion)
    if (roi == "HVR" && adj == "PRP") {
      table_subtitle <- sprintf("%s, %s", sex, side_label)
    } else {
      table_subtitle <- sprintf("%s, %s, %s", sex, ADJS[adj], side_label)
    }

    # Create HTML table
    html_table <- create_normative_table(
      norm_table.dt,
      title = table_title,
      subtitle = table_subtitle,
      roi_type = roi,
      adj_type = adj,
      family = TRUE,
      sample_size = TRUE,
      add_clinical_labels = TRUE,
      format = "html"
    )

    # Save HTML
    filename_base <- sprintf(
      "normative_%s_%s_%s_%s",
      roi, adj, side, tolower(sex)
    )
    gt::gtsave(
      html_table,
      file.path(table_dir, paste0(filename_base, ".html"))
    )

    # Create LaTeX table
    latex_table <- create_normative_table(
      norm_table.dt,
      title = table_title,
      subtitle = table_subtitle,
      roi_type = roi,
      adj_type = adj,
      family = TRUE,
      sample_size = TRUE,
      add_clinical_labels = TRUE,
      format = "latex"
    )

    # Save LaTeX and fix source notes formatting
    tex_file <- file.path(table_dir, paste0(filename_base, ".tex"))
    gt::gtsave(latex_table, tex_file)
    fix_latex_source_notes(tex_file)
  }

  log_info("Normative tables saved to:")
  log_info("  - Output directory: %s", table_dir)

  # Compile all normative tables into single organized document
  log_info("Compiling comprehensive normative tables document")
  compile_path <- compile_normative_tables(
    table_dir,
    output_path = file.path(get_output_path("tables"), "gamlss_report"),
    title = "UK Biobank Normative Centile Tables"
  )
  log_info("  - Compiled document: %s", compile_path)
}

# --- Create manuscript summary tables ---
summary_table_dir <- file.path(get_output_path("tables"), "gamlss_summary")
needs_summary <- REDO_TABLES || !dir.exists(summary_table_dir)
if (needs_summary) {
  log_section("Creating manuscript summary tables (year-bins)")
  ensure_directory(summary_table_dir)

  # Configuration for summary tables
  # summary_centiles <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  summary_centiles <- get_parameter("gamlss", "summary_centiles")
  age_bin_width <- get_parameter("gamlss", "age_bin")

  # Use normative combinations (same as detailed tables)
  sort.subdt <- norm_sort.dt[, unique(.SD), .SDcols = c("ROI", "ADJ")]

  pb <- progress_bar$new(
    format = "Creating summary tables | :what [:bar] :current/:total\n",
    total = sort.subdt[, .N], clear = FALSE, width = 100, show_after = 0
  )
  pb$tick(0)

  # Generate one summary table per ROI and adjustment method
  # Try to combine both sides and sexes if it fits
  for (i in seq_len(nrow(sort.subdt))) {
    row.dt <- sort.subdt[i]
    roi <- row.dt$ROI
    adj <- row.dt$ADJ

    pb$tick(tokens = list(what = sprintf("%s × %s", ROIS[roi], ADJS[adj])))

    # Collect data for all sex/side combinations
    summary_tables.lst <- list()
    for (sex in sort.dt[, unique(SEX)]) {
      for (side in sort.dt[!"LR", on = "SIDE", unique(SIDE)]) {
        key <- paste(sex, side, sep = "_")
        table_data.subdt <- norm_tables.lst[[sex]][[roi]][[adj]][[side]]
        if (!is.null(table_data.subdt)) {
          summary_tables.lst[[key]] <- table_data.subdt
        }
      }
    }

    if (length(summary_tables.lst) == 0) next

    # # Determine combine mode based on number of combinations
    # # If we have all 4 (Male_L, Male_R, Female_L, Female_R), try "both"
    # # For LaTeX, might be too wide, so we'll create multiple versions
    # n_combos <- length(summary_tables.lst)

    # Simplified title (no redundant info)
    table_title <- sprintf("%s (%s)", ROIS[roi], ADJS[adj])

    # HTML version: hierarchical structure with side spanners
    html_table <- create_summary_normative_table(
      summary_tables.lst,
      title = table_title,
      subtitle = NULL,
      format = "html",
      bin_width = age_bin_width,
      centiles = summary_centiles,
      roi_type = roi,
      adj_type = adj,
      sample_size = TRUE,
      add_clinical_labels = TRUE
    )
    # Save HTML
    filename_base <- sprintf("summary_%s_%s", roi, adj)
    gt::gtsave(
      html_table,
      file.path(summary_table_dir, paste0(filename_base, ".html"))
    )

    # LaTeX version: minimal formatting with side spanners
    latex_table <- create_summary_normative_table(
      summary_tables.lst,
      title = table_title,
      subtitle = NULL,
      format = "latex",
      bin_width = age_bin_width,
      centiles = summary_centiles,
      roi_type = roi,
      adj_type = adj,
      sample_size = TRUE,
      add_clinical_labels = TRUE
    )

    tex_path <- file.path(summary_table_dir, paste0(filename_base, ".tex"))
    gt::gtsave(latex_table, tex_path)
    fix_latex_source_notes(tex_path)

    # Create standalone LaTeX document - no extra title
    wrap_latex_table(
      tex_path,
      output_path = file.path(
        summary_table_dir,
        paste0(filename_base, "_standalone.tex")
      )
    )

    # --- Bilateral summary tables (side == "LR") ---
    # Collect bilateral data
    bilateral_tables.lst <- list()
    for (sex in sort.dt[, unique(SEX)]) {
      key <- paste(sex, "LR", sep = "_")
      table_data.subdt <- norm_tables.lst[[sex]][[roi]][[adj]][["LR"]]
      if (!is.null(table_data.subdt)) {
        bilateral_tables.lst[[key]] <- table_data.subdt
      }
    }

    if (length(bilateral_tables.lst) > 0) {
      # Bilateral tables - simpler structure (no L/R columns)
      bilat_title <- sprintf("%s (%s) - Bilateral", ROIS[roi], ADJS[adj])

      # For bilateral, we just need Sex grouping without side columns
      # Process bilateral data differently - no need for side spanners
      bilat_processed.lst <- lapply(names(bilateral_tables.lst), function(key) {
        dt <- copy(bilateral_tables.lst[[key]])
        sex <- strsplit(key, "_")[[1]][1]
        dt[, AGE_BIN := bin_ages(AGE, age_bin_width)]
        cent_cols <- paste0("p", summary_centiles * 100)

        # Scale proportions if needed
        if (adj == "PRP" && roi %in% c("HC", "LV")) {
          dt[, (cent_cols) := lapply(.SD, `*`, 1000), .SDcols = cent_cols]
        }

        binned <- merge(
          dt[, .(N = sum(N)), "AGE_BIN"],
          dt[, lapply(.SD, mean), AGE_BIN, .SDcols = cent_cols]
        )
        binned[, SEX := sex]
        binned
      })

      bilat_combined.dt <- rbindlist(bilat_processed.lst)
      setorder(bilat_combined.dt, SEX, AGE_BIN)

      # Create sex labels with N
      bilat_sex_n <- bilat_combined.dt[, .(N_total = sum(N)), SEX]
      bilat_combined.dt[bilat_sex_n, SEX_LABEL := sprintf("%ss (N = %s)", SEX, format(i.N_total, big.mark = ",")), on = "SEX"]

      # Create Age label with N
      bilat_combined.dt[, AGE_LABEL := sprintf("%s (%d)", AGE_BIN, N)]
      bilat_combined.dt[, c("AGE_BIN", "N", "SEX") := NULL]

      cent_cols <- paste0("p", summary_centiles * 100)

      bilat_gt <- bilat_combined.dt |>
        gt(groupname_col = "SEX_LABEL", rowname_col = "AGE_LABEL") |>
        tab_stubhead(label = "Age (N)") |>
        tab_header(title = bilat_title) |>
        fmt_number(columns = all_of(cent_cols), decimals = 2)

      # Add centile labels
      cent_labels <- setNames(
        paste0(as.numeric(sub("p", "", cent_cols)), "th"),
        cent_cols
      )
      bilat_gt <- bilat_gt |> cols_label(.list = cent_labels)

      # Source notes
      bilat_notes <- c(
        "Age in years",
        sprintf("Values represent means within %d-year age bins", age_bin_width),
        "Bilateral = sum of left and right hemispheres"
      )
      if (adj == "PRP" && roi %in% c("HC", "LV")) {
        bilat_notes <- c(bilat_notes[1],
          "Values \u00d7 10\u207b\u00b3 (multiply by 0.001 for proportions)",
          bilat_notes[-1])
      }
      for (note in bilat_notes) {
        bilat_gt <- bilat_gt |> tab_source_note(source_note = note)
      }

      # Save bilateral HTML
      bilat_filename <- sprintf("summary_%s_%s_bilateral", roi, adj)
      gt::gtsave(bilat_gt, file.path(summary_table_dir, paste0(bilat_filename, ".html")))

      # LaTeX version with styling
      bilat_gt_tex <- bilat_gt |>
        tab_options(
          latex.use_longtable = TRUE,
          table.font.size = px(9),
          heading.title.font.size = px(11)
        )
      bilat_tex_path <- file.path(summary_table_dir, paste0(bilat_filename, ".tex"))
      gt::gtsave(bilat_gt_tex, bilat_tex_path)
      fix_latex_source_notes(bilat_tex_path)
    }
  }
  log_info("Manuscript summary tables saved to: %s", summary_table_dir)
}

# Save normative tables (R objects)
log_section("Saving outputs")

# Validate normative tables before saving
n_tables <- sum(sapply(norm_tables.lst, function(sex) {
  sum(sapply(sex, function(roi) {
    sum(sapply(roi, function(adj) {
      sum(!sapply(adj, is.null))
    }))
  }))
}))
if (n_tables == 0) {
  log_error("No normative tables generated")
  stop("Failed to generate normative tables", call. = FALSE)
}
log_info("Validated %d normative tables for saving", n_tables)

write_rds_safe(
  norm_tables.lst,
  get_data_path("processed", "norm_tables"),
  description = "Normative centile tables"
)

# ----- Calculate Z-scores -----
log_section("Calculating Z-scores from GAMLSS models")

# Create new data.table for Z-scores (don't overwrite original)
zscores_data.dt <- copy(gamlss.lst$DATA$CRS)
zscores_data.dt[, ZSCORE := NA_real_]

pb <- progress_bar$new(
  format = "Z-scores | :what [:bar] :current/:total\n",
  total = sort.dt[, .N], clear = FALSE, width = 100, show_after = 0
)

for (i in seq_len(nrow(sort.dt))) {
  row.dt <- sort.dt[i]
  sex <- row.dt$SEX
  roi <- row.dt$ROI
  adj <- row.dt$ADJ
  side <- row.dt$SIDE

  pb$tick(tokens = list(what = sprintf(
    "%s — %s (%s; %s)", sex, roi, side, adj
  )))

  # Get final model
  if (is.null(gamlss.lst$FINAL[[sex]][[roi]][[adj]][[side]])) next

  final_mod.lst <- gamlss.lst$FINAL[[sex]][[roi]][[adj]][[side]]
  gamlss_final.subdt <- final_mod.lst$DATA
  fam <- final_mod.lst$FAMILY
  mod.fit <- final_mod.lst$FIT

  # Get data subset
  data.subdt <- zscores_data.dt[row.dt][VAL > 0][!is.na(VAL)]
  if (!nrow(data.subdt)) next

  # --- Diagnostic: Check for potential extrapolation ---
  train_data <- final_mod.lst$DATA
  train_age_range <- range(train_data$AGE)
  data_age_range <- range(data.subdt$AGE)
  train_icc_range <- range(train_data$ICC)
  data_icc_range <- range(data.subdt$ICC)

  zscore_age_extrap <- data_age_range[1] < train_age_range[1] ||
    data_age_range[2] > train_age_range[2]
  zscore_icc_extrap <- data_icc_range[1] < train_icc_range[1] ||
    data_icc_range[2] > train_icc_range[2]

  if ((zscore_age_extrap || zscore_icc_extrap) && fam %in% c("BCCG", "BE")) {
    log_debug(
      paste(
        "%s %s (%s; %s): Z-score data may extrapolate",
        "(Age: [%.1f-%.1f] vs [%.1f-%.1f], ICC: [%.0f-%.0f] vs [%.0f-%.0f])"
      ),
      sex, roi, side, adj,
      data_age_range[1], data_age_range[2],
      train_age_range[1], train_age_range[2],
      data_icc_range[1], data_icc_range[2],
      train_icc_range[1], train_icc_range[2]
    )
  }

  # Calculate Z-scores using fitted model with warning capture
  # Use .SD approach like temporal validation to ensure proper column selection
  # Transfer models don't have SITE, so adjust column selection accordingly
  is_transfer <- isTRUE(final_mod.lst$TRANSFER)
  pred_cols_base <- c("AGE", "EDUC_num", "ICC")
  pred_cols <- if (is_transfer) pred_cols_base else c(pred_cols_base, "SITE")

  zscore_refit_warning <- FALSE
  withCallingHandlers(
    {
      if (fam == "NO") {
        # Gaussian distribution - calculate percentile then convert to z-score
        data.subdt[
          ,
          ZSCORE := qnorm(
            pNO(
              VAL,
              mu = predict(mod.fit, new = .SD, type = "response"),
              sigma = predict(mod.fit, "sigma", new = .SD, type = "response")
            )
          ),
          .SDcols = c(pred_cols, "VAL")
        ]
      } else if (fam == "L_NO") {
        # Logit-Normal distribution - transform then calculate
        eps <- 1e-6
        data.subdt[, VAL_logit := VAL |> pmin(1 - eps) |> pmax(eps) |> qlogis()]
        data.subdt[
          ,
          ZSCORE := qnorm(
            pNO(
              VAL_logit,
              mu = predict(mod.fit, new = .SD, type = "response"),
              sigma = predict(mod.fit, "sigma", new = .SD, type = "response")
            )
          ),
          .SDcols = c(pred_cols, "VAL_logit")
        ]
      } else if (fam == "BE") {
        # Beta distribution
        data.subdt[
          ,
          ZSCORE := qnorm(
            pBE(
              VAL,
              mu = predict(mod.fit, new = .SD, type = "response"),
              sigma = predict(mod.fit, "sigma", new = .SD, type = "response")
            )
          ),
          .SDcols = c(pred_cols, "VAL")
        ]
      } else if (fam == "BCCG") {
        # Box-Cox Cole & Green distribution
        data.subdt[
          ,
          ZSCORE := qnorm(
            pBCCG(
              VAL,
              mu = predict(mod.fit, new = .SD, type = "response"),
              sigma = predict(mod.fit, "sigma", new = .SD, type = "response"),
              nu = predict(mod.fit, "nu", new = .SD, type = "response")
            )
          ),
          .SDcols = c(pred_cols, "VAL")
        ]
      }
    },
    warning = function(w) {
      if (grepl("discrepancy.*re-fit", w$message)) {
        zscore_refit_warning <<- TRUE
        log_debug(
          "%s %s (%s; %s): GAMLSS refit warning (z-scores)",
          sex, roi, side, adj
        )
      }
      invokeRestart("muffleWarning")
    }
  )

  # Store diagnostics
  if (!exists("zscore_diag.lst")) zscore_diag.lst <- list()
  zscore_diag.lst[[sex]][[roi]][[adj]][[side]] <- list(
    REFIT_WARNING = zscore_refit_warning,
    AGE_EXTRAP = zscore_age_extrap,
    ICC_EXTRAP = zscore_icc_extrap,
    N = nrow(data.subdt)
  )

  # Assign Z-scores back to z-scores data
  zscores_data.dt[
    data.subdt,
    on = c("EID", "INST", "SEX", "SIDE", "ADJ", "ROI"),
    ZSCORE := i.ZSCORE
  ]
}

# --- Diagnostic summary ---
log_section("Z-score calculation diagnostics")

# Extract diagnostics
zscore_diag.dt <- sort.dt[, {
  diag.lst <- zscore_diag.lst[[SEX]][[ROI]][[ADJ]][[SIDE]]
  if (is.null(diag.lst)) NULL else list(
    REFIT_WARNING = diag.lst$REFIT_WARNING,
    AGE_EXTRAP = diag.lst$AGE_EXTRAP,
    ICC_EXTRAP = diag.lst$ICC_EXTRAP,
    N = diag.lst$N
  )
}, by = .(SEX, ROI, ADJ, SIDE)]

# Summary statistics
n_zscore_models <- nrow(zscore_diag.dt)
n_zscore_refit <- sum(zscore_diag.dt$REFIT_WARNING, na.rm = TRUE)
n_zscore_age_extrap <- sum(zscore_diag.dt$AGE_EXTRAP, na.rm = TRUE)
n_zscore_icc_extrap <- sum(zscore_diag.dt$ICC_EXTRAP, na.rm = TRUE)
total_observations <- sum(zscore_diag.dt$N, na.rm = TRUE)

log_info("Total z-score calculations: %d models", n_zscore_models)
log_info("Total observations scored: %d", total_observations)
log_info(
  "GAMLSS refit warnings: %d (%.1f%%)",
  n_zscore_refit, 100 * n_zscore_refit / n_zscore_models
)
log_info(
  "Potential extrapolation: Age=%d (%.1f%%), ICC=%d (%.1f%%)",
  n_zscore_age_extrap, 100 * n_zscore_age_extrap / n_zscore_models,
  n_zscore_icc_extrap, 100 * n_zscore_icc_extrap / n_zscore_models
)

# Save z-scores
write_rds_safe(
  zscores_data.dt,
  get_data_path("processed", "z_scores"),
  description = "Z-scores from GAMLSS models"
)

# Save additional diagnostic summaries
log_section("Saving diagnostic summaries")
gamlss_diagnostics.lst <- list(
  NORMATIVE_DIAGNOSTICS = if (exists("norm_diag.dt")) norm_diag.dt else NULL,
  ZSCORE_DIAGNOSTICS = if (exists("zscore_diag.dt")) zscore_diag.dt else NULL
)
# Note: Model comparison, validation, and stability summaries
# are already saved in gamlss.lst$SUMMARIES
log_info("Saving normative table and z-score diagnostic summaries")
write_rds_safe(
  gamlss_diagnostics.lst,
  get_data_path("models", "diagnostics", "gamlss"),
  description = "Additional GAMLSS diagnostics"
)

# ----- Create Diagnostic Summary Tables -----
log_section("Creating diagnostic summary tables")

# Get table directory from config
summary_table_dir <- file.path(get_output_path("tables"), "gamlss_summary")
ensure_directory(summary_table_dir)

log_info("Preparing diagnostic tables for publication")

# --- Model Validation Tables (Reorganized) ---
log_info("Creating model validation summary tables")

# Get abbreviation helpers
adj_abbrev <- get_adjustment_abbrev()
fam_info <- get_family_info()
interp_symbols <- get_interpretation_symbols()

if (!is.null(gamlss.lst$SUMMARIES$VALIDATION_RESULTS)) {
  log_info("  - Extracting validation results from %d model combinations", nrow(sort.dt))
  validation_data.dt <- sort.dt[, {
    test.lst <- gamlss.lst$TEST[[SEX]][[ROI]][[ADJ]][[SIDE]]
    if (is.null(test.lst)) NULL else list(
      MAE = round(test.lst$MAE, 2),
      INTERPRETATION = test.lst$INTERPRETATION
    )
  }, by = .(SEX, ROI, ADJ, SIDE)]

  # Unlist the columns created by the by operation
  if (nrow(validation_data.dt) > 0) {
    validation_data.dt[, MAE := unlist(MAE)]
    validation_data.dt[, INTERPRETATION := unlist(INTERPRETATION)]
  }

  if (nrow(validation_data.dt) > 0) {
    # --- UNILATERAL TABLES (L/R only) ---
    val_unilat.dt <- validation_data.dt[SIDE != "LR"]

    if (nrow(val_unilat.dt) > 0) {
      # Ensure all columns are atomic (not list) before processing
      for (col in c("SEX", "ROI", "ADJ", "SIDE", "MAE", "INTERPRETATION")) {
        if (col %in% names(val_unilat.dt) && is.list(val_unilat.dt[[col]])) {
          val_unilat.dt[, (col) := unlist(get(col))]
        }
      }

      # Add interpretation symbols to MAE values
      val_unilat.dt[, MAE_SYMBOL := paste0(MAE, interp_symbols[INTERPRETATION])]

      # Use ROI labels and adjustment abbreviations - force to character
      val_unilat.dt[, ROI_LABEL := as.character(ROIS[ROI])]
      val_unilat.dt[, ADJ_LABEL := as.character(adj_abbrev[ADJ])]

      # Restructure: ROI grouping, ADJ rows, Sex×Side columns (MAE only, with symbols)
      val_unilat.dt[, COL_ID := as.character(paste(SEX, SIDE, sep = "_"))]

      validation_wide.dt <- dcast(
        val_unilat.dt,
        ROI_LABEL + ADJ_LABEL ~ COL_ID,
        value.var = "MAE_SYMBOL"
      )

      # Order columns
      col_order <- c("Female_L", "Female_R", "Male_L", "Male_R")
      col_order <- col_order[col_order %in% names(validation_wide.dt)]
      setcolorder(validation_wide.dt, c("ROI_LABEL", "ADJ_LABEL", col_order))

      # Create gt table
      validation_gt <- validation_wide.dt |>
        gt(groupname_col = "ROI_LABEL", rowname_col = "ADJ_LABEL") |>
        tab_stubhead(label = "Adj.") |>
        tab_header(
          title = "GAMLSS Model Validation: Hold-out Test Set Performance",
          subtitle = "Mean Absolute Error (MAE) with quality indicators"
        )

      # Add sex spanners
      for (sex in c("Female", "Male")) {
        sex_cols <- grep(paste0("^", sex, "_"), names(validation_wide.dt), value = TRUE)
        if (length(sex_cols) > 0) {
          validation_gt <- validation_gt |>
            tab_spanner(label = sex, columns = all_of(sex_cols), level = 2)
        }
      }

      # Add side spanners (use unique IDs to avoid conflict with column names)
      for (sex in c("Female", "Male")) {
        for (side in c("L", "R")) {
          side_col <- paste0(sex, "_", side)
          if (side_col %in% names(validation_wide.dt)) {
            validation_gt <- validation_gt |>
              tab_spanner(
                label = ifelse(side == "L", "L", "R"),
                columns = all_of(side_col),
                level = 1,
                id = paste0("spanner_", sex, "_", side)
              )
          }
        }
      }

      # Column labels (just show side)
      col_labels <- setNames(rep("MAE", length(col_order)), col_order)
      validation_gt <- validation_gt |> cols_label(.list = col_labels)

      # Source notes (each on separate line)
      val_notes <- c(
        get_abbrev_legend("interpretation"),
        get_abbrev_legend("adjustment"),
        "L = Left hemisphere; R = Right hemisphere"
      )
      for (note in val_notes) {
        validation_gt <- validation_gt |> tab_source_note(source_note = note)
      }

      save_table(validation_gt, "gamlss_validation_unilateral", summary_table_dir, formats = c("html", "tex"))
      log_info("  - Model validation table created (unilateral)")
    }

    # --- BILATERAL TABLES (LR only) ---
    val_bilat.dt <- validation_data.dt[SIDE == "LR"]

    if (nrow(val_bilat.dt) > 0) {
      # Ensure all columns are atomic (not list) before processing
      for (col in c("SEX", "ROI", "ADJ", "SIDE", "MAE", "INTERPRETATION")) {
        if (col %in% names(val_bilat.dt) && is.list(val_bilat.dt[[col]])) {
          val_bilat.dt[, (col) := unlist(get(col))]
        }
      }

      val_bilat.dt[, MAE_SYMBOL := paste0(MAE, interp_symbols[INTERPRETATION])]
      val_bilat.dt[, ROI_LABEL := as.character(ROIS[ROI])]
      val_bilat.dt[, ADJ_LABEL := as.character(adj_abbrev[ADJ])]

      val_bilat_wide.dt <- dcast(
        val_bilat.dt,
        ROI_LABEL + ADJ_LABEL ~ SEX,
        value.var = "MAE_SYMBOL"
      )

      bilat_val_gt <- val_bilat_wide.dt |>
        gt(groupname_col = "ROI_LABEL", rowname_col = "ADJ_LABEL") |>
        tab_stubhead(label = "Adj.") |>
        tab_header(
          title = "GAMLSS Model Validation: Bilateral (Sum)",
          subtitle = "Mean Absolute Error (MAE) with quality indicators"
        ) |>
        cols_label(Female = "Female", Male = "Male")

      bilat_val_notes <- c(
        get_abbrev_legend("interpretation"),
        get_abbrev_legend("adjustment"),
        "Bilateral = sum of left and right hemispheres"
      )
      for (note in bilat_val_notes) {
        bilat_val_gt <- bilat_val_gt |> tab_source_note(source_note = note)
      }

      save_table(bilat_val_gt, "gamlss_validation_bilateral", summary_table_dir, formats = c("html", "tex"))
      log_info("  - Model validation table created (bilateral)")
    }
  }
}

# --- Model Summary Tables (Reorganized) ---
log_info("Creating model summary tables")

# Main summary table: Final selected models only
log_info("  - Extracting final model summaries from %d combinations", nrow(sort.dt))
summary_data.dt <- sort.dt[, {
  mod.lst <- gamlss.lst$FINAL[[SEX]][[ROI]][[ADJ]][[SIDE]]
  if (is.null(mod.lst)) NULL else list(
    FAMILY = mod.lst$FAMILY,  # Keep code, convert to abbrev later
    AIC = round(mod.lst$AIC, 1),
    BIC = round(mod.lst$BIC, 1)
  )
}, by = .(SEX, ROI, ADJ, SIDE)]

if (nrow(summary_data.dt) > 0) {
  # Ensure all columns are atomic (not list) before processing
  for (col in c("SEX", "ROI", "ADJ", "SIDE", "FAMILY", "AIC", "BIC")) {
    if (col %in% names(summary_data.dt) && is.list(summary_data.dt[[col]])) {
      summary_data.dt[, (col) := unlist(get(col))]
    }
  }

  # --- UNILATERAL TABLES (L/R only) ---
  summ_unilat.dt <- summary_data.dt[SIDE != "LR"]

  if (nrow(summ_unilat.dt) > 0) {
    # Use abbreviations - force to character
    summ_unilat.dt[, FAMILY_ABBREV := as.character(fam_info$abbrev[FAMILY])]
    summ_unilat.dt[, ROI_LABEL := as.character(ROIS[ROI])]
    summ_unilat.dt[, ADJ_LABEL := as.character(adj_abbrev[ADJ])]

    # Restructure: ROI grouping, ADJ rows, Sex×Side columns
    summ_long.dt <- melt(
      summ_unilat.dt,
      id.vars = c("SEX", "ROI_LABEL", "ADJ_LABEL", "SIDE"),
      measure.vars = c("FAMILY_ABBREV", "AIC", "BIC"),
      variable.name = "METRIC",
      value.name = "VALUE"
    )

    summ_long.dt[, COL_ID := paste(METRIC, SEX, SIDE, sep = "_")]

    summ_wide.dt <- dcast(
      summ_long.dt,
      ROI_LABEL + ADJ_LABEL ~ COL_ID,
      value.var = "VALUE"
    )

    # Order columns by sex then side then metric
    col_order <- character(0)
    for (sex in c("Female", "Male")) {
      for (side in c("L", "R")) {
        for (metric in c("FAMILY_ABBREV", "AIC", "BIC")) {
          col <- paste(metric, sex, side, sep = "_")
          if (col %in% names(summ_wide.dt)) col_order <- c(col_order, col)
        }
      }
    }
    setcolorder(summ_wide.dt, c("ROI_LABEL", "ADJ_LABEL", col_order))

    # Create gt table
    summary_gt <- summ_wide.dt |>
      gt(groupname_col = "ROI_LABEL", rowname_col = "ADJ_LABEL") |>
      tab_stubhead(label = "Adj.") |>
      tab_header(
        title = "GAMLSS Model Summary: Final Selected Models",
        subtitle = "Family and fit statistics"
      )

    # Add sex spanners (level 2)
    for (sex in c("Female", "Male")) {
      sex_cols <- grep(paste0("_", sex, "_"), col_order, value = TRUE)
      if (length(sex_cols) > 0) {
        summary_gt <- summary_gt |>
          tab_spanner(label = sex, columns = all_of(sex_cols), level = 2)
      }
    }

    # Add side spanners (level 1)
    for (sex in c("Female", "Male")) {
      for (side in c("L", "R")) {
        side_cols <- grep(paste0("_", sex, "_", side, "$"), col_order, value = TRUE)
        if (length(side_cols) > 0) {
          summary_gt <- summary_gt |>
            tab_spanner(
              label = ifelse(side == "L", "L", "R"),
              columns = all_of(side_cols),
              level = 1,
              id = paste0(sex, "_", side)
            )
        }
      }
    }

    # Column labels
    col_labels_list <- list()
    for (col in col_order) {
      if (grepl("^FAMILY_ABBREV_", col)) col_labels_list[[col]] <- "Dist."
      else if (grepl("^AIC_", col)) col_labels_list[[col]] <- "AIC"
      else if (grepl("^BIC_", col)) col_labels_list[[col]] <- "BIC"
    }
    summary_gt <- summary_gt |> cols_label(.list = col_labels_list)

    # Source notes
    summ_notes <- c(
      get_abbrev_legend("family"),
      get_abbrev_legend("adjustment"),
      "L = Left hemisphere; R = Right hemisphere"
    )
    for (note in summ_notes) {
      summary_gt <- summary_gt |> tab_source_note(source_note = note)
    }

    save_table(summary_gt, "gamlss_summary_final_unilateral", summary_table_dir, formats = c("html", "tex"))
    log_info("  - Model summary table created (unilateral)")
  }

  # --- BILATERAL TABLES (LR only) ---
  summ_bilat.dt <- summary_data.dt[SIDE == "LR"]

  if (nrow(summ_bilat.dt) > 0) {
    summ_bilat.dt[, FAMILY_ABBREV := as.character(fam_info$abbrev[FAMILY])]
    summ_bilat.dt[, ROI_LABEL := as.character(ROIS[ROI])]
    summ_bilat.dt[, ADJ_LABEL := as.character(adj_abbrev[ADJ])]

    summ_bilat_long.dt <- melt(
      summ_bilat.dt,
      id.vars = c("SEX", "ROI_LABEL", "ADJ_LABEL"),
      measure.vars = c("FAMILY_ABBREV", "AIC", "BIC"),
      variable.name = "METRIC",
      value.name = "VALUE"
    )

    summ_bilat_long.dt[, COL_ID := paste(METRIC, SEX, sep = "_")]

    summ_bilat_wide.dt <- dcast(
      summ_bilat_long.dt,
      ROI_LABEL + ADJ_LABEL ~ COL_ID,
      value.var = "VALUE"
    )

    bilat_col_order <- c("FAMILY_ABBREV_Female", "AIC_Female", "BIC_Female",
                         "FAMILY_ABBREV_Male", "AIC_Male", "BIC_Male")
    bilat_col_order <- bilat_col_order[bilat_col_order %in% names(summ_bilat_wide.dt)]
    setcolorder(summ_bilat_wide.dt, c("ROI_LABEL", "ADJ_LABEL", bilat_col_order))

    bilat_summ_gt <- summ_bilat_wide.dt |>
      gt(groupname_col = "ROI_LABEL", rowname_col = "ADJ_LABEL") |>
      tab_stubhead(label = "Adj.") |>
      tab_header(
        title = "GAMLSS Model Summary: Bilateral (Sum)",
        subtitle = "Family and fit statistics"
      )

    for (sex in c("Female", "Male")) {
      sex_cols <- grep(paste0("_", sex, "$"), bilat_col_order, value = TRUE)
      if (length(sex_cols) > 0) {
        bilat_summ_gt <- bilat_summ_gt |>
          tab_spanner(label = sex, columns = all_of(sex_cols))
      }
    }

    bilat_col_labels <- list()
    for (col in bilat_col_order) {
      if (grepl("^FAMILY_ABBREV_", col)) bilat_col_labels[[col]] <- "Dist."
      else if (grepl("^AIC_", col)) bilat_col_labels[[col]] <- "AIC"
      else if (grepl("^BIC_", col)) bilat_col_labels[[col]] <- "BIC"
    }
    bilat_summ_gt <- bilat_summ_gt |> cols_label(.list = bilat_col_labels)

    bilat_summ_notes <- c(
      get_abbrev_legend("family"),
      get_abbrev_legend("adjustment"),
      "Bilateral = sum of left and right hemispheres"
    )
    for (note in bilat_summ_notes) {
      bilat_summ_gt <- bilat_summ_gt |> tab_source_note(source_note = note)
    }

    save_table(bilat_summ_gt, "gamlss_summary_final_bilateral", summary_table_dir, formats = c("html", "tex"))
    log_info("  - Model summary table created (bilateral)")
  }
}

# Supplementary tables: All fitted models by ROI
if (!is.null(gamlss.lst$SUMMARIES$MODEL_COMPARISON)) {
  model_comparison.dt <- gamlss.lst$SUMMARIES$MODEL_COMPARISON

  # Get side info
  side_info <- get_side_info()

  for (roi in c("HC", "LV", "HVR")) {
    roi_models.dt <- copy(model_comparison.dt[ROI == roi])

    if (nrow(roi_models.dt) > 0) {
      # Find selected model for each combination
      roi_models.dt[, SELECTED_FAMILY := NA_character_]
      for (i in seq_len(nrow(roi_models.dt))) {
        row <- roi_models.dt[i]
        final_mod <- gamlss.lst$FINAL[[row$SEX]][[row$ROI]][[row$ADJ]][[row$SIDE]]
        if (!is.null(final_mod)) {
          roi_models.dt[i, SELECTED_FAMILY := final_mod$FAMILY]
        }
      }

      # Mark if this row is the best for its combination
      roi_models.dt[, IS_BEST_AIC := DELTA_AIC == 0]
      roi_models.dt[, IS_BEST_BIC := DELTA_BIC == 0]

      # Apply abbreviations (remove ROI column since it's redundant)
      roi_models.dt[, ADJ_LABEL := adj_abbrev[ADJ]]
      roi_models.dt[, SIDE_LABEL := side_info$abbrev[SIDE]]
      roi_models.dt[, MOD_LABEL := fam_info$abbrev[MOD]]

      # Selected column: show family abbreviation or empty
      roi_models.dt[, SELECTED_LABEL := ifelse(
        MOD == SELECTED_FAMILY,
        fam_info$abbrev[MOD],
        ""
      )]

      # Format AIC/BIC - round to integers for cleaner display
      roi_models.dt[, AIC := round(AIC, 0)]
      roi_models.dt[, BIC := round(BIC, 0)]
      roi_models.dt[, DELTA_AIC := round(DELTA_AIC, 1)]
      roi_models.dt[, DELTA_BIC := round(DELTA_BIC, 1)]

      # Select and order columns for display
      display_cols <- c("ADJ_LABEL", "SIDE_LABEL", "MOD_LABEL", "AIC", "BIC", "DELTA_AIC", "DELTA_BIC", "SELECTED_LABEL")
      supp_display.dt <- roi_models.dt[, c("SEX", display_cols, "IS_BEST_AIC", "IS_BEST_BIC"), with = FALSE]

      # Create table
      supp_gt <- supp_display.dt |>
        gt(groupname_col = "SEX") |>
        tab_header(
          title = sprintf("%s: Family Comparison", ROIS[roi]),
          subtitle = "AIC/BIC comparison across fitted families"
        ) |>
        cols_label(
          ADJ_LABEL = "Adj.",
          SIDE_LABEL = "Side",
          MOD_LABEL = "Fam.",
          AIC = "AIC",
          BIC = "BIC",
          DELTA_AIC = "\u0394AIC",
          DELTA_BIC = "\u0394BIC",
          SELECTED_LABEL = "Selected"
        ) |>
        cols_hide(columns = c("IS_BEST_AIC", "IS_BEST_BIC"))

      # Add model/fit statistics spanner
      supp_gt <- supp_gt |>
        tab_spanner(
          label = "Fit Statistics",
          columns = c("AIC", "BIC", "DELTA_AIC", "DELTA_BIC")
        )

      # Bold the best AIC/BIC values
      supp_gt <- supp_gt |>
        tab_style(
          style = cell_text(weight = "bold"),
          locations = cells_body(columns = "AIC", rows = IS_BEST_AIC == TRUE)
        ) |>
        tab_style(
          style = cell_text(weight = "bold"),
          locations = cells_body(columns = "BIC", rows = IS_BEST_BIC == TRUE)
        )

      # Source notes
      comp_notes <- c(
        get_abbrev_legend("family"),
        get_abbrev_legend("adjustment"),
        get_abbrev_legend("side"),
        "Bold values indicate best fit (lowest AIC/BIC); Selected = final model chosen"
      )
      for (note in comp_notes) {
        supp_gt <- supp_gt |> tab_source_note(source_note = note)
      }

      save_table(
        supp_gt,
        sprintf("gamlss_comparison_%s", tolower(roi)),
        summary_table_dir,
        formats = c("html", "tex")
      )
    }
  }

  log_info("  - Supplementary model comparison tables created (all models)")
}

# --- Longitudinal Stability Tables (Reorganized) ---
log_info("Creating longitudinal stability tables")

log_info("  - Extracting temporal validation results from %d combinations", nrow(sort.dt))
stability_data.dt <- sort.dt[, {
  final.lst <- gamlss.lst$FINAL[[SEX]][[ROI]][[ADJ]][[SIDE]]
  temp.lst <- final.lst$TEMP_VAL
  if (is.null(temp.lst)) NULL else list(
    FAMILY = final.lst$FAMILY,  # Get selected family
    CORR = round(temp.lst$CORR, 3),
    MEAN_CHANGE = round(temp.lst$MEAN_CHANGE, 3),
    SD_CHANGE = round(temp.lst$SD_CHANGE, 3),
    CROSS_PROP = round(temp.lst$CROSS_PROP, 3),
    REFIT_WARNING = temp.lst$REFIT_WARNING,
    AGE_EXTRAP = temp.lst$AGE_EXTRAP,
    ICC_EXTRAP = temp.lst$ICC_EXTRAP
  )
}, by = .(SEX, ROI, ADJ, SIDE)]

if (nrow(stability_data.dt) > 0) {
  # Ensure all columns are atomic (not list) before processing
  for (col in c("SEX", "ROI", "ADJ", "SIDE", "FAMILY", "CORR", "MEAN_CHANGE", "SD_CHANGE",
                "CROSS_PROP", "REFIT_WARNING", "AGE_EXTRAP", "ICC_EXTRAP")) {
    if (col %in% names(stability_data.dt) && is.list(stability_data.dt[[col]])) {
      stability_data.dt[, (col) := unlist(get(col))]
    }
  }

  # --- UNILATERAL TABLES (L/R only) ---
  stab_unilat.dt <- stability_data.dt[SIDE != "LR"]

  if (nrow(stab_unilat.dt) > 0) {
    # Create combined adjustment label with family: "Unadj. (Norm)"
    stab_unilat.dt[, ADJ_FAMILY_LABEL := as.character(sprintf("%s (%s)", adj_abbrev[ADJ], fam_info$abbrev[FAMILY]))]
    stab_unilat.dt[, ROI_LABEL := as.character(ROIS[ROI])]

    # Restructure
    stab_long.dt <- melt(
      stab_unilat.dt,
      id.vars = c("SEX", "ROI_LABEL", "ADJ_FAMILY_LABEL", "SIDE", "REFIT_WARNING", "AGE_EXTRAP", "ICC_EXTRAP"),
      measure.vars = c("CORR", "MEAN_CHANGE", "SD_CHANGE", "CROSS_PROP"),
      variable.name = "METRIC",
      value.name = "VALUE"
    )

    stab_long.dt[, COL_ID := paste(METRIC, SEX, SIDE, sep = "_")]

    # Keep warning columns for footnote generation
    warning_data.dt <- unique(stab_unilat.dt[, .(SEX, ROI_LABEL, ADJ_FAMILY_LABEL, SIDE, REFIT_WARNING, AGE_EXTRAP, ICC_EXTRAP)])

    stab_wide.dt <- dcast(
      stab_long.dt,
      ROI_LABEL + ADJ_FAMILY_LABEL ~ COL_ID,
      value.var = "VALUE"
    )

    # Order columns
    col_order <- character(0)
    for (sex in c("Female", "Male")) {
      for (side in c("L", "R")) {
        for (metric in c("CORR", "MEAN_CHANGE", "SD_CHANGE", "CROSS_PROP")) {
          col <- paste(metric, sex, side, sep = "_")
          if (col %in% names(stab_wide.dt)) col_order <- c(col_order, col)
        }
      }
    }
    setcolorder(stab_wide.dt, c("ROI_LABEL", "ADJ_FAMILY_LABEL", col_order))

    # Create gt table
    stability_gt <- stab_wide.dt |>
      gt(groupname_col = "ROI_LABEL", rowname_col = "ADJ_FAMILY_LABEL") |>
      tab_stubhead(label = "Adj. (Dist.)") |>
      tab_header(
        title = "Longitudinal Stability: Temporal Validation",
        subtitle = "Correlation and centile crossing analysis"
      )

    # Add sex spanners (level 2)
    for (sex in c("Female", "Male")) {
      sex_cols <- grep(paste0("_", sex, "_"), col_order, value = TRUE)
      if (length(sex_cols) > 0) {
        stability_gt <- stability_gt |>
          tab_spanner(label = sex, columns = all_of(sex_cols), level = 2)
      }
    }

    # Add side spanners (level 1)
    for (sex in c("Female", "Male")) {
      for (side in c("L", "R")) {
        side_cols <- grep(paste0("_", sex, "_", side, "$"), col_order, value = TRUE)
        if (length(side_cols) > 0) {
          stability_gt <- stability_gt |>
            tab_spanner(
              label = ifelse(side == "L", "L", "R"),
              columns = all_of(side_cols),
              level = 1,
              id = paste0(sex, "_", side)
            )
        }
      }
    }

    # Column labels
    col_labels_list <- list()
    for (col in col_order) {
      if (grepl("^CORR_", col)) col_labels_list[[col]] <- "r"
      else if (grepl("^MEAN_CHANGE_", col)) col_labels_list[[col]] <- "\u0394\u03bc"
      else if (grepl("^SD_CHANGE_", col)) col_labels_list[[col]] <- "\u0394\u03c3"
      else if (grepl("^CROSS_PROP_", col)) col_labels_list[[col]] <- "X%"
    }
    stability_gt <- stability_gt |> cols_label(.list = col_labels_list)

    # Create footnotes for warnings
    n_refit <- sum(warning_data.dt$REFIT_WARNING, na.rm = TRUE)
    n_age_extrap <- sum(warning_data.dt$AGE_EXTRAP, na.rm = TRUE)
    n_icc_extrap <- sum(warning_data.dt$ICC_EXTRAP, na.rm = TRUE)
    total_combinations <- nrow(warning_data.dt)

    stab_notes <- c(
      "r = Pearson correlation; \u0394\u03bc = Mean centile change; \u0394\u03c3 = SD of centile change; X% = Proportion crossing clinical thresholds",
      get_abbrev_legend("family"),
      get_abbrev_legend("adjustment"),
      "L = Left hemisphere; R = Right hemisphere",
      sprintf("Warnings: %d/%d models required refitting, %d/%d had age extrapolation, %d/%d had TIV extrapolation",
        n_refit, total_combinations, n_age_extrap, total_combinations, n_icc_extrap, total_combinations)
    )
    for (note in stab_notes) {
      stability_gt <- stability_gt |> tab_source_note(source_note = note)
    }

    save_table(stability_gt, "gamlss_stability_unilateral", summary_table_dir, formats = c("html", "tex"))
    log_info("  - Longitudinal stability table created (unilateral)")
  }

  # --- BILATERAL TABLES (LR only) ---
  stab_bilat.dt <- stability_data.dt[SIDE == "LR"]

  if (nrow(stab_bilat.dt) > 0) {
    stab_bilat.dt[, ADJ_FAMILY_LABEL := as.character(sprintf("%s (%s)", adj_abbrev[ADJ], fam_info$abbrev[FAMILY]))]
    stab_bilat.dt[, ROI_LABEL := as.character(ROIS[ROI])]

    stab_bilat_long.dt <- melt(
      stab_bilat.dt,
      id.vars = c("SEX", "ROI_LABEL", "ADJ_FAMILY_LABEL", "REFIT_WARNING", "AGE_EXTRAP", "ICC_EXTRAP"),
      measure.vars = c("CORR", "MEAN_CHANGE", "SD_CHANGE", "CROSS_PROP"),
      variable.name = "METRIC",
      value.name = "VALUE"
    )

    stab_bilat_long.dt[, COL_ID := paste(METRIC, SEX, sep = "_")]

    stab_bilat_wide.dt <- dcast(
      stab_bilat_long.dt,
      ROI_LABEL + ADJ_FAMILY_LABEL ~ COL_ID,
      value.var = "VALUE"
    )

    bilat_col_order <- c("CORR_Female", "MEAN_CHANGE_Female", "SD_CHANGE_Female", "CROSS_PROP_Female",
                         "CORR_Male", "MEAN_CHANGE_Male", "SD_CHANGE_Male", "CROSS_PROP_Male")
    bilat_col_order <- bilat_col_order[bilat_col_order %in% names(stab_bilat_wide.dt)]
    setcolorder(stab_bilat_wide.dt, c("ROI_LABEL", "ADJ_FAMILY_LABEL", bilat_col_order))

    bilat_stab_gt <- stab_bilat_wide.dt |>
      gt(groupname_col = "ROI_LABEL", rowname_col = "ADJ_FAMILY_LABEL") |>
      tab_stubhead(label = "Adj. (Dist.)") |>
      tab_header(
        title = "Longitudinal Stability: Bilateral (Sum)",
        subtitle = "Correlation and centile crossing analysis"
      )

    for (sex in c("Female", "Male")) {
      sex_cols <- grep(paste0("_", sex, "$"), bilat_col_order, value = TRUE)
      if (length(sex_cols) > 0) {
        bilat_stab_gt <- bilat_stab_gt |>
          tab_spanner(label = sex, columns = all_of(sex_cols))
      }
    }

    bilat_col_labels <- list()
    for (col in bilat_col_order) {
      if (grepl("^CORR_", col)) bilat_col_labels[[col]] <- "r"
      else if (grepl("^MEAN_CHANGE_", col)) bilat_col_labels[[col]] <- "\u0394\u03bc"
      else if (grepl("^SD_CHANGE_", col)) bilat_col_labels[[col]] <- "\u0394\u03c3"
      else if (grepl("^CROSS_PROP_", col)) bilat_col_labels[[col]] <- "X%"
    }
    bilat_stab_gt <- bilat_stab_gt |> cols_label(.list = bilat_col_labels)

    bilat_stab_notes <- c(
      "r = Pearson correlation; \u0394\u03bc = Mean centile change; \u0394\u03c3 = SD of centile change; X% = Proportion crossing clinical thresholds",
      get_abbrev_legend("family"),
      get_abbrev_legend("adjustment"),
      "Bilateral = sum of left and right hemispheres"
    )
    for (note in bilat_stab_notes) {
      bilat_stab_gt <- bilat_stab_gt |> tab_source_note(source_note = note)
    }

    save_table(bilat_stab_gt, "gamlss_stability_bilateral", summary_table_dir, formats = c("html", "tex"))
    log_info("  - Longitudinal stability table created (bilateral)")
  }
}

log_info("Diagnostic summary tables saved to: %s", summary_table_dir)

log_script_end("08_normative_tables.R", success = TRUE)

# ----- Generate Diagnostic Plots -----
if (REDO_PLOTS) {
  log_section("Generating diagnostic plots")

  # Get directories from consolidated output paths
  fig_dir <- get_output_path("figures")
  centiles_dir <- file.path(fig_dir, "gamlss_centiles")
  validation_dir <- file.path(fig_dir, "gamlss_validation")
  stability_dir <- file.path(fig_dir, "gamlss_stability")
  summary_table_dir <- file.path(get_output_path("tables"), "gamlss_summary")

  ensure_directory(centiles_dir)
  ensure_directory(validation_dir)
  ensure_directory(stability_dir)

  # -------------------------------------------------------------------
  # Faceted Centile Plots by Sex
  # -------------------------------------------------------------------
  log_info("Creating faceted centile plots (by ROI and adjustment)")

  for (sex in sexes.v) {
    for (roi in c("HC", "LV", "HVR")) {
      log_info("  - %s: %s", sex, ROIS[roi])

      # Prepare data across all adjustments
      norm_list <- list()
      obs_list <- list()

      for (adj in c("NON", "PRP", "STX", "RES")) {
        if (is.null(norm_tables.lst[[sex]][[roi]][[adj]][["LR"]])) {
          next
        }

        # Normative table
        norm.dt <- copy(norm_tables.lst[[sex]][[roi]][[adj]][["LR"]])
        norm.dt[, ADJ := ADJS[[adj]]]
        norm_list[[adj]] <- norm.dt

        # Observed data
        obs.dt <- gamlss.lst$DATA$CRS[
          .(sex, roi, adj, "LR")
        ][
          !is.na(VAL), .(AGE, VAL)
        ]
        obs.dt[, ADJ := ADJS[[adj]]]
        obs_list[[adj]] <- obs.dt
      }

      if (length(norm_list) == 0) next

      # Create plot using utility function
      p <- plot_gamlss_faceted_by_adjustment(
        norm_tables_list = norm_list,
        obs_data_list = obs_list,
        roi_label = ROIS[roi],
        sex = sex,
        sex_color = get_palette("sex")[sex],
        y_label = if (roi == "HVR") "Ratio" else "Volume (cc)"
      )

      # Save plot
      save_plot(
        p,
        file.path(
          centiles_dir,
          sprintf("centiles_faceted_%s_%s.png", roi, tolower(sex))
        ),
        width = 12,
        height = 10
      )
    }
  }

  # -------------------------------------------------------------------
  # Sex Comparison Plots - Faceted by Adjustment
  # -------------------------------------------------------------------
  log_info("Creating sex-comparison faceted plots")

  for (roi in c("HC", "LV", "HVR")) {
    log_info("  - Sex comparison: %s", ROIS[roi])

    # Prepare data across sexes and adjustments
    norm_list <- list()
    obs_list <- list()

    for (adj in c("NON", "PRP", "STX", "RES")) {
      for (sex in sexes.v) {
        key <- paste(adj, sex, sep = "_")

        if (is.null(norm_tables.lst[[sex]][[roi]][[adj]][["LR"]])) {
          next
        }

        # Normative table
        norm.dt <- copy(norm_tables.lst[[sex]][[roi]][[adj]][["LR"]])
        norm.dt[, let(ADJ = ADJS[[adj]], SEX = sex)]
        norm_list[[key]] <- norm.dt

        # Observed data
        obs.dt <- gamlss.lst$DATA$CRS[
          .(sex, roi, adj, "LR")
        ][
          !is.na(VAL), .(AGE, VAL, SEX)
        ]
        obs.dt[, ADJ := ADJS[[adj]]]
        obs_list[[key]] <- obs.dt
      }
    }

    if (length(norm_list) == 0) next

    # Create plot using utility function
    p <- plot_gamlss_sexcomp_faceted(
      norm_tables_list = norm_list,
      obs_data_list = obs_list,
      roi_label = ROIS[roi],
      sex_colors = get_palette("sex"),
      y_label = if (roi == "HVR") "Ratio" else "Volume (cc)"
    )

    # Save plot
    save_plot(
      p,
      file.path(
        centiles_dir,
        sprintf("centiles_sexcomp_faceted_%s.png", roi)
      ),
      width = 12,
      height = 10
    )
  }

  # -------------------------------------------------------------------
  # Validation & Stability Plots - Summarized
  # -------------------------------------------------------------------
  log_info("Creating validation/stability summary plots")

  # Validation plots
  for (roi in c("HC", "LV", "HVR")) {
    for (sex in sexes.v) {
      # Collect validation data across adjustments
      valid_list <- list()

      for (adj in c("NON", "PRP", "STX", "RES")) {
        test.lst <- gamlss.lst$TEST[[sex]][[roi]][[adj]][["LR"]]
        if (is.null(test.lst)) next

        valid.dt <- copy(test.lst$CENT_VALID)
        valid.dt[, ADJ := ADJS[[adj]]]
        valid_list[[adj]] <- valid.dt
      }

      if (length(valid_list) == 0) next

      # Create plot using utility function
      p <- plot_validation_faceted(
        valid_data_list = valid_list,
        roi_label = ROIS[roi],
        sex = sex
      )

      # Save plot
      save_plot(
        p,
        file.path(
          validation_dir,
          sprintf("validation_faceted_%s_%s.png", roi, tolower(sex))
        ),
        width = 10,
        height = 10
      )
    }
  }

  # Stability plots (wrap in tryCatch since data structure may vary)
  tryCatch({
    for (roi in c("HC", "LV", "HVR")) {
      for (sex in sexes.v) {
        # Collect stability data across adjustments
        stab_list <- list()

        for (adj in c("NON", "PRP", "STX", "RES")) {
          temp.lst <- gamlss.lst$FINAL[[sex]][[roi]][[adj]][[
            "LR"
          ]]$TEMP_VAL
          if (is.null(temp.lst)) next

          stab.dt <- copy(temp.lst$CENTILES)
          stab.dt[, ADJ := ADJS[[adj]]]
          stab_list[[adj]] <- stab.dt
        }

        if (length(stab_list) == 0) next

        # Create plot using utility function
        p <- plot_stability_faceted(
          stability_data_list = stab_list,
          roi_label = ROIS[roi],
          sex_colors = SEX_COLORS
        )

        # Save plot
        save_plot(
          p,
          file.path(
            stability_dir,
            sprintf("stability_faceted_%s_%s.png", roi, tolower(sex))
          ),
          width = 10,
          height = 10
        )
      }
    }
  }, error = function(e) {
    log_warn("Stability plots could not be generated: %s", conditionMessage(e))
  })

  log_info("GAMLSS diagnostic plots saved:")
  log_info("  - Centile plots: %s (%d plots)", centiles_dir, 3 * 2 + 3)
  log_info("  - Validation plots: %s (%d plots)", validation_dir, 3 * 2)
  log_info("  - Stability plots: %s (%d plots)", stability_dir, 3 * 2)

  # -------------------------------------------------------------------
  # Z-Score Standardized Plots for Cross-Method Comparison
  # -------------------------------------------------------------------
  # These plots show HC across all adjustment methods on standardized scale
  # This solves the problem of proportions having different scale than raw volumes
  log_info("Creating z-score standardized comparison plots")

  for (sex in sexes.v) {
    log_info("  - Z-score HC trajectories: %s", sex)

    # Collect observed data across all adjustment methods
    obs_list <- list()
    for (adj in c("NON", "PRP", "STX", "RES")) {
      obs.dt <- gamlss.lst$DATA$CRS[
        .(sex, "HC", adj, "LR")
      ][
        !is.na(VAL), .(AGE, VAL)
      ]
      obs.dt[, ADJ := ADJS[[adj]]]
      obs.dt[, SEX := sex]  # Add SEX column for the plotting function
      obs_list[[adj]] <- obs.dt
    }

    if (length(obs_list) > 0) {
      p <- plot_zscore_trajectories(
        obs_data_list = obs_list,
        roi_label = "Hippocampal Volume",
        sex_colors = get_palette("sex"),
        show_points = FALSE
      )

      save_plot(
        p,
        file.path(centiles_dir, sprintf("zscore_hc_%s.png", tolower(sex))),
        width = 10,
        height = 8
      )
    }
  }

  # -------------------------------------------------------------------
  # Simple HVR Centile Plot (Self-Normalizing)
  # -------------------------------------------------------------------
  # HVR = HC/(HC+LV) is inherently self-normalizing - the ratio naturally
  # cancels out head-size effects, so no explicit adjustment is needed.
  # Note: We use PRP (proportions) as the data key since HVR has no adjustments.
  log_info("Creating HVR centile plot (self-normalizing, both sexes)")

  norm_hvr_list <- list()
  obs_hvr_list <- list()

  # Try PRP first (proportions), then RES as fallback
  hvr_adj <- "PRP"

  for (sex in sexes.v) {
    if (!is.null(norm_tables.lst[[sex]][["HVR"]][[hvr_adj]][["LR"]])) {
      norm.dt <- copy(norm_tables.lst[[sex]][["HVR"]][[hvr_adj]][["LR"]])
      norm.dt[, SEX := sex]
      norm_hvr_list[[sex]] <- norm.dt

      obs.dt <- gamlss.lst$DATA$CRS[
        .(sex, "HVR", hvr_adj, "LR")
      ][
        !is.na(VAL), .(AGE, VAL)
      ]
      obs.dt[, SEX := sex]  # Add SEX column explicitly
      obs_hvr_list[[sex]] <- obs.dt
    }
  }

  if (length(norm_hvr_list) > 0) {
    p <- plot_hvr_simple(
      norm_tables_list = norm_hvr_list,
      obs_data_list = obs_hvr_list,
      sex_colors = get_palette("sex")
    )

    save_plot(
      p,
      file.path(centiles_dir, "centiles_hvr_self_normalizing.png"),
      width = 10,
      height = 8
    )
  }

  log_info("  - Z-score plots: 2 (Female, Male)")
  log_info("  - HVR self-normalizing plot: 1")
}

log_info("GAMLSS normative tables and z-scores generated successfully")

# =============================================================================
# TRANSFER MODELS (for cross-study application)
# =============================================================================
# The site-controlled models include random(SITE) for within-UKB accuracy.
# For external cohorts (e.g., ADNI), we need models without SITE since:
# 1. gamlss.add's random() returns NA for levels not in training data
# 2. Using a UKB site for external data is methodologically inappropriate
#
# References:
# - PCNtoolkit: https://pcntoolkit.readthedocs.io/en/latest/pages/tutorials/06_transfer_extend.html
# - Kia et al. (2022): Accommodating site variation in neuroimaging data
#
# Strategy: Refit FINAL models without SITE for transfer/external use
# =============================================================================

# Check if transfer models already exist (have TRANSFER = TRUE flag)
transfer_models_exist <- FALSE
if (check_files_exist(gamlss.path, stop_on_missing = FALSE)) {
  existing.lst <- read_rds_safe(gamlss.path, "existing GAMLSS file")
  # Check if any FINAL model has TRANSFER flag
  if (!is.null(existing.lst$FINAL)) {
    first_model <- existing.lst$FINAL[[1]][[1]][[1]][[1]]
    transfer_models_exist <- isTRUE(first_model$TRANSFER)
  }
  rm(existing.lst)
}

# Skip if transfer models already exist and not force regenerating
FORCE_TRANSFER <- get_script_setting("gamlss", "force_transfer", default = FALSE)
if (transfer_models_exist && !FORCE_TRANSFER) {
  log_info("Transfer models already exist in %s", basename(gamlss.path))
  log_info("Set force_regenerate: transfer: yes to refit")
} else {
  log_section("Generating Transfer Models (no SITE random effect)")

# Save site-controlled models to new path
gamlss_site.path <- sub("gamlss\\.rds$", "gamlss_site.rds", gamlss.path)
log_info("Saving site-controlled models to: %s", gamlss_site.path)
write_rds_safe(gamlss.lst, gamlss_site.path, "GAMLSS site-controlled models")

# Create transfer model structure with train/test validation
log_info("Refitting models without SITE for external application (with validation)")
transfer.lst <- list(
  DATA = gamlss.lst$DATA,  # Include DATA for table/plot generation
  TRAIN = list(),
  TEST = list(),
  FINAL = list(),
  SUMMARIES = NULL
)

# Get the full data for refitting (includes SPLIT column from site-controlled models)
gamlss.dt <- gamlss.lst$DATA$CRS

# Verify SPLIT column exists
if (!"SPLIT" %in% names(gamlss.dt)) {
  log_error("SPLIT column not found in data - cannot perform transfer model validation")
  stop("Missing SPLIT column for transfer model validation", call. = FALSE)
}

split_n.dt <- gamlss.dt[!duplicated(EID), .N, keyby = .(SPLIT, SEX)]
log_info(
  "Transfer model validation split: Train=%d, Test=%d (same as site-controlled)",
  split_n.dt[SPLIT == "train", sum(N)],
  split_n.dt[SPLIT == "test", sum(N)]
)

# ---- Phase 1: Training transfer models ----
log_info("Phase 1/3: Fitting transfer models on training data")
pb <- progress_bar$new(
  format = "Transfer Training | :what [:bar] :current/:total\n",
  total = nrow(sort.dt), clear = FALSE, width = 100, show_after = 0
)
pb$tick(0)
setkeyv(gamlss.dt, group_cols.v)

for (i in seq_len(nrow(sort.dt))) {
  row.dt <- sort.dt[i]
  sex <- row.dt$SEX
  roi <- row.dt$ROI
  adj <- row.dt$ADJ
  side <- row.dt$SIDE

  # --- Data preparation (without SITE, with SPLIT) ---
  gamlss.subdt <- gamlss.dt[row.dt][
    VAL > 0, .(EID, AGE, EDUC_num, ICC, VAL, SPLIT)
  ] |> na.omit()

  gamlss_train.subdt <- gamlss.subdt[SPLIT == "train", -"SPLIT"]

  if (!nrow(gamlss_train.subdt)) {
    pb$tick()
    next
  }

  # Use same model family as site-controlled version
  site_mod.lst <- gamlss.lst$FINAL[[sex]][[roi]][[adj]][[side]]
  if (is.null(site_mod.lst)) {
    pb$tick()
    next
  }
  best_mod <- site_mod.lst$FAMILY

  pb$tick(tokens = list(what = sprintf(
    "%s — %s (%s; %s): %s",
    sex, roi, side, adj, fam_names.v[[best_mod]]
  )))

  # Formula WITHOUT random(SITE)
  formula_base <- if (roi != "HVR" && adj == "NON") {
    VAL ~ cs(AGE) + ICC + EDUC_num
  } else {
    VAL ~ cs(AGE) + EDUC_num
  }

  # Initialize nested list structure for TRAIN
  if (is.null(transfer.lst$TRAIN[[sex]])) transfer.lst$TRAIN[[sex]] <- list()
  if (is.null(transfer.lst$TRAIN[[sex]][[roi]])) transfer.lst$TRAIN[[sex]][[roi]] <- list()
  if (is.null(transfer.lst$TRAIN[[sex]][[roi]][[adj]])) transfer.lst$TRAIN[[sex]][[roi]][[adj]] <- list()

  # Fit model without SITE on training data
  tryCatch({
    if (best_mod == "NO") {
      output.fit <- gamlss(
        formula = formula_base,
        sigma.formula = ~ cs(AGE),
        data = gamlss_train.subdt,
        family = NO()
      )
    } else if (best_mod == "L_NO") {
      eps <- 1e-6
      gamlss_train.subdt[
        , VAL_logit := VAL |> pmin(1 - eps) |> pmax(eps) |> qlogis()
      ]
      formula_logit <- if (roi != "HVR" && adj == "NON") {
        VAL_logit ~ cs(AGE) + ICC + EDUC_num
      } else {
        VAL_logit ~ cs(AGE) + EDUC_num
      }
      output.fit <- gamlss(
        formula = formula_logit,
        sigma.formula = ~ cs(AGE),
        data = gamlss_train.subdt,
        family = NO()
      )
    } else if (best_mod == "BCCG") {
      output.fit <- gamlss(
        formula = formula_base,
        sigma.formula = ~ cs(AGE),
        nu.formula = ~ cs(AGE),
        data = gamlss_train.subdt,
        family = BCCG()
      )
    } else if (best_mod == "BE") {
      output.fit <- gamlss(
        formula = formula_base,
        sigma.formula = ~ cs(AGE),
        nu.formula = ~ cs(AGE),
        data = gamlss_train.subdt,
        family = BE()
      )
    }

    # Store training model
    transfer.lst$TRAIN[[sex]][[roi]][[adj]][[side]] <- list(
      FAMILY = best_mod,
      DATA = gamlss_train.subdt,
      FIT = output.fit,
      AIC = output.fit$aic,
      BIC = BIC(output.fit)
    )

    rm(output.fit)
  }, error = function(e) {
    log_warn(
      "%s — %s (%s; %s): Transfer training model fit failed: %s",
      sex, roi, side, adj, e$message
    )
  })
}

# ---- Phase 2: Validating transfer models on test data ----
log_info("Phase 2/3: Validating transfer models on hold-out test data")
cent_thresholds.v <- c(1, 5, 10, 25, 50, 75, 90, 95, 99)

pb <- progress_bar$new(
  format = "Transfer Validation | :what [:bar] :current/:total\n",
  total = nrow(sort.dt), clear = FALSE, width = 100, show_after = 0
)
pb$tick(0)

for (i in seq_len(nrow(sort.dt))) {
  row.dt <- sort.dt[i]
  sex <- row.dt$SEX
  roi <- row.dt$ROI
  adj <- row.dt$ADJ
  side <- row.dt$SIDE

  pb$tick(tokens = list(what = sprintf("%s — %s (%s; %s)", sex, roi, side, adj)))

  # Check if training model exists
  train_mod.lst <- transfer.lst$TRAIN[[sex]][[roi]][[adj]][[side]]
  if (is.null(train_mod.lst)) next

  best_mod <- train_mod.lst$FAMILY

  # Get test data
  gamlss.subdt <- gamlss.dt[row.dt][
    VAL > 0, .(EID, AGE, EDUC_num, ICC, VAL, SPLIT)
  ] |> na.omit()

  gamlss_test.subdt <- gamlss.subdt[SPLIT == "test", -"SPLIT"]

  if (!nrow(gamlss_test.subdt)) next

  # GAMLSS requires training data in environment for prediction
  # Assign it to the variable name the model expects
  gamlss_train.subdt <- train_mod.lst$DATA

  # Initialize nested list structure for TEST
  if (is.null(transfer.lst$TEST[[sex]])) transfer.lst$TEST[[sex]] <- list()
  if (is.null(transfer.lst$TEST[[sex]][[roi]])) transfer.lst$TEST[[sex]][[roi]] <- list()
  if (is.null(transfer.lst$TEST[[sex]][[roi]][[adj]])) transfer.lst$TEST[[sex]][[roi]][[adj]] <- list()

  # --- Predicted centiles ---
  # Prepare prediction data with exact columns from training
  # Transfer models use: AGE, EDUC_num, (ICC for NON non-HVR)
  pred_cols <- if (roi != "HVR" && adj == "NON") {
    c("AGE", "EDUC_num", "ICC")
  } else {
    c("AGE", "EDUC_num")
  }
  pred_data <- gamlss_test.subdt[, ..pred_cols]

  refit_warning <- FALSE
  tryCatch({
    if (best_mod == "NO") {
      gamlss_test.subdt[
        , CENT_pred := pNO(
          VAL,
          mu = predict(train_mod.lst$FIT, newdata = pred_data, type = "response"),
          sigma = predict(train_mod.lst$FIT, "sigma", newdata = pred_data, type = "response")
        ) * 100
      ] |> invisible()
    } else if (best_mod == "L_NO") {
      eps <- 1e-6
      gamlss_test.subdt[
        , VAL_logit := VAL |> pmin(1 - eps) |> pmax(eps) |> qlogis()
      ]
      gamlss_test.subdt[
        , CENT_pred := pNO(
          VAL_logit,
          mu = predict(train_mod.lst$FIT, newdata = pred_data, type = "response"),
          sigma = predict(train_mod.lst$FIT, "sigma", newdata = pred_data, type = "response")
        ) * 100
      ] |> invisible()
    } else if (best_mod == "BE") {
      withCallingHandlers(
        {
          gamlss_test.subdt[
            , CENT_pred := pBE(
              VAL,
              mu = predict(train_mod.lst$FIT, newdata = pred_data, type = "response"),
              sigma = predict(train_mod.lst$FIT, "sigma", newdata = pred_data, type = "response")
            ) * 100
          ] |> invisible()
        },
        warning = function(w) {
          if (grepl("discrepancy.*re-fit", w$message)) refit_warning <<- TRUE
          invokeRestart("muffleWarning")
        }
      )
    } else if (best_mod == "BCCG") {
      withCallingHandlers(
        {
          gamlss_test.subdt[
            , CENT_pred := pBCCG(
              VAL,
              mu = predict(train_mod.lst$FIT, newdata = pred_data, type = "response"),
              sigma = predict(train_mod.lst$FIT, "sigma", newdata = pred_data, type = "response"),
              nu = predict(train_mod.lst$FIT, "nu", newdata = pred_data, type = "response")
            ) * 100
          ] |> invisible()
        },
        warning = function(w) {
          if (grepl("discrepancy.*re-fit", w$message)) refit_warning <<- TRUE
          invokeRestart("muffleWarning")
        }
      )
    }

    # Calculate calibration metrics
    cent_calibration.v <- sapply(
      cent_thresholds.v,
      \(t) mean(gamlss_test.subdt$CENT_pred < t, na.rm = TRUE) * 100
    )
    mae <- mean(abs(cent_calibration.v - cent_thresholds.v))

    calibration.dt <- data.table(
      CENT_exp = cent_thresholds.v,
      CENT_obs = round(cent_calibration.v, 2),
      DIFF = round(cent_calibration.v - cent_thresholds.v, 2)
    )

    transfer.lst$TEST[[sex]][[roi]][[adj]][[side]] <- list(
      MOD = best_mod,
      FIT = train_mod.lst$FIT,
      CENT_VALID = calibration.dt,
      MAE = mae,
      INTERPRETATION = fcase(
        mae < 2, "Excellent",
        mae < 5, "Good",
        default = "Poor"
      ),
      REFIT_WARNING = refit_warning
    )
  }, error = function(e) {
    log_warn(
      "%s — %s (%s; %s): Transfer validation failed: %s",
      sex, roi, side, adj, e$message
    )
  })
}

# --- Transfer validation diagnostics ---
log_section("Transfer model validation diagnostics")
transfer_test_results.dt <- sort.dt[, {
  res <- transfer.lst$TEST[[SEX]][[ROI]][[ADJ]][[SIDE]]
  if (is.null(res)) NULL else list(
    MOD = res$MOD,
    MAE = res$MAE,
    INTERPRETATION = res$INTERPRETATION,
    REFIT_WARNING = res$REFIT_WARNING
  )
}, by = .(SEX, ROI, ADJ, SIDE)]

n_transfer_models <- nrow(transfer_test_results.dt)
n_transfer_excellent <- sum(transfer_test_results.dt$INTERPRETATION == "Excellent", na.rm = TRUE)
n_transfer_good <- sum(transfer_test_results.dt$INTERPRETATION == "Good", na.rm = TRUE)
n_transfer_poor <- sum(transfer_test_results.dt$INTERPRETATION == "Poor", na.rm = TRUE)

log_info("Transfer models validated: %d", n_transfer_models)
log_info(
  "Transfer model quality: %d Excellent (%.1f%%), %d Good (%.1f%%), %d Poor (%.1f%%)",
  n_transfer_excellent, 100 * n_transfer_excellent / n_transfer_models,
  n_transfer_good, 100 * n_transfer_good / n_transfer_models,
  n_transfer_poor, 100 * n_transfer_poor / n_transfer_models
)

# Store transfer validation summary
transfer.lst$SUMMARIES <- list(
  VALIDATION_RESULTS = transfer_test_results.dt,
  CALIBRATION_DETAILS = sort.dt[, {
    res <- transfer.lst$TEST[[SEX]][[ROI]][[ADJ]][[SIDE]]
    if (is.null(res) || is.null(res$CENT_VALID)) NULL else res$CENT_VALID
  }, by = .(SEX, ROI, ADJ, SIDE)]
)

# ---- Phase 3: Final transfer models on full data ----
log_info("Phase 3/3: Fitting final transfer models on full data")
pb <- progress_bar$new(
  format = "Transfer Final | :what [:bar] :current/:total\n",
  total = nrow(sort.dt), clear = FALSE, width = 100, show_after = 0
)
pb$tick(0)

for (i in seq_len(nrow(sort.dt))) {
  row.dt <- sort.dt[i]
  sex <- row.dt$SEX
  roi <- row.dt$ROI
  adj <- row.dt$ADJ
  side <- row.dt$SIDE

  # --- Data preparation (without SITE, full data) ---
  gamlss_final.subdt <- gamlss.dt[row.dt][
    VAL > 0, .(EID, AGE, EDUC_num, ICC, VAL)
  ] |> na.omit()

  if (!nrow(gamlss_final.subdt)) {
    pb$tick()
    next
  }

  # Use same model family as site-controlled version
  site_mod.lst <- gamlss.lst$FINAL[[sex]][[roi]][[adj]][[side]]
  if (is.null(site_mod.lst)) {
    pb$tick()
    next
  }
  best_mod <- site_mod.lst$FAMILY

  pb$tick(tokens = list(what = sprintf(
    "%s — %s (%s; %s): %s",
    sex, roi, side, adj, fam_names.v[[best_mod]]
  )))

  # Formula WITHOUT random(SITE)
  formula_base <- if (roi != "HVR" && adj == "NON") {
    VAL ~ cs(AGE) + ICC + EDUC_num
  } else {
    VAL ~ cs(AGE) + EDUC_num
  }

  # Initialize nested list structure for FINAL
  if (is.null(transfer.lst$FINAL[[sex]])) transfer.lst$FINAL[[sex]] <- list()
  if (is.null(transfer.lst$FINAL[[sex]][[roi]])) transfer.lst$FINAL[[sex]][[roi]] <- list()
  if (is.null(transfer.lst$FINAL[[sex]][[roi]][[adj]])) transfer.lst$FINAL[[sex]][[roi]][[adj]] <- list()

  # Fit model without SITE on full data
  tryCatch({
    if (best_mod == "NO") {
      output.fit <- gamlss(
        formula = formula_base,
        sigma.formula = ~ cs(AGE),
        data = gamlss_final.subdt,
        family = NO()
      )
    } else if (best_mod == "L_NO") {
      eps <- 1e-6
      gamlss_final.subdt[
        , VAL_logit := VAL |> pmin(1 - eps) |> pmax(eps) |> qlogis()
      ]
      formula_logit <- if (roi != "HVR" && adj == "NON") {
        VAL_logit ~ cs(AGE) + ICC + EDUC_num
      } else {
        VAL_logit ~ cs(AGE) + EDUC_num
      }
      output.fit <- gamlss(
        formula = formula_logit,
        sigma.formula = ~ cs(AGE),
        data = gamlss_final.subdt,
        family = NO()
      )
    } else if (best_mod == "BCCG") {
      output.fit <- gamlss(
        formula = formula_base,
        sigma.formula = ~ cs(AGE),
        nu.formula = ~ cs(AGE),
        data = gamlss_final.subdt,
        family = BCCG()
      )
    } else if (best_mod == "BE") {
      output.fit <- gamlss(
        formula = formula_base,
        sigma.formula = ~ cs(AGE),
        nu.formula = ~ cs(AGE),
        data = gamlss_final.subdt,
        family = BE()
      )
    }

    # Store final transfer model
    transfer.lst$FINAL[[sex]][[roi]][[adj]][[side]] <- list(
      FAMILY = best_mod,
      DATA = gamlss_final.subdt,
      FIT = output.fit,
      AIC = output.fit$aic,
      BIC = BIC(output.fit),
      TRANSFER = TRUE  # Flag indicating this is a transfer model
    )

    rm(output.fit)
  }, error = function(e) {
    log_warn(
      "%s — %s (%s; %s): Transfer final model fit failed: %s",
      sex, roi, side, adj, e$message
    )
  })
}

# Add metadata
transfer.lst$METADATA <- list(
  created = Sys.time(),
  description = "Transfer models for cross-study application (no SITE random effect)",
  source_models = gamlss_site.path,
  reference = "For external cohorts, use these models instead of site-controlled versions"
)

# Save transfer models to original path (backward compatible)
log_info("Saving transfer models to: %s", gamlss.path)
write_rds_safe(transfer.lst, gamlss.path, "GAMLSS transfer models (no SITE)")

log_info("Transfer model generation complete")
log_info("  - Site-controlled models: %s", basename(gamlss_site.path))
log_info("  - Transfer models: %s", basename(gamlss.path))
}  # End of transfer model generation block

# ===========================================================================
# Export Manuscript Figures
# ===========================================================================
log_section("Exporting manuscript figures")

export_data <- list(
  norm_tables = norm_tables.lst,
  brain_data = hc_hvr.lst[[NETWORK]]$CRS$ALL,
  gamlss_calib = gamlss.lst$TEST
)

export_manuscript_figures("gamlss", export_data)
log_info("Manuscript assets exported to: %s", get_output_path("figures"))

log_script_end("08_normative_tables.R", success = TRUE)
