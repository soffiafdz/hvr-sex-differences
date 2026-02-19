#!/usr/bin/env Rscript

# =============================================================================
# Structural Equation Modeling - Brain-Cognition Associations
# =============================================================================
# Tests brain-cognition associations using SEM with bifactor cognitive structure.
#
# Model 1 (HC_COG): HC → Cognition
#   - Latent HC from subfield indicators (head, body, tail)
#   - Controls: Age, ICV, SES (latent)
#
# Model 2 (ROI_COG): HC + LV → Cognition [MAIN ANALYSIS]
#   - Latent HC and LV from subfield indicators
#   - Parallel pathways represent hippocampal-ventricular components
#   - Controls: Age, ICV, SES (latent)
#
# Model 3 (HVR_COG): HVR → Cognition [SUPPLEMENTARY]
#   - Observed HVR ratio variable
#   - Model comparison: test if ratio fits better than components
#   - Controls: Age, ICV, SES (latent)
#
# Cognitive structure: Bifactor model (g + MEM + PRSP)
# - g: General cognitive factor (all 12 tests)
# - MEM: Memory specific factor (orthogonal to g)
# - PRSP: Processing speed specific factor (orthogonal to g)
#
# Key features:
# - Unadjusted brain volumes + ICV covariate (preserves sex differences)
# - Multi-group SEM by sex (tests sex-specific pathways)
# - Bootstrap inference (2000 samples for accurate indirect effect CIs)
# - Subfield measurement models (accounts for measurement error)
#
# Inputs:
#   - data/derivatives/hc-hvr_adj.rds: Brain volumes (unadjusted used)
#   - data/derivatives/lat-cog_vals.rds: Cognitive factor scores
#   - data/derivatives/ses_data.rds: SES indicator variables
#   - models/lavaan/sem_*.txt: SEM model definitions
#
# Outputs:
#   - models/fits/sem_analysis.rds: Fitted SEM models
#   - models/results/sem_params.rds: Parameter estimates
#   - models/results/sem_fit_measures.rds: Fit statistics
#   - outputs/tables/sem/*.{html,tex}: Publication tables
#   - outputs/figures/sem/*.png: Forest plots and path diagrams
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
library(lavaan)
library(semhelpinghands)  # For bootstrap CIs on standardized estimates
library(ggplot2)
library(ggtext)
library(gt)
library(progress)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))
source(here("R/utils/tables_core.R"))
source(here("R/utils/plotting_core.R"))
source(here("R/utils/plotting_pipeline.R"))
source(here("R/utils/export.R"))

# ----- Initialize -----
init_logger(log_level = "INFO")
log_script_start("11_sem_analysis.R")

# Load configuration
config <- load_config()
set_seed()

# ----- Constants -----
REFIT_MODELS <- get_script_setting("sem", "refit_models", default = FALSE)
REDO_TABLES <- get_script_setting("sem", "redo_tables", default = TRUE)
REDO_PLOTS <- get_script_setting("sem", "redo_plots", default = TRUE)
BOOTSTRAP <- get_script_setting("sem", "bootstrap", default = TRUE)
BOOTSTRAP_N <- get_script_setting("sem", "bootstrap_samples", default = 2000)
NETWORK <- get_script_setting("sem", "network", default = "LPP_CNN")

# Cognitive factors (from bifactor model)
# Note: 'g' is saved as 'COG' in the factor scores file
COG_FACS <- c("COG", "MEM", "PRSP")

# Define model specifications
# These are COMPETING models compared via AIC/BIC/R² (not combined in same equation)
MODELS.lst <- list(
  HC_COG          = "sem_hc_cog",           # HC → Cognition (ICV as covariate)
  HVR_COG         = "sem_hvr_cog",          # HVR → Cognition (self-normalizing)
  HC_RES_COG      = "sem_hc_res_cog"        # HC-residualized → Cognition (r=0 with ICV)
)

log_info("Network: %s", NETWORK)
log_info("Refit models: %s, Bootstrap: %s (%d samples)",
         REFIT_MODELS, BOOTSTRAP, BOOTSTRAP_N)
log_info("Models to fit: %s", paste(names(MODELS.lst), collapse = ", "))
log_info("Cognitive factors: %s", paste(COG_FACS, collapse = ", "))

# ===========================================================================
# Load Data
# ===========================================================================
log_section("Loading input data")

# Brain volumes (will use UNADJUSTED with ICV covariate)
hc_hvr_path <- get_data_path("processed", "hc_hvr_adjusted")
if (!check_files_exist(hc_hvr_path, stop_on_missing = FALSE)) {
  log_warn("Head-size adjusted data not found, running 05_adjust_headsize.R")
  source(here("R", "05_adjust_headsize.R"))
}
hc_hvr.lst <- read_rds_safe(hc_hvr_path, "Brain volume data")
validate_not_empty(hc_hvr.lst, "Brain volume data")

# Individual cognitive test scores (for joint measurement-structural SEM)
cog_clean_path <- here("data/derivatives/cog_clean-data.rds")
if (!check_files_exist(cog_clean_path, stop_on_missing = FALSE)) {
  log_warn("Cognitive test data not found, running 06_cognitive_factors.R")
  source(here("R", "06_cognitive_factors.R"))
}
cog_clean.lst <- read_rds_safe(cog_clean_path, "Cognitive test data")
validate_not_empty(cog_clean.lst, "Cognitive test data")

# IMDP (Index of Multiple Deprivation for England)
ses_path <- get_data_path("processed", "ses_data")
if (!check_files_exist(ses_path, stop_on_missing = FALSE)) {
  log_warn("SES data not found, running 04_clean_ses_cognitive.R")
  source(here("R", "04_clean_ses_cognitive.R"))
}
ses_data.dt <- read_rds_safe(ses_path, "SES data")
validate_not_empty(ses_data.dt, "SES data")
validate_columns(ses_data.dt, c("EID", "COUNTRY", "IMDP"), "IMDP data")

# Filter to England only (IMDP not comparable across countries)
ses_data.dt <- ses_data.dt[COUNTRY == "ENG", .(EID, IMDP)]
log_info("SES data filtered to England: %d subjects", ses_data.dt[, .N])

# Education and Scanner Site data (from covariates file)
log_info("Loading education and scanner site data")
covars_path <- get_data_path("processed", "covars_fst")
if (check_files_exist(covars_path, stop_on_missing = FALSE)) {
  covars.dt <- fst::read_fst(covars_path, as.data.table = TRUE)

  # Compute years of education from education level
  edu_years_map <- get_education_years_map()
  educ.dt <- covars.dt[, .(EDUC_char = as.character(EDUC_lv_ses2_0)), EID]
  educ.dt[, EDUC := edu_years_map[EDUC_char]]
  educ.dt[, EDUC_char := NULL]
  setkey(educ.dt, EID)

  log_info("Education data computed: %d subjects with valid education",
           educ.dt[!is.na(EDUC), .N])

  # Extract scanner site (Field 54 Assessment Centre) for clustered SEs
  site.dt <- covars.dt[, .(SITE = as.factor(SITE_ses2)), EID]
  site.dt <- site.dt[!is.na(SITE)]
  setkey(site.dt, EID)

  log_info("Scanner site data extracted: %d subjects with valid SITE",
           site.dt[, .N])
  log_info("  Sites: %s", paste(levels(site.dt$SITE), collapse = ", "))
} else {
  log_warn("Covariates file not found, education/site covariates will not be used")
  educ.dt <- NULL
  site.dt <- NULL
}

# ===========================================================================
# Data Preparation
# ===========================================================================
log_section("Preparing SEM data")

# Extract cross-sectional brain volumes
brain.dt <- hc_hvr.lst[[NETWORK]]$CRS$ALL

# Filter to unadjusted total volumes (bilateral sum computed below)
# Reference: Mak et al. (2017) - total hippocampal volume = sum of L + R
log_info("Filtering to unadjusted total volumes")

brain_total.dt <- brain.dt[

  ADJ == "NON" &
  SUBFIELD == "total" &
  SIDE %in% c("L", "R")
]

log_info("Brain data filtered: %d observations", nrow(brain_total.dt))

# Compute bilateral sums (L + R) for HC and LV
# Note: Sum, not mean - total bilateral volume is L + R
brain_wide.dt <- brain_total.dt[,
  .(
    HC = sum(HC, na.rm = TRUE),
    LV = sum(LV, na.rm = TRUE),
    HVR = sum(HC, na.rm = TRUE) / (sum(HC, na.rm = TRUE) + sum(LV, na.rm = TRUE))
  ),
  by = .(EID, INST, SEX, AGE, ICC)
]

log_info("Bilateral sums computed: %d subjects", brain_wide.dt[, uniqueN(EID)])

# Compute HC_RES for HC_RES_COG sensitivity analysis
# Residualize HC against ICV (ICC) within sex
# Reference: O'Brien et al. (2011) - residual method for ICV correction
# This removes ICV dependence completely (r = 0 by definition)
brain_wide.dt[, HC_RES := resid(lm(HC ~ ICC)), by = SEX]
log_info("HC_RES computed as within-sex residualized HC (r with ICV = 0)")

# Rename ICC to ICV for clarity in models
setnames(brain_wide.dt, "ICC", "ICV")

# Note: ICV remains in mL here; z-standardization below handles scale differences

# Merge individual cognitive test scores (for joint SEM measurement model)
log_info("Merging cognitive test scores")

# Each domain has duplicates from online/offline, so we need to handle carefully
# Merge domains keeping only unique EID + SESSION combinations
memory.dt <- cog_clean.lst$Memory[, .(EID, SESSION, PRS_mean_time, PRS_mean_inc, NUM, PRMEM_res_n)]
procspeed.dt <- cog_clean.lst$Proc_speed[, .(EID, SESSION, REACT, TRLS_alnum_time, TRLS_num_time, SYM_corr, SYM_try)]
reasexec.dt <- cog_clean.lst$Reas_Exec[, .(EID, SESSION, FLINT, MATS_corr, TOWER_corr)]

# Merge on EID and SESSION
cog_tests.dt <- memory.dt[procspeed.dt, on = c("EID", "SESSION"), nomatch = 0]
cog_tests.dt <- cog_tests.dt[reasexec.dt, on = c("EID", "SESSION"), nomatch = 0]

# Rename SESSION to INST for consistency and convert format
setnames(cog_tests.dt, "SESSION", "INST")

# Convert INST format from "ses2" to "ses-2"
cog_tests.dt[, INST := gsub("ses([0-9])", "ses-\\1", INST)]

# Select only cross-sectional data (ses-2)
cog_tests.dt <- cog_tests.dt[INST == "ses-2"]

log_info("Cognitive tests prepared: %d subjects", cog_tests.dt[, uniqueN(EID)])

# Merge with brain data
data.dt <- cog_tests.dt[
  brain_wide.dt,
  on = c("EID", "INST"),
  nomatch = 0
]

# Merge IMDP (deprivation index for England)
data.dt <- ses_data.dt[
  data.dt,
  on = "EID",
  nomatch = NA
]

# Merge education (years of schooling)
if (!is.null(educ.dt)) {
  data.dt <- educ.dt[
    data.dt,
    on = "EID",
    nomatch = NA
  ]
  log_info("Education merged: %d subjects with valid EDUC", data.dt[!is.na(EDUC), .N])
}

# Merge scanner site (for clustered SEs)
if (!is.null(site.dt)) {
  data.dt <- site.dt[
    data.dt,
    on = "EID",
    nomatch = NA
  ]
  log_info("Scanner site merged: %d subjects with valid SITE", data.dt[!is.na(SITE), .N])
}

# Recode SEX to factor for multi-group (BP11: explicit level ordering)
# Female = group 1, Male = group 2 (alphabetical order, explicitly set for clarity)
# This ensures consistent group ordering across all lavaan analyses
sex_levels <- get_factor_levels("SEX")  # Returns c("Female", "Male")
data.dt[, SEX := factor(SEX, levels = sex_levels)]
log_info("SEX factor levels set: %s (group 1 = %s, group 2 = %s)",
         paste(sex_levels, collapse = ", "), sex_levels[1], sex_levels[2])

# Remove rows with missing critical variables
data.dt <- data.dt[!is.na(EID) & !is.na(SEX) & !is.na(AGE) & !is.na(ICV)]

log_info("Data merged: N = %d", data.dt[, uniqueN(EID)])
log_info("  Females: %d, Males: %d",
         data.dt[SEX == "Female", uniqueN(EID)],
         data.dt[SEX == "Male", uniqueN(EID)])

# Validate final data
validate_columns(
  data.dt,
  c("EID", "SEX", "AGE", "ICV", "HC", "LV", "HVR", "HC_RES"),
  "Brain structure variables"
)
validate_columns(
  data.dt,
  c("PRS_mean_time", "PRS_mean_inc", "NUM", "PRMEM_res_n", "FLINT", "MATS_corr",
    "TRLS_alnum_time", "TOWER_corr", "REACT", "TRLS_num_time", "SYM_corr", "SYM_try"),
  "Cognitive test variables"
)
validate_columns(data.dt, "IMDP", "Deprivation index")

# Check for sufficient data in each sex
log_info("Data availability by sex:")
log_info("  Female: HC=%d, LV=%d, HVR=%d, PRS_mean_time=%d",
         data.dt[SEX == "Female" & !is.na(HC), .N],
         data.dt[SEX == "Female" & !is.na(LV), .N],
         data.dt[SEX == "Female" & !is.na(HVR), .N],
         data.dt[SEX == "Female" & !is.na(PRS_mean_time), .N])
log_info("  Male: HC=%d, LV=%d, HVR=%d, PRS_mean_time=%d",
         data.dt[SEX == "Male" & !is.na(HC), .N],
         data.dt[SEX == "Male" & !is.na(LV), .N],
         data.dt[SEX == "Male" & !is.na(HVR), .N],
         data.dt[SEX == "Male" & !is.na(PRS_mean_time), .N])

# ---------------------------------------------------------------------------
# Z-standardization within sex groups
# ---------------------------------------------------------------------------
# Purpose: Fix numerical scale issues (ICV ~1500 mL vs HVR ~0.75)
#
# Why within-sex (not pooled)?
# 1. Multi-group SEM compares path SLOPES across groups, not means.
#    Within-group standardization preserves variance structure for slopes.
# 2. For pooled models, SEX is a control covariate. The SEX→brain path
#    becomes ~0 by construction, but this is irrelevant—we care about
#    brain→cognition paths, which are unaffected (slopes depend on
#    covariance, not means).
#
# Limitation: Do not interpret SEX coefficients on brain measures in pooled
# models; these are artificially zeroed by within-group standardization.
#
# NOTE: Cognitive test indicators are NOT standardized to preserve the
# measurement model structure (factor loadings).
# ---------------------------------------------------------------------------
log_info("Standardizing structural variables within sex groups")

# Create quadratic age term from CENTERED age (before standardization)
# Rationale: Brain-cognition relationships and brain atrophy are non-linear
# across the lifespan; AGE_sq captures accelerating decline at older ages
# IMPORTANT: Center AGE first to avoid collinearity (r(AGE, AGE^2) ≈ 0.997 without centering)
data.dt[, AGE_c := AGE - mean(AGE), by = SEX]  # Center within sex
data.dt[, AGE_sq := AGE_c^2]                    # Square centered age
data.dt[, AGE_c := NULL]                        # Remove temporary column
log_info("Created AGE_sq (quadratic age term, centered to avoid collinearity)")

structural_vars <- c(
  # Brain volumes
  "HC", "LV", "HVR", "HC_RES",
  # Covariates
  "AGE", "AGE_sq", "ICV", "IMDP", "EDUC"
)

# Standardize within each sex group
for (var in structural_vars) {
  data.dt[, (var) := scale(.SD[[1]]), by = SEX, .SDcols = var]
}

log_info("Structural variables standardized: %d variables", length(structural_vars))

# ===========================================================================
# SEM Model Fitting (Multi-Group by Sex)
# ===========================================================================
# MISSING DATA HANDLING (CL6):
# Models use FIML (Full Information Maximum Likelihood) via missing = "fiml"
#
# FIML assumes data is Missing At Random (MAR), meaning:
# - Missingness can depend on OBSERVED variables in the model
# - Missingness should NOT depend on the missing values themselves
#
# Potential violations (cognitive outcomes):
# - Participants with cognitive decline may be more likely to miss tests
# - Those with severe impairment may have more missing data
#
# Mitigation:
# - We include age and ICV as covariates (observed predictors of missingness)
# - Sample is population-based (UK Biobank), not clinical
# - Severe cognitive impairment is rare in this age range
#
# Users should interpret results with awareness of this assumption.
# For sensitivity, listwise deletion can be tested by setting missing = "listwise"

# ===========================================================================
# Fit Pooled SEM Models (Main Effects)
# ===========================================================================
# Pooled models estimate overall brain-cognition effects on the full sample.
# These provide main effect estimates before testing sex differences.
# SEX is included as a covariate in pooled models.
#
# Path coefficient comparisons (β_HVR→g vs β_HC→g vs β_HC_RES→g) are
# performed in manuscript.qmd using bootstrap CIs.
# ===========================================================================

log_section("Fitting pooled SEM models (main effects)")

# Output files
sem_fits_path <- get_data_path("models", "fit", "sem_analysis")
sem_pooled_path <- get_data_path("models", "fit", "sem_pooled")

# Define pooled model files (separate from multi-group models)
POOLED_MODELS.lst <- list(
  HC_COG_POOLED     = "sem_hc_cog_pooled",
  HVR_COG_POOLED    = "sem_hvr_cog_pooled",
  HC_RES_COG_POOLED = "sem_hc_res_cog_pooled"
)

# Define use_cluster (needed for both pooled and multi-group)
use_cluster <- !is.null(site.dt) && "SITE" %in% names(data.dt) && !BOOTSTRAP

# Load existing pooled fits if available
if (!REFIT_MODELS && check_files_exist(sem_pooled_path, stop_on_missing = FALSE)) {
  log_info("Loading existing pooled SEM fits")
  pooled_fits.lst <- read_rds_safe(sem_pooled_path, "Pooled SEM fits")
} else {
  log_info("Fitting new pooled SEM models")
  pooled_fits.lst <- list()
}

# Fit pooled models
for (mod_name in names(POOLED_MODELS.lst)) {

  # Skip if already fitted
  if (!REFIT_MODELS &&
      !is.null(pooled_fits.lst[[mod_name]]) &&
      !"ERROR" %in% names(pooled_fits.lst[[mod_name]])) {
    log_info("  Skipping %s (already fitted)", mod_name)
    next
  }

  # Load pooled model syntax
  mod_key <- POOLED_MODELS.lst[[mod_name]]
  modfile <- get_data_path("models", "definition", mod_key)

  if (!file.exists(modfile)) {
    log_error("Model file not found: %s", modfile)
    next
  }

  model_syntax <- paste(readLines(modfile), collapse = "\n")

  log_info("  Fitting %s (pooled, N = %d)", mod_name, nrow(data.dt))

  pooled_fits.lst[[mod_name]] <- tryCatch(
    {
      if (use_cluster) {
        log_info("    Using robust SE (MLR) with site-clustered SEs")
        sem(
          model = model_syntax,
          data = data.dt,
          missing = "fiml",
          estimator = "MLR",
          se = "robust",
          cluster = "SITE"
        )
      } else {
        log_info("    Using robust SE (MLR estimator)")
        sem(
          model = model_syntax,
          data = data.dt,
          missing = "fiml",
          estimator = "MLR",
          se = "robust"
        )
      }
    },
    error = function(e) {
      log_error("Error fitting %s: %s", mod_name, e$message)
      list(ERROR = e$message)
    }
  )

  # Save after each fit
  ensure_directory(dirname(sem_pooled_path))
  write_rds_safe(pooled_fits.lst, sem_pooled_path, "Pooled SEM fits")
  log_info("  Saved: %s", sem_pooled_path)
}

log_info("Pooled SEM model fitting complete")

# ===========================================================================
# Fit Multi-Group SEM Models (Sex Differences)
# ===========================================================================
# Multi-group models estimate separate paths for males and females.
# Defined parameters (diff_*) test whether paths differ by sex.
# Bootstrap CIs recommended for final inference on sex differences.
# ===========================================================================

log_section("Fitting multi-group SEM models (sex differences)")

# Load existing fits if available
if (!REFIT_MODELS && check_files_exist(sem_fits_path, stop_on_missing = FALSE)) {
  log_info("Loading existing multi-group SEM fits")
  fits.lst <- read_rds_safe(sem_fits_path, "SEM model fits")
} else {
  log_info("Fitting new multi-group SEM models")
  fits.lst <- list()
}

# Progress bar
pb <- progress_bar$new(
  format = "SEM | :what [:bar] :current/:total\n",
  total = length(MODELS.lst),
  clear = FALSE,
  width = 75
)

# Fit multi-group models
for (mod_name in names(MODELS.lst)) {
  pb$tick(tokens = list(what = sprintf("Fitting: %s", mod_name)))

  # Skip if already fitted
  if (!REFIT_MODELS &&
      !is.null(fits.lst[[mod_name]]) &&
      !"ERROR" %in% names(fits.lst[[mod_name]])) {
    log_info("  Skipping %s (already fitted)", mod_name)
    next
  }

  # Load model definition
  mod_key <- MODELS.lst[[mod_name]]
  modfile <- get_data_path("models", "definition", mod_key)

  if (!file.exists(modfile)) {
    log_error("Model file not found: %s", modfile)
    next
  }

  model_syntax <- paste(readLines(modfile), collapse = "\n")

  log_info("  Fitting %s with multi-group by sex", mod_name)

  fits.lst[[mod_name]] <- tryCatch(
    {
      if (BOOTSTRAP) {
        # Bootstrap: cluster not supported, but bootstrap accounts for sampling variability
        log_info("    Using bootstrap SE (n = %d)", BOOTSTRAP_N)
        sem(
          model = model_syntax,
          data = data.dt,
          group = "SEX",
          missing = "fiml",
          se = "bootstrap",
          bootstrap = BOOTSTRAP_N,
          parallel = "snow",
          ncpus = 18,
          estimator = "ML"
        )
      } else {
        # Robust SEs (MLR) - note: cluster argument incompatible with multigroup in lavaan
        log_info("    Using robust SE (MLR estimator)")
        sem(
          model = model_syntax,
          data = data.dt,
          group = "SEX",
          missing = "fiml",
          estimator = "MLR",
          se = "robust"
        )
      }
    },
    error = function(e) {
      log_error("Error fitting %s: %s", mod_name, e$message)
      list(ERROR = e$message)
    }
  )

  # Save after each fit
  ensure_directory(dirname(sem_fits_path))
  write_rds_safe(fits.lst, sem_fits_path, "SEM model fits")
  log_info("  Saved: %s", sem_fits_path)
}

log_info("Multi-group SEM model fitting complete")

# ===========================================================================
# Extract Parameters
# ===========================================================================
log_section("Extracting model parameters")

params.lst <- list()
fit_measures.lst <- list()

# Combine pooled and multi-group fits for parameter extraction
all_fits.lst <- c(pooled_fits.lst, fits.lst)

for (mod_name in names(all_fits.lst)) {
  fit_obj <- all_fits.lst[[mod_name]]

  if (is.null(fit_obj) || "ERROR" %in% names(fit_obj)) {
    log_warn("Skipping %s (fit failed or not available)", mod_name)
    next
  }

  # Check if model converged
  if (!lavInspect(fit_obj, "converged")) {
    log_warn("Model did not converge: %s. Skipping parameter extraction.", mod_name)
    next
  }

  # Determine CI method used
  ci_method <- if (BOOTSTRAP) "bootstrap percentile" else "robust Wald"

  # Extract standardized parameters with proper bootstrap CIs
  # Using semhelpinghands::standardizedSolution_boot_ci() for bootstrap models
  # This computes bootstrap CIs directly on standardized estimates (not delta-method)
  if (BOOTSTRAP) {
    log_info("    Computing bootstrap CIs for standardized estimates")
    params_dt <- tryCatch({
      standardizedSolution_boot_ci(fit_obj, level = 0.95, boot_ci_type = "perc") |>
        as.data.table() |>
        setkey(lhs, op, rhs)
    }, error = function(e) {
      log_warn("    Failed to compute bootstrap CIs: %s. Using delta-method.", e$message)
      standardizedSolution(fit_obj) |>
        as.data.table() |>
        setkey(lhs, op, rhs)
    })
  } else {
    params_dt <- standardizedSolution(fit_obj) |>
      as.data.table() |>
      setkey(lhs, op, rhs)
  }

  # Add metadata
  params_dt[, `:=`(
    MODEL = mod_name,
    CI_METHOD = ci_method
  )]

  params.lst[[mod_name]] <- params_dt

  # Extract fit measures
  fit_indices <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea",
                   "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic")

  # For robust estimation, also get scaled indices
  if (!BOOTSTRAP) {
    fit_indices_robust <- c("chisq.scaled", "cfi.scaled", "tli.scaled",
                           "rmsea.scaled")
    fit_indices <- c(fit_indices, fit_indices_robust)
  }

  fit_measures.lst[[mod_name]] <- tryCatch({
    fitMeasures(fit_obj, fit_indices) |>
      as.list() |>
      as.data.table() |>
      (\(DT) {
        DT[, `:=`(
          MODEL = mod_name,
          CI_METHOD = ci_method,
          BOOTSTRAP = BOOTSTRAP,
          BOOTSTRAP_N = if (BOOTSTRAP) BOOTSTRAP_N else NA
        )]
        setcolorder(DT, "MODEL")
        DT
      })()
  }, error = function(e) {
    log_warn("Error extracting fit measures for %s: %s", mod_name, e$message)
    NULL
  })
}

# Save extracted parameters
log_info("Saving parameter estimates and fit measures")

params_path <- get_data_path("models", "results", "sem_params")
ensure_directory(dirname(params_path))
write_rds_safe(params.lst, params_path, "SEM model parameters")

fit_measures_path <- get_data_path("models", "results", "sem_fit_measures")
ensure_directory(dirname(fit_measures_path))
write_rds_safe(fit_measures.lst, fit_measures_path, "SEM fit measures")

# ===========================================================================
# Test Measurement Invariance
# ===========================================================================
invariance_path <- get_data_path("models", "results", "sem_invariance")

# Skip if invariance results already exist and not refitting
if (!REFIT_MODELS && check_files_exist(invariance_path, stop_on_missing = FALSE)) {
  log_info("Loading existing measurement invariance results")
  invariance.lst <- read_rds_safe(invariance_path, "Measurement invariance tests")
} else {
  log_section("Testing measurement invariance across sex")

  invariance.lst <- list()

  # Helper: check if model has valid (non-NA) SEs
  has_valid_se <- function(fit) {
    if (is.null(fit)) return(FALSE)
    se_vals <- tryCatch(lavInspect(fit, "se"), error = function(e) NULL)
    if (is.null(se_vals)) return(FALSE)
    !any(is.na(unlist(se_vals)))
  }

  # Helper: safely extract fit measure with fallback
  safe_fit <- function(fit, measure) {
    if (is.null(fit)) return(NA_real_)
    tryCatch(fitMeasures(fit)[[measure]], error = function(e) NA_real_)
  }

  for (mod_name in names(MODELS.lst)) {
    # Skip pooled models (no multi-group structure)
    if (mod_name == "HC_HVR_COMPARE") {
      log_info("Skipping invariance test for %s (pooled model, no multi-group)", mod_name)
      next
    }

    log_info("Testing invariance for %s", mod_name)

    mod_key <- MODELS.lst[[mod_name]]
    modfile <- get_data_path("models", "definition", mod_key)
    model_syntax <- paste(readLines(modfile), collapse = "\n")

    # Note: cluster argument incompatible with multigroup in lavaan
    # Using robust MLR without site clustering for invariance tests

    # Configural invariance (baseline - different parameters by group)
    fit_configural <- tryCatch(
      sem(model_syntax, data = data.dt, group = "SEX", missing = "fiml",
          estimator = "MLR"),
      error = function(e) NULL
    )

    # Metric invariance (equal loadings)
    fit_metric <- tryCatch(
      sem(model_syntax, data = data.dt, group = "SEX", group.equal = "loadings",
          missing = "fiml", estimator = "MLR"),
      error = function(e) NULL
    )

    # Scalar invariance (equal loadings + intercepts)
    fit_scalar <- tryCatch(
      sem(model_syntax, data = data.dt, group = "SEX",
          group.equal = c("loadings", "intercepts"),
          missing = "fiml", estimator = "MLR"),
      error = function(e) NULL
    )

    # Store fits
    invariance.lst[[mod_name]] <- list(
      configural = fit_configural,
      metric = fit_metric,
      scalar = fit_scalar
    )

    # Report fit indices (even if LRT fails)
    if (!is.null(fit_configural)) {
      # Check SE validity for each model
      se_valid_conf <- has_valid_se(fit_configural)
      se_valid_metr <- has_valid_se(fit_metric)
      se_valid_scal <- has_valid_se(fit_scalar)

      if (!se_valid_conf) {
        log_warn("  %s configural: SE computation failed (information matrix not invertible)", mod_name)
        log_warn("    This typically indicates near-collinearity or boundary estimates")
        log_warn("    Reporting naive fit indices (robust CFI/RMSEA may be NA)")
      }

      # Extract fit indices (use naive if robust fails)
      cfi_conf <- safe_fit(fit_configural, "cfi")
      rmsea_conf <- safe_fit(fit_configural, "rmsea")

      log_info("  Invariance test results:")
      log_info("    Configural: CFI=%.3f, RMSEA=%.3f%s",
               cfi_conf, rmsea_conf,
               if (!se_valid_conf) " (naive indices)" else "")

      if (!is.null(fit_metric)) {
        cfi_metr <- safe_fit(fit_metric, "cfi")
        delta_cfi_metr <- cfi_metr - cfi_conf
        log_info("    Metric: CFI=%.3f (ΔCFI=%.3f)%s",
                 cfi_metr, delta_cfi_metr,
                 if (!se_valid_metr) " (naive indices)" else "")
      }

      if (!is.null(fit_scalar) && !is.null(fit_metric)) {
        cfi_scal <- safe_fit(fit_scalar, "cfi")
        cfi_metr <- safe_fit(fit_metric, "cfi")
        delta_cfi_scal <- cfi_scal - cfi_metr
        log_info("    Scalar: CFI=%.3f (ΔCFI=%.3f)%s",
                 cfi_scal, delta_cfi_scal,
                 if (!se_valid_scal) " (naive indices)" else "")
      }

      # Attempt LRT only if all models have valid SEs
      if (se_valid_conf && se_valid_metr) {
        comparison <- tryCatch(
          lavTestLRT(fit_configural, fit_metric, fit_scalar),
          error = function(e) {
            log_warn("  lavTestLRT failed for %s: %s", mod_name, e$message)
            NULL
          }
        )
        if (!is.null(comparison)) {
          log_info("  Chi-square difference tests computed successfully")
          invariance.lst[[mod_name]]$lrt <- comparison
        }
      } else {
        log_info("  Skipping chi-square difference test (SE issues in one or more models)")
      }
    }
  }

  # Save invariance tests
  ensure_directory(dirname(invariance_path))
  write_rds_safe(invariance.lst, invariance_path, "Measurement invariance tests")
}  # End of invariance testing block

# ===========================================================================
# Generate Summary Tables
# ===========================================================================
if (REDO_TABLES) {
  log_section("Generating SEM tables")

  tables_dir <- get_output_path("tables")
  ensure_directory(tables_dir)

  # --- Table 1: Model Fit Comparison ---
  fit_measures.dt <- rbindlist(fit_measures.lst, fill = TRUE)

  if (nrow(fit_measures.dt) > 0) {
    # Select key fit indices
    fit_cols <- c("MODEL", "chisq", "df", "cfi", "tli", "rmsea",
                  "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic")
    fit_summary.dt <- fit_measures.dt[, .SD, .SDcols = intersect(names(fit_measures.dt), fit_cols)]

    # Create gt table
    fit_summary_gt <- fit_summary.dt |>
      gt() |>
      tab_header(
        title = "SEM Model Fit Statistics",
        subtitle = "Multi-group models (by sex) for brain-cognition relationships"
      ) |>
      cols_label(
        MODEL = "Model",
        chisq = md("*χ*²"),
        df = "df",
        cfi = "CFI",
        tli = "TLI",
        rmsea = "RMSEA",
        rmsea.ci.lower = "RMSEA 90% CI Lower",
        rmsea.ci.upper = "RMSEA 90% CI Upper",
        srmr = "SRMR",
        aic = "AIC",
        bic = "BIC"
      ) |>
      fmt_number(
        columns = c(chisq, cfi, tli, rmsea, rmsea.ci.lower, rmsea.ci.upper, srmr),
        decimals = 3
      ) |>
      fmt_number(
        columns = c(aic, bic),
        decimals = 1
      ) |>
      tab_source_note("Good fit: CFI/TLI > 0.95, RMSEA < 0.06, SRMR < 0.08") |>
      tab_source_note("Model comparison: Lower AIC/BIC indicates better fit")

    # Save
    gt::gtsave(fit_summary_gt, file.path(tables_dir, "model_fit_comparison.html"))
    gt::gtsave(fit_summary_gt, file.path(tables_dir, "model_fit_comparison.tex"))

    log_info("Model fit table saved")
  }

  # --- Table 2: Key Structural Paths by Sex ---
  # Extract structural regressions from each model
  structural.dt <- rbindlist(
    lapply(names(params.lst), function(mod) {
      dt <- params.lst[[mod]]

      # Filter to structural regressions (brain → cognition)
      struct_dt <- dt[
        op == "~" &
        lhs %in% c("g", "MEM_s", "PRSP_s") &
        rhs %in% c("HC", "LV", "HVR")
      ]

      if (nrow(struct_dt) == 0) return(NULL)

      # Handle pooled models (no group column)
      has_group <- "group" %in% names(struct_dt)

      struct_dt[, .(
        MODEL = mod,
        GROUP = if (has_group) group else NA_integer_,
        OUTCOME = lhs,
        PREDICTOR = rhs,
        EST_STD = est.std,
        SE = se,
        CI_LOWER = ci.lower,
        CI_UPPER = ci.upper,
        PVALUE = pvalue
      )]
    }),
    fill = TRUE
  )

  if (!is.null(structural.dt) && nrow(structural.dt) > 0) {
    # Format for table
    structural.dt[, `:=`(
      PATH = sprintf("%s → %s", PREDICTOR, OUTCOME),
      EST_CI = sprintf("%.3f [%.3f, %.3f]%s",
                       EST_STD, CI_LOWER, CI_UPPER,
                       ifelse(PVALUE < 0.001, "***",
                       ifelse(PVALUE < 0.01, "**",
                       ifelse(PVALUE < 0.05, "*", ""))))
    )]

    # Pivot to show Male vs Female
    structural_wide.dt <- dcast(
      structural.dt,
      MODEL + PATH ~ GROUP,
      value.var = "EST_CI"
    )

    # Rename group columns: 1=Female, 2=Male, NA=Pooled
    old_names <- names(structural_wide.dt)
    new_names <- old_names
    new_names[old_names == "1"] <- "Female"
    new_names[old_names == "2"] <- "Male"
    new_names[old_names == "NA"] <- "Pooled"
    setnames(structural_wide.dt, old_names, new_names)

    structural_gt <- structural_wide.dt |>
      gt(groupname_col = "MODEL") |>
      tab_header(
        title = "Structural Paths: Brain → Cognition",
        subtitle = "Standardized coefficients [95% CI] by sex"
      ) |>
      tab_source_note("* p<.05, ** p<.01, *** p<.001")

    gt::gtsave(structural_gt, file.path(tables_dir, "structural_paths_by_sex.html"))
    gt::gtsave(structural_gt, file.path(tables_dir, "structural_paths_by_sex.tex"))

    log_info("Structural paths table saved")
  }

  log_info("SEM tables saved to: %s", tables_dir)
}

# ===========================================================================
# Generate Forest Plots
# ===========================================================================
if (REDO_PLOTS) {
  log_section("Generating SEM forest plots")

  plots_dir <- file.path(get_output_path("figures"), "sem")
  ensure_directory(plots_dir)

  # --- Plot 1: Structural Paths Forest Plot ---
  if (!is.null(structural.dt) && nrow(structural.dt) > 0) {
    # Convert GROUP to factor with labels
    structural.dt[, GROUP := factor(GROUP, levels = c(1, 2), labels = c("Female", "Male"))]

    # Add significance indicator
    structural.dt[, SIGN := PVALUE < 0.05]

    p1 <- ggplot(
      structural.dt,
      aes(x = PATH, y = EST_STD, ymin = CI_LOWER, ymax = CI_UPPER,
          color = GROUP, shape = SIGN)
    ) +
      geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      facet_wrap(~ MODEL, scales = "free_y", ncol = 1) +
      coord_flip() +
      scale_color_manual(values = c("Female" = "darkred", "Male" = "midnightblue")) +
      scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1)) +
      labs(
        title = "Brain → Cognition Structural Paths",
        subtitle = sprintf("Standardized coefficients with 95%% %s CI",
                          if (BOOTSTRAP) "bootstrap" else "robust"),
        x = "Structural Path",
        y = "Standardized Effect",
        color = "Sex",
        shape = "p < 0.05"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )

    save_plot(
      p1,
      file.path(plots_dir, "forest_structural_paths.png"),
      width = 7,
      height = 5
    )

    log_info("Forest plot saved")
  }

  # --- Plot 2: Model Comparison (AIC/BIC) ---
  if (nrow(fit_measures.dt) > 0 && all(c("aic", "bic") %in% names(fit_measures.dt))) {
    comparison.dt <- melt(
      fit_measures.dt[, .(MODEL, aic, bic)],
      id.vars = "MODEL",
      variable.name = "METRIC",
      value.name = "VALUE"
    )

    p2 <- ggplot(comparison.dt, aes(x = MODEL, y = VALUE, fill = METRIC)) +
      geom_col(position = "dodge") +
      scale_fill_manual(
        values = c("aic" = "#00A087", "bic" = "#3C5488"),
        labels = c("AIC", "BIC")
      ) +
      labs(
        title = "Model Comparison: Information Criteria",
        subtitle = "Lower values indicate better fit",
        x = "Model",
        y = "Information Criterion Value",
        fill = "Metric"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )

    save_plot(
      p2,
      file.path(plots_dir, "model_comparison_ic.png"),
      width = 8,
      height = 6
    )

    log_info("Model comparison plot saved")
  }

  log_info("SEM plots saved to: %s", plots_dir)
}

# ===========================================================================
# Pooled Model Beta Comparison (CI Overlap Method)
# ===========================================================================
# Compare path coefficients across pooled models to test:
#   - Is HVR a better predictor of cognition than HC?
#   - Does residualization (HC_RES) hurt prediction compared to HVR?
# Method: Non-overlapping 95% CIs indicate significant difference (p < 0.05)
# Reference: Cumming & Finch (2005) - CI-based inference approach
# ===========================================================================
log_section("Comparing pooled model path coefficients")

# Extract brain→cognition paths from pooled models
pooled_paths.dt <- rbindlist(

  lapply(names(pooled_fits.lst), function(mod_name) {
    fit_obj <- pooled_fits.lst[[mod_name]]
    if (is.null(fit_obj) || "ERROR" %in% names(fit_obj)) return(NULL)

    # Get standardized solution with CIs
    if (BOOTSTRAP) {
      params_dt <- tryCatch({
        standardizedSolution_boot_ci(fit_obj, level = 0.95, boot_ci_type = "perc") |>
          as.data.table()
      }, error = function(e) {
        standardizedSolution(fit_obj) |> as.data.table()
      })
    } else {
      params_dt <- standardizedSolution(fit_obj) |> as.data.table()
    }

    # Filter to structural paths (brain → cognition)
    brain_vars <- c("HC", "HVR", "HC_RES")
    cog_vars <- c("g", "MEM_s", "PRSP_s")

    struct_dt <- params_dt[
      op == "~" &
      lhs %in% cog_vars &
      rhs %in% brain_vars
    ]

    if (nrow(struct_dt) == 0) return(NULL)

    struct_dt[, .(
      MODEL = mod_name,
      OUTCOME = lhs,
      PREDICTOR = rhs,
      BETA = est.std,
      SE = se,
      CI_LOWER = ci.lower,
      CI_UPPER = ci.upper,
      PVALUE = pvalue
    )]
  }),
  fill = TRUE
)

if (!is.null(pooled_paths.dt) && nrow(pooled_paths.dt) > 0) {

  # Function to test CI overlap (Cumming & Finch, 2005)
  # Non-overlapping CIs at 95% level implies p < 0.05
  ci_overlap_test <- function(ci1_lower, ci1_upper, ci2_lower, ci2_upper) {
    # Handle NA values (e.g., from Heywood cases)
    if (any(is.na(c(ci1_lower, ci1_upper, ci2_lower, ci2_upper)))) {
      return(list(overlap = NA, overlap_prop = NA))
    }
    # Check if intervals overlap
    overlap <- !(ci1_upper < ci2_lower || ci2_upper < ci1_lower)
    # Calculate overlap proportion (for reporting)
    if (overlap) {
      overlap_start <- max(ci1_lower, ci2_lower)
      overlap_end <- min(ci1_upper, ci2_upper)
      overlap_width <- overlap_end - overlap_start
      avg_width <- ((ci1_upper - ci1_lower) + (ci2_upper - ci2_lower)) / 2
      overlap_prop <- overlap_width / avg_width
    } else {
      overlap_prop <- 0
    }
    list(overlap = overlap, overlap_prop = overlap_prop)
  }

  # Compare models for each cognitive outcome
  log_info("Comparing brain→cognition paths across pooled models:")

  comparison_results <- list()

  for (outcome in c("g", "MEM_s", "PRSP_s")) {
    outcome_dt <- pooled_paths.dt[OUTCOME == outcome]

    if (nrow(outcome_dt) < 2) next

    # Get paths for each model
    hc_path <- outcome_dt[PREDICTOR == "HC"]
    hvr_path <- outcome_dt[PREDICTOR == "HVR"]
    hc_res_path <- outcome_dt[PREDICTOR == "HC_RES"]

    log_info("\n  %s (general cognition):", outcome)

    # Compare HVR vs HC
    if (nrow(hvr_path) > 0 && nrow(hc_path) > 0) {
      test_hvr_hc <- ci_overlap_test(
        hvr_path$CI_LOWER, hvr_path$CI_UPPER,
        hc_path$CI_LOWER, hc_path$CI_UPPER
      )
      diff_hvr_hc <- hvr_path$BETA - hc_path$BETA
      sig_hvr_hc <- if (is.na(test_hvr_hc$overlap)) "NA (missing CI)" else if (!test_hvr_hc$overlap) "SIGNIFICANT (CIs don't overlap)" else "not significant"

      log_info("    HVR vs HC: β_HVR=%.3f [%.3f, %.3f], β_HC=%.3f [%.3f, %.3f]",
               hvr_path$BETA, hvr_path$CI_LOWER, hvr_path$CI_UPPER,
               hc_path$BETA, hc_path$CI_LOWER, hc_path$CI_UPPER)
      log_info("      Δβ = %.3f, %s", diff_hvr_hc, sig_hvr_hc)

      comparison_results[[paste0(outcome, "_HVR_vs_HC")]] <- data.table(
        OUTCOME = outcome,
        COMPARISON = "HVR vs HC",
        BETA_1 = hvr_path$BETA,
        CI_1 = sprintf("[%.3f, %.3f]", hvr_path$CI_LOWER, hvr_path$CI_UPPER),
        BETA_2 = hc_path$BETA,
        CI_2 = sprintf("[%.3f, %.3f]", hc_path$CI_LOWER, hc_path$CI_UPPER),
        DELTA_BETA = diff_hvr_hc,
        CI_OVERLAP = test_hvr_hc$overlap,
        SIGNIFICANT = if (is.na(test_hvr_hc$overlap)) NA else !test_hvr_hc$overlap
      )
    }

    # Compare HVR vs HC_RES
    if (nrow(hvr_path) > 0 && nrow(hc_res_path) > 0) {
      test_hvr_res <- ci_overlap_test(
        hvr_path$CI_LOWER, hvr_path$CI_UPPER,
        hc_res_path$CI_LOWER, hc_res_path$CI_UPPER
      )
      diff_hvr_res <- hvr_path$BETA - hc_res_path$BETA
      sig_hvr_res <- if (is.na(test_hvr_res$overlap)) "NA (missing CI)" else if (!test_hvr_res$overlap) "SIGNIFICANT (CIs don't overlap)" else "not significant"

      log_info("    HVR vs HC_RES: β_HVR=%.3f [%.3f, %.3f], β_HC_RES=%.3f [%.3f, %.3f]",
               hvr_path$BETA, hvr_path$CI_LOWER, hvr_path$CI_UPPER,
               hc_res_path$BETA, hc_res_path$CI_LOWER, hc_res_path$CI_UPPER)
      log_info("      Δβ = %.3f, %s", diff_hvr_res, sig_hvr_res)

      comparison_results[[paste0(outcome, "_HVR_vs_HC_RES")]] <- data.table(
        OUTCOME = outcome,
        COMPARISON = "HVR vs HC_RES",
        BETA_1 = hvr_path$BETA,
        CI_1 = sprintf("[%.3f, %.3f]", hvr_path$CI_LOWER, hvr_path$CI_UPPER),
        BETA_2 = hc_res_path$BETA,
        CI_2 = sprintf("[%.3f, %.3f]", hc_res_path$CI_LOWER, hc_res_path$CI_UPPER),
        DELTA_BETA = diff_hvr_res,
        CI_OVERLAP = test_hvr_res$overlap,
        SIGNIFICANT = if (is.na(test_hvr_res$overlap)) NA else !test_hvr_res$overlap
      )
    }

    # Compare HC vs HC_RES (sanity check)
    if (nrow(hc_path) > 0 && nrow(hc_res_path) > 0) {
      test_hc_res <- ci_overlap_test(
        hc_path$CI_LOWER, hc_path$CI_UPPER,
        hc_res_path$CI_LOWER, hc_res_path$CI_UPPER
      )
      diff_hc_res <- hc_path$BETA - hc_res_path$BETA
      sig_hc_res <- if (is.na(test_hc_res$overlap)) "NA (missing CI)" else if (!test_hc_res$overlap) "SIGNIFICANT (CIs don't overlap)" else "not significant"

      log_info("    HC vs HC_RES: β_HC=%.3f [%.3f, %.3f], β_HC_RES=%.3f [%.3f, %.3f]",
               hc_path$BETA, hc_path$CI_LOWER, hc_path$CI_UPPER,
               hc_res_path$BETA, hc_res_path$CI_LOWER, hc_res_path$CI_UPPER)
      log_info("      Δβ = %.3f, %s", diff_hc_res, sig_hc_res)

      comparison_results[[paste0(outcome, "_HC_vs_HC_RES")]] <- data.table(
        OUTCOME = outcome,
        COMPARISON = "HC vs HC_RES",
        BETA_1 = hc_path$BETA,
        CI_1 = sprintf("[%.3f, %.3f]", hc_path$CI_LOWER, hc_path$CI_UPPER),
        BETA_2 = hc_res_path$BETA,
        CI_2 = sprintf("[%.3f, %.3f]", hc_res_path$CI_LOWER, hc_res_path$CI_UPPER),
        DELTA_BETA = diff_hc_res,
        CI_OVERLAP = test_hc_res$overlap,
        SIGNIFICANT = if (is.na(test_hc_res$overlap)) NA else !test_hc_res$overlap
      )
    }
  }

  # Create comparison summary table
  if (length(comparison_results) > 0) {
    beta_comparison.dt <- rbindlist(comparison_results)

    # Save comparison table
    if (REDO_TABLES) {
      comparison_gt <- beta_comparison.dt |>
        gt() |>
        tab_header(
          title = "Pooled Model Path Coefficient Comparison",
          subtitle = sprintf("95%% %s CI overlap test (Cumming & Finch, 2005)",
                            if (BOOTSTRAP) "bootstrap" else "robust")
        ) |>
        cols_label(
          OUTCOME = "Outcome",
          COMPARISON = "Comparison",
          BETA_1 = "β₁",
          CI_1 = "95% CI₁",
          BETA_2 = "β₂",
          CI_2 = "95% CI₂",
          DELTA_BETA = "Δβ",
          CI_OVERLAP = "CIs Overlap?",
          SIGNIFICANT = "Sig. Different?"
        ) |>
        fmt_number(columns = c(BETA_1, BETA_2, DELTA_BETA), decimals = 3) |>
        tab_style(
          style = cell_fill(color = "#E8F5E9"),
          locations = cells_body(columns = SIGNIFICANT, rows = SIGNIFICANT == TRUE)
        ) |>
        tab_style(
          style = cell_fill(color = "#FFEBEE"),
          locations = cells_body(columns = SIGNIFICANT, rows = SIGNIFICANT == FALSE)
        ) |>
        tab_source_note("Non-overlapping 95% CIs indicate significant difference (p < 0.05)") |>
        tab_source_note("β₁ = first model in comparison, β₂ = second model")

      gt::gtsave(comparison_gt, file.path(tables_dir, "pooled_beta_comparison.html"))
      gt::gtsave(comparison_gt, file.path(tables_dir, "pooled_beta_comparison.tex"))

      log_info("Beta comparison table saved")
    }

    # Create forest plot comparing pooled model paths
    if (REDO_PLOTS) {
      # Reshape for plotting
      pooled_plot.dt <- pooled_paths.dt[OUTCOME == "g"]  # Focus on general cognition
      pooled_plot.dt[, MODEL_LABEL := fcase(
        PREDICTOR == "HC", "HC (raw)",
        PREDICTOR == "HVR", "HVR (ratio)",
        PREDICTOR == "HC_RES", "HC (residualized)"
      )]

      p_compare <- ggplot(
        pooled_plot.dt,
        aes(x = reorder(MODEL_LABEL, BETA), y = BETA,
            ymin = CI_LOWER, ymax = CI_UPPER)
      ) +
        geom_pointrange(size = 1, color = "#3C5488") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        coord_flip() +
        labs(
          title = "Brain Measure → General Cognition (g)",
          subtitle = sprintf("Pooled model comparison | 95%% %s CI",
                            if (BOOTSTRAP) "bootstrap" else "robust"),
          x = "Brain Measure",
          y = "Standardized Path Coefficient (β)"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11),
          panel.grid.minor = element_blank()
        )

      save_plot(
        p_compare,
        file.path(plots_dir, "pooled_beta_comparison.png"),
        width = 7,
        height = 4
      )

      log_info("Beta comparison plot saved")
    }
  }

  log_info("Pooled model comparison complete")
} else {
  log_warn("No pooled model paths available for comparison")
}

# ===========================================================================
# Summary
# ===========================================================================
log_section("SEM Analysis Summary")

log_info("Models fitted: %d", length(fits.lst))
log_info("Model types: %s", paste(names(MODELS.lst), collapse = ", "))
log_info("Multi-group by sex: Female (N=%d), Male (N=%d)",
         data.dt[SEX == "Female", uniqueN(EID)],
         data.dt[SEX == "Male", uniqueN(EID)])

if (REDO_TABLES) {
  log_info("Tables saved: %s", tables_dir)
}

if (REDO_PLOTS) {
  log_info("Plots saved: %s", plots_dir)
}

# ===========================================================================
# Export Manuscript Figures
# ===========================================================================
log_section("Exporting manuscript figures")

export_data <- list(
  sem_params = params.lst
)

export_manuscript_figures("sem", export_data)
log_info("Manuscript assets exported to: %s", get_output_path("figures"))

log_script_end("11_sem_analysis.R", success = TRUE)
