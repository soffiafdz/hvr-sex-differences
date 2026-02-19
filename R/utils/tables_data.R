# =============================================================================
# Tables Data Extraction Utilities
# =============================================================================
# Functions for extracting and transforming data for table display.
# Part of the R/utils/ module system.
# =============================================================================

# Source dependencies
if (!exists("format_p")) {
  source(here::here("R/utils/formatting.R"))
}

# =============================================================================
# Global Ordering Constants
# =============================================================================
# Use these throughout to ensure consistent table ordering

ROI_ORDER <- c("HVR", "HC", "LV")
ROI_LABELS <- c(
  "HVR" = "Hippocampal-to-Ventricle Ratio",
  "HC" = "Hippocampus",
  "LV" = "Lateral Ventricles"
)

ADJ_ORDER <- c("NON", "PRP", "STX", "RES")
ADJ_LABELS <- c(
  "NON" = "Unadjusted",
  "PRP" = "Proportions",
  "STX" = "Stereotaxic",
  "RES" = "Residualized"
)

#' Apply ROI ordering to a data.table
#' @param dt A data.table with ROI column
#' @param col Column name containing ROI codes (default: "ROI")
#' @return data.table with ROI as ordered factor
#' @export
apply_roi_order <- function(dt, col = "ROI") {
  if (!is.data.table(dt)) dt <- as.data.table(dt)
  if (col %in% names(dt)) {
    dt[, (col) := factor(get(col), levels = ROI_ORDER)]
    setorderv(dt, col)
  }
  dt
}

#' Apply adjustment method ordering to a data.table
#' @param dt A data.table with ADJ column
#' @param col Column name containing adjustment codes (default: "ADJ")
#' @return data.table with ADJ as ordered factor
#' @export
apply_adj_order <- function(dt, col = "ADJ") {
  if (!is.data.table(dt)) dt <- as.data.table(dt)
  if (col %in% names(dt)) {
    dt[, (col) := factor(get(col), levels = ADJ_ORDER)]
    setorderv(dt, col)
  }
  dt
}

#' Apply both ROI and ADJ ordering
#' @param dt A data.table
#' @param roi_col ROI column name
#' @param adj_col ADJ column name
#' @return Ordered data.table
#' @export
apply_table_order <- function(dt, roi_col = "ROI", adj_col = "ADJ") {
  dt <- apply_roi_order(dt, roi_col)
  dt <- apply_adj_order(dt, adj_col)
  dt
}

# ----- SEM Data Extraction -----

#' Get SEM fit summary from lavaan objects
#' @param sem_fits Named list of lavaan fit objects
#' @return data.table with fit indices
#' @export
get_sem_fit_summary <- function(sem_fits) {
  if (is.null(sem_fits)) {
    return(NULL)
  }

  fit_list <- list()
  for (model_name in names(sem_fits)) {
    fit <- sem_fits[[model_name]]
    if (inherits(fit, "lavaan")) {
      fm <- lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli",
                                       "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr"))
      fit_list[[model_name]] <- data.table(
        Model = model_name,
        Chi_sq = fm["chisq"],
        df = fm["df"],
        p = fm["pvalue"],
        CFI = fm["cfi"],
        TLI = fm["tli"],
        RMSEA = fm["rmsea"],
        RMSEA_CI_L = fm["rmsea.ci.lower"],
        RMSEA_CI_U = fm["rmsea.ci.upper"],
        SRMR = fm["srmr"]
      )
    }
  }

  if (length(fit_list) == 0) return(NULL)
  rbindlist(fit_list)
}

#' Extract structural paths from SEM parameters
#' @param sem_params SEM parameters data.table
#' @return data.table with structural path estimates
#' @export
get_sem_structural_paths <- function(sem_params) {
  if (is.null(sem_params)) {
    return(NULL)
  }

  # Filter to regression paths (op == "~")
  paths <- sem_params[op == "~"]

  if (nrow(paths) == 0) return(NULL)

  # Add sex labels
  paths[, Sex := fifelse(group == 1, "Female", "Male")]

  # Filter to key paths: brain -> cog and covariates -> cog
  cog_outcomes <- c("g", "MEM_s", "PRSP_s")
  key_paths <- paths[lhs %in% cog_outcomes | rhs %in% c("HC", "HVR", "HC_RES")]

  # Use est.std for standardized estimates
  key_paths[, .(
    Outcome = lhs,
    Predictor = rhs,
    Sex,
    Estimate = est.std,
    SE = se,
    CI_Lower = ci.lower,
    CI_Upper = ci.upper,
    Std_Est = est.std,
    P_Value = pvalue
  )]
}

#' Get HVR vs HC comparison results from SEM
#' @param sem_params SEM parameters from load_analysis_data()
#' @return data.table with comparison results
#' @export
get_hvr_hc_sem_comparison <- function(sem_params) {
  if (is.null(sem_params) || !"HC_HVR_COMPARE" %in% names(sem_params)) {
    return(NULL)
  }

  params <- as.data.table(sem_params$HC_HVR_COMPARE)

  # Extract defined parameters (differences and individual paths)
  defined <- params[op == ":="]

  if (nrow(defined) == 0) {
    return(NULL)
  }

  # Extract key results
  result <- defined[, .(
    Parameter = lhs,
    Estimate = est,
    SE = se,
    CI_Lower = ci.lower,
    CI_Upper = ci.upper,
    P_Value = pvalue,
    Estimate_Std = est.std
  )]

  result
}

#' Get individual path coefficients from HVR vs HC comparison
#' @param comparison_dt data.table from get_hvr_hc_sem_comparison()
#' @return data.table with path coefficients
#' @export
get_hvr_hc_paths <- function(comparison_dt) {
  if (is.null(comparison_dt)) return(NULL)

  paths <- comparison_dt[Parameter %like% "^path_"]
  if (nrow(paths) == 0) return(NULL)

  paths[, `:=`(
    Predictor = fifelse(Parameter %like% "_hc_", "HC", "HVR"),
    Outcome = fcase(
      Parameter %like% "_g$", "g",
      Parameter %like% "_mem$", "MEM_s",
      Parameter %like% "_prs$", "PRSP_s"
    )
  )]

  paths
}

#' Extract cognitive domain-specific paths from SEM
#' @param sem_params SEM parameters from load_analysis_data()
#' @param model_name Which model to use (default: "HC_COG")
#' @return data.table with paths to g, MEM_s, and PRSP_s
#' @export
get_cognitive_paths <- function(sem_params, model_name = "HC_COG") {
  if (is.null(sem_params) || !model_name %in% names(sem_params)) {
    return(NULL)
  }

  params <- sem_params[[model_name]]

  # Extract brain -> cognitive factor paths
  # Determine brain variable based on model name
  brain_var <- if (grepl("HVR", model_name)) "HVR" else if (grepl("RES", model_name)) "HC_RES" else "HC"

  paths <- params[rhs == brain_var & lhs %in% c("g", "MEM_s", "PRSP_s"), .(
    Outcome = lhs,
    group,
    beta = est.std,
    ci_lower = ci.lower,
    ci_upper = ci.upper,
    pvalue
  )]

  # Add group labels
  paths[, Sex := fifelse(group == 1, "Female", "Male")]

  paths
}

# ----- Sex Differences Data Extraction -----

#' Get hemisphere comparison data (L vs R)
#' @param sex_diff Sex differences data from load_analysis_data()
#' @return data.table with hemisphere effect sizes
#' @export
get_hemisphere_comparison <- function(sex_diff) {
  if (is.null(sex_diff) || !"OVERALL" %in% names(sex_diff)) {
    return(NULL)
  }

  # Filter to L and R only: HC/LV residualized, HVR self-normalizing
  hc_lv <- sex_diff$OVERALL[SIDE %in% c("L", "R") & ADJ == "RES" & ROI %in% c("HC", "LV"), .(
    ROI, ROI_LABEL, SIDE, ESTIMATE, CI_LOWER, CI_UPPER
  )]
  hvr <- sex_diff$OVERALL[SIDE %in% c("L", "R") & ADJ == "NON" & ROI == "HVR", .(
    ROI, ROI_LABEL, SIDE, ESTIMATE, CI_LOWER, CI_UPPER
  )]
  rbind(hc_lv, hvr)
}

#' Get age-stratified effect sizes
#' @param sex_diff Sex differences data from load_analysis_data()
#' @return data.table with effect sizes by age bin
#' @export
get_age_stratified <- function(sex_diff) {
  if (is.null(sex_diff) || !"AGE_STRATIFIED" %in% names(sex_diff)) {
    return(NULL)
  }

  sex_diff$AGE_STRATIFIED[SIDE == "LR" & ADJ == "RES"]
}

#' Get sensitivity analysis comparison
#' @param sex_diff Sex differences data from load_analysis_data()
#' @return data.table with primary vs sensitivity effect sizes
#' @export
get_sensitivity_comparison <- function(sex_diff) {
  if (is.null(sex_diff) || !"SENS_COMPARISON" %in% names(sex_diff)) {
    return(NULL)
  }

  # SENS_COMPARISON may not have SIDE column (bilateral only)
  dt <- sex_diff$SENS_COMPARISON
  if ("SIDE" %in% names(dt)) {
    return(dt[SIDE == "LR"])
  }
  return(dt)
}

#' Get HVR-ICV correlation validation data
#' @param sex_diff Sex differences data from load_analysis_data()
#' @return data.table with HVR, HC, LV correlations with ICV
#' @export
get_hvr_icv_validation <- function(sex_diff) {
  if (is.null(sex_diff) || !"HVR_ICV_VALIDATION" %in% names(sex_diff)) {
    return(NULL)
  }
  sex_diff$HVR_ICV_VALIDATION
}
