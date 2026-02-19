# ===========================================================================
# GAMLSS Utility Functions
# ===========================================================================
# Helper functions for GAMLSS modeling, prediction, and validation
# ===========================================================================

library(gamlss)
library(gamlss.dist)
library(data.table)

# ---------------------------------------------------------------------------
# Family Information
# ---------------------------------------------------------------------------

#' Get GAMLSS family display names
#' @return Named character vector
get_gamlss_family_names <- function() {
 c(

    BE = "Beta",
    BCCG = "Box-Cox Cole & Green",
    L_NO = "Gaussian (logit-transformed)",
    NO = "Gaussian"
  )
}

#' Get parameter names for each GAMLSS family
#' @param family Family code (NO, L_NO, BE, BCCG)
#' @return Character vector of parameter names
get_family_params <- function(family) {
  switch(family,
    NO = c("mu", "sigma"),
    L_NO = c("mu", "sigma"),
    BE = c("mu", "sigma"),
    BCCG = c("mu", "sigma", "nu"),
    stop("Unknown family: ", family, call. = FALSE)
  )
}

# ---------------------------------------------------------------------------
# Prediction Helpers
# ---------------------------------------------------------------------------

#' Predict GAMLSS parameters for new data
#' @param fit GAMLSS fit object
#' @param newdata Data for prediction
#' @param family Family code
#' @return List of predicted parameters
predict_gamlss_params <- function(fit, newdata, family) {
  params <- list(
    mu = predict(fit, newdata = newdata, type = "response"),
    sigma = predict(fit, "sigma", newdata = newdata, type = "response")
  )
  if (family == "BCCG") {
    params$nu <- predict(fit, "nu", newdata = newdata, type = "response")
  }
  params
}

#' Compute centiles (CDF) for observed values
#' @param values Observed values
#' @param fit GAMLSS fit object
#' @param newdata Data for prediction (must match values length)
#' @param family Family code (NO, L_NO, BE, BCCG)
#' @param scale Multiply result by this (default 100 for percentiles)
#' @return Numeric vector of centiles
gamlss_compute_centile <- function(
    values, fit, newdata, family, scale = 100) {
  params <- predict_gamlss_params(fit, newdata, family)

  centiles <- switch(family,
    NO = pNO(values, mu = params$mu, sigma = params$sigma),
    L_NO = {
      eps <- 1e-6
      val_logit <- values |> pmin(1 - eps) |> pmax(eps) |> qlogis()
      pNO(val_logit, mu = params$mu, sigma = params$sigma)
    },
    BE = pBE(values, mu = params$mu, sigma = params$sigma),
    BCCG = pBCCG(
      values,
      mu = params$mu, sigma = params$sigma, nu = params$nu
    ),
    stop("Unknown family: ", family, call. = FALSE)
  )

  centiles * scale
}

#' Compute quantiles for given centile probabilities
#' @param probs Centile probabilities (0-1)
#' @param fit GAMLSS fit object
#' @param newdata Data for prediction
#' @param family Family code (NO, L_NO, BE, BCCG)
#' @return Matrix with rows = newdata rows, cols = probs
gamlss_compute_quantile <- function(probs, fit, newdata, family) {
  params <- predict_gamlss_params(fit, newdata, family)

  # Compute quantiles for each probability
  result <- sapply(probs, function(p) {
    switch(family,
      NO = qNO(p, mu = params$mu, sigma = params$sigma),
      L_NO = plogis(qNO(p, mu = params$mu, sigma = params$sigma)),
      BE = qBE(p, mu = params$mu, sigma = params$sigma),
      BCCG = qBCCG(
        p,
        mu = params$mu, sigma = params$sigma, nu = params$nu
      ),
      stop("Unknown family: ", family, call. = FALSE)
    )
  })

  # Name columns by centile
  if (is.matrix(result)) {
    colnames(result) <- paste0("p", probs * 100)
  } else {
    names(result) <- paste0("p", probs * 100)
  }

  result
}

# ---------------------------------------------------------------------------
# Warning Handling
# ---------------------------------------------------------------------------
#' Execute expression with GAMLSS refit warning capture
#' @param expr Expression to evaluate
#' @param context Optional context string for logging
#' @return List with result and refit_warning flag
with_gamlss_warnings <- function(expr, context = NULL) {
  refit_warning <- FALSE

  result <- withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("discrepancy.*re-fit", w$message)) {
        refit_warning <<- TRUE
        if (!is.null(context)) {
          log_debug("%s: GAMLSS refit warning", context)
        }
      }
      invokeRestart("muffleWarning")
    }
  )

  list(result = result, refit_warning = refit_warning)
}

# ---------------------------------------------------------------------------
# Centile Prediction for Data Tables
# ---------------------------------------------------------------------------

#' Add centile predictions to data.table
#' @param dt Data table with VAL column
#' @param fit GAMLSS fit object
#' @param family Family code
#' @param pred_cols Columns to use for prediction
#' @param context Context string for warnings
#' @return List with modified dt and refit_warning
dt_add_centiles <- function(
    dt, fit, family, pred_cols = c("AGE", "EDUC_num", "ICC"),
    context = NULL) {

  # Prepare newdata
  newdata <- dt[, ..pred_cols]

  # Handle L_NO transformation
  if (family == "L_NO") {
    eps <- 1e-6
    dt[, VAL_logit := VAL |> pmin(1 - eps) |> pmax(eps) |> qlogis()]
  }

  # Compute centiles with warning capture
  if (family %in% c("BE", "BCCG")) {
    res <- with_gamlss_warnings({
      values <- if (family == "L_NO") dt$VAL_logit else dt$VAL
      gamlss_compute_centile(values, fit, newdata, family)
    }, context)
    dt[, CENT_pred := res$result]
    refit_warning <- res$refit_warning
  } else {
    values <- if (family == "L_NO") dt$VAL_logit else dt$VAL
    dt[, CENT_pred := gamlss_compute_centile(values, fit, newdata, family)]
    refit_warning <- FALSE
  }

  list(dt = dt, refit_warning = refit_warning)
}

#' Add quantile columns to data.table
#' @param dt Data table (will be modified in place)
#' @param fit GAMLSS fit object
#' @param family Family code
#' @param centiles Centile probabilities
#' @param newdata Data for prediction (if NULL, uses dt)
#' @param pred_cols Columns to use for prediction
#' @param context Context string for warnings
#' @return List with refit_warning flag
dt_add_quantiles <- function(
    dt, fit, family, centiles,
    newdata = NULL, pred_cols = c("AGE", "EDUC_num", "ICC"),
    context = NULL) {

  if (is.null(newdata)) {
    newdata <- dt[, ..pred_cols]
  }

  col_names <- paste0("p", centiles * 100)

  # Compute with warning capture for BE/BCCG
  if (family %in% c("BE", "BCCG")) {
    res <- with_gamlss_warnings({
      gamlss_compute_quantile(centiles, fit, newdata, family)
    }, context)
    quantiles <- res$result
    refit_warning <- res$refit_warning
  } else {
    quantiles <- gamlss_compute_quantile(centiles, fit, newdata, family)
    refit_warning <- FALSE
  }

  # Add columns to dt
  for (i in seq_along(col_names)) {
    col <- col_names[i]
    if (is.matrix(quantiles)) {
      dt[, (col) := quantiles[, i]]
    } else {
      dt[, (col) := quantiles[i]]
    }
  }

  list(refit_warning = refit_warning)
}

# ---------------------------------------------------------------------------
# Model Comparison Helpers
# ---------------------------------------------------------------------------

#' Save model fit to comparison list
#' @param mod_lst List to store results
#' @param mod_name Model name
#' @param mod_fit GAMLSS fit object
#' @return Updated mod_lst
save_gamlss_fit <- function(mod_lst, mod_name, mod_fit) {
  comp_row <- data.table(
    MOD = mod_name,
    AIC = mod_fit$aic,
    BIC = BIC(mod_fit),
    DF = mod_fit$df.fit,
    GDEV = deviance(mod_fit)
  )

  if (is.null(mod_lst$COMP)) {
    mod_lst$COMP <- comp_row
  } else {
    mod_lst$COMP <- rbind(mod_lst$COMP, comp_row)
  }
  setkey(mod_lst$COMP, MOD)
  mod_lst[[mod_name]] <- mod_fit
  mod_lst
}

#' Pick best model using AIC with BIC check
#' @param comp_dt Comparison data.table with MOD, AIC, BIC columns
#' @param bic_crit BIC difference threshold (default 10)
#' @return List with MOD and WARN flag
pick_best_gamlss <- function(comp_dt, bic_crit = 10) {
  # Heuristics:
  # 1) Use AIC to select
  # 2) If BIC strongly disagrees (diff > bic_crit), warn
  # Reference: Burnham & Anderson (2002)
  best_mod <- comp_dt[which.min(AIC), MOD]
  bic_diff <- comp_dt[
    , .SD[best_mod, on = "MOD", BIC] - min(BIC)
  ]

  if (bic_diff > bic_crit) {
    return(list(MOD = comp_dt[which.min(BIC), MOD], WARN = TRUE))
  }
  list(MOD = best_mod, WARN = FALSE)
}

# ---------------------------------------------------------------------------
# Validation Helpers
# ---------------------------------------------------------------------------

#' Compute centile calibration metrics
#' @param predicted_centiles Vector of predicted centiles (0-100)
#' @param thresholds Threshold values to check calibration at
#' @return Named list with calibration metrics
compute_centile_calibration <- function(
    predicted_centiles,
    thresholds = c(5, 10, 25, 50, 75, 90, 95)) {

  observed_props <- sapply(thresholds, function(t) {
    mean(predicted_centiles < t) * 100
  })

  mae <- mean(abs(observed_props - thresholds))

  interpretation <- if (mae < 5) {
    "Excellent"
  } else if (mae < 10) {
    "Good"
  } else if (mae < 15) {
    "Acceptable"
  } else {
    "Poor"
  }

  list(
    CALIBRATION = setNames(observed_props, paste0("obs_", thresholds)),
    MAE = mae,
    INTERPRETATION = interpretation
  )
}

#' Compute longitudinal stability metrics
#' @param centiles_t1 Centiles at time 1
#' @param centiles_t2 Centiles at time 2
#' @param crossing_thresholds Thresholds for clinical crossing
#' @return Named list with stability metrics
compute_stability_metrics <- function(
    centiles_t1, centiles_t2,
    crossing_thresholds = c(5, 95)) {

  valid_idx <- !is.na(centiles_t1) & !is.na(centiles_t2)
  c1 <- centiles_t1[valid_idx]
  c2 <- centiles_t2[valid_idx]

  if (length(c1) < 10) {
    return(list(
      CORR = NA_real_,
      MEAN_CHANGE = NA_real_,
      SD_CHANGE = NA_real_,
      CROSS_PROP = NA_real_,
      N = length(c1)
    ))
  }

  change <- c2 - c1

  # Count clinical threshold crossings
  crossings <- sapply(crossing_thresholds, function(t) {
    sum((c1 < t & c2 >= t) | (c1 >= t & c2 < t))
  })

  list(
    CORR = cor(c1, c2),
    MEAN_CHANGE = mean(change),
    SD_CHANGE = sd(change),
    CROSS_PROP = sum(crossings) / length(c1),
    N = length(c1)
  )
}

# ---------------------------------------------------------------------------
# Z-Score Calculation
# ---------------------------------------------------------------------------

#' Compute Z-scores from GAMLSS model
#' @param values Observed values
#' @param fit GAMLSS fit object
#' @param newdata Data for prediction
#' @param family Family code
#' @return Numeric vector of Z-scores
gamlss_compute_zscore <- function(values, fit, newdata, family) {
  # Get centiles (0-1 scale)
  centiles <- gamlss_compute_centile(values, fit, newdata, family, scale = 1)
  # Convert to Z-scores

  qnorm(centiles)
}
