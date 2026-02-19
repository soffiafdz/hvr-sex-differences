# =============================================================================
# Formatting Utilities
# =============================================================================
# Common formatting functions for numeric values, p-values, effect sizes, etc.
# These are small helpers used across tables, plots, and inline reporting.
# =============================================================================

#' Format Cohen's d effect size with confidence interval
#' @param d Effect size estimate
#' @param ci_lo Lower confidence interval
#' @param ci_hi Upper confidence interval
#' @return Formatted string: "d = X.XX [X.XX, X.XX]"
format_d <- function(d, ci_lo = NULL, ci_hi = NULL) {
  if (is.null(ci_lo) || is.null(ci_hi)) {
    sprintf("d = %.2f", d)
  } else {
    sprintf("d = %.2f [%.2f, %.2f]", d, ci_lo, ci_hi)
  }
}

#' Format large numbers with comma separators
#' @param x Numeric value
#' @return Formatted string
format_n <- function(x) {
  format(x, big.mark = ",", scientific = FALSE)
}

#' Format percentage
#' @param x Numeric value (proportion)
#' @param digits Number of decimal places
#' @return Formatted string with % symbol
format_pct <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f%%"), x * 100)
}

#' Format p-value with appropriate precision
#' @param p P-value
#' @return Formatted string (scientific notation for very small values)
format_p <- function(p) {
  if (p < 0.001) {
    # Use scientific notation for very small p-values
    sprintf("p = %.2e", p)
  } else if (p < 0.01) {
    sprintf("p = %.3f", p)
  } else {
    sprintf("p = %.2f", p)
  }
}

#' Format standardized beta coefficient
#' @param beta Standardized coefficient
#' @param se Standard error (optional)
#' @param p P-value (optional)
#' @return Formatted string
format_beta <- function(beta, se = NULL, p = NULL) {
  result <- sprintf("\u03B2 = %.3f", beta)
  if (!is.null(se)) {
    result <- paste0(result, sprintf(" (SE = %.3f)", se))
  }
  if (!is.null(p)) {
    result <- paste0(result, ", ", format_p(p))
  }
  result
}

#' Format Williams t-test result
#' @param ratio HVR/HC ratio
#' @param p P-value
#' @return Formatted string
format_williams <- function(ratio, p) {
  if (p < 2.2e-16) {
    sprintf("%.2fx (p < 10\u207b\u00b9\u2076)", ratio)
  } else {
    sprintf("%.2fx (%s)", ratio, format_p(p))
  }
}

#' Format slope ratio for age interactions
#' @param slope_ratio Male/Female slope ratio
#' @param p P-value
#' @return Formatted string
format_slope_ratio <- function(slope_ratio, p) {
  if (p < 2.2e-16) {
    sprintf("%.1fx faster (p < 10\u207b\u00b9\u2076)", abs(slope_ratio))
  } else {
    sprintf("%.1fx faster (%s)", abs(slope_ratio), format_p(p))
  }
}

#' Safely format a numeric value for display
#' Handles NA, NaN, Inf, and -Inf gracefully
#' @param x Numeric value to format
#' @param digits Number of decimal places (default 3)
#' @param na_string String to display for NA/NaN/Inf values (default "—")
#' @return Formatted string
safe_round <- function(x, digits = 3, na_string = "\u2014") {
  if (is.null(x) || length(x) == 0) return(na_string)
  if (is.na(x) || is.nan(x) || is.infinite(x)) return(na_string)
  round(x, digits)
}

#' Safely format a p-value for display
#' @param p P-value to format
#' @param na_string String to display for NA values (default "—")
#' @return Formatted string (scientific notation for very small values)
safe_pvalue <- function(p, na_string = "\u2014") {
  if (is.null(p) || length(p) == 0) return(na_string)
  if (is.na(p) || is.nan(p) || is.infinite(p)) return(na_string)
  if (p == 0) return("< 2.2e-16")
  if (p < 0.001) return(sprintf("%.2e", p))
  sprintf("%.3f", p)
}

#' Safely format a confidence interval for display
#' @param lower Lower CI bound
#' @param upper Upper CI bound
#' @param digits Number of decimal places (default 3)
#' @param na_string String to display for NA values (default "—")
#' @return Formatted CI string "[lower, upper]" or na_string
safe_ci <- function(lower, upper, digits = 3, na_string = "\u2014") {
  if (is.null(lower) || is.null(upper)) return(na_string)
  if (is.na(lower) || is.na(upper) || is.infinite(lower) || is.infinite(upper)) {
    return(na_string)
  }
  sprintf("[%.3f, %.3f]", round(lower, digits), round(upper, digits))
}
