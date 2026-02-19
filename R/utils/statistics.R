# =============================================================================
# Statistical Utility Functions
# =============================================================================
# Helper functions for statistical tests and effect size calculations
# =============================================================================

#' Williams t-test for comparing dependent correlations
#'
#' Tests whether two correlations sharing a common variable differ significantly.
#' Used to compare r(X,Y) vs r(X,Z) where Y and Z are correlated.
#'
#' @param r_xy Correlation between X and Y
#' @param r_xz Correlation between X and Z
#' @param r_yz Correlation between Y and Z
#' @param n Sample size
#' @return List with t statistic, p-value, and correlation difference
#' @references Williams, E. J. (1959). The comparison of regression variables.
#'   Journal of the Royal Statistical Society, Series B, 21, 396-399.
#' @export
williams_t_test <- function(r_xy, r_xz, r_yz, n) {
  # Validate inputs
  if (!is.numeric(r_xy) || !is.numeric(r_xz) || !is.numeric(r_yz)) {
    stop("All correlations must be numeric")
  }
  if (any(abs(c(r_xy, r_xz, r_yz)) > 1)) {
    stop("Correlations must be between -1 and 1")
  }
  if (!is.numeric(n) || n < 4) {
    stop("Sample size must be numeric and >= 4")
  }


  # Determinant of correlation matrix

det_R <- 1 - r_xy^2 - r_xz^2 - r_yz^2 + 2 * r_xy * r_xz * r_yz

  # Williams t-statistic
  t_stat <- (r_xy - r_xz) * sqrt(
    (n - 1) * (1 + r_yz) /
    (2 * det_R + (r_xy - r_xz)^2 * (1 - r_yz)^3 / (4 * (n - 1)))
  )

  # Two-tailed p-value
  p_val <- 2 * pt(-abs(t_stat), df = n - 3)

  list(
    t = t_stat,
    df = n - 3,
    p = p_val,
    r_diff = r_xy - r_xz
  )
}

#' Compute partial correlation
#'
#' Calculates the correlation between two variables after controlling for
#' one or more covariates.
#'
#' @param x Numeric vector for first variable
#' @param y Numeric vector for second variable
#' @param z Data frame or matrix of covariates to control for
#' @return Partial correlation coefficient
#' @export
partial_cor <- function(x, y, z) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length")
  }

  # Residualize x and y on z
  if (is.vector(z)) {
    z <- data.frame(z = z)
  }

  df <- data.frame(x = x, y = y, z)
  df <- df[complete.cases(df), ]

  resid_x <- residuals(lm(x ~ ., data = df[, -1, drop = FALSE]))
  resid_y <- residuals(lm(y ~ ., data = df[, -2, drop = FALSE]))

 cor(resid_x, resid_y)
}

#' Extract interaction term from linear model
#'
#' Extracts coefficient, standard error, t-value and p-value for an
#' interaction term in a linear model.
#'
#' @param model lm object
#' @param interaction_pattern Regex pattern to match interaction term
#' @return Data frame with interaction statistics
#' @export
extract_interaction <- function(model, interaction_pattern = ":") {
  summ <- summary(model)$coefficients
  int_row <- grep(interaction_pattern, rownames(summ))

  if (length(int_row) == 0) {
    warning("No interaction term found matching pattern: ", interaction_pattern
)
    return(NULL)
  }

  data.frame(
    term = rownames(summ)[int_row],
    estimate = summ[int_row, "Estimate"],
    std_error = summ[int_row, "Std. Error"],
    t_value = summ[int_row, "t value"],
    p_value = summ[int_row, "Pr(>|t|)"]
  )
}
