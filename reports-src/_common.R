# =============================================================================
# Common Setup for Quarto Documents
# =============================================================================
# This file provides standardized setup for all Quarto documents in reports-src/
# Source this file at the beginning of each document's setup chunk.
#
# Usage in Quarto document:
#   ```{r setup}
#   library(here)
#   here::i_am("reports-src/[document_name].qmd")
#   source(here("reports-src/_common.R"))
#   ```
#
# Note: All plotting, table, and data loading functions have been consolidated
# into R/utils/ modules. This file is a thin loader that sources those modules.
# =============================================================================

# ----- Core Libraries -----
library(data.table)
library(gt)
library(yaml)
library(ggplot2)
library(ggtext)
library(patchwork)
library(gamlss)  # For GAMLSS coefficient extraction

# ----- Source Pipeline Utilities -----
# These modules contain all analysis functions:
# - config.R: Configuration and path resolution
# - data_io.R: Data loading functions (load_analysis_data, load_normative_tables, etc.)
# - formatting.R: Display formatting helpers (format_d, format_p, safe_round, etc.)
# - plotting.R: Figure functions (fig_effect_sizes_by_adjustment, fig_hvr_centiles, etc.)
# - tables.R: GT table functions (format_sem_fit_gt, format_centile_gt, etc.)
# - validation.R: Data validation utilities
# - logging.R: Logging utilities

source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/validation.R"))
source(here("R/utils/formatting.R"))
source(here("R/utils/data_io.R"))
# Table utilities (split into focused modules)
source(here("R/utils/tables_core.R"))
source(here("R/utils/tables_gt.R"))
source(here("R/utils/tables_data.R"))
source(here("R/utils/tables_normative.R"))

# Plotting utilities (split into focused modules)
source(here("R/utils/plotting_core.R"))
source(here("R/utils/plotting_figures.R"))
source(here("R/utils/plotting_pipeline.R"))

# ----- Load Configuration -----
config <- load_config()

# ----- Global Options -----
options(scipen = 999, digits = 3)

# ----- LaTeX Table Rendering Override -----
# Override gt's knit_print to produce clean LaTeX output when rendering to PDF.
# Post-processing fixes: full-width tables, HTML tags, Unicode symbols, alignment.
# For HTML output, gt's default rendering is used unchanged.
registerS3method("knit_print", "gt_tbl", function(x, ...) {
  if (knitr::is_latex_output()) {
    latex_code <- as.character(gt::as_latex(x))

    # 0. Strip gt's outer \begin{table}...\end{table} float wrapper.
    #    Quarto provides its own float via tbl-cap; keeping gt's wrapper creates
    #    nested floats which can cause LaTeX to silently discard content.
    #    Strategy: keep only from the fontsize+tabular onward, drop table/caption.
    strip_pattern <- "(?s)^.*?(\\\\fontsize\\{[^}]+\\}\\{[^}]+\\}\\\\selectfont\\s*\n?\\\\begin\\{(?:tabular|longtable))"
    if (grepl(strip_pattern, latex_code, perl = TRUE)) {
      latex_code <- sub(strip_pattern, "\\1", latex_code, perl = TRUE)
    } else {
      # Fallback: strip to tabular/longtable start directly
      latex_code <- sub("(?s)^.*?(\\\\begin\\{(?:tabular|longtable))", "\\1", latex_code, perl = TRUE)
    }
    latex_code <- sub("\\\\end\\{table\\}\\s*$", "", latex_code)

    # 1. Fix table width: tabular*{\linewidth}{@{\extracolsep{\fill}}...} → tabular{...}
    latex_code <- gsub(
      "\\\\begin\\{tabular\\*\\}\\{\\\\linewidth\\}\\{@\\{\\\\extracolsep\\{\\\\fill\\}\\}",
      "\\\\begin{tabular}{",
      latex_code
    )
    latex_code <- gsub("\\\\end\\{tabular\\*\\}", "\\\\end{tabular}", latex_code)

    # 2. Fix HTML tags that fmt_markdown() produces (should be LaTeX commands)
    latex_code <- gsub("<em>([^<]+)</em>", "\\\\textit{\\1}", latex_code)
    latex_code <- gsub("<strong>([^<]+)</strong>", "\\\\textbf{\\1}", latex_code)
    latex_code <- gsub("<sub>([^<]+)</sub>", "\\\\textsubscript{\\1}", latex_code)
    latex_code <- gsub("<sup>([^<]+)</sup>", "\\\\textsuperscript{\\1}", latex_code)

    # 3. Fix Unicode symbols not in default LaTeX font
    latex_code <- gsub("\u0394", "$\\\\Delta$", latex_code)  # Δ → $\Delta$
    latex_code <- gsub("\u03B2", "$\\\\beta$", latex_code)   # β → $\beta$
    latex_code <- gsub("\u03BB", "$\\\\lambda$", latex_code)  # λ → $\lambda$
    latex_code <- gsub("\u00B2", "\\\\textsuperscript{2}", latex_code)  # ² → ^2
    latex_code <- gsub("\u2264", "$\\\\leq$", latex_code)    # ≤
    latex_code <- gsub("\u2265", "$\\\\geq$", latex_code)    # ≥
    latex_code <- gsub("\u2192", "$\\\\rightarrow$", latex_code)  # → → $\rightarrow$
    latex_code <- gsub("\u209b", "\\\\textsubscript{s}", latex_code)  # ₛ → \textsubscript{s}
    latex_code <- gsub("\u2014", "---", latex_code)          # — → ---

    # 4. Fix footnote/source note minipage width (constrain to table width)
    #    gt wraps footnotes in \minipage{\linewidth} which spans full page.
    #    Replace with narrower width and add raggedright to prevent overflow.
    latex_code <- gsub(
      "\\\\begin\\{minipage\\}\\{\\\\linewidth\\}",
      "\\\\begin{minipage}[t]{0.65\\\\textwidth}\\\\raggedright\\\\small",
      latex_code
    )

    # 5. Enforce consistent column alignment
    #    Standard format: first column left-aligned (labels), rest right-aligned (data)
    #    Also handle longtable environments
    fix_alignment <- function(code) {
      # Pattern matches \begin{tabular}{...} or \begin{longtable}{...}
      pattern <- "(\\\\begin\\{(?:tabular|longtable)\\}\\{)([lcr]+)(\\})"
      m <- gregexpr(pattern, code, perl = TRUE)
      if (m[[1]][1] == -1) return(code)

      matches <- regmatches(code, m)[[1]]
      for (match in matches) {
        # Extract column spec
        col_spec <- gsub(".*\\{([lcr]+)\\}$", "\\1", match)
        n_cols <- nchar(col_spec)
        if (n_cols > 1) {
          # First column left, rest right
          new_spec <- paste0("l", paste(rep("r", n_cols - 1), collapse = ""))
          new_match <- gsub("\\{[lcr]+\\}$", paste0("{", new_spec, "}"), match)
          code <- sub(match, new_match, code, fixed = TRUE)
        }
      }
      code
    }
    latex_code <- fix_alignment(latex_code)

    # 6. Add consistent table spacing
    latex_code <- gsub(
      "\\\\begin\\{tabular\\}",
      "\\\\renewcommand{\\\\arraystretch}{1.15}\\\\begin{tabular}",
      latex_code
    )

    # 7. Center each tabular inside the float
    latex_code <- gsub(
      "(\\\\renewcommand\\{\\\\arraystretch\\}\\{1\\.15\\})(\\\\begin\\{tabular\\})",
      "\\\\begin{center}\\1\\2",
      latex_code
    )
    latex_code <- gsub(
      "\\\\end\\{tabular\\}",
      "\\\\end{tabular}\\\\end{center}",
      latex_code
    )

    knitr::asis_output(latex_code)
  } else {
    gt:::knit_print.gt_tbl(x, ...)
  }
}, envir = asNamespace("knitr"))

# ----- Quarto-Specific Helpers -----

#' Generate inline reference to a figure
#' @param fig_num Figure number (e.g., 1, 2, "S1")
#' @return Formatted figure reference string
fig_ref <- function(fig_num) {
  sprintf("Figure %s", fig_num)
}

#' Generate inline reference to a table
#' @param tab_num Table number (e.g., 1, 2, "S1")
#' @return Formatted table reference string
tab_ref <- function(tab_num) {
  sprintf("Table %s", tab_num)
}

#' Generate inline reference to supplementary material
#' @param type Type of supplement ("Figure", "Table", "Text")
#' @param num Supplement number
#' @return Formatted supplement reference string
supp_ref <- function(type = "Figure", num) {
  sprintf("%s S%s", type, num)
}

# ----- Supplementary Reference System -----
# Order must match the supplementary file. Update vectors here if items are reordered.
supp_tables <- c(
  "tbl-cog-battery",
  "tbl-sem-covariates-fields",
  "tbl-temporal-stability",
  "tbl-hvr-centiles-female",
  "tbl-hvr-centiles-male",
  "tbl-hc-centiles-female",
  "tbl-hc-centiles-male",
  "tbl-lv-centiles-female",
  "tbl-lv-centiles-male",
  "tbl-gamlss-hvr-coef",
  "tbl-gamlss-hc-coef",
  "tbl-measurement-invariance",
  "tbl-sem-paths-by-sex",
  "tbl-sem-fit",
  "tbl-sem-loadings",
  "tbl-sem-covariate-effects",
  "tbl-site-adjusted",
  "tbl-hemisphere",
  "tbl-sensitivity"
)

supp_figures <- c(
  "fig-distributions",
  "fig-gamlss-calibration",
  "fig-hemisphere",
  "fig-sensitivity"
)

#' Generate reference to a supplementary item by its original chunk label
#' @param label Original chunk label (e.g., "tbl-cog-battery", "fig-distributions")
#' @return Formatted reference (e.g., "Table S1", "Figure S3")
S <- function(label) {
  if (label %in% supp_tables) {
    sprintf("Table S%d", which(supp_tables == label))
  } else if (label %in% supp_figures) {
    sprintf("Figure S%d", which(supp_figures == label))
  } else {
    warning(sprintf("Unknown supplementary label: %s", label))
    sprintf("??%s", label)
  }
}

#' Generate a range reference for consecutive supplementary items
#' @param label_first First chunk label
#' @param label_last Last chunk label
#' @return Formatted range (e.g., "Tables S4--S9")
S_range <- function(label_first, label_last) {
  if (label_first %in% supp_tables && label_last %in% supp_tables) {
    n1 <- which(supp_tables == label_first)
    n2 <- which(supp_tables == label_last)
    sprintf("Tables S%d--S%d", n1, n2)
  } else if (label_first %in% supp_figures && label_last %in% supp_figures) {
    n1 <- which(supp_figures == label_first)
    n2 <- which(supp_figures == label_last)
    sprintf("Figures S%d--S%d", n1, n2)
  } else {
    warning("Labels must be of the same type (both tables or both figures)")
    "??range"
  }
}

message("=== _common.R loaded (formatting + gt override only) ===")
