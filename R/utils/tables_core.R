# =============================================================================
# Table Core Utilities
# =============================================================================
# Core table functions: save, helpers, labels, and formatting utilities.
# Part of the R/utils/ module system.
# =============================================================================

# ----- Helper Functions -----

#' Wrap LaTeX table fragment in standalone document
#' @param tex_fragment_path Path to .tex file with table code
#' @param output_path Path for standalone .tex document (optional)
#' @param title Document title (optional, if NULL table title is used)
#' @return Path to standalone document
wrap_latex_table <- function(
    tex_fragment_path, output_path = NULL, title = NULL) {
  if (is.null(output_path)) {
    output_path <- sub("\\.tex$", "_standalone.tex", tex_fragment_path)
  }

  # Read the table fragment
  table_code <- readLines(tex_fragment_path)

  # Create standalone document - minimal header, table has its own title
  doc_header <- c(
    "\\documentclass[11pt]{article}",
    "\\usepackage{booktabs}",
    "\\usepackage{longtable}",
    "\\usepackage{geometry}",
    "\\geometry{a4paper, margin=0.75in}",
    "\\usepackage{caption}",
    "\\begin{document}",
    "\\pagestyle{empty}" # No page numbers
  )

  # Only add section title if explicitly provided
  # Otherwise the table's own title/subtitle will be used
  if (!is.null(title)) {
    doc_header <- c(doc_header, "", paste0("\\section*{", title, "}"), "")
  }

  doc_footer <- c("", "\\end{document}")

  # Combine
  standalone_doc <- c(doc_header, "", table_code, doc_footer)

  # Write
  writeLines(standalone_doc, output_path)

  return(output_path)
}

#' Bin ages into intervals
#' @param ages Vector of ages
#' @param bin_width Width of age bins in years
#' @return Character vector of age bin labels
bin_ages <- function(ages, bin_width = 5) {
  breaks <- seq(
    floor(min(ages) / bin_width) * bin_width,
    ceiling(max(ages) / bin_width) * bin_width,
    by = bin_width
  )
  bins <- cut(ages, breaks = breaks, include.lowest = TRUE, right = FALSE)
  # Format as "45-50", "50-55", etc.
  levels(bins) <- sapply(seq_along(breaks[-length(breaks)]), \(i) {
    sprintf("%d-%d", breaks[i], breaks[i + 1])
  })
  return(bins)
}

# ----- Table Styling -----

#' Standard gt table theme
#' @param gt_table gt table object
#' @return Styled gt table
theme_publication_table <- function(gt_table) {
  gt_table |>
    gt::tab_options(
      footnotes.multiline = FALSE,
      latex.tbl.pos = "h"
    ) |>
    gt::cols_align("center", columns = gt::everything())
}

#' Apply significance styling to table
#' @param gt_table gt table object
#' @param sig_col Column indicating significance
#' @param style_cols Columns to style
#' @return Styled gt table
apply_significance_style <- function(
    gt_table, sig_col = "SIGN", style_cols = NULL) {
  if (is.null(style_cols)) {
    style_cols <- gt::everything()
  }

  gt_table |>
    gt::tab_style(
      style = gt::cell_text(style = "italic"),
      locations = gt::cells_body(
        columns = style_cols,
        rows = .data[[sig_col]] == FALSE
      )
    ) |>
    gt::cols_hide(columns = gt::contains(sig_col))
}

#' Format p-values consistently
#' @param pval P-value
#' @param threshold Threshold for "<" notation
#' @return Formatted p-value string
format_pval <- function(pval, threshold = 0.001, stars = FALSE) {
  if (stars) {
    # Return significance stars
    ifelse(pval < 0.001, "***",
      ifelse(pval < 0.01, "**",
        ifelse(pval < 0.05, "*", "")
      )
    )
  } else {
    # Return formatted p-value
    ifelse(pval < threshold,
      sprintf("<%s", threshold),
      sprintf("%.3f", pval)
    )
  }
}

#' Format effect size with CI
#' @param estimate Point estimate
#' @param ci_lower Lower CI bound
#' @param ci_upper Upper CI bound
#' @param digits Number of decimal places
#' @return Formatted string
format_effect_ci <- function(estimate, ci_lower, ci_upper, digits = 2) {
  fmt <- sprintf("%%.%df [%%.%df, %%.%df]", digits, digits, digits)
  sprintf(fmt, estimate, ci_lower, ci_upper)
}

#' Format mean and SD
#' @param values Vector of values
#' @param is_icc Whether this is ICC (use comma separator)
#' @param use_cog_format Whether to use cognitive format (3 decimals)
#' @return Formatted string
format_mean_sd <- function(values, is_icc = FALSE, use_cog_format = FALSE) {
  if (use_cog_format) {
    sprintf("%.3f (%.2f)", mean(values, na.rm = TRUE), sd(values, na.rm = TRUE))
  } else if (is_icc) {
    sprintf(
      "%s (%s)",
      format(mean(values, na.rm = TRUE), big.mark = ",", digits = 4),
      format(sd(values, na.rm = TRUE), big.mark = ",", digits = 4)
    )
  } else {
    sprintf("%.1f (%.1f)", mean(values, na.rm = TRUE), sd(values, na.rm = TRUE))
  }
}

#' Format standardized coefficient
#' @param value Coefficient value
#' @param threshold Threshold for scientific notation
#' @return Formatted string
format_beta_std <- function(value, threshold = 0.01) {
  ifelse(abs(value) < threshold,
    sprintf("%.0e", value),
    sprintf("%.2f", value)
  )
}

#' Format confidence interval
#' @param lower Lower bound
#' @param upper Upper bound
#' @param threshold Threshold for scientific notation
#' @return Formatted string
format_ci <- function(lower, upper, threshold = 0.01) {
  lower_fmt <- ifelse(abs(lower) < threshold, "%.0e", "%.2f")
  upper_fmt <- ifelse(abs(upper) < threshold, "%.0e", "%.2f")
  sprintf(paste(lower_fmt, upper_fmt, sep = ", "), lower, upper)
}

#' Calculate t-test and format results
#' @param data Data table
#' @param variable Variable name
#' @param group_var Grouping variable
#' @return List with test statistics
ttest_summary <- function(data, variable, group_var = "SEX") {
  tt <- data[, t.test(get(variable) ~ get(group_var), na.rm = TRUE)]
  list(
    X = variable,
    Tstat = tt$statistic,
    DF = tt$parameter,
    Pval = tt$p.value,
    Pval_fmt = format_pval(tt$p.value)
  )
}

# ----- Abbreviations and Labels -----

#' Get family/distribution abbreviations and labels
#' @return Named list with abbrev and label vectors
get_family_info <- function() {
  list(
    abbrev = c(
      NO = "Norm",
      L_NO = "L-Norm",
      BCCG = "BCCG",
      BE = "Beta"
    ),
    label = c(
      NO = "Normal (Gaussian)",
      L_NO = "Logit-Normal",
      BCCG = "Box-Cox Cole & Green",
      BE = "Beta"
    )
  )
}

#' Get adjustment method abbreviations (lowercase)
#' @return Named vector of abbreviations
get_adjustment_abbrev <- function() {
  c(
    NON = "Unadj.",
    PRP = "Prop.",
    STX = "Stereo.",
    RES = "Resid."
  )
}

#' Get side labels and abbreviations
#' @return Named list with abbrev and label vectors
get_side_info <- function() {
  list(
    abbrev = c(L = "L", R = "R", LR = "Bilat."),
    label = c(L = "Left", R = "Right", LR = "Bilateral")
  )
}

#' Get interpretation symbols for validation quality
#' @return Named vector mapping interpretation to symbol
get_interpretation_symbols <- function() {
  c(
    Excellent = "***",
    Good = "**",
    Acceptable = "*",
    Poor = "\u2020" # dagger symbol
  )
}

#' Create abbreviation legend for source notes
#' @param type One of "adjustment", "family", "side", "interpretation"
#' @return Character string with abbreviation definitions
get_abbrev_legend <- function(
    type = c("adjustment", "family", "side", "interpretation"),
    filter_values = NULL) {
  type <- match.arg(type)

  switch(type,
    adjustment = paste(
      "Adjustment methods:",
      "Unadj. = Unadjusted (raw volumes);",
      "Prop. = Proportions (volume/TIV);",
      "Stereo. = Stereotaxic (atlas-normalized);",
      "Resid. = Residuals (TIV regressed out)"
    ),
    family = {
      all_fams <- c(
        "Norm" = "Norm = Normal (Gaussian)",
        "L-Norm" = "L-Norm = Logit-Normal",
        "BCCG" = "BCCG = Box-Cox Cole & Green",
        "Beta" = "Beta = Beta distribution"
      )
      if (!is.null(filter_values)) {
        all_fams <- all_fams[names(all_fams) %in% filter_values]
      }
      if (length(all_fams) == 0) {
        "Distributions: None"
      } else {
        paste(
          "Distributions:", paste(all_fams, collapse = "; ")
        )
      }
    },
    side = paste(
      "L = Left hemisphere; R = Right hemisphere;",
      "Bilat. = Bilateral (average)"
    ),
    interpretation = paste(
      "*** = Excellent (MAE < 5); ** = Good (5-10);",
      "* = Acceptable (10-15); \u2020 = Poor (> 15)"
    )
  )
}

#' Get adjustment method labels from config
#' @return Named vector
get_adjustment_labels <- function() {
  tryCatch(
    {
      params <- get_parameter("adjustment_methods")
      setNames(unlist(params), names(params))
    },
    error = function(e) {
      # Fallback if config not available (standalone mode)
      c(
        NON = "Unadjusted",
        PRP = "Proportions",
        STX = "Stereotaxic",
        RES = "Residuals"
      )
    }
  )
}

#' Get ROI labels from config
#' @return Named vector
get_roi_labels <- function() {
  tryCatch(
    {
      params <- get_parameter("rois")
      setNames(unlist(params), names(params))
    },
    error = function(e) {
      # Fallback if config not available (standalone mode)
      c(
        HC  = "Hippocampus",
        LV  = "Lateral Ventricles",
        HVR = "Hippocampus-to-Ventricle ratio",
        ICC = "Intracranial volume"
      )
    }
  )
}

#' Create variable labels for demographics
#' @param use_markdown Whether to use LaTeX/markdown formatting
#' @return Named vector of labels
get_demographic_labels <- function(use_markdown = TRUE) {
  if (use_markdown) {
    c(
      AGE          = "$\\text{Age (years)}$",
      TIME_diff    = "$\\text{Time (months)}$",
      EDUC         = "$\\text{Education (years)}$",
      EDUC_NA      = "$\\text{Missing}$",
      ICC          = "$\\text{Head-size (TIV)}$",
      HC           = "$\\text{Hippocampus (cc)}$",
      VC           = "$\\text{Lat. ventricles (cc)}$",
      MEM_EXEC_z   = "$\\text{Memory \\& exec. func.}$",
      MEM_EXEC_NA  = "$\\text{Missing}$",
      PRSP_z       = "$\\text{Processing speed}$",
      PRSP_NA      = "$\\text{Missing}$"
    )
  } else {
    c(
      AGE          = "Age (years)",
      TIME_diff    = "Time (months)",
      EDUC         = "Education (years)",
      EDUC_NA      = "Missing",
      ICC          = "Head-size (TIV)",
      HC           = "Hippocampus (cc)",
      VC           = "Lat. ventricles (cc)",
      MEM_EXEC_z   = "Memory & exec. func.",
      MEM_EXEC_NA  = "Missing",
      PRSP_z       = "Processing speed",
      PRSP_NA      = "Missing"
    )
  }
}

#' Create path labels for SEM models
#' @param use_markdown Whether to use LaTeX formatting
#' @return Named vector of path labels
get_sem_path_labels <- function(use_markdown = TRUE) {
  if (use_markdown) {
    list(
      tot_SES_COG = "$\\text{SES} \\to \\text{Cog (Total)}$",
      dir_SES_COG = "$\\text{SES} \\to \\text{Cog (Direct)}$",
      SES_ICC_HC_COG = paste0(
        "$\\text{SES} \\to \\text{TIV} \\to ",
        "\\text{HC} \\to \\text{Cog}$"
      ),
      SES_ICC_VC_COG = paste0(
        "$\\text{SES} \\to \\text{TIV} \\to ",
        "\\text{VC} \\to \\text{Cog}$"
      ),
      SES_HC_COG = "$\\text{SES} \\to \\text{HC} \\to \\text{Cog}$",
      SES_VC_COG = "$\\text{SES} \\to \\text{VC} \\to \\text{Cog}$",
      SES_ICC_COG = "$\\text{SES} \\to \\text{TIV} \\to \\text{Cog}$",
      tot_SEX_COG = "$\\text{Sex} \\to \\text{Cog (Total)}$",
      dir_SEX_COG = "$\\text{Sex} \\to \\text{Cog (Direct)}$",
      SEX_ICC_HC_COG = paste0(
        "$\\text{Sex} \\to \\text{TIV} \\to ",
        "\\text{HC} \\to \\text{Cog}$"
      ),
      SEX_ICC_VC_COG = paste0(
        "$\\text{Sex} \\to \\text{TIV} \\to ",
        "\\text{VC} \\to \\text{Cog}$"
      ),
      SEX_HC_COG = "$\\text{Sex} \\to \\text{HC} \\to \\text{Cog}$",
      SEX_VC_COG = "$\\text{Sex} \\to \\text{VC} \\to \\text{Cog}$",
      SEX_ICC_COG = "$\\text{Sex} \\to \\text{TIV} \\to \\text{Cog}$"
    )
  } else {
    list(
      tot_SES_COG     = "SES → Cog (Total)",
      dir_SES_COG     = "SES → Cog (Direct)",
      SES_ICC_HC_COG  = "SES → TIV → HC → Cog",
      SES_ICC_VC_COG  = "SES → TIV → VC → Cog",
      SES_HC_COG      = "SES → HC → Cog",
      SES_VC_COG      = "SES → VC → Cog",
      SES_ICC_COG     = "SES → TIV → Cog",
      tot_SEX_COG     = "Sex → Cog (Total)",
      dir_SEX_COG     = "Sex → Cog (Direct)",
      SEX_ICC_HC_COG  = "Sex → TIV → HC → Cog",
      SEX_ICC_VC_COG  = "Sex → TIV → VC → Cog",
      SEX_HC_COG      = "Sex → HC → Cog",
      SEX_VC_COG      = "Sex → VC → Cog",
      SEX_ICC_COG     = "Sex → TIV → Cog"
    )
  }
}

# ----- LaTeX Post-processing -----

#' Fix LaTeX source notes formatting
#' @param tex_path Path to .tex file
#' @param table_type Type of table: "normative" or "summary"
#' @return Invisibly returns the path
fix_latex_source_notes <- function(
    tex_path, table_type = "normative", bilateral = FALSE) {
  if (!file.exists(tex_path)) {
    return(invisible(tex_path))
  }

  lines <- readLines(tex_path, warn = FALSE)

  # Title/subtitle sizes - consistent across all tables
  title_size <- "16"
  title_leading <- "20"
  subtitle_size <- "14"
  subtitle_leading <- "18"

  # Define table body and source notes font sizes based on table type
  if (table_type == "normative" || table_type == "diagnostic") {
    # Normative and validation diagnostic tables: same formatting
    table_size <- "11"
    table_leading <- "14"
    source_notes_size <- "9"
  } else if (table_type == "diagnostic_stability") {
    # Stability diagnostic tables: slightly smaller
    table_size <- "9"
    table_leading <- "11.5"
    source_notes_size <- "8"
  } else if (table_type == "diagnostic_summary") {
    # Model summary diagnostic tables: smallest (wide tables)
    table_size <- "7.5"
    table_leading <- "9.5"
    source_notes_size <- "7"
  } else if (table_type == "summary") {
    # Summary tables: smaller fonts to fit wide tables
    if (bilateral) {
      table_size <- "9.5"
      table_leading <- "10"
      source_notes_size <- "6"
    } else {
      table_size <- "6"
      table_leading <- "7.5"
      source_notes_size <- "6"
    }
  } else {
    stop("table_type must be 'normative', 'summary', 'diagnostic', 'diagnostic_stability', or 'diagnostic_summary'", call. = FALSE)
  }

  # Fix table font size
  fontsize_line <- grep("^\\\\fontsize\\{", lines)[1]
  if (!is.na(fontsize_line)) {
    lines[fontsize_line] <- sprintf(
      "\\fontsize{%spt}{%spt}\\selectfont",
      table_size, table_leading
    )
  }

  # Fix title/subtitle font sizes in caption
  for (i in seq_along(lines)) {
    if (grepl("\\\\caption\\*\\{", lines[i])) {
      # Next few lines should contain title and subtitle
      # Title line
      if (
        i + 1 <= length(lines) &&
          grepl("\\{\\\\fontsize\\{", lines[i + 1])
      ) {
        lines[i + 1] <- gsub(
          "\\{\\\\fontsize\\{[0-9.]+\\}\\{[0-9.]+\\}",
          sprintf("{\\\\fontsize{%s}{%s}", title_size, title_leading),
          lines[i + 1]
        )
      }
      # Subtitle line
      if (
        i + 2 <= length(lines) &&
          grepl("\\{\\\\fontsize\\{", lines[i + 2])
      ) {
        lines[i + 2] <- gsub(
          "\\{\\\\fontsize\\{[0-9.]+\\}\\{[0-9.]+\\}",
          sprintf(
            "{\\\\fontsize{%s}{%s}",
            subtitle_size, subtitle_leading
          ),
          lines[i + 2]
        )
      }
      # Add vertical space between title and subtitle
      # Replace \\ at end of title line with \\[6pt]
      if (i + 1 <= length(lines)) {
        lines[i + 1] <- gsub(
          "\\\\\\\\\\s*$",
          "\\\\\\\\[6pt]",
          lines[i + 1]
        )
      }
      break
    }
  }

  # Add longtable centering if not present
  ltpost_line <- grep("\\\\setlength\\{\\\\LTpost\\}", lines)
  if (length(ltpost_line) > 0) {
    # Check if LTleft/LTright already exist
    has_ltleft <- any(grepl("\\\\setlength\\{\\\\LTleft\\}", lines))
    if (!has_ltleft) {
      # Insert centering settings after LTpost
      lines <- c(
        lines[1:ltpost_line[1]],
        "\\setlength{\\LTleft}{0pt plus 1fill}",
        "\\setlength{\\LTright}{0pt plus 1fill}",
        lines[(ltpost_line[1] + 1):length(lines)]
      )
    }
  }

  # Find minipage or center with source notes
  minipage_start <- grep("\\\\begin\\{minipage\\}", lines)
  center_start <- grep("\\\\begin\\{center\\}", lines)

  # Process minipage if found
  if (length(minipage_start) > 0) {
    for (start_idx in minipage_start) {
      # Find end of minipage
      end_idx <- start_idx
      while (
        end_idx <= length(lines) &&
          !grepl("\\\\end\\{minipage\\}", lines[end_idx])
      ) {
        end_idx <- end_idx + 1
      }

      if (end_idx <= length(lines)) {
        # Extract source note lines (between minipage tags)
        note_lines <- lines[(start_idx + 1):(end_idx - 1)]

        # Remove gt's formatting commands
        note_lines <- note_lines[
          !grepl("^\\\\centering", note_lines) &
            !grepl("^\\\\fontsize", note_lines)
        ]

        # Clean up trailing \\ but keep separate lines if multiple
        note_lines <- gsub("\\\\\\\\$", "", note_lines)

        # Remove empty or whitespace-only lines
        note_lines <- note_lines[trimws(note_lines) != ""]

        # Only consolidate if there are many short fragments (old behavior)
        # Keep separate if they're already organized as intended
        if (length(note_lines) > 3) {
          # Many fragments - consolidate
          combined_note <- paste(note_lines, collapse = " ")
          note_lines <- combined_note
        } else if (length(note_lines) > 1) {
          # 2-3 lines: add LaTeX line breaks between them
          note_lines <- paste(note_lines, collapse = " \\\\\n")
        }

        # Replace minipage - remove width spec and use centering
        # Add negative vspace to reduce gap between table and notes
        lines <- c(
          lines[1:(start_idx - 1)],
          "\\vspace{-12pt}",
          "\\begin{center}",
          sprintf(
            "\\fontsize{%s}{%s}\\selectfont", source_notes_size,
            as.character(as.numeric(source_notes_size) * 1.2)
          ),
          note_lines,
          "\\end{center}",
          lines[(end_idx + 1):length(lines)]
        )
      }
    }
  } else if (length(center_start) > 0) {
    # If center already exists, add vspace before it if missing
    for (start_idx in center_start) {
      # Check if vspace already exists before this center
      has_vspace <- start_idx > 1 &&
        grepl("\\\\vspace", lines[start_idx - 1])

      if (!has_vspace) {
        # Insert vspace before center
        lines <- c(
          lines[1:(start_idx - 1)],
          "\\vspace{-12pt}",
          lines[start_idx:length(lines)]
        )
      }
    }
  }

  # Fix escaped dollar signs (gt converts $ to \$)
  lines <- gsub("\\\\\\$", "$", lines)

  # Fix \times that gt escaped as \textbackslash{}times
  lines <- gsub(
    "\\\\textbackslash\\{\\}times", "\\\\times", lines
  )

  # Fix superscripts: \textasciicircum{}\{-3\} -> ^{-3}
  lines <- gsub(
    "\\\\textasciicircum\\{\\}\\\\\\{(-?[0-9]+)\\\\\\}",
    "^{\\1}",
    lines
  )

  # Fix percent signs - multiple patterns
  lines <- gsub("\\\\textbackslash\\{\\}\\\\%", "\\\\%", lines)
  lines <- gsub("\\\\textbackslash\\{\\}%", "\\\\%", lines)
  lines <- gsub("textbackslash\\{\\}%", "\\\\%", lines)

  # Fix Unicode symbols
  lines <- gsub("\u2020", "\\\\dag", lines)
  lines <- gsub("\u00d7", "$\\\\times$", lines)
  lines <- gsub("\u207b\u00b3", "$^{-3}$", lines)

  writeLines(lines, tex_path)

  invisible(tex_path)
}

# ----- Table Saving -----

# Note: Use ensure_directory() from data_io.R instead
# Kept for backward compatibility
#' @export
ensure_table_dir <- function(table_dir) {
  ensure_directory(table_dir)
}

#' Save table in multiple formats
#' @param gt_table gt table object
#' @param filename Base filename (without extension)
#' @param output_dir Output directory
#' @param formats Vector of formats ("html", "tex", "rtf")
#' @param fix_latex Whether to apply LaTeX post-processing fixes
#'   (default TRUE)
#' @param create_standalone Whether to create standalone tex documents
#'   (default TRUE)
save_table <- function(
    gt_table, filename, output_dir, formats = c("html", "tex"),
    fix_latex = TRUE, create_standalone = TRUE) {
  ensure_table_dir(output_dir)

  for (fmt in formats) {
    fpath <- file.path(output_dir, paste0(filename, ".", fmt))
    gt::gtsave(gt_table, fpath)

    # Apply LaTeX fixes and create standalone document
    if (fmt == "tex") {
      if (fix_latex) {
        fix_latex_source_notes(fpath)
      }
      if (create_standalone) {
        wrap_latex_table(fpath)
      }
    }
  }

  invisible(filename)
}

#' Save table to both HTML and LaTeX
#' @param gt_table gt table object
#' @param filename Base filename
#' @param html_dir HTML output directory
#' @param tex_dir LaTeX output directory
save_table_dual <- function(gt_table, filename, html_dir, tex_dir) {
  ensure_table_dir(html_dir)
  ensure_table_dir(tex_dir)

  gt::gtsave(gt_table, file.path(html_dir, paste0(filename, ".html")))
  gt::gtsave(gt_table, file.path(tex_dir, paste0(filename, ".tex")))

  invisible(filename)
}

#' Add clinical Z-score spanners to gt table
#' @param gt_table gt table object
#' @param col_names Column names (excluding SEX and AGE_BIN)
#' @param roi_type Type of ROI ("HC", "LV", or "HVR") to determine thresholds
#' @return Modified gt table with Z-score spanners and clinical footnotes
#' @details Hierarchy:
#'   Level 1 (percentiles) <- Level 2 (clinical) <- Level 3 (side)
add_clinical_spanners <- function(gt_table, col_names, roi_type = "HC") {
  # --- 1. Define Z-Score Map (Universal) ---
  # Maps percentile numeric value -> Z-score label
  # Note: 50th percentile is "0 SD" (Average)
  z_map <- c(
    "1"    = "-2.3 SD",
    "2.5"  = "-2.0 SD",
    "5"    = "-1.6 SD",
    "10"   = "-1.3 SD",
    "25"   = "-0.7 SD",
    "50"   = "0 SD", # Median
    "75"   = "+0.7 SD",
    "90"   = "+1.3 SD",
    "95"   = "+1.6 SD",
    "97.5" = "+2.0 SD",
    "99"   = "+2.3 SD"
  )

  # Detect whether columns have:
  # - "SIDE_" prefix (Summary Table)
  # - just "pXX" (Main Table)
  has_side_prefix <- any(grepl("_[p][0-9]", col_names))

  # Extract numeric values
  if (has_side_prefix) {
    # Format: L_p1, R_p99
    p_vals <- unique(sub(".*_p", "", col_names))
    sides <- unique(sub("_p.*", "", col_names))
  } else {
    # Format: p1, p99 (Main Table)
    p_vals <- unique(sub("^p", "", col_names))
    sides <- "dummy" # No side prefix to span
  }

  # Filter map to only what exists in data
  z_map <- z_map[names(z_map) %in% p_vals]

  # Apply Spanners
  for (pz in names(z_map)) {
    lbl <- z_map[[pz]]

    if (has_side_prefix) {
      for (side in sides) {
        target <- paste0(side, "_p", pz)
        if (target %in% col_names) {
          gt_table <- gt_table |>
            gt::tab_spanner(
              label = lbl,
              columns = all_of(target),
              level = 2,
              id = paste0(side, pz)
            )
        }
      }
    } else {
      # Main table (no side prefix)
      target <- paste0("p", pz)
      if (target %in% col_names) {
        gt_table <- gt_table |>
          gt::tab_spanner(
            label = lbl,
            columns = all_of(target),
            level = 2,
            id = paste0("z", pz)
          )
      }
    }
  }

  # --- Add Clinical Interpretation Footnotes ---
  # This keeps the header clean but provides the clinical context requested.

  if (roi_type %in% c("HC", "HVR")) {
    # ATROPHY (Low values are bad)
    note_text <- gt::md(paste(
      "**Clinical Interpretation (Atrophy):**",
      "&lt;1% (-2.3 SD): Severe;",
      "&lt;2.5% (-2.0 SD): Abnormal;",
      "&lt;5% (-1.6 SD): Borderline."
    ))
  } else {
    # ENLARGEMENT (High values are bad)
    note_text <- gt::md(paste(
      "**Clinical Interpretation (Enlargement):**",
      "&gt;99% (+2.3 SD): Severe;",
      "&gt;97.5% (+2.0 SD): Abnormal;",
      "&gt;95% (+1.6 SD): Borderline."
    ))
  }

  gt_table <- gt_table |> gt::tab_source_note(source_note = note_text)

  return(gt_table)
}
