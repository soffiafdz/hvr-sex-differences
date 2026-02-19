# =============================================================================
# Tables Normative Builders
# =============================================================================
# Functions for creating GAMLSS normative centile tables.
# Part of the R/utils/ module system.
# =============================================================================

# ----- Normative Table Compilation -----

#' Compile individual normative tables into a single LaTeX document
#' @param table_dir Directory containing normative table .tex files
#' @param output_path Output file path (default: table_dir/all_normative_tables.tex)
#' @param pattern Regex pattern to match table files
#' @param title Document title
#' @return Path to compiled document
#' @export
compile_normative_tables <- function(
    table_dir,
    output_path = NULL,
    pattern = "normative_.*\\.tex$",
    title = "Normative Centile Tables") {
  if (is.null(output_path)) {
    output_path <- file.path(table_dir, "all_normative_tables.tex")
  }

  # Find all matching table files
  table_files <- list.files(table_dir, pattern = pattern, full.names = TRUE)

  if (length(table_files) == 0) {
    stop("No table files found matching pattern: ", pattern, call. = FALSE)
  }

  # Parse filenames to organize by ROI -> ADJ -> SEX -> SIDE
  # Format: normative_ROI_ADJ_SIDE_sex.tex
  file_info <- lapply(basename(table_files), \(fname) {
    parts <- strsplit(sub("\\.tex$", "", fname), "_")[[1]]
    list(
      file = fname,
      roi = parts[2],
      adj = parts[3],
      side = parts[4],
      sex = parts[5]
    )
  })

  # Sort by ROI, ADJ, SEX, SIDE
  file_info <- file_info[order(
    sapply(file_info, `[[`, "roi"),
    sapply(file_info, `[[`, "adj"),
    sapply(file_info, `[[`, "sex"),
    factor(sapply(file_info, `[[`, "side"), levels = c("L", "R", "LR"))
  )]

  # Create document header with TOC
  doc_lines <- c(
    "\\documentclass[11pt]{article}",
    "\\usepackage{booktabs}",
    "\\usepackage{longtable}",
    "\\usepackage{geometry}",
    "\\geometry{a4paper, margin=1in}",
    "\\usepackage{caption}",
    "\\usepackage{titlesec}",
    "% Make section headers larger",
    paste0(
      "\\titleformat{\\section}{\\Large\\bfseries}",
      "{\\thesection}{1em}{}"
    ),
    paste0(
      "\\titleformat{\\subsection}{\\large\\bfseries}",
      "{\\thesubsection}{1em}{}"
    ),
    "\\begin{document}",
    "\\pagestyle{plain}",
    "",
    sprintf("\\begin{center}\\Large\\textbf{%s}\\end{center}", title),
    "\\vspace{1em}",
    "",
    "\\tableofcontents",
    "\\clearpage",
    ""
  )

  # Get labels from config
  roi_labels <- get_roi_labels()
  adj_labels <- get_adjustment_labels()

  # Add tables organized by sections
  current_roi <- NULL
  current_adj <- NULL

  for (info in file_info) {
    # ROI section
    if (is.null(current_roi) || current_roi != info$roi) {
      roi_label <- if (info$roi %in% names(roi_labels)) {
        roi_labels[info$roi]
      } else {
        info$roi
      }
      doc_lines <- c(
        doc_lines,
        "",
        sprintf("\\section{%s}", roi_label),
        ""
      )
      current_roi <- info$roi
      current_adj <- NULL
    }

    # Adjustment subsection
    if (is.null(current_adj) || current_adj != info$adj) {
      adj_label <- if (info$adj %in% names(adj_labels)) {
        adj_labels[info$adj]
      } else {
        info$adj
      }
      doc_lines <- c(
        doc_lines,
        sprintf("\\subsection{%s}", adj_label),
        ""
      )
      current_adj <- info$adj
    }

    # Include table - sex and side are in table title
    table_path <- file.path(table_dir, info$file)
    doc_lines <- c(
      doc_lines,
      readLines(table_path),
      "\\clearpage",
      ""
    )
  }

  # Document footer
  doc_lines <- c(doc_lines, "\\end{document}")

  # Write
  writeLines(doc_lines, output_path)

  return(output_path)
}

# ----- Main Normative Table -----

#' Create normative centile table
#' @param norm_data Normative data.table with AGE and centile columns
#' @param title Table title
#' @param subtitle Table subtitle
#' @param roi_type ROI type ("HC", "LV", "HVR")
#' @param adj_type Adjustment type ("NON", "PRP", "STX", "RES")
#' @param family Whether to include distribution family in source notes
#' @param sample_size Whether to include sample size
#' @param add_clinical_labels Whether to add clinical threshold notes
#' @param format Output format ("html" or "latex")
#' @return gt table object
#' @export
create_normative_table <- function(
    norm_data,
    title = "Normative Centiles",
    subtitle = NULL,
    roi_type = "HC",
    adj_type = "NON",
    family = FALSE,
    sample_size = FALSE,
    add_clinical_labels = FALSE,
    format = "html") {
  # Work with a copy to avoid modifying original data

  norm_data <- copy(norm_data)

  # Get centile columns
  cent_cols <- grep("^p[0-9.]+$", names(norm_data), value = TRUE)

  # For Proportions (PRP) on HC/LV, values are very small - multiply by 1000

  is_scaled <- FALSE
  if (adj_type == "PRP" && roi_type %in% c("HC", "LV")) {
    norm_data[, (cent_cols) := lapply(.SD, `*`, 1000), .SDcols = cent_cols]
    is_scaled <- TRUE
  }

  # Create labels - Age without (y), will add note in source notes
  cent_labels <- as.list(
    setNames(
      paste0(as.numeric(sub("p", "", cent_cols)), "th"),
      cent_cols
    )
  )
  cent_labels$AGE <- "Age"

  # Create base table
  gt_table <- norm_data[, c("AGE", cent_cols), with = FALSE] |>
    gt() |>
    tab_header(
      title = title,
      subtitle = subtitle
    ) |>
    cols_label(.list = cent_labels) |>
    fmt_number(columns = all_of(cent_cols), decimals = 2)

  gt_table <- gt_table |>
    cols_align(align = "left", columns = "AGE") |>
    tab_style(
      style = cell_text(align = "left"),
      locations = cells_column_labels(columns = "AGE")
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    )

  # Source notes - each on separate line
  source_notes.lst <- list()
  source_notes.lst[[1]] <- character(0)

  if (sample_size && "N" %in% names(norm_data)) {
    source_notes.lst[[1]] <- c(
      source_notes.lst[[1]],
      sprintf("N = %s;", norm_data[, format(sum(N), big.mark = ",")])
    )
  }

  if (family && "FAMILY" %in% names(norm_data)) {
    fam_info <- get_family_info()
    fam_code <- unique(norm_data$FAMILY)
    fam_label <- ifelse(
      fam_code %in% names(fam_info$label),
      fam_info$label[fam_code],
      fam_code
    )
    source_notes.lst[[1]] <- c(
      source_notes.lst[[1]],
      sprintf("%s distribution.", fam_label)
    )
  }

  # Add Age note
  source_notes.lst[[1]] <- c(
    source_notes.lst[[1]], "Age in years."
  )

  # Add scaling note for proportions
  if (is_scaled) {
    if (format == "latex") {
      source_notes.lst[[2]] <- paste0(
        "Values $\\times 10^{-3}$ ",
        "(multiply by 0.001 to obtain proportions)."
      )
    } else {
      source_notes.lst[[2]] <- paste0(
        "Values \u00d710\u207b\u00b3 ",
        "(multiply by 0.001 to obtain proportions)."
      )
    }
  }

  # Add clinical threshold information to source notes
  if (add_clinical_labels) {
    list_n <- length(source_notes.lst) + 1
    if (roi_type %in% c("HC", "HVR")) {
      if (format == "latex") {
        source_notes.lst[[list_n]] <- paste0(
          "Clinical thresholds: $<$2.5\\% = Severe atrophy (-2SD), ",
          "2.5-5\\% = Borderline (-1.6SD)"
        )
      } else {
        source_notes.lst[[list_n]] <- paste0(
          "Clinical thresholds: <2.5% = Severe atrophy (-2SD), ",
          "2.5-5% = Borderline (-1.6SD)"
        )
      }
    } else if (roi_type == "LV") {
      if (format == "latex") {
        source_notes.lst[[list_n]] <- paste0(
          "Clinical thresholds: $>$97.5\\% = ",
          "Abnormal enlargement (+2SD), ",
          "95-97.5\\% = Borderline (+1.6SD)"
        )
      } else {
        source_notes.lst[[list_n]] <- paste0(
          "Clinical thresholds: >97.5% = ",
          "Abnormal enlargement (+2SD), ",
          "95-97.5% = Borderline (+1.6SD)"
        )
      }
    }
  }

  # Add all source notes as single line
  for (i in seq_along(source_notes.lst)) {
    gt_table <- gt_table |>
      tab_source_note(source_note = paste(
        source_notes.lst[[i]],
        collapse = " "
      ))
  }

  # Format-specific styling
  if (format == "latex") {
    # Formal LaTeX: no colors, clean borders, centered table
    # Note: Font sizes are set by fix_latex_source_notes()
    gt_table <- gt_table |>
      tab_options(
        latex.use_longtable = TRUE,
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.top.style = "solid",
        column_labels.border.bottom.style = "solid"
      )
  } else {
    # HTML: with colors
    gt_table <- gt_table |>
      tab_style(
        style = cell_fill(color = "#E8F4F8"),
        locations = cells_body(columns = "AGE")
      ) |>
      tab_options(
        table.font.size = px(12),
        heading.background.color = "#4A90E2",
        heading.title.font.size = px(16),
        heading.subtitle.font.size = px(14),
        column_labels.font.weight = "bold",
        column_labels.background.color = "#F0F0F0",
        data_row.padding = px(3),
        table.border.top.style = "solid",
        table.border.bottom.style = "solid"
      )
  }

  gt_table
}

# NOTE: add_clinical_spanners() is defined in tables_core.R — do not duplicate here.

# ----- Summary Normative Table -----

#' Create summary normative table with hierarchical year-bins structure
#' @param norm_data_list Named list of normative tables in "Sex_Side" format
#'   (e.g., "Male_L", "Female_R")
#' @param title Table title
#' @param subtitle Table subtitle (optional, NULL for cleaner tables)
#' @param format Output format ("html" or "latex")
#' @param bin_width Age bin width in years (default 5)
#' @param centiles Vector of centiles to display
#'   (default c(.05, .25, .5, .75, .95))
#' @param roi_type Type of ROI for clinical thresholds ("HC", "LV", "HVR")
#' @param adj_type Adjustment type ("NON", "PRP", "STX", "RES")
#' @param add_clinical_labels Whether to add clinical threshold spanners
#' @return gt table object with hierarchical structure:
#'   - Row groups by Sex (Female/Male) with N counts
#'   - Column spanners by Side (Left/Right hemisphere)
#'   - Optional clinical threshold spanners
#'   - Columns for centiles within each side
#' @details Creates compact tables suitable for portrait pages by using
#'   vertical grouping for sex and horizontal spanners for sides
#' @export
create_summary_normative_table <- function(
    norm_data_list,
    title = "Normative Centiles Summary",
    subtitle = NULL,
    format = "html",
    bin_width = 5,
    centiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
    roi_type = "HC",
    adj_type = "NON",
    family = FALSE,
    sample_size = FALSE,
    add_clinical_labels = FALSE) {
  # Process and combine data with proper hierarchical structure
  # Structure: Sex (row groups) -> Side (column spanners) -> Centiles (columns)
  # For bilateral: Sex (row groups) -> Centiles (columns only)

  # Determine if we need to scale for proportions
  is_scaled <- adj_type == "PRP" && roi_type %in% c("HC", "LV")

  # Detect if this is bilateral data (all sides are "LR")
  all_sides <- unique(sapply(strsplit(names(norm_data_list), "_"), `[`, 2))
  is_bilateral <- length(all_sides) == 1 && all_sides[1] == "LR"

  # Extract family information if available (check first table)
  family_label <- NULL
  if (family && length(norm_data_list) > 0) {
    first_dt <- norm_data_list[[1]]
    if ("FAMILY" %in% names(first_dt)) {
      fam_info <- get_family_info()
      fam_code <- unique(first_dt$FAMILY)
      family_label <- ifelse(
        fam_code %in% names(fam_info$label),
        fam_info$label[fam_code],
        fam_code
      )
    }
  }

  processed_list <- lapply(names(norm_data_list), function(key) {
    dt <- copy(norm_data_list[[key]])

    # Parse key to extract sex and side (format: "Sex_Side")
    parts <- strsplit(key, "_")[[1]]
    sex <- parts[1]
    side <- parts[2]

    # Add age bins
    dt[, AGE_BIN := bin_ages(AGE, bin_width)]

    # Get centile column names
    cent_cols <- paste0("p", centiles * 100)

    # Scale proportions if needed (multiply by 1000)
    if (is_scaled) {
      dt[, (cent_cols) := lapply(.SD, `*`, 1000), .SDcols = cent_cols]
    }

    # Average within bins
    binned <- merge(
      dt[, .(N = sum(N)), "AGE_BIN"],
      dt[, lapply(.SD, mean), AGE_BIN, .SDcols = cent_cols]
    )

    # Add grouping variables for hierarchical structure
    binned[, `:=`(SEX = sex, SIDE = side)]

    binned
  })

  # Combine all data
  combined_dt <- rbindlist(processed_list)
  id_cols <- c("SEX", "AGE_BIN")
  sex_labels <- total_n <- NULL
  # Create sex labels with N counts (pluralized)
  if (sample_size && "N" %in% names(combined_dt)) {
    id_cols <- c(id_cols, "N")
    sex_n <- combined_dt[, .(N = sum(N)), .(SEX, SIDE)]
    if (any(sex_n[, any(uniqueN(N) != 1), SEX]$V1)) {
      stop("Left/Right sample sizes are unequal within at least one SEX", call. = FALSE)
    }
    total_n <- sex_n[, unique(.SD[, -"SIDE"])][, sum(N)]
    sex_labels <- sex_n[, unique(.SD[, -"SIDE"])][, setNames(
      sprintf("%ss (N = %s)", SEX, format(N, big.mark = ",")),
      SEX
    )]
  }

  # Reshape data based on whether it's bilateral or unilateral
  cent_cols <- paste0("p", centiles * 100)

  if (is_bilateral) {
    # Bilateral: Simple structure - Sex | Age | p5 | p25 | p50 | ...
    # No SIDE dimension needed
    display_dt <- copy(combined_dt)
    display_dt[, SIDE := NULL]
    setorder(display_dt, SEX, AGE_BIN)
    all_cols <- cent_cols
  } else {
    # Unilateral: Age | 5th_L | 5th_R | 25th_L | 25th_R | ...
    # Centile spanners on top, L/R spanners below

    # Melt to long
    long_dt <- melt(
      combined_dt,
      measure.vars = cent_cols,
      variable.name = "CENTILE",
      value.name = "VALUE"
    )

    # Create centile-side identifier for columns
    long_dt[, CENT_SIDE := paste(CENTILE, SIDE, sep = "_")]

    # Cast to wide: one row per Sex + Age_Bin (+ N) combination
    display_dt <- dcast(
      long_dt[, -c("SIDE", "CENTILE")],
      ... ~ CENT_SIDE,
      value.var = "VALUE"
    )

    # Sort by sex and age
    setorder(display_dt, SEX, AGE_BIN)

    # Reorder columns to ensure proper grouping
    # Expected order: AGE_BIN, p5_L, p5_R, p25_L, p25_R, ...
    data_cols <- setdiff(names(display_dt), id_cols)

    # Separate by centile and ensure L/R order within each centile
    centiles_in_data <- unique(
      sapply(data_cols, \(x) strsplit(x, "_")[[1]][1])
    )
    # Sort centiles numerically
    cent_values <- as.numeric(sub("p", "", centiles_in_data))
    centiles_in_data <- centiles_in_data[order(cent_values)]

    ordered_cols <- unlist(lapply(centiles_in_data, \(c) {
      cent_cols_side <- grep(paste0("^", c, "_"), data_cols, value = TRUE)
      # Ensure L comes before R
      cent_cols_side[order(sub(".*_", "", cent_cols_side))]
    }))

    # Reorder display_dt columns
    setcolorder(display_dt, c(id_cols, ordered_cols))
    all_cols <- ordered_cols
  }

  # Replace SEX values with labeled versions (with N counts)
  if (!is.null(sex_labels)) {
    display_dt[, SEX := sex_labels[SEX]]
  } else {
    display_dt[, SEX := paste0(SEX, "s")]
  }

  # Create the gt table - use AGE_BIN as row name, hide N if present
  # Remove N column and AGE_BIN column from data columns for display
  if (sample_size && "N" %in% names(display_dt)) {
    # Create combined Age(N) label for row names
    display_dt[, AGE_LABEL := sprintf("%s (%d)", AGE_BIN, N)]
    display_dt[, c("AGE_BIN", "N") := NULL]

    gt_table <- display_dt |>
      gt(groupname_col = "SEX", rowname_col = "AGE_LABEL") |>
      tab_stubhead(label = "Age (N)") |>
      fmt_number(columns = all_of(all_cols), decimals = 2)
  } else {
    gt_table <- display_dt |>
      gt(groupname_col = "SEX", rowname_col = "AGE_BIN") |>
      tab_stubhead(label = "Age") |>
      fmt_number(columns = all_of(all_cols), decimals = 2)
  }

  # Add title if provided
  if (!is.null(title)) {
    gt_table <- gt_table |>
      tab_header(
        title = title,
        subtitle = subtitle
      )
  }

  # Create column labels and spanners based on data type
  if (is_bilateral) {
    # Bilateral: Simple centile labels (1th, 2.5th, 5th, etc.)
    col_labels <- setNames(
      paste0(as.numeric(sub("p", "", all_cols)), "th"),
      all_cols
    )
    gt_table <- gt_table |> cols_label(.list = col_labels)
    # No spanners needed for bilateral
  } else {
    # Unilateral: Show side (L or R) as column headers
    col_labels <- setNames(
      sub("^.*_", "", all_cols),
      all_cols
    )
    gt_table <- gt_table |> cols_label(.list = col_labels)

    # Add centile spanners (groups L/R columns by centile)
    # Extract unique centiles from column names
    centiles_in_data <- unique(
      sapply(all_cols, \(x) strsplit(x, "_")[[1]][1])
    )
    for (centile in centiles_in_data) {
      cent_label <- paste0(as.numeric(sub("p", "", centile)), "th")
      cent_pattern_cols <- grep(
        paste0("^", centile, "_"), all_cols,
        value = TRUE
      )

      gt_table <- gt_table |>
        tab_spanner(
          label = cent_label,
          id = paste0("centile_", centile),
          columns = all_of(cent_pattern_cols),
          level = 1
        )
    }
  }

  # Source notes - each on separate line
  source_notes.lst <- list()
  source_notes.lst[[1]] <- character(0)

  if (!is.null(total_n)) {
    source_notes.lst[[1]] <- c(
      source_notes.lst[[1]],
      sprintf("Total N = %s.", format(total_n, big.mark = ","))
    )
  }

  # Add family distribution info if available
  if (!is.null(family_label)) {
    source_notes.lst[[1]] <- c(
      source_notes.lst[[1]],
      sprintf("%s distribution.", family_label)
    )
  }

  # Side legend - different for bilateral vs unilateral
  if (!is_bilateral) {
    source_notes.lst[[1]] <- c(
      source_notes.lst[[1]],
      ifelse(
        is_scaled,
        "L: Left hemisphere & R: Right hemisphere.",
        "L/R: Left/Right hemisphere."
      )
    )
  }

  if (is_scaled) {
    if (is_bilateral || format == "latex") i <- 2 else i <- 1
    if (i == 2) source_notes.lst[[2]] <- character(0)
    source_notes.lst[[i]] <- c(
      source_notes.lst[[i]],
      sprintf(
        "%s %s averaged within %d-year age bins.",
        ifelse(is_bilateral, "Sum of hemispheres", "Values"),
        ifelse(format == "latex", "$\\times 10^{-3}$", "\u00d710\u207b\u00b3"),
        bin_width
      )
    )
  } else {
    if (format == "latex" && is_bilateral) i <- 2 else i <- 1
    if (i == 2) source_notes.lst[[2]] <- character(0)
    source_notes.lst[[i]] <- c(
      source_notes.lst[[i]],
      sprintf(
        "%s within %d-year age bins.",
        ifelse(is_bilateral, "Sum of hemispheres averaged", "Means"),
        bin_width
      )
    )
  }


  # Add clinical threshold information with proper LaTeX escaping
  i <- length(source_notes.lst) + 1
  source_notes.lst[[i]] <- character(0)
  if (add_clinical_labels) {
    if (roi_type %in% c("HC", "HVR")) {
      if (format == "latex") {
        source_notes.lst[[i]] <- paste(
          "Clinical thresholds: $<$2.5\\% = Severe atrophy (-2SD),",
          "2.5-5\\% = Borderline (-1.6SD)."
        )
      } else {
        source_notes.lst[[i]] <- paste(
          "Clinical thresholds: <2.5% = Severe atrophy (-2SD),",
          "2.5-5% = Borderline (-1.6SD)."
        )
      }
    } else if (roi_type == "LV") {
      if (format == "latex") {
        source_notes.lst[[i]] <- paste(
          "Clinical thresholds: $>$97.5\\% =",
          "Abnormal enlargement (+2SD),",
          "95-97.5\\% = Borderline (+1.6SD)."
        )
      } else {
        source_notes.lst[[i]] <- paste(
          "Clinical thresholds: >97.5% =",
          "Abnormal enlargement (+2SD),",
          "95-97.5% = Borderline (+1.6SD)."
        )
      }
    }
  }

  # Add source notes - each as a separate line
  for (i in seq_along(source_notes.lst)) {
    gt_table <- gt_table |>
      tab_source_note(source_note = paste(
        source_notes.lst[[i]],
        collapse = " "
      ))
  }

  # Format-specific styling
  if (format == "latex") {
    # Formal LaTeX style: smaller font to fit wide tables, no colors,
    # clean typography
    # Note: Font sizes are set by fix_latex_source_notes()
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        latex.use_longtable = TRUE,
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.top.style = "solid",
        column_labels.border.bottom.style = "solid"
      )
  } else {
    # HTML style with colors
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        table.font.size = px(10),
        heading.background.color = "#4A90E2",
        heading.title.font.size = px(13),
        heading.subtitle.font.size = px(10),
        column_labels.font.weight = "bold",
        column_labels.background.color = "#F0F0F0",
        data_row.padding = px(3),
        table.border.top.style = "solid",
        table.border.bottom.style = "solid"
      )
  }

  gt_table
}

# ----- Validation and Summary Tables -----

#' Create GAMLSS model validation table
#' @param validation_data Data table with columns: SEX, ROI, ADJ, SIDE,
#'   MAE, INTERPRETATION
#' @param roi_code ROI code to filter for (e.g., "HC", "HVR", "LV")
#' @param format Output format ("html" or "latex")
#' @return gt table object with Sex as row groups
#' @export
create_validation_table <- function(
    validation_data,
    roi_code,
    format = "html") {
  # Get labels
  roi_labels <- get_roi_labels()
  adj_labels <- get_adjustment_labels()
  interp_symbols <- get_interpretation_symbols()

  # Filter for this ROI
  val.dt <- copy(validation_data)[ROI == roi_code]
  if (nrow(val.dt) == 0) return(NULL)

  # Check if all interpretations are the same (no need for indicators)
  use_indicators <- length(unique(val.dt$INTERPRETATION)) > 1

  # Add symbols if needed
  if (use_indicators) {
    val.dt[, MAE_DISPLAY := paste0(MAE, interp_symbols[INTERPRETATION])]
  } else {
    val.dt[, MAE_DISPLAY := as.character(MAE)]
  }

  # Add labels
  val.dt[, ADJ_LABEL := as.character(adj_labels[ADJ])]
  val.dt[, SIDE_LABEL := fcase(
    SIDE == "L", "Left",
    SIDE == "R", "Right",
    SIDE == "LR", "Bilateral"
  )]

  # Reshape: rows = SEX + ADJ, columns = SIDE
  val_wide.dt <- dcast(
    val.dt, SEX + ADJ_LABEL ~ SIDE_LABEL,
    value.var = "MAE_DISPLAY"
  )

  # Column order
  side_order <- c("Left", "Right", "Bilateral")
  side_cols <- side_order[side_order %in% names(val_wide.dt)]
  setcolorder(val_wide.dt, c("SEX", "ADJ_LABEL", side_cols))

  # Create table with Sex grouping
  title <- paste(roi_labels[[roi_code]], "- Model Validation (MAE)")
  gt_table <- val_wide.dt |>
    gt(groupname_col = "SEX", rowname_col = "ADJ_LABEL") |>
    tab_stubhead(label = "Adjustment") |>
    tab_header(title = title)

  # Column labels (no "MAE" label, just the side names)
  col_labels <- setNames(side_cols, side_cols)
  gt_table <- gt_table |> cols_label(.list = col_labels)

  # Source notes
  source_note <- "Mean absolute error on hold-out test set."
  if (use_indicators) {
    source_note <- paste(
      source_note,
      get_abbrev_legend("interpretation"),
      sep = " "
    )
  }
  gt_table <- gt_table |> tab_source_note(source_note = source_note)

  # Format-specific styling
  if (format == "latex") {
    # Note: Font sizes are set by fix_latex_source_notes()
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        latex.use_longtable = TRUE,
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.top.style = "solid",
        column_labels.border.bottom.style = "solid",
        heading.align = "center"
      )
  } else {
    # HTML with colors
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        table.font.size = px(10),
        heading.background.color = "#4A90E2",
        heading.title.font.size = px(13),
        heading.align = "center",
        column_labels.font.weight = "bold",
        column_labels.background.color = "#F0F0F0",
        data_row.padding = px(3),
        table.border.top.style = "solid",
        table.border.bottom.style = "solid"
      )
  }

  gt_table
}

#' Create longitudinal stability table
#' @param stability_data Data table with columns: SEX, ROI, ADJ, SIDE,
#'   CORR, MEAN_CHG, SD_CHG, CROSS
#' @param roi_code ROI code to filter for (e.g., "HC", "HVR", "LV")
#' @param format Output format ("html" or "latex")
#' @return gt table object with Sex as row groups
#' @export
create_stability_table <- function(
    stability_data,
    roi_code,
    format = "html") {
  # Get labels
  roi_labels <- get_roi_labels()
  adj_labels <- get_adjustment_labels()

  # Filter for this ROI
  stab.dt <- copy(stability_data)[ROI == roi_code]
  if (nrow(stab.dt) == 0) return(NULL)

  # Add labels
  stab.dt[, ADJ_LABEL := as.character(adj_labels[ADJ])]
  stab.dt[, SIDE_LABEL := fcase(
    SIDE == "L", "Left",
    SIDE == "R", "Right",
    SIDE == "LR", "Bilateral"
  )]

  # Reshape: rows = SEX + ADJ, columns = SIDE x METRIC
  stab_long.dt <- melt(
    stab.dt,
    id.vars = c("SEX", "ADJ_LABEL", "SIDE_LABEL"),
    measure.vars = c("CORR", "MEAN_CHG", "SD_CHG", "CROSS"),
    variable.name = "METRIC", value.name = "VALUE"
  )
  stab_long.dt[, COL := paste(SIDE_LABEL, METRIC, sep = "_")]

  stab_wide.dt <- dcast(
    stab_long.dt, SEX + ADJ_LABEL ~ COL,
    value.var = "VALUE"
  )

  # Build column order: Left (r, Mean, SD, Cross), Right (...), Bilateral (...)
  cols <- character(0)
  for (side in c("Left", "Right", "Bilateral")) {
    for (metric in c("CORR", "MEAN_CHG", "SD_CHG", "CROSS")) {
      col <- paste(side, metric, sep = "_")
      if (col %in% names(stab_wide.dt)) cols <- c(cols, col)
    }
  }
  setcolorder(stab_wide.dt, c("SEX", "ADJ_LABEL", cols))

  # Create table with Sex grouping
  title <- paste(roi_labels[[roi_code]], "- Longitudinal Stability")
  gt_table <- stab_wide.dt |>
    gt(groupname_col = "SEX", rowname_col = "ADJ_LABEL") |>
    tab_stubhead(label = "Adjustment") |>
    tab_header(title = title)

  # Add side spanners
  for (side in c("Left", "Right", "Bilateral")) {
    side_cols <- grep(paste0("^", side, "_"), cols, value = TRUE)
    if (length(side_cols) > 0) {
      gt_table <- gt_table |>
        tab_spanner(
          label = side, columns = all_of(side_cols),
          id = paste0("spanner_", side)
        )
    }
  }

  # Column labels
  labels <- list()
  for (col in cols) {
    labels[[col]] <- fcase(
      grepl("_CORR$", col), "r",
      grepl("_MEAN_CHG$", col), "Mean",
      grepl("_SD_CHG$", col), "SD",
      grepl("_CROSS$", col), "Cross"
    )
  }
  gt_table <- gt_table |> cols_label(.list = labels)

  # Source notes
  source_note <- paste(
    "r = correlation between timepoints;",
    "Mean/SD = centile change (mean and standard deviation);",
    "Cross = proportion crossing clinical thresholds.",
    "Bilateral = sum of left and right hemispheres."
  )
  gt_table <- gt_table |> tab_source_note(source_note = source_note)

  # Format-specific styling
  if (format == "latex") {
    # Note: Font sizes are set by fix_latex_source_notes()
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        latex.use_longtable = TRUE,
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.top.style = "solid",
        column_labels.border.bottom.style = "solid"
      )
  } else {
    # HTML with colors (smaller font for stability)
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        table.font.size = px(9),
        heading.background.color = "#4A90E2",
        heading.title.font.size = px(12),
        column_labels.font.weight = "bold",
        column_labels.background.color = "#F0F0F0",
        data_row.padding = px(3),
        table.border.top.style = "solid",
        table.border.bottom.style = "solid"
      )
  }

  gt_table
}

#' Create GAMLSS model summary table
#' @param summary_data Data table with columns: SEX, ROI, ADJ, SIDE,
#'   FAMILY, AIC, BIC
#' @param roi_code ROI code to filter for (e.g., "HC", "HVR", "LV")
#' @param format Output format ("html" or "latex")
#' @return gt table object with Sex as row groups
#' @export
create_model_summary_table <- function(
    summary_data,
    roi_code,
    format = "html") {
  # Get labels
  roi_labels <- get_roi_labels()
  adj_labels <- get_adjustment_labels()
  fam_info <- get_family_info()

  # Filter for this ROI
  summ.dt <- copy(summary_data)[ROI == roi_code]
  if (nrow(summ.dt) == 0) return(NULL)

  # Add labels
  summ.dt[, FAM_ABBREV := as.character(fam_info$abbrev[FAMILY])]
  summ.dt[, ADJ_LABEL := as.character(adj_labels[ADJ])]
  summ.dt[, SIDE_LABEL := fcase(
    SIDE == "L", "Left",
    SIDE == "R", "Right",
    SIDE == "LR", "Bilateral"
  )]

  # Convert numeric columns to character for consistent melting
  summ.dt[, AIC := as.character(AIC)]
  summ.dt[, BIC := as.character(BIC)]

  # Reshape: rows = SEX + ADJ, columns = SIDE x METRIC
  summ_long.dt <- melt(
    summ.dt,
    id.vars = c("SEX", "ADJ_LABEL", "SIDE_LABEL"),
    measure.vars = c("FAM_ABBREV", "AIC", "BIC"),
    variable.name = "METRIC", value.name = "VALUE"
  )
  summ_long.dt[, COL := paste(SIDE_LABEL, METRIC, sep = "_")]

  summ_wide.dt <- dcast(
    summ_long.dt, SEX + ADJ_LABEL ~ COL,
    value.var = "VALUE"
  )

  # Build column order: Left (Family, AIC, BIC), Right (...), Bilateral (...)
  cols <- character(0)
  for (side in c("Left", "Right", "Bilateral")) {
    for (metric in c("FAM_ABBREV", "AIC", "BIC")) {
      col <- paste(side, metric, sep = "_")
      if (col %in% names(summ_wide.dt)) cols <- c(cols, col)
    }
  }
  setcolorder(summ_wide.dt, c("SEX", "ADJ_LABEL", cols))

  # Create table with Sex grouping
  title <- paste(roi_labels[[roi_code]], "- Model Summary")
  gt_table <- summ_wide.dt |>
    gt(groupname_col = "SEX", rowname_col = "ADJ_LABEL") |>
    tab_stubhead(label = "Adjustment") |>
    tab_header(title = title)

  # Add side spanners
  for (side in c("Left", "Right", "Bilateral")) {
    side_cols <- grep(paste0("^", side, "_"), cols, value = TRUE)
    if (length(side_cols) > 0) {
      gt_table <- gt_table |>
        tab_spanner(
          label = side, columns = all_of(side_cols),
          id = paste0("spanner_", side)
        )
    }
  }

  # Column labels
  labels <- list()
  for (col in cols) {
    labels[[col]] <- fcase(
      grepl("_FAM_ABBREV$", col), "Family",
      grepl("_AIC$", col), "AIC",
      grepl("_BIC$", col), "BIC"
    )
  }
  gt_table <- gt_table |> cols_label(.list = labels)

  # Source notes
  used_families <- unique(summ.dt$FAM_ABBREV)
  source_note <- paste(
    get_abbrev_legend("family", filter_values = used_families),
    "AIC = Akaike Information Criterion; BIC = Bayesian Information Criterion.",
    sep = " "
  )
  gt_table <- gt_table |> tab_source_note(source_note = source_note)

  # Format-specific styling
  if (format == "latex") {
    # Note: Font sizes are set by fix_latex_source_notes()
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        latex.use_longtable = TRUE,
        table.border.top.style = "solid",
        table.border.bottom.style = "solid",
        column_labels.border.top.style = "solid",
        column_labels.border.bottom.style = "solid"
      )
  } else {
    # HTML with colors (smaller font for wide model summary)
    gt_table <- gt_table |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) |>
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_spanners()
      ) |>
      tab_options(
        table.font.size = px(7.5),
        heading.background.color = "#4A90E2",
        heading.title.font.size = px(11),
        column_labels.font.weight = "bold",
        column_labels.background.color = "#F0F0F0",
        data_row.padding = px(2),
        table.border.top.style = "solid",
        table.border.bottom.style = "solid"
      )
  }

  gt_table
}
