# =============================================================================
# Plotting Core Utilities
# =============================================================================
# Core plotting functions: themes, palettes, and save helpers.
# Part of the R/utils/ module system.
# =============================================================================

# ----- Color Palettes -----

#' Get standard color palette
#' @param type Color palette type ("default", "sex", "roi")
#' @return Named vector of colors
get_palette <- function(type = "default") {
  palettes <- list(
    default = c(
      "#999999", "#E69F00", "#56B4E9", "#009E73",
      "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
    ),
    sex = c(
      Female = "darkred",
      Male = "midnightblue"
    ),
    sex_md = c(
      Female = "<span style='color: darkred;'>Females</span>",
      Male = "<span style='color: midnightblue;'>Males</span>"
    ),
    roi = c(
      HC = "#0072B2",
      LV = "#D55E00",  # Lateral Ventricles (was incorrectly labeled VC)
      HVR = "#009E73",
      ICC = "#999999"
    ),
    # Named adjustment method colors for consistent plotting
    # HVR is self-normalizing (NOT "Unadjusted") - it's a ratio measure
    adjustment = c(
      "Unadjusted" = "#E64B35",
      "Proportions" = "#FA8072",  # salmon (distinct from teal STX)
      "Stereotaxic" = "#00A087",
      "Residualized" = "#3C5488",
      "HVR (Self-Normalizing)" = "#9467BD"
    )
  )

  palettes[[type]]
}

# ----- Plot Themes -----

#' Standard ggplot theme for publications
#' @param base_size Base font size
#' @param use_markdown Whether to use ggtext markdown
#' @return ggplot2 theme object
theme_publication <- function(base_size = 11, use_markdown = FALSE) {
  base_theme <- ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1),
      plot.caption = ggplot2::element_text(size = base_size - 3),
      legend.position = "bottom"
    )

  if (use_markdown) {
    base_theme <- base_theme +
      ggplot2::theme(
        axis.text.x = ggtext::element_markdown(),
        axis.text.y = ggtext::element_markdown(),
        strip.text = ggtext::element_markdown(),
        plot.title = ggtext::element_markdown()
      )
  }

  base_theme
}

# ----- Plot Helpers -----

# Note: Use ensure_directory() from data_io.R instead
# Kept for backward compatibility
#' @export
ensure_plot_dir <- function(plot_dir) {
  ensure_directory(plot_dir)
}

#' Save plot with standard settings
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution
save_plot <- function(plot, filename, width = 7, height = 7, dpi = 300) {
  ensure_plot_dir(dirname(filename))

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )

  invisible(filename)
}

#' Wrap text to specified width
#' @param text Text to wrap
#' @param width Maximum line width
#' @param indent Indentation for wrapped lines
#' @return Wrapped text
wrap_text <- function(text, width = 100, indent = 0) {
  stringr::str_wrap(text, width = width, exdent = indent)
}
