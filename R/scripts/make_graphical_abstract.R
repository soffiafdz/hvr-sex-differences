#!/usr/bin/env Rscript
# Graphical abstract for the HBM submission.
#
# Single-panel dot plot on one Cohen's d axis: hippocampal volume flips sign
# across head-size adjustment methods, while HVR gives one consistent answer.
# All effect sizes are read from the manuscript environment so the figure cannot
# drift from the manuscript's reported values.
#
# Run: Rscript R/scripts/make_graphical_abstract.R
# Out: submission/graphical-abstract/graphical-abstract.{pdf,tiff}

suppressMessages({
  renv::load(here::here())
  library(here)
  library(data.table)
  library(ggplot2)
})

here::i_am("R/scripts/make_graphical_abstract.R")
outdir <- here("submission/graphical-abstract")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ms <- readRDS(here("data/derivatives/manuscript_env.rds"))
bil <- as.data.table(ms$data$bilateral)
cmp <- as.data.table(ms$data$sex_diff$COMPARISON)
n_primary <- ms$scalars$n_primary

# Categorical slots 1 and 2 from the validated reference palette
# (validate_palette.js, light mode on white: all checks pass, CVD dE 24.7).
col_hc  <- "#eb6834"
col_hvr <- "#2a78d6"
ink     <- "#0b0b0b"
ink2    <- "#52514e"
rule    <- "#b8b7b2"

hc <- bil[ROI == "HC" & SIDE == "LR",
          .(x = ESTIMATE, method = ADJ_LABEL, y = 2)]
hvr_full <- bil[ROI == "HVR" & SIDE == "LR", ESTIMATE]
hvr_mtch <- cmp[ROI == "HVR" & ADJ == "NON", ESTIMATE_MTCH]

hvr <- data.table(
  x      = c(hvr_full, hvr_mtch),
  method = c("Full sample", "Head-size matched"),
  y      = 1,
  fill   = c(col_hvr, "white")
)

hc_span <- data.table(y = 2, xmin = min(hc$x), xmax = max(hc$x))
hc_range <- diff(range(hc$x))

# Proportions (.58) and Stereotaxic (.59) are visually coincident at this scale:
# two labels on one dot reads as an error. Collapse them into a single honest
# label placed at their midpoint rather than pretending they resolve.
hc_lab <- rbind(
  hc[!method %in% c("Proportions", "Stereotaxic"), .(x, label = method)],
  hc[method %in% c("Proportions", "Stereotaxic"),
     .(x = mean(x), label = "Proportions &\nStereotaxic")]
)

p <- ggplot() +
  # zero reference
  geom_vline(xintercept = 0, linewidth = 0.4, colour = rule) +
  # hippocampus: span bar showing the spread across methods
  geom_segment(data = hc_span,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               linewidth = 3.2, colour = col_hc, alpha = 0.18,
               lineend = "round") +
  geom_point(data = hc, aes(x = x, y = y),
             size = 3.6, colour = col_hc, stroke = 0) +
  geom_text(data = hc_lab, aes(x = x, y = 2, label = label),
            size = 2.4, colour = ink2, vjust = 1, lineheight = 0.95,
            nudge_y = -0.13) +
  # HVR: full sample and matched subsample
  geom_segment(aes(x = hvr_mtch, xend = hvr_full, y = 1, yend = 1),
               linewidth = 0.5, colour = col_hvr, linetype = "22") +
  geom_point(data = hvr, aes(x = x, y = y), shape = 21,
             size = 3.6, colour = col_hvr, fill = hvr$fill, stroke = 1.2) +
  geom_text(data = hvr, aes(x = x, y = y, label = method),
            size = 2.4, colour = ink2, vjust = 1, nudge_y = -0.13) +
  # row labels
  annotate("text", x = -1.12, y = 2.30, hjust = 0, size = 3.0,
           fontface = "bold", colour = ink, label = "Hippocampal volume") +
  annotate("text", x = -1.12, y = 2.15, hjust = 0, size = 2.5, colour = ink2,
           label = sprintf("4 head-size corrections — sign flips, %.1f SD range",
                           hc_range)) +
  annotate("text", x = -1.12, y = 1.30, hjust = 0, size = 3.0,
           fontface = "bold", colour = ink,
           label = "Hippocampal-to-ventricle ratio (HVR)") +
  annotate("text", x = -1.12, y = 1.15, hjust = 0, size = 2.5, colour = ink2,
           label = "No correction needed — one consistent answer") +
  # direction cues, sitting just above the axis
  annotate("text", x = -0.04, y = 0.60, hjust = 1, size = 2.5, colour = ink2,
           label = "◀  Males larger") +
  annotate("text", x = 0.04, y = 0.60, hjust = 0, size = 2.5, colour = ink2,
           label = "Females larger  ▶") +
  scale_x_continuous(
    name   = "Sex difference in brain structure (Cohen's d)",
    limits = c(-1.15, 0.92),
    breaks = seq(-1.0, 0.5, 0.5)
  ) +
  scale_y_continuous(limits = c(0.45, 2.44)) +
  labs(
    title    = "Hippocampal sex differences depend on the head-size correction",
    subtitle = sprintf("%s cognitively healthy UK Biobank adults",
                       format(n_primary, big.mark = ","))
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title       = element_text(face = "bold", size = 10.5, colour = ink,
                                    margin = margin(b = 2)),
    plot.subtitle    = element_text(size = 8.5, colour = ink2,
                                    margin = margin(b = 8)),
    axis.title.x     = element_text(size = 8.5, colour = ink2,
                                    margin = margin(t = 4)),
    axis.text.x      = element_text(size = 8, colour = ink2),
    axis.title.y     = element_blank(),
    axis.text.y      = element_blank(),
    panel.grid       = element_blank(),
    axis.line.x      = element_line(linewidth = 0.4, colour = rule),
    axis.ticks.x     = element_line(linewidth = 0.4, colour = rule),
    plot.margin      = margin(10, 12, 8, 12)
  )

W <- 6.6; H <- 3.1
ggsave(file.path(outdir, "graphical-abstract.pdf"), p,
       width = W, height = H, device = cairo_pdf)
ggsave(file.path(outdir, "graphical-abstract.tiff"), p,
       width = W, height = H, dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(outdir, "graphical-abstract.png"), p,
       width = W, height = H, dpi = 300, bg = "white")

cat(sprintf("HC range %.2f SD (%.2f to %.2f); HVR full %.2f, matched %.2f\n",
            hc_range, min(hc$x), max(hc$x), hvr_full, hvr_mtch))
cat("Wrote graphical-abstract.{pdf,tiff,png} to", outdir, "\n")
