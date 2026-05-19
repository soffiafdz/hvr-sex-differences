#!/usr/bin/env Rscript

# =============================================================================
# Bifactor Anchor Audit
# =============================================================================
# Diagnostic for the bifactor specification in models/lavaan/cfa_cog.txt and
# the joint SEM measurement models. The current spec fixes THREE indicators
# of the g factor to a unit loading:
#
#   g =~ 1*PRS_mean_time + PRS_mean_inc + NUM + PRMEM_res_n +
#        1*FLINT        + MATS_corr   + TRLS_alnum_time + TOWER_corr +
#        1*REACT        + TRLS_num_time + SYM_corr + SYM_try
#
# Only one indicator is needed for scale identification. The extra two `1*`
# anchors impose raw-scale equality constraints on the g-loadings of those
# three indicators. This script fits three competing specifications on the
# same data and compares them:
#
#   M1  current   : 3 anchors on g (PRS_mean_time, FLINT, REACT = 1)
#   M2  one_anchor: 1 anchor on g (PRS_mean_time = 1)
#   M3  std_lv    : no anchors; var(g) = var(MEM_s) = var(PRSP_s) = 1
#
# The script is READ-ONLY with respect to pipeline artefacts:
#   - reads  data/derivatives/cog_clean-data.rds          (pipeline input)
#   - writes outputs/audit/bifactor_anchors.rds           (scratch only)
#   - writes outputs/audit/bifactor_anchors_summary.csv   (scratch only)
#
# Decision rule:
#   - factor-score correlations > 0.98 across specs  -> downstream 07/10 stable
#   - factor-score correlations < 0.90 in any factor -> expect headline shifts
# =============================================================================

# ----- Libraries -----
library(here)
library(data.table)
library(lavaan)

# ----- Load Utilities -----
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))

init_logger(log_level = "INFO")
log_script_start("audit_bifactor_anchors.R")

config <- load_config()
set_seed()

# ----- Output Locations -----
out_dir <- here("outputs/audit")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_rds <- file.path(out_dir, "bifactor_anchors.rds")
out_csv <- file.path(out_dir, "bifactor_anchors_summary.csv")

# ----- Input -----
log_section("Loading cognitive test data")
cog.lst <- read_rds_safe(
  here("data/derivatives/cog_clean-data.rds"),
  "Cleaned cognitive data"
)

# Reproduce the 12-indicator wide table from 11_sem_analysis.R / 06.
memory.dt <- cog.lst$Memory[,     .(EID, SESSION, PRS_mean_time, PRS_mean_inc, NUM, PRMEM_res_n)]
procsp.dt <- cog.lst$Proc_speed[, .(EID, SESSION, REACT, TRLS_alnum_time, TRLS_num_time, SYM_corr, SYM_try)]
exec.dt   <- cog.lst$Reas_Exec[,  .(EID, SESSION, FLINT, MATS_corr, TOWER_corr)]

cog_tests.dt <- memory.dt[procsp.dt, on = c("EID", "SESSION"), nomatch = 0]
cog_tests.dt <- cog_tests.dt[exec.dt, on = c("EID", "SESSION"), nomatch = 0]
cog_tests.dt <- cog_tests.dt[SESSION %in% c("ses2", "ses-2")]

indicators <- c(
  "PRS_mean_time", "PRS_mean_inc", "NUM", "PRMEM_res_n",
  "FLINT", "MATS_corr", "TRLS_alnum_time", "TOWER_corr",
  "REACT", "TRLS_num_time", "SYM_corr", "SYM_try"
)

# z-score indicators to mirror script 06 (cfa.scaled_scores: yes).
# Comment the next line out to re-run on raw scales (mirrors script 11's SEM).
cog_tests.dt[, (indicators) := lapply(.SD, \(x) as.numeric(scale(x))),
             .SDcols = indicators]

log_info("N rows = %d (12 indicators, ses-2 only, complete-case after merge)",
         nrow(cog_tests.dt))

# ----- Model Specifications -----

mod_current <- '
  g =~ 1*PRS_mean_time + PRS_mean_inc + NUM + PRMEM_res_n +
       1*FLINT + MATS_corr + TRLS_alnum_time + TOWER_corr +
       1*REACT + TRLS_num_time + SYM_corr + SYM_try

  MEM_s  =~ PRS_mean_time + PRS_mean_inc + NUM + PRMEM_res_n
  PRSP_s =~ REACT + TRLS_num_time + SYM_corr + SYM_try

  g      ~~ 0*MEM_s + 0*PRSP_s
  MEM_s  ~~ 0*PRSP_s

  PRS_mean_time ~~ PRS_mean_inc
  SYM_corr      ~~ SYM_try
  TRLS_num_time ~~ TRLS_alnum_time
'

# M2 keeps only the default first-indicator marker for each factor.
mod_one_anchor <- '
  g =~ PRS_mean_time + PRS_mean_inc + NUM + PRMEM_res_n +
       FLINT + MATS_corr + TRLS_alnum_time + TOWER_corr +
       REACT + TRLS_num_time + SYM_corr + SYM_try

  MEM_s  =~ PRS_mean_time + PRS_mean_inc + NUM + PRMEM_res_n
  PRSP_s =~ REACT + TRLS_num_time + SYM_corr + SYM_try

  g      ~~ 0*MEM_s + 0*PRSP_s
  MEM_s  ~~ 0*PRSP_s

  PRS_mean_time ~~ PRS_mean_inc
  SYM_corr      ~~ SYM_try
  TRLS_num_time ~~ TRLS_alnum_time
'

# M3 fixes the latent variances to 1 instead (std.lv=TRUE behaviour).
mod_std_lv <- mod_one_anchor  # same syntax; std.lv flag handled at fit time

# ----- Fit -----
log_section("Fitting three bifactor specifications")

fit_one <- function(model, std_lv = FALSE, label) {
  log_info("  Fitting: %s", label)
  tryCatch(
    cfa(
      model     = model,
      data      = cog_tests.dt,
      estimator = "MLR",
      missing   = "fiml",
      std.lv    = std_lv
    ),
    error = function(e) {
      log_error("    Failed: %s", e$message)
      NULL
    }
  )
}

fits <- list(
  current    = fit_one(mod_current,    std_lv = FALSE, label = "M1 (current, 3 anchors)"),
  one_anchor = fit_one(mod_one_anchor, std_lv = FALSE, label = "M2 (1 anchor)"),
  std_lv     = fit_one(mod_std_lv,     std_lv = TRUE,  label = "M3 (std.lv = TRUE)")
)

# Drop any that failed
fits <- Filter(\(f) !is.null(f) && lavInspect(f, "converged"), fits)
if (length(fits) < 2) {
  stop("Need at least two converged fits to compare; got ", length(fits))
}

# ----- Fit indices side by side -----
log_section("Fit indices")

measures <- c("chisq.scaled", "df.scaled", "pvalue.scaled",
              "cfi.robust", "tli.robust", "rmsea.robust",
              "rmsea.ci.lower.robust", "rmsea.ci.upper.robust",
              "srmr", "aic", "bic", "npar")

fit_dt <- rbindlist(
  lapply(names(fits), \(nm) {
    fm <- tryCatch(fitMeasures(fits[[nm]], measures), error = \(e) NULL)
    if (is.null(fm)) {
      fm <- fitMeasures(fits[[nm]], intersect(measures, names(fitMeasures(fits[[nm]]))))
    }
    cbind(SPEC = nm, as.data.table(as.list(fm)))
  }),
  fill = TRUE
)
print(fit_dt)

# ----- Standardized loadings side by side -----
log_section("Standardized loadings on g (std.all)")

load_dt <- rbindlist(
  lapply(names(fits), \(nm) {
    std <- standardizedSolution(fits[[nm]])
    setDT(std)
    std[op == "=~", .(SPEC = nm, FACTOR = lhs, INDICATOR = rhs,
                      LOADING = round(est.std, 3),
                      SE = round(se, 3))]
  })
)

load_wide <- dcast(load_dt, FACTOR + INDICATOR ~ SPEC, value.var = "LOADING")
log_info("Standardized loadings (rounded):")
print(load_wide)

# ----- Factor score correlations across specs -----
log_section("Factor-score correlations across specifications")

scores.lst <- lapply(fits, \(f) {
  as.data.table(lavPredict(f, type = "lv"))
})

# Match on row order (lavPredict returns rows aligned to input)
factors_present <- Reduce(intersect, lapply(scores.lst, names))
log_info("Factors recovered in all specs: %s", paste(factors_present, collapse = ", "))

cor_dt <- rbindlist(lapply(factors_present, \(fac) {
  combos <- utils::combn(names(scores.lst), 2, simplify = FALSE)
  rbindlist(lapply(combos, \(pair) {
    a <- scores.lst[[pair[1]]][[fac]]
    b <- scores.lst[[pair[2]]][[fac]]
    data.table(
      FACTOR = fac,
      COMPARISON = sprintf("%s vs %s", pair[1], pair[2]),
      r_pearson  = round(cor(a, b, use = "complete.obs"), 4),
      r_spearman = round(cor(a, b, use = "complete.obs", method = "spearman"), 4)
    )
  }))
}))

log_info("Factor-score correlations:")
print(cor_dt)

# ----- Decision summary -----
log_section("Decision summary")

min_r <- cor_dt[, min(r_pearson, na.rm = TRUE)]
log_info("Minimum factor-score correlation across all spec pairs: r = %.4f", min_r)

verdict <- if (min_r >= 0.98) {
  "STABLE: factor scores are near-identical across specs; downstream 07/10 unlikely to shift."
} else if (min_r >= 0.90) {
  "MILD: factor scores correlate >0.90 but <0.98; expect small shifts in r(brain, COG) magnitudes."
} else {
  "MATERIAL: factor scores diverge across specs; expect headline conclusions in 10 to change."
}
log_info("Verdict: %s", verdict)

# ----- Save artefacts -----
log_section("Saving diagnostic outputs")

saveRDS(
  list(
    fits        = fits,
    fit_indices = fit_dt,
    loadings    = load_wide,
    score_cor   = cor_dt,
    min_r       = min_r,
    verdict     = verdict,
    n_rows      = nrow(cog_tests.dt),
    indicators  = indicators
  ),
  out_rds
)
log_info("Saved: %s", out_rds)

fwrite(fit_dt, out_csv)
log_info("Saved: %s", out_csv)

log_script_end("audit_bifactor_anchors.R", success = TRUE)
