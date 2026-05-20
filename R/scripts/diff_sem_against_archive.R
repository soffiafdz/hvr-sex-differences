#!/usr/bin/env Rscript

# =============================================================================
# Diff SEM results: new vs archived (pre-bifactor-fix)
# =============================================================================
# Compares brain -> cognition path coefficients between the archived snapshot
# (pre-fix, 3 anchors on g) and the current results (post-fix, 1 anchor).
#
# Read-only. Prints side-by-side tables of:
#   - pooled brain -> g paths
#   - sex-stratified brain -> g paths
#   - fit indices
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

new_root <- here()
old_root <- here("archive/2026-05-19_pre-bifactor-fix")

# ---- helpers ----
get_brain_cog <- function(params.lst) {
  rbindlist(lapply(names(params.lst), \(m) {
    dt <- params.lst[[m]]
    if (is.null(dt)) return(NULL)
    setDT(dt)
    out <- dt[op == "~" &
              lhs %in% c("g", "MEM_s", "PRSP_s") &
              rhs %in% c("HC", "HVR", "HC_RES")]
    if (nrow(out) == 0) return(NULL)
    out[, MODEL := m]
    out[, GROUP := if ("group" %in% names(out)) group else NA_integer_]
    out[, .(MODEL, GROUP, lhs, rhs, est.std, se, ci.lower, ci.upper, pvalue)]
  }), fill = TRUE)
}

format_side_by_side <- function(old, new, key_cols) {
  setnames(old, c("est.std","se","ci.lower","ci.upper","pvalue"),
                c("beta_OLD","se_OLD","lo_OLD","hi_OLD","p_OLD"))
  setnames(new, c("est.std","se","ci.lower","ci.upper","pvalue"),
                c("beta_NEW","se_NEW","lo_NEW","hi_NEW","p_NEW"))
  merged <- merge(old, new, by = key_cols, all = TRUE)
  merged[, dBeta := beta_NEW - beta_OLD]
  merged[, .SD, .SDcols = c(key_cols,
                            "beta_OLD","lo_OLD","hi_OLD","p_OLD",
                            "beta_NEW","lo_NEW","hi_NEW","p_NEW",
                            "dBeta")]
}

# ---- load ----
old_params <- readRDS(file.path(old_root, "models/results/sem_params.rds"))
new_params <- readRDS(file.path(new_root, "models/results/sem_params.rds"))
old_fit    <- readRDS(file.path(old_root, "models/results/sem_fit_measures.rds"))
new_fit    <- readRDS(file.path(new_root, "models/results/sem_fit_measures.rds"))

old_paths <- get_brain_cog(old_params)
new_paths <- get_brain_cog(new_params)

# ---- pooled diff (GROUP == NA) ----
cat("\n====================================================================\n")
cat("POOLED brain -> cognition paths (old vs new)\n")
cat("====================================================================\n")
pooled_old <- old_paths[is.na(GROUP)]
pooled_new <- new_paths[is.na(GROUP)]
pooled_diff <- format_side_by_side(pooled_old, pooled_new,
                                   key_cols = c("MODEL","lhs","rhs"))
setnames(pooled_diff, c("lhs","rhs"), c("OUTCOME","PREDICTOR"))
options(width = 200)
print(pooled_diff[order(MODEL, OUTCOME, PREDICTOR)],
      digits = 3, nrows = 200)

# ---- multigroup diff (GROUP 1 = Female, 2 = Male) ----
cat("\n====================================================================\n")
cat("MULTI-GROUP brain -> cognition paths by sex (old vs new)\n")
cat("  GROUP 1 = Female, GROUP 2 = Male\n")
cat("====================================================================\n")
mg_old <- old_paths[!is.na(GROUP)]
mg_new <- new_paths[!is.na(GROUP)]
mg_diff <- format_side_by_side(mg_old, mg_new,
                               key_cols = c("MODEL","GROUP","lhs","rhs"))
setnames(mg_diff, c("lhs","rhs"), c("OUTCOME","PREDICTOR"))
print(mg_diff[order(MODEL, OUTCOME, PREDICTOR, GROUP)],
      digits = 3, nrows = 200)

# ---- fit indices ----
cat("\n====================================================================\n")
cat("FIT INDICES (old vs new)\n")
cat("====================================================================\n")
fit_compare <- rbindlist(lapply(union(names(old_fit), names(new_fit)), \(m) {
  o <- if (!is.null(old_fit[[m]])) as.data.table(old_fit[[m]]) else NULL
  n <- if (!is.null(new_fit[[m]])) as.data.table(new_fit[[m]]) else NULL
  if (!is.null(o)) o[, VER := "OLD"]
  if (!is.null(n)) n[, VER := "NEW"]
  rbindlist(list(o, n), fill = TRUE)
}), fill = TRUE)

if (nrow(fit_compare) > 0) {
  keep <- intersect(names(fit_compare),
                    c("MODEL","VER","chisq","df","cfi","tli","rmsea","srmr","aic","bic"))
  print(fit_compare[, ..keep][order(MODEL, VER)], digits = 4, nrows = 200)
}

cat("\n====================================================================\n")
cat("Summary: max |dBeta| across all brain->cognition paths\n")
cat("====================================================================\n")
all_diff <- rbind(pooled_diff[, .(MODEL, dBeta)],
                  mg_diff[,     .(MODEL, dBeta)])
cat(sprintf("Pooled:      max|dBeta| = %.4f\n", max(abs(pooled_diff$dBeta), na.rm=TRUE)))
cat(sprintf("Multigroup:  max|dBeta| = %.4f\n", max(abs(mg_diff$dBeta),     na.rm=TRUE)))
cat(sprintf("Overall:     max|dBeta| = %.4f\n", max(abs(all_diff$dBeta),    na.rm=TRUE)))
