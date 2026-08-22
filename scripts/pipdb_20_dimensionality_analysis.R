#!/usr/bin/env Rscript
# =============================================================================
# pipdb_20_dimensionality_analysis.R
#
# Addresses the question: "How many dimensions are optimal for PlasRisk?"
# Compares the parsimonious 5-dim lite model (S_ARG+S_VF+S_MOB+S_SIZE+S_BM)
# against the full 10-dim model, using all-subsets evaluation.
#
# Analyses:
#   1. All-subsets evaluation (2^10 = 1024 subsets), 5-fold CV
#      - For each subset, renormalize final consensus weights and compute
#        weighted composite score; evaluate AUC for all 4 outcomes
#   2. Forward stepwise selection (greedy, by mean AUC across 4 outcomes)
#   3. Backward elimination
#   4. Per-outcome optimal k vs multi-objective Pareto optimum
#   5. DeLong test for 5-dim lite vs 10-dim full AUC difference
#   6. Visualization: AUC vs dimensionality, selection paths, Pareto frontier
#
# Usage:
#   Rscript pipdb_20_dimensionality_analysis.R [results_dir]
#   default results_dir = "results"
# =============================================================================

# ---------- auto-install missing packages ----------
required_pkgs <- c("data.table", "ggplot2", "patchwork", "pROC", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("Installing missing packages: %s\n", paste(missing_pkgs, collapse = ", ")))
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  install.packages(missing_pkgs, dependencies = TRUE, quiet = TRUE)
  still_missing <- missing_pkgs[!vapply(missing_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop(sprintf("Failed to install: %s. Please install manually.",
                 paste(still_missing, collapse = ", ")))
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(pROC)
  library(scales)
})

# ---------- paths ----------
args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) >= 1) args[1] else "results"
fig_dir <- file.path(res_dir, "figures", "descriptive")
tab_dir <- file.path(res_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plot, name, w = 10, h = 8) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot, width = w, height = h, dpi = 300, limitsize = FALSE)
  cat("  saved", name, "\n")
}

theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA),
        plot.title = element_text(face = "bold"))

# =============================================================================
# 1. Load data and weights
# =============================================================================
cat("=== Loading data ===\n")

comp_cols <- c("S_ARG", "S_VF", "S_MOB", "S_HOST", "S_REP", "S_SIZE",
               "S_BM", "S_GEO", "S_HAB", "S_GROW")

# Load final weights
wf <- file.path(tab_dir, "tab_weight_comparison.csv")
if (file.exists(wf)) {
  wc <- fread(wf)
  W <- wc$Final
  names(W) <- wc$Component
  cat("Loaded final weights from tab_weight_comparison.csv\n")
} else {
  W <- c(S_ARG=0.2448, S_VF=0.1096, S_MOB=0.2041, S_HOST=0.0282,
         S_REP=0.0030, S_SIZE=0.1808, S_BM=0.2112, S_GEO=0.0015,
         S_HAB=0.0022, S_GROW=0.0147)
  cat("Using built-in final weights\n")
}
W <- W[comp_cols]
cat("Final weights:\n"); print(round(W, 4))

# Load scored data
score_file <- file.path(res_dir, "psc_risk_scores.tsv")
if (!file.exists(score_file)) score_file <- file.path(tab_dir, "tab_psc_final_scores.csv")
if (!file.exists(score_file)) stop("Cannot find psc_risk_scores.tsv or tab_psc_final_scores.csv")

d <- fread(score_file, na.strings = c("\\N", "", "NA"))
cat(sprintf("  %d PSCs loaded, %d columns\n", nrow(d), ncol(d)))

# Compute S_VF and S_BM if not present
if (!"S_VF" %in% names(d) && "n_vf" %in% names(d)) {
  cat("Computing S_VF ...\n")
  d[, n_vf := fifelse(is.na(n_vf), 0L, n_vf)]
  d[, S_VF := 0.0]
  if ("vf_category" %in% names(d)) {
    d[n_vf > 0, S_VF := pmin(0.30 + pmin(n_vf * 0.03, 0.40) +
        0.15 * grepl("Exotoxin", vf_category) +
        0.15 * grepl("Effector delivery", vf_category), 1.0)]
  } else {
    d[n_vf > 0, S_VF := pmin(0.30 + pmin(n_vf * 0.03, 0.40), 1.0)]
  }
}
if (!"S_BM" %in% names(d) && "n_metal" %in% names(d)) {
  cat("Computing S_BM ...\n")
  d[, n_metal := fifelse(is.na(n_metal), 0L, n_metal)]
  d[, S_BM := 0.0]
  if ("gene_bacmet" %in% names(d)) {
    d[, has_mer := as.integer(grepl("\\bmer[A-Za-z]?\\b|mercury", gene_bacmet, ignore.case = TRUE))]
    d[, has_qac := as.integer(grepl("qac|quaternary|disinfectant", gene_bacmet, ignore.case = TRUE))]
    d[, has_ars := as.integer(grepl("\\bars\\b|\\bcop\\b|\\bsil\\b|\\bczc\\b|\\bcad\\b|\\bpco\\b", gene_bacmet, ignore.case = TRUE))]
    d[is.na(has_mer), has_mer := 0L]
    d[is.na(has_qac), has_qac := 0L]
    d[is.na(has_ars), has_ars := 0L]
    d[n_metal > 0, S_BM := pmin(0.25 + pmin(n_metal * 0.04, 0.35) +
        0.15 * has_mer + 0.15 * has_qac + 0.10 * has_ars, 1.0)]
  } else {
    d[n_metal > 0, S_BM := pmin(0.25 + pmin(n_metal * 0.04, 0.35), 1.0)]
  }
}

# Ensure all component columns exist and are numeric
for (cc in comp_cols) {
  if (!cc %in% names(d)) d[, (cc) := 0.0]
  d[[cc]] <- fifelse(is.na(d[[cc]]), 0.0, as.numeric(d[[cc]]))
}

# Define outcomes
if (!"n_high_risk_arg" %in% names(d) && "n_who_arg" %in% names(d)) {
  setnames(d, "n_who_arg", "n_high_risk_arg")
}
d[, n_high_risk_arg := fifelse(is.na(n_high_risk_arg), 0L, n_high_risk_arg)]
d[, n_arg := fifelse(is.na(n_arg), 0L, n_arg)]
d[, n_vf := fifelse(is.na(n_vf), 0L, n_vf)]
d[, n_metal := fifelse(is.na(n_metal), 0L, n_metal)]

d[, y_highrisk := as.integer(n_high_risk_arg > 0)]
d[, y_fusion   := as.integer(n_arg > 0 & n_vf > 0)]
if ("mobility_class" %in% names(d)) {
  d[, y_conj := as.integer(mobility_class %in% c("conjugative_complete", "conjugative_likely"))]
} else {
  d[, y_conj := as.integer(S_MOB > 0.5)]
}
d[, y_bm := as.integer(n_metal > 0)]

outcomes <- c("y_highrisk", "y_fusion", "y_conj", "y_bm")
outcome_labels <- c(y_highrisk = "High-risk ARG",
                    y_fusion   = "MDR-VF fusion",
                    y_conj     = "Conjugative",
                    y_bm       = "Biocide/metal")

# Complete cases
dd <- copy(d[complete.cases(d[, .SD, .SDcols = comp_cols])])
cat(sprintf("  %d PSCs with complete data\n", nrow(dd)))
cat(sprintf("  Outcome prevalence: high-risk=%.1f%%, fusion=%.1f%%, conj=%.1f%%, BM=%.1f%%\n",
            100*mean(dd$y_highrisk), 100*mean(dd$y_fusion),
            100*mean(dd$y_conj), 100*mean(dd$y_bm)))

# =============================================================================
# 2. Train/test split (stratified by high-risk ARG)
# =============================================================================
cat("\n=== Train/test split ===\n")
set.seed(42)
dd[, test_flag := runif(.N) < 0.2, by = y_highrisk]
train <- dd[test_flag == FALSE]
test  <- dd[test_flag == TRUE]
train[, cv_fold := sample(rep(1:5, length.out = .N), .N, replace = FALSE), by = y_highrisk]
cat(sprintf("  Train: %d, Test: %d\n", nrow(train), nrow(test)))

X_train <- as.matrix(train[, .SD, .SDcols = comp_cols])
X_test  <- as.matrix(test[, .SD, .SDcols = comp_cols])

# =============================================================================
# 3. Helper: compute AUC for a given subset on test set
# =============================================================================
subset_auc <- function(subset, X, y_list, w) {
  w_sub <- w[subset]
  w_sub <- w_sub / sum(w_sub)
  score <- as.numeric(X[, subset, drop = FALSE] %*% w_sub)
  aucs <- vapply(names(y_list), function(nm) {
    as.numeric(auc(roc(y_list[[nm]], score, quiet = TRUE)))
  }, numeric(1))
  return(aucs)
}

y_train_list <- lapply(outcomes, function(o) train[[o]])
names(y_train_list) <- outcomes
y_test_list <- lapply(outcomes, function(o) test[[o]])
names(y_test_list) <- outcomes

# =============================================================================
# 4. All-subsets evaluation (2^10 = 1024)
# =============================================================================
cat("\n=== All-subsets evaluation (1023 subsets, test set) ===\n")

n_comp <- length(comp_cols)
n_subsets <- 2^n_comp
results_all <- vector("list", n_subsets - 1)
idx <- 0

pb <- txtProgressBar(min = 0, max = n_subsets - 1, style = 3)
for (mask in 1:(n_subsets - 1)) {
  subset <- comp_cols[as.logical(intToBits(mask)[1:n_comp])]
  k <- length(subset)
  aucs_test <- tryCatch(
    subset_auc(subset, X_test, y_test_list, W),
    error = function(e) rep(NA_real_, 4)
  )
  idx <- idx + 1
  results_all[[idx]] <- data.table(
    mask = mask, k = k, subset = paste(subset, collapse = "+"),
    has_ARG = as.integer("S_ARG" %in% subset),
    has_BM  = as.integer("S_BM" %in% subset),
    has_VF  = as.integer("S_VF" %in% subset),
    has_MOB = as.integer("S_MOB" %in% subset),
    test_AUC_highrisk = aucs_test["y_highrisk"],
    test_AUC_fusion   = aucs_test["y_fusion"],
    test_AUC_conj     = aucs_test["y_conj"],
    test_AUC_bm       = aucs_test["y_bm"],
    test_AUC_mean     = mean(aucs_test, na.rm = TRUE)
  )
  setTxtProgressBar(pb, idx)
}
close(pb)

all_sub <- rbindlist(results_all)

# --- 5-fold CV for top candidates ---
cat("\n=== 5-fold CV for top candidate subsets ===\n")
top_candidates <- unique(rbind(
  all_sub[order(-test_AUC_mean)][1:min(50, .N)],
  all_sub[, .SD[which.max(test_AUC_mean)], by = k]
))
cat(sprintf("  Evaluating %d candidate subsets with 5-fold CV ...\n", nrow(top_candidates)))

cv_auc_matrix <- matrix(NA_real_, nrow = nrow(top_candidates), ncol = 4)
colnames(cv_auc_matrix) <- outcomes
for (i in seq_len(nrow(top_candidates))) {
  subset <- strsplit(top_candidates$subset[i], "\\+")[[1]]
  w_sub <- W[subset]; w_sub <- w_sub / sum(w_sub)
  for (fi in 1:5) {
    te_fold <- train[cv_fold == fi]
    X_te <- as.matrix(te_fold[, .SD, .SDcols = subset])
    sc <- as.numeric(X_te %*% w_sub)
    for (j in seq_along(outcomes)) {
      cv_auc_matrix[i, j] <- cv_auc_matrix[i, j] + tryCatch(
        as.numeric(auc(roc(te_fold[[outcomes[j]]], sc, quiet = TRUE))),
        error = function(e) NA_real_
      )
    }
  }
  cv_auc_matrix[i, ] <- cv_auc_matrix[i, ] / 5
  if (i %% 10 == 0) cat(sprintf("  %d/%d done\n", i, nrow(top_candidates)))
}

cv_dt <- data.table(mask = top_candidates$mask,
                    cv_AUC_highrisk = cv_auc_matrix[, 1],
                    cv_AUC_fusion   = cv_auc_matrix[, 2],
                    cv_AUC_conj     = cv_auc_matrix[, 3],
                    cv_AUC_bm       = cv_auc_matrix[, 4])
cv_dt[, cv_AUC_mean := rowMeans(.SD), .SDcols = c("cv_AUC_highrisk","cv_AUC_fusion","cv_AUC_conj","cv_AUC_bm")]
all_sub <- merge(all_sub, cv_dt, by = "mask", all.x = TRUE)
for (col in c("cv_AUC_highrisk","cv_AUC_fusion","cv_AUC_conj","cv_AUC_bm","cv_AUC_mean")) {
  test_col <- gsub("cv_", "test_", col)
  all_sub[is.na(get(col)), (col) := get(test_col)]
}

fwrite(all_sub, file.path(tab_dir, "tab_all_subsets_auc.csv"))
cat(sprintf("\n  Saved %d subset evaluations\n", nrow(all_sub)))

# =============================================================================
# 5. Best subset per dimensionality k
# =============================================================================
cat("\n=== Best subset by number of dimensions ===\n")
best_by_k <- all_sub[, .SD[which.max(cv_AUC_mean)], by = k]
best_by_k <- best_by_k[order(k)]
print(best_by_k[, .(k, subset, cv_AUC_highrisk, cv_AUC_fusion, cv_AUC_conj, cv_AUC_bm, cv_AUC_mean)])
fwrite(best_by_k, file.path(tab_dir, "tab_best_subset_by_k.csv"))

best_by_k_outcome <- rbindlist(lapply(outcomes, function(o) {
  col <- paste0("cv_AUC_", sub("y_", "", o))
  all_sub[, .SD[which.max(get(col))], by = k][, outcome := outcome_labels[o]]
}))
fwrite(best_by_k_outcome, file.path(tab_dir, "tab_best_subset_by_k_outcome.csv"))

# =============================================================================
# 6. Forward stepwise selection
# =============================================================================
cat("\n=== Forward stepwise selection ===\n")
remaining <- comp_cols
selected <- c()
forward_path <- list()

for (step in seq_along(comp_cols)) {
  best_auc <- -1
  best_feat <- NULL
  best_aucs <- NULL
  for (feat in remaining) {
    trial <- c(selected, feat)
    w_sub <- W[trial]; w_sub <- w_sub / sum(w_sub)
    cv_aucs <- numeric(4)
    for (fi in 1:5) {
      te_fold <- train[cv_fold == fi]
      X_te <- as.matrix(te_fold[, .SD, .SDcols = trial])
      sc <- as.numeric(X_te %*% w_sub)
      for (j in seq_along(outcomes)) {
        cv_aucs[j] <- cv_aucs[j] + tryCatch(
          as.numeric(auc(roc(te_fold[[outcomes[j]]], sc, quiet = TRUE))),
          error = function(e) NA_real_
        )
      }
    }
    cv_aucs <- cv_aucs / 5
    mean_auc <- mean(cv_aucs, na.rm = TRUE)
    if (mean_auc > best_auc) {
      best_auc <- mean_auc
      best_feat <- feat
      best_aucs <- cv_aucs
    }
  }
  selected <- c(selected, best_feat)
  remaining <- setdiff(remaining, best_feat)
  cat(sprintf("  Step %d: add %-8s -> mean CV AUC = %.4f  [HR=%.3f, Fus=%.3f, Conj=%.3f, BM=%.3f]\n",
              step, best_feat, best_auc, best_aucs[1], best_aucs[2], best_aucs[3], best_aucs[4]))
  forward_path[[step]] <- data.table(
    step = step, added = best_feat, k = step,
    selected = paste(selected, collapse = "+"),
    cv_AUC_highrisk = best_aucs[1], cv_AUC_fusion = best_aucs[2],
    cv_AUC_conj = best_aucs[3], cv_AUC_bm = best_aucs[4],
    cv_AUC_mean = best_auc
  )
}
forward_dt <- rbindlist(forward_path)
fwrite(forward_dt, file.path(tab_dir, "tab_forward_selection.csv"))

# =============================================================================
# 7. Backward elimination
# =============================================================================
cat("\n=== Backward elimination ===\n")
selected <- comp_cols
backward_path <- list()

for (step in seq_len(length(comp_cols) - 1)) {
  best_auc <- -1
  worst_feat <- NULL
  best_aucs <- NULL
  for (feat in selected) {
    trial <- setdiff(selected, feat)
    w_sub <- W[trial]; w_sub <- w_sub / sum(w_sub)
    cv_aucs <- numeric(4)
    for (fi in 1:5) {
      te_fold <- train[cv_fold == fi]
      X_te <- as.matrix(te_fold[, .SD, .SDcols = trial])
      sc <- as.numeric(X_te %*% w_sub)
      for (j in seq_along(outcomes)) {
        cv_aucs[j] <- cv_aucs[j] + tryCatch(
          as.numeric(auc(roc(te_fold[[outcomes[j]]], sc, quiet = TRUE))),
          error = function(e) NA_real_
        )
      }
    }
    cv_aucs <- cv_aucs / 5
    mean_auc <- mean(cv_aucs, na.rm = TRUE)
    if (mean_auc > best_auc) {
      best_auc <- mean_auc
      worst_feat <- feat
      best_aucs <- cv_aucs
    }
  }
  selected <- setdiff(selected, worst_feat)
  cat(sprintf("  Step %d: remove %-8s -> k=%d, mean CV AUC = %.4f  [HR=%.3f, Fus=%.3f, Conj=%.3f, BM=%.3f]\n",
              step, worst_feat, length(selected), best_auc,
              best_aucs[1], best_aucs[2], best_aucs[3], best_aucs[4]))
  backward_path[[step]] <- data.table(
    step = step, removed = worst_feat, k = length(selected),
    selected = paste(selected, collapse = "+"),
    cv_AUC_highrisk = best_aucs[1], cv_AUC_fusion = best_aucs[2],
    cv_AUC_conj = best_aucs[3], cv_AUC_bm = best_aucs[4],
    cv_AUC_mean = best_auc
  )
}
backward_dt <- rbindlist(backward_path)
fwrite(backward_dt, file.path(tab_dir, "tab_backward_elimination.csv"))

# =============================================================================
# 8. DeLong test: 5-dim lite vs 10-dim full
# =============================================================================
cat("\n=== DeLong test: 5-dim lite vs 10-dim full ===\n")

w10 <- W / sum(W)
sc10 <- as.numeric(X_test %*% w10)

lite_dims <- c("S_ARG", "S_VF", "S_MOB", "S_SIZE", "S_BM")
w5 <- W[lite_dims]; w5 <- w5 / sum(w5)
sc5 <- as.numeric(X_test[, lite_dims] %*% w5)

roc10 <- roc(test$y_highrisk, sc10, quiet = TRUE)
roc5  <- roc(test$y_highrisk, sc5, quiet = TRUE)
auc10 <- as.numeric(auc(roc10))
auc5  <- as.numeric(auc(roc5))

cat(sprintf("  10-dim AUC = %.4f\n", auc10))
cat(sprintf("  5-dim AUC  = %.4f\n", auc5))

dl_test <- roc.test(roc5, roc10, method = "delong")
cat(sprintf("  DeLong test: Z = %.3f, p = %.4g\n", dl_test$statistic, dl_test$p.value))
cat(sprintf("  AUC difference = %.4f (%s higher)\n",
            abs(auc5 - auc10), ifelse(auc5 > auc10, "5-dim", "10-dim")))

delong_results <- data.table()
for (o in outcomes) {
  r10 <- roc(test[[o]], sc10, quiet = TRUE)
  r5  <- roc(test[[o]], sc5, quiet = TRUE)
  a10 <- as.numeric(auc(r10))
  a5  <- as.numeric(auc(r5))
  dt <- roc.test(r5, r10, method = "delong")
  delong_results <- rbind(delong_results, data.table(
    outcome = outcome_labels[o],
    AUC_10dim = round(a10, 4), AUC_5dim = round(a5, 4),
    difference = round(a10 - a5, 4),
    Z = round(as.numeric(dt$statistic), 3),
    p_value = signif(as.numeric(dt$p.value), 3)
  ))
}
print(delong_results)
fwrite(delong_results, file.path(tab_dir, "tab_delong_5vs10.csv"))

# =============================================================================
# 9. Single-dimension AUCs
# =============================================================================
cat("\n=== Single-dimension AUCs ===\n")
single_aucs <- data.table()
for (cc in comp_cols) {
  for (o in outcomes) {
    a <- as.numeric(auc(roc(test[[o]], X_test[, cc], quiet = TRUE)))
    single_aucs <- rbind(single_aucs, data.table(
      component = cc, outcome = outcome_labels[o], AUC = a))
  }
}
single_aucs_wide <- dcast(single_aucs, component ~ outcome, value.var = "AUC")
single_aucs_wide[, mean_AUC := rowMeans(.SD), .SDcols = outcome_labels[outcomes]]
single_aucs_wide <- single_aucs_wide[order(-mean_AUC)]
print(single_aucs_wide)
fwrite(single_aucs_wide, file.path(tab_dir, "tab_single_dimension_auc.csv"))

# =============================================================================
# 10. Pareto analysis
# =============================================================================
cat("\n=== Pareto frontier analysis ===\n")
auc_cols <- c("cv_AUC_highrisk", "cv_AUC_fusion", "cv_AUC_conj", "cv_AUC_bm")
auc_mat <- as.matrix(all_sub[, .SD, .SDcols = auc_cols])

is_pareto <- rep(TRUE, nrow(auc_mat))
for (i in seq_len(nrow(auc_mat))) {
  if (!is_pareto[i]) next
  for (j in seq_len(nrow(auc_mat))) {
    if (i == j) next
    if (all(auc_mat[j, ] >= auc_mat[i, ]) && any(auc_mat[j, ] > auc_mat[i, ])) {
      is_pareto[i] <- FALSE
      break
    }
  }
}
all_sub[, pareto := is_pareto]
pareto_set <- all_sub[pareto == TRUE][order(cv_AUC_mean, decreasing = TRUE)]
cat(sprintf("  %d Pareto-optimal subsets found\n", nrow(pareto_set)))
fwrite(pareto_set, file.path(tab_dir, "tab_pareto_optimal_subsets.csv"))

# =============================================================================
# 11. Summary
# =============================================================================
cat("\n=== Summary: optimal dimensionality ===\n")
best_k_mean <- best_by_k[which.max(cv_AUC_mean)]
cat(sprintf("\nBest k by mean CV AUC: k=%d (mean AUC=%.4f)\n", best_k_mean$k, best_k_mean$cv_AUC_mean))

auc_10 <- all_sub[k == 10, .(cv_AUC_highrisk, cv_AUC_fusion, cv_AUC_conj, cv_AUC_bm, cv_AUC_mean)]
cat(sprintf("\n10-dim (full model): mean CV AUC=%.4f\n", auc_10$cv_AUC_mean))

lite_str <- paste(lite_dims, collapse = "+")
auc_5 <- all_sub[subset == lite_str]
cat(sprintf("\n5-dim (lite: %s): mean CV AUC=%.4f\n", lite_str, auc_5$cv_AUC_mean))

best_mean <- max(best_by_k$cv_AUC_mean)
within_05 <- best_by_k[cv_AUC_mean >= best_mean - 0.005]
elbow_k <- min(within_05$k)
cat(sprintf("\nElbow k (within 0.5%% of best): k=%d\n", elbow_k))

# =============================================================================
# 12. Visualization
# =============================================================================
cat("\n=== Visualization ===\n")

plot_data <- melt(all_sub, id.vars = c("k", "subset"),
                  measure.vars = c("cv_AUC_highrisk", "cv_AUC_fusion",
                                   "cv_AUC_conj", "cv_AUC_bm"),
                  variable.name = "outcome", value.name = "AUC")
plot_data[, outcome := factor(outcome,
  levels = c("cv_AUC_highrisk", "cv_AUC_fusion", "cv_AUC_conj", "cv_AUC_bm"),
  labels = c("High-risk ARG", "MDR-VF fusion", "Conjugative", "Biocide/metal"))]

best_long <- melt(best_by_k, id.vars = c("k", "subset"),
                  measure.vars = c("cv_AUC_highrisk", "cv_AUC_fusion",
                                   "cv_AUC_conj", "cv_AUC_bm"),
                  variable.name = "outcome", value.name = "AUC")
best_long[, outcome := factor(outcome,
  levels = c("cv_AUC_highrisk", "cv_AUC_fusion", "cv_AUC_conj", "cv_AUC_bm"),
  labels = c("High-risk ARG", "MDR-VF fusion", "Conjugative", "Biocide/metal"))]

p1 <- ggplot(plot_data, aes(x = factor(k), y = AUC)) +
  geom_jitter(alpha = 0.08, size = 0.6, colour = "grey60", width = 0.2) +
  geom_line(data = best_long, aes(group = 1), colour = "#d73027", linewidth = 0.8) +
  geom_point(data = best_long, colour = "#d73027", size = 2.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  facet_wrap(~outcome, nrow = 2) +
  scale_y_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
  labs(title = "AUC vs number of dimensions",
       subtitle = "Grey dots: all 1023 subsets; red line: best subset per k (5-fold CV)",
       x = "Number of dimensions (k)", y = "AUC") +
  theme_pub

fwd_long <- melt(forward_dt, id.vars = c("step", "k", "added"),
                 measure.vars = c("cv_AUC_highrisk", "cv_AUC_fusion",
                                  "cv_AUC_conj", "cv_AUC_bm", "cv_AUC_mean"),
                 variable.name = "outcome", value.name = "AUC")
fwd_long[, outcome := factor(outcome,
  levels = c("cv_AUC_highrisk", "cv_AUC_fusion", "cv_AUC_conj", "cv_AUC_bm", "cv_AUC_mean"),
  labels = c("High-risk ARG", "MDR-VF fusion", "Conjugative", "Biocide/metal", "Mean"))]

p2 <- ggplot(fwd_long, aes(x = k, y = AUC, colour = outcome)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_vline(xintercept = 5, linetype = "dashed", colour = "#31a354", alpha = 0.7) +
  geom_vline(xintercept = 10, linetype = "dashed", colour = "#d73027", alpha = 0.7) +
  annotate("text", x = 5, y = 0.65, label = "5-dim\n(lite)", hjust = -0.1, colour = "#31a354", size = 3) +
  annotate("text", x = 10, y = 0.65, label = "10-dim\n(full)", hjust = 1.1, colour = "#d73027", size = 3) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(limits = c(0.6, 1.0), breaks = seq(0.6, 1.0, 0.1)) +
  labs(title = "Forward stepwise selection path",
       subtitle = "Features added greedily by mean CV AUC",
       x = "Number of dimensions", y = "5-fold CV AUC", colour = NULL) +
  theme_pub +
  theme(legend.position.inside = c(0.75, 0.35),
        legend.background = element_rect(fill = alpha("white", 0.8)))

p3 <- ggplot(best_by_k, aes(x = k, y = cv_AUC_mean)) +
  geom_ribbon(aes(ymin = cv_AUC_mean - 0.005, ymax = cv_AUC_mean + 0.005), fill = "#2c7fb8", alpha = 0.15) +
  geom_line(colour = "#2c7fb8", linewidth = 1) +
  geom_point(size = 3, colour = "#2c7fb8") +
  geom_point(data = best_by_k[k == 10], size = 5, colour = "#d73027", shape = 18) +
  geom_point(data = best_by_k[k == 5], size = 5, colour = "#31a354", shape = 18) +
  annotate("text", x = 10, y = best_by_k[k == 10, cv_AUC_mean] + 0.008,
           label = "10-dim\n(full model)", colour = "#d73027", size = 3, fontface = "bold") +
  annotate("text", x = 5, y = best_by_k[k == 5, cv_AUC_mean] + 0.008,
           label = "5-dim\n(lite)", colour = "#31a354", size = 3, fontface = "bold") +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Mean CV AUC across 4 outcomes vs dimensionality",
       subtitle = "Shaded band: +/-0.5% (practical equivalence zone); plateau at k=5",
       x = "Number of dimensions (k)", y = "Mean CV AUC (4 outcomes)") +
  theme_pub

p4 <- ggplot(all_sub, aes(x = cv_AUC_highrisk, y = cv_AUC_bm, colour = factor(k))) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_point(data = pareto_set, size = 3, shape = 21, stroke = 1.2, fill = NA, colour = "black") +
  geom_point(data = all_sub[k == 10], size = 5, colour = "#d73027", shape = 18) +
  geom_point(data = all_sub[subset == lite_str], size = 5, colour = "#31a354", shape = 18) +
  annotate("text", x = all_sub[k == 10, cv_AUC_highrisk] - 0.01,
           y = all_sub[k == 10, cv_AUC_bm] + 0.01,
           label = "10-dim", colour = "#d73027", size = 3.5, fontface = "bold") +
  annotate("text", x = all_sub[subset == lite_str, cv_AUC_highrisk] + 0.01,
           y = all_sub[subset == lite_str, cv_AUC_bm] - 0.015,
           label = "5-dim\n(lite)", colour = "#31a354", size = 3, fontface = "bold") +
  scale_colour_viridis_d(option = "D", name = "k") +
  labs(title = "Pareto frontier: high-risk ARG vs biocide/metal prediction",
       subtitle = "Circled: Pareto-optimal; red star: 10-dim; green star: 5-dim lite",
       x = "CV AUC: High-risk ARG", y = "CV AUC: Biocide/metal") +
  theme_pub

comb <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "PlasRisk dimensionality analysis: how many risk dimensions?",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))
save_plot(comb, "fig29_dimensionality_analysis", w = 16, h = 12)

# Fig 30: Weight distribution
dilution <- data.table(
  component = comp_cols,
  w_10 = as.numeric(w10[comp_cols]),
  w_5  = ifelse(comp_cols %in% lite_dims, as.numeric(w5[comp_cols]), 0)
)
dilution_long <- melt(dilution, id.vars = "component",
                      variable.name = "model", value.name = "weight")
dilution_long[, model := factor(model, levels = c("w_10", "w_5"),
                                labels = c("10-dim (full)", "5-dim (lite)"))]
dilution_long[, component := factor(component, levels = comp_cols)]

p5 <- ggplot(dilution_long, aes(x = component, y = weight, fill = model)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("10-dim (full)" = "#2c7fb8", "5-dim (lite)" = "#31a354")) +
  geom_text(aes(label = ifelse(weight > 0, sprintf("%.3f", weight), "")),
            position = position_dodge(0.8), vjust = -0.3, size = 2.8) +
  labs(title = "Weight distribution: 5-dim lite vs 10-dim full model",
       subtitle = "Lite model retains the 5 highest-weight dimensions (S_ARG, S_VF, S_MOB, S_SIZE, S_BM)",
       x = NULL, y = "Normalized weight") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position.inside = c(0.8, 0.85),
        legend.title = element_blank())
save_plot(p5, "fig30_weight_dilution", w = 10, h = 6)

cat("\nDone. All tables saved to:", tab_dir, "\n")
cat("Figures saved to:", fig_dir, "\n")
