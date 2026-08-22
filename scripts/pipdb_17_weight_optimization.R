#!/usr/bin/env Rscript
# =============================================================================
# pipdb_17_weight_optimization.R
# Data-driven weight determination for the PlasRisk 10-dimension model.
# Methods: RF MDG, LASSO, Entropy, Grid-search AUC optimization, LORO-CV.
# =============================================================================
required_pkgs <- c("data.table", "randomForest", "pROC", "ggplot2", "patchwork", "glmnet")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("Installing missing packages: %s\n", paste(missing_pkgs, collapse=", ")))
  install.packages(missing_pkgs, repos="https://cloud.r-project.org", dependencies=TRUE, quiet=TRUE)
}
suppressPackageStartupMessages({
  library(data.table); library(randomForest); library(pROC); library(ggplot2); library(patchwork); library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) >= 1) args[1] else "results"
fig_dir <- file.path(res_dir, "figures", "descriptive")
tab_dir <- file.path(res_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading data ...\n")
master <- fread(file.path(res_dir, "psc_master.tsv"), na.strings = c("\\N", "", "NA"))
risk   <- fread(file.path(res_dir, "psc_risk_scores.tsv"), na.strings = c("\\N", "", "NA"))
d <- copy(merge(master, risk, by = "id", suffixes = c("", ".r")))
cat(sprintf("  %d PSCs loaded\n", nrow(d)))

d[, n_vf := fifelse(is.na(n_vf), 0L, n_vf)]
d[, n_metal := fifelse(is.na(n_metal), 0L, n_metal)]
d[, n_arg := fifelse(is.na(n_arg), 0L, n_arg)]

# S_VF
d[, S_VF := 0.0]
d[n_vf > 0, S_VF := pmin(0.30 + pmin(n_vf * 0.03, 0.40) +
  0.15 * grepl("Exotoxin", vf_category) + 0.15 * grepl("Effector delivery", vf_category), 1.0)]

# S_BM
d[, has_mer := as.integer(grepl("mer|mercury", gene_bacmet, ignore.case = TRUE))]
d[, has_qac := as.integer(grepl("qac|quaternary|disinfectant", gene_bacmet, ignore.case = TRUE))]
d[, has_ars := as.integer(grepl("ars|cop|sil|czc|cad|pco", gene_bacmet, ignore.case = TRUE))]
d[is.na(has_mer), has_mer := 0L]; d[is.na(has_qac), has_qac := 0L]; d[is.na(has_ars), has_ars := 0L]
d[, S_BM := 0.0]
d[n_metal > 0, S_BM := pmin(0.25 + pmin(n_metal * 0.04, 0.35) + 0.15 * has_mer + 0.15 * has_qac + 0.10 * has_ars, 1.0)]

comp_cols <- c("S_ARG","S_VF","S_MOB","S_HOST","S_REP","S_SIZE","S_BM","S_GEO","S_HAB","S_GROW")
for (cc in comp_cols) d[[cc]] <- fifelse(is.na(d[[cc]]), 0.0, as.numeric(d[[cc]]))

d[, y_highrisk := as.integer(n_high_risk_arg > 0)]
d[, y_fusion   := as.integer(n_arg > 0 & n_vf > 0)]
d[, y_conj     := as.integer(mobility_class %in% c("conjugative_complete","conjugative_likely"))]
d[, y_bm       := as.integer(n_metal > 0)]

dd <- copy(d[!is.na(S_ARG) & !is.na(S_VF) & !is.na(S_MOB) & !is.na(S_HOST) &
  !is.na(S_REP) & !is.na(S_SIZE) & !is.na(S_GEO) & !is.na(S_HAB) & !is.na(S_GROW)])
cat(sprintf("  %d PSCs with complete component data\n", nrow(dd)))

expert_w <- c(S_ARG=0.30, S_VF=0.30, S_MOB=0.20, S_HOST=0.15, S_REP=0.10, S_SIZE=0.10, S_BM=0.10, S_GEO=0.05, S_HAB=0.05, S_GROW=0.05)
expert_w_norm <- expert_w / sum(expert_w)

# 1. RF MDG
cat("\n=== 1. Random Forest MDG-derived weights ===\n")
set.seed(42)
train_idx <- sample(seq_len(nrow(dd)), min(100000, nrow(dd)))
dtrain <- dd[train_idx]
outcomes <- list(highrisk="y_highrisk", fusion="y_fusion", conj="y_conj", bm="y_bm")
rf_mdg_list <- list(); rf_auc <- c()
for (nm in names(outcomes)) {
  ycol <- outcomes[[nm]]
  cat(sprintf("  Training RF for %s ...\n", nm))
  X <- as.data.frame(dtrain[, .SD, .SDcols = comp_cols])
  y <- as.factor(dtrain[[ycol]])
  rf <- randomForest(X, y, ntree = 300, importance = TRUE)
  imp <- importance(rf, type = 2)[, 1]
  rf_mdg_list[[nm]] <- imp
  pred_oob <- predict(rf, type = "prob")[, 2]
  roc_obj <- roc(dtrain[[ycol]], pred_oob, quiet = TRUE)
  rf_auc[nm] <- as.numeric(auc(roc_obj))
  cat(sprintf("    AUC = %.3f\n", rf_auc[nm]))
}
mdg_mat <- do.call(cbind, rf_mdg_list)
mdg_norm <- apply(mdg_mat, 2, function(x) x / sum(x))
mdg_mean <- rowMeans(mdg_norm)
rf_weights <- mdg_mean / sum(mdg_mean)

# 2. LASSO
cat("\n=== 2. LASSO logistic regression ===\n")
X_mat <- as.matrix(dtrain[, .SD, .SDcols = comp_cols])
lasso_weights_list <- list()
for (nm in names(outcomes)) {
  ycol <- outcomes[[nm]]
  cat(sprintf("  Fitting LASSO for %s ...\n", nm))
  y <- dtrain[[ycol]]
  cvfit <- cv.glmnet(X_mat, y, family = "binomial", alpha = 1, nfolds = 5, type.measure = "auc", parallel = FALSE)
  coefs <- as.matrix(coef(cvfit, s = "lambda.min"))[-1, 1]
  coefs_pos <- pmax(coefs, 0)
  if (sum(coefs_pos) > 0) lasso_weights_list[[nm]] <- coefs_pos / sum(coefs_pos) else {
    lasso_weights_list[[nm]] <- rep(1/length(comp_cols), length(comp_cols)); names(lasso_weights_list[[nm]]) <- comp_cols
  }
}
lasso_mat <- do.call(cbind, lasso_weights_list)
lasso_weights <- rowMeans(lasso_mat); lasso_weights <- lasso_weights / sum(lasso_weights)

# 3. Entropy
cat("\n=== 3. Entropy weight method ===\n")
entropy_weights <- function(X) {
  X <- as.matrix(X)
  for (j in seq_len(ncol(X))) {
    rng <- range(X[, j], na.rm = TRUE)
    if (rng[2] > rng[1]) X[, j] <- (X[, j] - rng[1]) / (rng[2] - rng[1]) else X[, j] <- 0
  }
  P <- sweep(X, 2, colSums(X), "/"); P[is.nan(P)] <- 0
  n <- nrow(X); k <- 1 / log(n)
  P_log <- P; P_log[P > 0] <- P[P > 0] * log(P[P > 0])
  e <- -k * colSums(P_log); e <- pmin(pmax(e, 0), 1)
  d <- 1 - e
  if (sum(d) == 0) return(rep(1 / ncol(X), ncol(X)))
  w <- d / sum(d); return(w)
}
set.seed(123)
ent_sample <- dd[sample(.N, min(200000, .N))]
ent_w <- entropy_weights(as.matrix(ent_sample[, .SD, .SDcols = comp_cols])); names(ent_w) <- comp_cols

# 4. Grid-search
cat("\n=== 4. Grid-search weight optimization ===\n")
optimize_weights_auc <- function(X, y_list, init_w, grid = seq(0, 0.5, 0.05), n_iter = 3, seed = 42) {
  set.seed(seed); w <- init_w / sum(init_w)
  mean_auc <- function(wv) {
    aucs <- vapply(names(y_list), function(nm) {
      sc <- as.numeric(X %*% wv); as.numeric(auc(roc(y_list[[nm]], sc, quiet = TRUE)))
    }, numeric(1)); mean(aucs)
  }
  best_auc <- mean_auc(w); cat(sprintf("  Initial mean AUC = %.4f\n", best_auc))
  for (iter in seq_len(n_iter)) {
    for (j in seq_along(w)) {
      best_wj <- w[j]
      for (wj in grid) { w_try <- w; w_try[j] <- wj; w_try <- w_try / sum(w_try); a <- mean_auc(w_try); if (a > best_auc) { best_auc <- a; best_wj <- wj } }
      w[j] <- best_wj; w <- w / sum(w)
    }
    cat(sprintf("  Iteration %d: mean AUC = %.4f\n", iter, best_auc))
  }
  return(list(weights = w, auc = best_auc))
}
set.seed(42)
opt_idx <- sample(seq_len(nrow(dd)), min(50000, nrow(dd)))
X_opt <- as.matrix(dd[opt_idx, .SD, .SDcols = comp_cols])
y_opt_list <- list(y_highrisk=dd[opt_idx]$y_highrisk, y_fusion=dd[opt_idx]$y_fusion, y_conj=dd[opt_idx]$y_conj, y_bm=dd[opt_idx]$y_bm)
opt_result <- optimize_weights_auc(X_opt, y_opt_list, expert_w, grid = seq(0, 0.5, 0.05), n_iter = 3)
opt_weights <- opt_result$weights; opt_auc <- opt_result$auc

# 5. Comparison
cat("\n=== 5. Weight scheme comparison ===\n")
weight_compare <- data.table(Component = comp_cols, Expert = round(expert_w_norm[comp_cols], 4),
  RF_MDG = round(rf_weights[comp_cols], 4), LASSO = round(lasso_weights[comp_cols], 4),
  Entropy = round(ent_w[comp_cols], 4), Optimized = round(opt_weights[comp_cols], 4))
weight_compare[, Final := round((RF_MDG + LASSO + Optimized) / 3, 4)]
weight_compare[, Final := round(Final / sum(Final), 4)]
print(weight_compare)
fwrite(weight_compare, file.path(tab_dir, "tab_weight_comparison.csv"))

set.seed(99)
test_idx <- sample(setdiff(seq_len(nrow(dd)), train_idx), min(50000, nrow(dd) - length(train_idx)))
X_test <- as.matrix(dd[test_idx, .SD, .SDcols = comp_cols])
y_test <- dd[test_idx]
auc_results <- data.table()
for (label in c("Expert","RF_MDG","LASSO","Entropy","Optimized","Final")) {
  wv <- weight_compare[[label]]
  for (outcome in c("y_highrisk","y_fusion","y_conj","y_bm")) {
    sc <- as.numeric(X_test %*% wv)
    a <- as.numeric(auc(roc(y_test[[outcome]], sc, quiet = TRUE)))
    auc_results <- rbind(auc_results, data.table(Scheme = label, Outcome = outcome, AUC = round(a, 4)))
  }
}
fwrite(auc_results, file.path(tab_dir, "tab_weight_auc_comparison.csv"))

# 6. LORO-CV
cat("\n=== 6. Leave-one-replicon-out cross-validation ===\n")
major_reps <- dd[!is.na(replicon_primary), .N, by = replicon_primary][N >= 500]$replicon_primary
loro_aucs <- c()
for (rep in major_reps) {
  te <- dd[replicon_primary == rep]
  if (nrow(te) < 50) next
  sc_te <- as.numeric(as.matrix(te[, .SD, .SDcols = comp_cols]) %*% weight_compare$Final)
  a_te <- tryCatch(as.numeric(auc(roc(te$y_highrisk, sc_te, quiet = TRUE))), error = function(e) NA_real_)
  loro_aucs <- c(loro_aucs, a_te)
}
cat(sprintf("  LORO CV AUC: mean = %.3f, range = [%.3f, %.3f] across %d replicons\n",
  mean(loro_aucs, na.rm = TRUE), min(loro_aucs, na.rm = TRUE), max(loro_aucs, na.rm = TRUE), sum(!is.na(loro_aucs))))

# 7. Visualization
cat("\n=== 7. Visualization ===\n")
wc_long <- melt(weight_compare, id.vars = "Component", variable.name = "Scheme", value.name = "Weight")
wc_long[, Component := factor(Component, levels = comp_cols)]
wc_long[, Scheme := factor(Scheme, levels = c("Expert","RF_MDG","LASSO","Entropy","Optimized","Final"))]
p1 <- ggplot(wc_long, aes(x = Component, y = Weight, fill = Scheme)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "white", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2") + labs(title = "Weight scheme comparison", x = NULL, y = "Weight") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))
p2 <- ggplot(auc_results, aes(x = Outcome, y = AUC, fill = Scheme)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "white", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2") + geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") + ylim(0, 1) +
  labs(title = "Predictive AUC by weight scheme", x = NULL, y = "AUC") + theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))
cat("  Computing ROC curves ...\n")
roc_data_list <- list()
for (label in c("Final","Expert","S_ARG only")) {
  if (label == "S_ARG only") sc <- as.numeric(X_test[, "S_ARG"]) else sc <- as.numeric(X_test %*% weight_compare[[label]])
  roc_obj <- roc(y_test$y_highrisk, sc, quiet = TRUE)
  roc_data_list[[label]] <- data.table(FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities, Scheme = label, AUC = as.numeric(auc(roc_obj)))
}
roc_df <- rbindlist(roc_data_list); roc_df[, Label := sprintf("%s (AUC=%.3f)", Scheme, AUC)]
p3 <- ggplot(roc_df, aes(x = FPR, y = TPR, colour = Label)) + geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") + scale_colour_brewer(palette = "Set1") +
  labs(title = "ROC: high-risk ARG prediction", x = "FPR", y = "TPR") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), legend.position.inside = c(0.75, 0.25), legend.title = element_blank())
p4 <- ggplot(data.table(AUC = loro_aucs), aes(x = AUC)) +
  geom_histogram(bins = 20, fill = "steelblue", colour = "white", alpha = 0.8) +
  geom_vline(xintercept = mean(loro_aucs, na.rm = TRUE), linetype = "dashed", colour = "red", linewidth = 0.8) +
  labs(title = "Leave-one-replicon-out CV AUC", x = "AUC", y = "Count") + theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))
comb <- (p1 | p2) / (p3 | p4) + plot_annotation(title = "PlasRisk weight optimization and validation", theme = theme(plot.title = element_text(face = "bold", size = 14)))
ggsave(file.path(fig_dir, "fig24_weight_optimization.png"), comb, width = 14, height = 10, dpi = 300)
ggsave(file.path(fig_dir, "fig24_weight_optimization.pdf"), comb, width = 14, height = 10)
cat("  saved fig24_weight_optimization\n")

# 8. Final weights
cat("\n=== 8. Final weights ===\n")
final_w <- weight_compare$Final; names(final_w) <- comp_cols
cat("  Final (consensus) weights:\n"); print(round(final_w, 4)); cat(sprintf("  Sum = %.4f\n", sum(final_w)))
dd[, S_final := as.numeric(as.matrix(dd[, .SD, .SDcols = comp_cols]) %*% final_w)]
out_cols <- c("id", "replicon_primary", comp_cols, "S_final", "y_highrisk", "y_fusion", "y_conj", "y_bm")
fwrite(dd[, .SD, .SDcols = out_cols], file.path(tab_dir, "tab_psc_final_scores.csv"))
rep_final <- dd[, lapply(.SD, mean), by = replicon_primary, .SDcols = c(comp_cols, "S_final")][order(-S_final)]
rep_final[, rank := .I]
fwrite(rep_final, file.path(tab_dir, "tab_replicon_final_ranking.csv"))
cat("\nDone. Weight optimization complete.\n")
