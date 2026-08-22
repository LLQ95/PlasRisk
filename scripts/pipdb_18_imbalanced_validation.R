#!/usr/bin/env Rscript
# =============================================================================
# pipdb_18_imbalanced_validation.R
# Imbalanced / natural-prevalence validation of PlasRisk.
# Addresses reviewer concern: balanced 20v20 external set may overstate AUC.
#   1. Natural-prevalence internal holdout (stratified by replicon, ~12% high-risk)
#   2. AUC, PR-AUC, precision/recall/F1 at multiple thresholds
#   3. Calibration analysis (observed vs predicted by decile)
#   4. Decision-curve analysis (net benefit across thresholds)
#   5. Helper shell script to score random NCBI RefSeq plasmids via PlasRisk CLI
# =============================================================================
required_pkgs <- c("data.table", "pROC", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("Installing missing packages: %s\n", paste(missing_pkgs, collapse=", ")))
  install.packages(missing_pkgs, repos="https://cloud.r-project.org", dependencies=TRUE)
}
has_PRROC <- requireNamespace("PRROC", quietly=TRUE)
if (!has_PRROC) {
  cat("Installing PRROC ...\n")
  tryCatch(install.packages("PRROC", repos="https://cloud.r-project.org", dependencies=TRUE),
    error=function(e) cat("PRROC install failed, will use pROC-based PR-AUC fallback\n"))
  has_PRROC <- requireNamespace("PRROC", quietly=TRUE)
}
suppressPackageStartupMessages({
  library(data.table); library(pROC); library(ggplot2)
  if (has_PRROC) library(PRROC)
})

# PR-AUC computation: use PRROC if available, otherwise pROC + trapezoid fallback
compute_pr_auc <- function(scores, labels) {
  if (has_PRROC) {
    pr <- pr.curve(scores.class0=scores[labels==1], scores.class1=scores[labels==0], curve=FALSE)
    return(as.numeric(pr$auc.integral))
  }
  # Fallback: compute precision-recall curve via pROC and trapezoidal integration
  roc_obj <- roc(labels, scores, quiet=TRUE)
  # Use pROC::coords to get sens/spec at all thresholds
  coords_list <- coords(roc_obj, "all", ret=c("threshold","sens","spec"), transpose=FALSE)
  sens <- coords_list$sens; spec <- coords_list$spec
  # precision = TP/(TP+FP) = sens*prev / (sens*prev + (1-spec)*(1-prev))
  prev <- mean(labels==1)
  precision <- ifelse(sens*prev + (1-spec)*(1-prev) > 0,
    sens*prev / (sens*prev + (1-spec)*(1-prev)), 1)
  recall <- sens
  ord <- order(recall)
  auc <- sum(diff(recall[ord]) * (precision[ord][-length(precision)] + precision[ord][-1]) / 2, na.rm=TRUE)
  return(abs(auc))
}

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results"
tab_dir <- file.path(results_dir, "tables")
fig_dir <- file.path(results_dir, "figures", "descriptive")
dir.create(tab_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

cat("Loading data ...\n")
score_file <- file.path(tab_dir, "tab_psc_final_scores.csv")
if (!file.exists(score_file)) score_file <- file.path(results_dir, "tab_psc_final_scores.csv")
if (!file.exists(score_file)) stop("Cannot find tab_psc_final_scores.csv. Run pipdb_15 first.")
d <- fread(score_file)
cat(sprintf("  %d PSCs loaded\n", nrow(d)))

# Normalize column names (handle both n_arg/n_ARG conventions)
if ("n_ARG" %in% names(d) & !"n_arg" %in% names(d)) setnames(d, "n_ARG", "n_arg")
if ("n_VF" %in% names(d) & !"n_vf" %in% names(d)) setnames(d, "n_VF", "n_vf")
if ("n_BM" %in% names(d) & !"n_metal" %in% names(d)) setnames(d, "n_BM", "n_metal")
if ("S_final" %in% names(d) & !"S_norm" %in% names(d)) setnames(d, "S_final", "S_norm")

d[, y_highrisk := as.integer(n_high_risk_arg > 0)]
d[, y_fusion := as.integer(n_arg > 0 & n_vf > 0)]
if ("mobility_class" %in% names(d)) d[, y_conj := as.integer(mobility_class %in% c("conjugative_complete","conjugative_likely"))] else d[, y_conj := 0L]
if ("n_metal" %in% names(d)) d[, y_bm := as.integer(n_metal > 0)] else d[, y_bm := 0L]
cat(sprintf("  High-risk ARG prevalence: %.1f%%\n", 100*mean(d$y_highrisk)))
cat(sprintf("  MDR-VF fusion prevalence: %.1f%%\n", 100*mean(d$y_fusion)))

# Natural-prevalence holdout
cat("\n=== Natural-prevalence holdout ===\n")
set.seed(42)
major_reps <- d[, .N, by=replicon_primary][N >= 100]$replicon_primary
d[, holdout := 0L]
for (rep in major_reps) {
  idx <- which(d$replicon_primary == rep)
  pos_idx <- idx[d$y_highrisk[idx]==1]; neg_idx <- idx[d$y_highrisk[idx]==0]
  n_pos_hold <- max(1, floor(length(pos_idx)*0.2)); n_neg_hold <- max(1, floor(length(neg_idx)*0.2))
  d[sample(pos_idx, min(n_pos_hold, length(pos_idx))), holdout := 1L]
  d[sample(neg_idx, min(n_neg_hold, length(neg_idx))), holdout := 1L]
}
train <- d[holdout==0]; test <- d[holdout==1]
cat(sprintf("  Training: %d PSCs (%.1f%% high-risk)\n", nrow(train), 100*mean(train$y_highrisk)))
cat(sprintf("  Holdout:  %d PSCs (%.1f%% high-risk)\n", nrow(test), 100*mean(test$y_highrisk)))

# Evaluate
cat("\n=== Performance at natural prevalence ===\n")
evaluate_model <- function(scores, labels, model_name) {
  roc_obj <- roc(labels, scores, quiet=TRUE); auc_roc <- as.numeric(auc(roc_obj))
  auc_pr <- compute_pr_auc(scores, labels)
  pred_ab <- as.integer(scores >= 0.45)
  tp <- sum(pred_ab==1 & labels==1); fp <- sum(pred_ab==1 & labels==0)
  fn <- sum(pred_ab==0 & labels==1); tn <- sum(pred_ab==0 & labels==0)
  precision <- ifelse(tp+fp>0, tp/(tp+fp), 0); recall <- ifelse(tp+fn>0, tp/(tp+fn), 0)
  f1 <- ifelse(precision+recall>0, 2*precision*recall/(precision+recall), 0)
  data.table(model=model_name, AUC_ROC=round(auc_roc,4), AUC_PR=round(auc_pr,4), threshold=0.45,
    TP=tp, FP=fp, FN=fn, TN=tn, precision=round(precision,4), recall=round(recall,4), F1=round(f1,4))
}
results_list <- list()
results_list[[1]] <- evaluate_model(test$S_norm, test$y_highrisk, "PlasRisk (10-dim)")
if ("S_BM" %in% names(test)) {
  w9 <- c(S_ARG=0.2448, S_VF=0.1096, S_MOB=0.2041, S_HOST=0.0282, S_REP=0.0030, S_SIZE=0.1808, S_GEO=0.0015, S_HAB=0.0022, S_GROW=0.0147)
  w9 <- w9/sum(w9); cc9 <- names(w9)
  test[, S9 := as.matrix(.SD) %*% w9, .SDcols=cc9]
  results_list[[2]] <- evaluate_model(test$S9, test$y_highrisk, "PlasRisk (9-dim, no S_BM)")
}
results_list[[3]] <- evaluate_model(test$S_ARG, test$y_highrisk, "S_ARG only")
arg_count <- test$n_arg / max(test$n_arg, na.rm=TRUE)
results_list[[4]] <- evaluate_model(arg_count, test$y_highrisk, "Raw ARG count")
eval_results <- rbindlist(results_list, fill=TRUE)
print(eval_results)
fwrite(eval_results, file.path(tab_dir, "tab_imbalanced_validation.csv"))

# Calibration
cat("\n=== Calibration analysis ===\n")
test[, score_decile := cut(S_norm, breaks=quantile(S_norm, probs=seq(0,1,0.1)), include.lowest=TRUE, labels=FALSE)]
calib <- test[, .(n=.N, mean_predicted=mean(S_norm), observed_highrisk=mean(y_highrisk),
  observed_fusion=mean(y_fusion), observed_conj=mean(y_conj), observed_bm=mean(y_bm)), by=score_decile][order(score_decile)]
print(calib); fwrite(calib, file.path(tab_dir, "tab_calibration_deciles.csv"))
p_calib <- ggplot(calib, aes(x=mean_predicted, y=observed_highrisk)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  geom_point(size=3, color="#2166ac") + geom_line(color="#2166ac", linewidth=0.8) +
  geom_errorbar(aes(ymin=observed_highrisk-1.96*sqrt(observed_highrisk*(1-observed_highrisk)/n),
    ymax=observed_highrisk+1.96*sqrt(observed_highrisk*(1-observed_highrisk)/n)), width=0.01, color="#2166ac") +
  labs(x="Mean predicted risk (S_norm)", y="Observed high-risk ARG rate", title="Calibration (natural-prevalence holdout)") +
  theme_bw(base_size=12) + theme(plot.title=element_text(face="bold"))
ggsave(file.path(fig_dir, "fig27_calibration.png"), p_calib, width=6, height=5, dpi=300)
ggsave(file.path(fig_dir, "fig27_calibration.pdf"), p_calib, width=6, height=5)
cat("  saved fig27_calibration\n")

# Decision curve
cat("\n=== Decision curve analysis ===\n")
thresholds <- seq(0.05, 0.80, by=0.05)
dca <- rbindlist(lapply(thresholds, function(pt) {
  pred <- as.integer(test$S_norm >= pt)
  tp <- sum(pred==1 & test$y_highrisk==1); fp <- sum(pred==1 & test$y_highrisk==0)
  n <- nrow(test); prevalence <- mean(test$y_highrisk)
  data.table(threshold=pt, net_benefit=tp/n - fp/n*(pt/(1-pt)),
    treat_all=prevalence-(1-prevalence)*(pt/(1-pt)), treat_none=0)
}))
p_dca <- ggplot(dca, aes(x=threshold)) +
  geom_line(aes(y=net_benefit, color="PlasRisk"), linewidth=1) +
  geom_line(aes(y=treat_all, color="Treat all"), linetype="dashed") +
  geom_line(aes(y=treat_none, color="Treat none"), linetype="dotted") +
  scale_color_manual(values=c("PlasRisk"="#2166ac","Treat all"="grey50","Treat none"="grey50"), name=NULL) +
  labs(x="Risk threshold", y="Net benefit", title="Decision curve analysis") +
  theme_bw(base_size=12) + theme(plot.title=element_text(face="bold"), legend.position=c(0.8,0.8))
ggsave(file.path(fig_dir, "fig28_decision_curve.png"), p_dca, width=6, height=5, dpi=300)
ggsave(file.path(fig_dir, "fig28_decision_curve.pdf"), p_dca, width=6, height=5)
cat("  saved fig28_decision_curve\n")

# Grade distribution
cat("\n=== Grade distribution at natural prevalence ===\n")
test[, grade := ifelse(S_norm>=0.60,"A", ifelse(S_norm>=0.45,"B", ifelse(S_norm>=0.30,"C", ifelse(S_norm>=0.15,"D","E"))))]
grade_dist <- test[, .N, by=.(grade, y_highrisk)]
grade_dist_wide <- dcast(grade_dist, grade ~ y_highrisk, value.var="N", fill=0)
if ("0" %in% names(grade_dist_wide)) setnames(grade_dist_wide, "0", "n_negative") else grade_dist_wide[, n_negative := 0L]
if ("1" %in% names(grade_dist_wide)) setnames(grade_dist_wide, "1", "n_positive") else grade_dist_wide[, n_positive := 0L]
grade_dist_wide[, total := n_negative+n_positive]
grade_dist_wide[, pct_positive := round(100*n_positive/total,1)]
grade_dist_wide[, pct_of_total := round(100*total/sum(total),1)]
print(grade_dist_wide); fwrite(grade_dist_wide, file.path(tab_dir, "tab_grade_distribution_natural.csv"))

# NCBI helper script
cat("\n=== Optional: NCBI random sample validation ===\n")
ncbi_script <- file.path(results_dir, "score_ncbi_sample.sh")
writeLines(c("#!/bin/bash", "# Download and score random NCBI RefSeq plasmids (natural-prevalence validation)", "set -e",
  "OUTDIR=ncbi_validation", "mkdir -p $OUTDIR",
  "if [ ! -f assembly_summary.txt ]; then wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/assembly_summary.txt; fi",
  "awk -F'\\t' '$12==\"Complete Genome\" && $11==\"latest\" {print $20}' assembly_summary.txt | shuf -n 200 > $OUTDIR/ftp_paths.txt",
  "while read ftp_path; do acc=$(basename $ftp_path); if [ ! -f $OUTDIR/${acc}.fna ]; then wget -q ${ftp_path}/${acc}_genomic.fna.gz -O $OUTDIR/${acc}.fna.gz; gunzip -f $OUTDIR/${acc}.fna.gz 2>/dev/null || true; fi; done < $OUTDIR/ftp_paths.txt",
  "plasrisk -o $OUTDIR/results --threads 8 $OUTDIR/*.fna",
  "echo 'Done.'"), ncbi_script)
Sys.chmod(ncbi_script, "0755")
cat(sprintf("  NCBI validation script: %s\n", ncbi_script))
cat("\nDone.\n")
