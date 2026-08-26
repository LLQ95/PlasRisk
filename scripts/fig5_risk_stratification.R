#!/usr/bin/env Rscript
# =============================================================================
# Figure 5: Risk stratification, benchmarking, external validation, case studies
# (A) Quartile characteristics
# (B) ROC comparison: PlasRisk vs PIPdb vs single-component
# (C) Precision-recall curves at natural prevalence
# (D) Grade distribution and high-risk ARG rates
# (E) Case study radar charts
# (F) External validation ROC (367 independent NCBI plasmids)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(pROC)
})

args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
tdir   <- file.path(resdir, "tables")
fdir   <- file.path(resdir, "figures")
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)

theme_pub <- theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.key.size  = unit(0.35, "cm"),
    legend.text      = element_text(size = 7.5),
    legend.title     = element_text(size = 8),
    axis.text        = element_text(size = 8, color = "#3E2723"),
    axis.title       = element_text(size = 9, color = "#3E2723"),
    plot.tag         = element_text(face = "bold", size = 12, color = "#3E2723"),
    plot.margin      = margin(3, 5, 3, 3)
  )

model_colors <- c("PlasRisk (10D)" = "#E76F51", "PIPdb ordinal" = "#F4A261",
                  "S_ARG only" = "#2A9D8F", "S_BM only" = "#6D4C41")
grade_colors <- c("A" = "#9C2C2C", "B" = "#E76F51", "C" = "#F4A261",
                  "D" = "#A8DADC", "E" = "#2A9D8F")

# ---- (A) Quartile characteristics ----
qg <- fread(file.path(tdir, "tab_risk_qgrade_metrics.csv"))
qg[, Q := factor(Q_grade, levels = c("Q1", "Q2", "Q3", "Q4"),
                 labels = c("Q1\n(Highest)", "Q2", "Q3", "Q4\n(Lowest)"))]
qg_long <- melt(qg, id.vars = "Q",
  measure.vars = c("pct_arg", "pct_highrisk", "pct_conj", "pct_integron", "pct_human"),
  variable.name = "Feature", value.name = "Percent")
qg_long[, Feature := factor(Feature,
  levels = c("pct_arg", "pct_highrisk", "pct_conj", "pct_integron", "pct_human"),
  labels = c("ARG+", "High-risk ARG", "Conjugative", "Integron", "Human host"))]

pA <- ggplot(qg_long, aes(Q, Percent, fill = Feature)) +
  geom_col(position = position_dodge(0.75), width = 0.7, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("#E76F51", "#9C2C2C", "#F4A261", "#2A9D8F", "#6D4C41")) +
  labs(x = "PlasRisk quartile", y = "Percent of PSCs (%)", fill = NULL, tag = "A") +
  theme_pub +
  guides(fill = guide_legend(nrow = 2))

# ---- (B) ROC comparison ----
psc <- fread(file.path(tdir, "tab_psc_final_scores.csv"),
             select = c("S_ARG", "S_VF", "S_MOB", "S_HOST", "S_HAB",
                        "S_SIZE", "S_BM", "S_GROW", "S_final",
                        "y_highrisk", "y_fusion"))
set.seed(42)
test_idx <- sample(seq_len(nrow(psc)), size = floor(0.2 * nrow(psc)))
pst <- psc[test_idx]
pst[, pipdb_raw := (S_HOST + S_HOST + S_HAB + S_ARG + S_VF +
                    2 * S_ARG + S_MOB + S_GROW) / 8 + 0.6]

roc_models <- list(
  "PlasRisk (10D)" = pst$S_final,
  "PIPdb ordinal"  = pst$pipdb_raw,
  "S_ARG only"     = pst$S_ARG,
  "S_BM only"      = pst$S_BM
)
roc_df <- rbindlist(lapply(names(roc_models), function(nm) {
  r <- roc(pst$y_fusion, roc_models[[nm]], quiet = TRUE)
  data.table(Model = nm, FPR = 1 - r$specificities, TPR = r$sensitivities,
             AUC = as.numeric(auc(r)))
}))
auc_lab <- roc_df[, .(AUC = AUC[1]), by = Model]
auc_lab[, lab := sprintf("%s (%.3f)", Model, AUC)]

pB <- ggplot(roc_df, aes(FPR, TPR, color = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#BCAAA4") +
  geom_path(linewidth = 0.6) +
  scale_color_manual(values = model_colors, labels = auc_lab$lab) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "1 - Specificity", y = "Sensitivity", color = NULL, tag = "B") +
  theme_pub +
  guides(color = guide_legend(nrow = 2))

# ---- (C) PR curves at natural prevalence ----
set.seed(123)
pos_idx <- which(pst$y_highrisk == 1)
neg_idx <- which(pst$y_highrisk == 0)
n_pos <- length(pos_idx)
n_neg_target <- floor(n_pos * (1 - 0.031) / 0.031)
neg_samp <- sample(neg_idx, min(n_neg_target, length(neg_idx)))
nat <- pst[c(pos_idx, neg_samp)]

pr_curve <- function(y, score) {
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  tp <- cumsum(y); fp <- cumsum(1 - y)
  prec <- tp / (tp + fp); rec <- tp / sum(y)
  data.table(recall = c(0, rec), precision = c(1, prec))
}
pr_models <- list("PlasRisk (10D)" = nat$S_final, "PIPdb ordinal" = nat$pipdb_raw,
                  "S_ARG only" = nat$S_ARG)
pr_df <- rbindlist(lapply(names(pr_models), function(nm) {
  cbind(Model = nm, pr_curve(nat$y_highrisk, pr_models[[nm]]))
}))
pr_auc <- function(y, score) {
  pr <- pr_curve(y, score)
  sum(diff(pr$recall) * (pr$precision[-1] + pr$precision[-nrow(pr)]) / 2)
}
pr_aucs <- data.table(Model = names(pr_models),
  AUC = sapply(names(pr_models), function(nm) pr_auc(nat$y_highrisk, pr_models[[nm]])))
pr_aucs[, lab := sprintf("%s (%.3f)", Model, AUC)]

pC <- ggplot(pr_df, aes(recall, precision, color = Model)) +
  geom_hline(yintercept = mean(nat$y_highrisk), linetype = "dashed", color = "#BCAAA4") +
  geom_path(linewidth = 0.6) +
  scale_color_manual(values = model_colors, labels = pr_aucs$lab) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Recall", y = "Precision", color = NULL, tag = "C") +
  theme_pub +
  guides(color = guide_legend(nrow = 2))

# ---- (D) Grade distribution ----
set.seed(123)
gd_pos <- pst[y_highrisk == 1]
gd_neg <- pst[y_highrisk == 0]
n_neg_gd <- floor(nrow(gd_pos) * (1 - 0.031) / 0.031)
gd_neg_s <- gd_neg[sample(.N, min(n_neg_gd, nrow(gd_neg)))]
gd_dat <- rbindlist(list(gd_pos, gd_neg_s))
gd_dat[, grade := cut(S_final, breaks = c(-Inf, 0.15, 0.30, 0.45, 0.60, Inf),
                      labels = c("E", "D", "C", "B", "A"))]
gd <- gd_dat[, .(n_positive = sum(y_highrisk), n_negative = sum(1 - y_highrisk),
                 total = .N), by = grade]
gd[, pct_positive := round(n_positive / total * 100, 1)]
gd[, pct_of_total := round(total / sum(total) * 100, 1)]
gd[, grade := factor(grade, levels = c("A", "B", "C", "D", "E"))]
gd <- gd[order(grade)]
gd_long <- melt(gd, id.vars = "grade",
  measure.vars = c("pct_of_total", "pct_positive"),
  variable.name = "Metric", value.name = "Percent")
gd_long[, Metric := factor(Metric,
  levels = c("pct_of_total", "pct_positive"),
  labels = c("% of total PSCs", "% high-risk ARG"))]

pD <- ggplot(gd_long, aes(grade, Percent, fill = grade, alpha = Metric)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = grade_colors, guide = "none") +
  scale_alpha_manual(values = c(0.9, 0.5)) +
  labs(x = "PlasRisk grade", y = "Percent (%)", alpha = NULL, tag = "D") +
  theme_pub

# ---- (E) Case study radar ----
cs <- fread(file.path(tdir, "tab_case_study_scores.csv"))
radar_dims <- c("S_ARG", "S_VF", "S_MOB", "S_HOST", "S_SIZE", "S_BM", "S_GROW")
cs_long <- melt(cs, id.vars = c("seq_id", "grade"), measure.vars = radar_dims,
                variable.name = "Dimension", value.name = "Score")
cs_long[, Dimension := gsub("S_", "", Dimension)]
cs_long[, Dimension := factor(Dimension,
  levels = c("ARG", "VF", "MOB", "HOST", "SIZE", "BM", "GROW"))]
cs_long[, seq_id := factor(seq_id, levels = cs$seq_id)]
case_colors <- c("#E76F51", "#F4A261", "#2A9D8F")

pE <- ggplot(cs_long, aes(Dimension, Score, group = seq_id, color = seq_id, fill = seq_id)) +
  geom_polygon(alpha = 0.1, linewidth = 0.5) +
  geom_point(size = 1.5) +
  coord_polar(start = -pi / length(radar_dims)) +
  scale_color_manual(values = case_colors) +
  scale_fill_manual(values = case_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(color = NULL, fill = NULL, tag = "E") +
  theme_pub +
  theme(axis.text.x = element_text(size = 8), axis.title = element_blank(),
        legend.text = element_text(size = 7))

# ---- (F) External validation ROC ----
ext <- fread(file.path(tdir, "tab_external_400_scores.csv"))
ext[, y := ifelse(label == "high_risk", 1, 0)]
r_ext <- roc(ext$y, ext$S_final, quiet = TRUE)
roc_ext <- data.table(FPR = 1 - r_ext$specificities, TPR = r_ext$sensitivities)
auc_ext <- as.numeric(auc(r_ext))
ci_ext  <- ci.auc(r_ext, quiet = TRUE)

pF <- ggplot(roc_ext, aes(FPR, TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#BCAAA4") +
  geom_path(color = "#E76F51", linewidth = 0.8) +
  annotate("text", x = 0.55, y = 0.12,
           label = sprintf("AUC = %.3f\n(95%% CI %.3f-%.3f)\nn = %d (%d HR / %d LR)",
                           auc_ext, ci_ext[1], ci_ext[3],
                           nrow(ext), sum(ext$y == 1), sum(ext$y == 0)),
           size = 2.8, hjust = 0, color = "#3E2723") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "1 - Specificity", y = "Sensitivity", tag = "F") +
  theme_pub

# ---- assemble ----
fig <- (pA | pB) / (pC | pD) / (pE | pF)

ggsave(file.path(fdir, "Figure5_risk_stratification.pdf"), fig,
       width = 7.2, height = 9, units = "in")
ggsave(file.path(fdir, "Figure5_risk_stratification.png"), fig,
       width = 7.2, height = 9, units = "in", dpi = 300)
cat("Saved Figure5_risk_stratification.pdf/png to", fdir, "\n")
