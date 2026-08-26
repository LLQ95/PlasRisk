#!/usr/bin/env Rscript
# =============================================================================
# Figure 2: Data-driven weight comparison and validation
# (A) Weight heatmap across methods
# (B) ROC curves for final weights across four outcomes
# (C) LORO cross-validation AUC distribution
# (D) Weight perturbation stability
# (E) Decile calibration at natural prevalence
# (F) Decision-curve analysis
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(pROC)
  library(ggrepel)
})

args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
tdir   <- file.path(resdir, "tables")
fdir   <- file.path(resdir, "figures")
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)

# ---- warm palette ----
outcome_colors <- c(
  "High-risk ARG" = "#E76F51",
  "MDR-VF fusion" = "#F4A261",
  "Conjugative"   = "#2A9D8F",
  "Biocide/metal" = "#6D4C41"
)
outcome_levels <- names(outcome_colors)

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
    plot.margin      = margin(4, 6, 4, 4)
  )

# ---- (A) Weight heatmap ----
wtab <- fread(file.path(tdir, "tab_weight_comparison.csv"))
wlong <- melt(wtab, id.vars = "Component", variable.name = "Method", value.name = "Weight")
wlong[, Component := factor(Component, levels = rev(wtab$Component))]
wlong[, Method := factor(Method,
  levels = c("Expert", "RF_MDG", "LASSO", "Entropy", "Optimized", "Final"),
  labels = c("Expert", "RF-MDG", "LASSO", "Entropy", "Grid", "Final"))]

pA <- ggplot(wlong, aes(Method, Component, fill = Weight)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Weight)), size = 2.6, color = "#3E2723") +
  scale_fill_gradient2(low = "#FFF8E7", mid = "#F4A261", high = "#9C2C2C",
                      midpoint = 0.18, limits = c(0, 0.36)) +
  labs(x = NULL, y = NULL, fill = "Weight", tag = "A") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))

# ---- (B) ROC curves ----
psc <- fread(file.path(tdir, "tab_psc_final_scores.csv"),
             select = c("S_final", "y_highrisk", "y_fusion", "y_conj", "y_bm"))
set.seed(42)
test_idx <- sample(seq_len(nrow(psc)), size = floor(0.2 * nrow(psc)))
pst <- psc[test_idx]

roc_list <- list(
  roc(pst$y_highrisk, pst$S_final, quiet = TRUE),
  roc(pst$y_fusion,   pst$S_final, quiet = TRUE),
  roc(pst$y_conj,     pst$S_final, quiet = TRUE),
  roc(pst$y_bm,       pst$S_final, quiet = TRUE)
)
roc_df <- rbindlist(lapply(seq_along(roc_list), function(i) {
  r <- roc_list[[i]]
  data.table(Outcome = outcome_levels[i],
             FPR = 1 - r$specificities, TPR = r$sensitivities,
             AUC = as.numeric(auc(r)))
}))
auc_labels <- roc_df[, .(AUC = AUC[1]), by = Outcome]
auc_labels[, lab := sprintf("%s (AUC=%.3f)", Outcome, AUC)]

pB <- ggplot(roc_df, aes(FPR, TPR, color = Outcome)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#BCAAA4") +
  geom_path(linewidth = 0.6) +
  scale_color_manual(values = outcome_colors, labels = auc_labels$lab) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "1 - Specificity", y = "Sensitivity", color = NULL, tag = "B") +
  theme_pub +
  guides(color = guide_legend(nrow = 2))

# ---- (C) LORO-CV AUC distribution ----
psc_all <- fread(file.path(tdir, "tab_psc_final_scores.csv"),
                 select = c("replicon_primary", "S_final",
                            "y_highrisk", "y_fusion", "y_conj", "y_bm"))
rep_counts <- psc_all[, .N, by = replicon_primary][N >= 500][order(-N)]
top_reps <- rep_counts$replicon_primary[1:min(40, nrow(rep_counts))]

loro_l <- list()
for (rep in top_reps) {
  d <- psc_all[replicon_primary == rep]
  for (oc in c("y_highrisk", "y_fusion", "y_conj", "y_bm")) {
    y <- d[[oc]]
    a <- if (length(unique(y)) == 2) as.numeric(auc(roc(y, d$S_final, quiet = TRUE))) else NA_real_
    loro_l[[length(loro_l) + 1]] <- data.table(replicon = rep, outcome = oc, AUC = a)
  }
}
loro <- rbindlist(loro_l)
loro[, outcome := factor(outcome,
  levels = c("y_highrisk", "y_fusion", "y_conj", "y_bm"),
  labels = outcome_levels)]
loro_sum <- loro[, .(mean_AUC = mean(AUC, na.rm = TRUE)), by = outcome]

pC <- ggplot(loro, aes(AUC, fill = outcome)) +
  geom_histogram(binwidth = 0.03, alpha = 0.75, color = "white", linewidth = 0.2,
                 position = "identity") +
  geom_vline(data = loro_sum, aes(xintercept = mean_AUC, color = outcome),
             linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = outcome_colors) +
  scale_color_manual(values = outcome_colors, guide = "none") +
  scale_x_continuous(limits = c(0.5, 1.02), breaks = seq(0.5, 1.0, 0.1)) +
  labs(x = "AUC (held-out replicon)", y = "Count", fill = NULL, tag = "C") +
  theme_pub +
  guides(fill = guide_legend(nrow = 2))

# ---- (D) Weight perturbation stability ----
pert <- fread(file.path(tdir, "tab_weight_sensitivity.csv"))
mean_rho <- mean(pert$spearman_rho)
mean_ovl <- mean(pert$top10_overlap)

pD <- ggplot(pert, aes(spearman_rho)) +
  geom_histogram(binwidth = 0.002, fill = "#E76F51", color = "white",
                 linewidth = 0.2, alpha = 0.8) +
  geom_vline(xintercept = mean_rho, linetype = "dashed", color = "#6D4C41",
             linewidth = 0.6) +
  annotate("text", x = mean_rho - 0.001, y = Inf, vjust = 1.5, hjust = 1,
           label = sprintf("mean rho = %.3f", mean_rho), size = 2.8,
           color = "#3E2723") +
  annotate("text", x = min(pert$spearman_rho), y = Inf, vjust = 1.5, hjust = 0,
           label = sprintf("mean top-10 overlap = %.1f/10", mean_ovl), size = 2.8,
           color = "#3E2723") +
  scale_x_continuous(limits = c(0.985, 1.0), breaks = seq(0.985, 1.0, 0.005)) +
  labs(x = "Spearman rank correlation", y = "Count", tag = "D") +
  theme_pub

# ---- (E) Decile calibration ----
set.seed(42)
cal_pos <- pst[y_highrisk == 1]
cal_neg <- pst[y_highrisk == 0]
n_neg_cal <- floor(nrow(cal_pos) * (1 - 0.031) / 0.031)
cal_neg_s <- cal_neg[sample(.N, min(n_neg_cal, nrow(cal_neg)))]
cal_dat <- rbindlist(list(cal_pos, cal_neg_s))
cal_glm <- glm(y_highrisk ~ S_final, data = cal_dat, family = binomial)
cal_dat[, pred := predict(cal_glm, type = "response")]
cal_dat[, decile := cut(pred, breaks = quantile(pred, probs = seq(0, 1, 0.1)),
                        include.lowest = TRUE, labels = FALSE)]
cal <- cal_dat[, .(mean_predicted = mean(pred),
                   observed_highrisk = mean(y_highrisk)), by = decile][order(decile)]

pE <- ggplot(cal, aes(mean_predicted, observed_highrisk)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#BCAAA4") +
  geom_point(size = 2.2, color = "#E76F51") +
  geom_line(color = "#E76F51", linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  scale_y_continuous(limits = c(0, 0.32), breaks = seq(0, 0.3, 0.1)) +
  labs(x = "Mean predicted risk", y = "Observed high-risk rate", tag = "E") +
  theme_pub

# ---- (F) Decision-curve analysis ----
dca <- fread(file.path(tdir, "tab_dca_all_outcomes.csv"))
dca[, outcome := factor(outcome,
  levels = c("High-risk ARG", "MDR-VF fusion", "Conjugative", "Biocide/metal"),
  labels = outcome_levels)]
dca_hr <- dca[outcome == "High-risk ARG"]
pF <- ggplot(dca_hr, aes(threshold)) +
  geom_line(aes(y = NB_model, color = "PlasRisk"), linewidth = 0.7) +
  geom_line(aes(y = NB_treat_all, color = "Treat all"), linetype = "dashed", linewidth = 0.5) +
  geom_line(aes(y = NB_treat_none, color = "Treat none"), linetype = "dotted", linewidth = 0.5) +
  scale_color_manual(values = c("PlasRisk" = "#E76F51", "Treat all" = "#BCAAA4",
                                "Treat none" = "#6D4C41")) +
  scale_x_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "Risk threshold", y = "Net benefit", color = NULL, tag = "F") +
  theme_pub +
  guides(color = guide_legend(nrow = 3))

# ---- assemble ----
fig <- (pA | pB) / (pC | pD) / (pE | pF)

ggsave(file.path(fdir, "Figure2_weight_validation.pdf"), fig,
       width = 7.2, height = 9, units = "in")
ggsave(file.path(fdir, "Figure2_weight_validation.png"), fig,
       width = 7.2, height = 9, units = "in", dpi = 300)
cat("Saved Figure2_weight_validation.pdf/png to", fdir, "\n")
