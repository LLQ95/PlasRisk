#!/usr/bin/env Rscript
# =============================================================================
# Figure 4: Dimensionality analysis and model parsimony
# (A) AUC vs k for all 1023 subsets + best per k
# (B) Forward stepwise selection path
# (C) Mean CV AUC vs k with full/lite highlights
# (D) Pareto frontier (high-risk ARG vs biocide/metal AUC)
# (E) Weight distribution: 5-dim lite vs 10-dim full
# (F) Backward elimination path
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
tdir   <- file.path(resdir, "tables")
fdir   <- file.path(resdir, "figures")
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)

outcome_colors <- c(
  "High-risk ARG" = "#E76F51",
  "MDR-VF fusion" = "#F4A261",
  "Conjugative"   = "#2A9D8F",
  "Biocide/metal" = "#6D4C41",
  "Mean"          = "#3E2723"
)

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

# ---- (A) All subsets jitter + best line ----
allsub <- fread(file.path(tdir, "tab_all_subsets_auc.csv"))
bestk  <- fread(file.path(tdir, "tab_best_subset_by_k.csv"))

all_long <- melt(allsub,
  id.vars = c("mask", "k", "subset"),
  measure.vars = c("test_AUC_highrisk", "test_AUC_fusion",
                   "test_AUC_conj", "test_AUC_bm"),
  variable.name = "outcome", value.name = "AUC")
all_long[, outcome := gsub("test_AUC_", "", outcome)]
all_long[, outcome := factor(outcome,
  levels = c("highrisk", "fusion", "conj", "bm"),
  labels = c("High-risk ARG", "MDR-VF fusion", "Conjugative", "Biocide/metal"))]

best_long <- melt(bestk, id.vars = "k",
  measure.vars = c("test_AUC_highrisk", "test_AUC_fusion",
                   "test_AUC_conj", "test_AUC_bm"),
  variable.name = "outcome", value.name = "AUC")
best_long[, outcome := gsub("test_AUC_", "", outcome)]
best_long[, outcome := factor(outcome,
  levels = c("highrisk", "fusion", "conj", "bm"),
  labels = c("High-risk ARG", "MDR-VF fusion", "Conjugative", "Biocide/metal"))]

pA <- ggplot(all_long, aes(factor(k), AUC)) +
  geom_jitter(aes(color = outcome), width = 0.25, alpha = 0.12, size = 0.4) +
  geom_line(data = best_long, aes(group = outcome, color = outcome), linewidth = 0.7) +
  geom_point(data = best_long, aes(color = outcome), size = 1.5) +
  scale_color_manual(values = outcome_colors) +
  coord_cartesian(ylim = c(0.5, 1.0)) +
  labs(x = "Number of dimensions (k)", y = "Test AUC", color = NULL, tag = "A") +
  theme_pub +
  guides(color = guide_legend(nrow = 2))

# ---- (B) Forward selection ----
fwd <- fread(file.path(tdir, "tab_forward_selection.csv"))
fwd_long <- melt(fwd, id.vars = c("step", "added", "k"),
  measure.vars = c("cv_AUC_highrisk", "cv_AUC_fusion",
                   "cv_AUC_conj", "cv_AUC_bm", "cv_AUC_mean"),
  variable.name = "outcome", value.name = "AUC")
fwd_long[, outcome := gsub("cv_AUC_", "", outcome)]
fwd_long[, outcome := factor(outcome,
  levels = c("highrisk", "fusion", "conj", "bm", "mean"),
  labels = c("High-risk ARG", "MDR-VF fusion", "Conjugative", "Biocide/metal", "Mean"))]

pB <- ggplot(fwd_long, aes(k, AUC, color = outcome)) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "#BCAAA4", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  scale_color_manual(values = outcome_colors) +
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim = c(0.55, 1.0)) +
  labs(x = "Step (k)", y = "CV AUC", color = NULL, tag = "B") +
  theme_pub +
  guides(color = guide_legend(nrow = 2))

# ---- (C) Mean CV AUC vs k ----
pC <- ggplot(bestk, aes(k, cv_AUC_mean)) +
  geom_ribbon(aes(ymin = cv_AUC_mean - 0.005, ymax = cv_AUC_mean + 0.005),
              fill = "#F4A261", alpha = 0.2) +
  geom_line(linewidth = 0.6, color = "#6D4C41") +
  geom_point(size = 1.5, color = "#6D4C41") +
  annotate("point", x = 10, y = bestk$cv_AUC_mean[bestk$k == 10],
           color = "#9C2C2C", size = 3.5, shape = 18) +
  annotate("point", x = 5, y = bestk$cv_AUC_mean[bestk$k == 5],
           color = "#2A9D8F", size = 3.5, shape = 18) +
  annotate("text", x = 9.2, y = bestk$cv_AUC_mean[bestk$k == 10] + 0.008,
           label = "Full (10D)", size = 2.8, color = "#9C2C2C") +
  annotate("text", x = 5.8, y = bestk$cv_AUC_mean[bestk$k == 5] - 0.012,
           label = "Lite (5D)", size = 2.8, color = "#2A9D8F") +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(limits = c(0.74, 0.94), breaks = seq(0.75, 0.95, 0.05)) +
  labs(x = "Number of dimensions (k)", y = "Mean CV AUC", tag = "C") +
  theme_pub

# ---- (D) Pareto frontier ----
pareto <- fread(file.path(tdir, "tab_pareto_optimal_subsets.csv"))
pareto[, k_fac := factor(k)]
full_k10 <- pareto[k == 10]
lite_k5  <- pareto[k == 5]
no_bm    <- pareto[k == 9 & has_BM == 0]
if (nrow(no_bm) == 0) {
  allsub9 <- allsub[k == 9 & has_BM == 0]
  if (nrow(allsub9) > 0) no_bm <- allsub9[which.max(test_AUC_mean)]
}

pD <- ggplot(pareto, aes(test_AUC_highrisk, test_AUC_bm)) +
  geom_point(aes(color = k), size = 1.5, alpha = 0.7) +
  geom_point(data = full_k10, color = "#9C2C2C", size = 4, shape = 18) +
  geom_point(data = lite_k5, color = "#2A9D8F", size = 4, shape = 18) +
  {if (nrow(no_bm) > 0) geom_point(data = no_bm, color = "#6D4C41", size = 2.5, shape = 1)} +
  annotate("text", x = full_k10$test_AUC_highrisk - 0.002,
           y = full_k10$test_AUC_bm + 0.006, label = "Full", size = 2.8, color = "#9C2C2C") +
  annotate("text", x = lite_k5$test_AUC_highrisk + 0.002,
           y = lite_k5$test_AUC_bm - 0.008, label = "Lite", size = 2.8, color = "#2A9D8F") +
  scale_color_gradient(low = "#FFE8D6", high = "#E76F51", name = "k") +
  coord_cartesian(xlim = c(0.93, 0.995), ylim = c(0.84, 0.93)) +
  labs(x = "AUC: High-risk ARG", y = "AUC: Biocide/metal", tag = "D") +
  theme_pub

# ---- (E) Weight distribution comparison ----
wfull <- fread(file.path(tdir, "tab_weight_comparison.csv"))
lite_dims <- c("S_ARG", "S_VF", "S_MOB", "S_SIZE", "S_BM")
wlite <- wfull[Component %in% lite_dims, .(Component, Weight = Final)]
wlite[, Weight := Weight / sum(Weight)]
wlite[, Model := "Lite (5D)"]
wfull_long <- wfull[, .(Component, Weight = Final, Model = "Full (10D)")]
wcomp <- rbindlist(list(wfull_long, wlite), fill = TRUE)
wcomp[, Component := factor(Component, levels = wfull$Component)]

pE <- ggplot(wcomp, aes(Component, Weight, fill = Model)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("Full (10D)" = "#9C2C2C", "Lite (5D)" = "#2A9D8F")) +
  labs(x = NULL, y = "Weight", fill = NULL, tag = "E") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))

# ---- (F) Backward elimination ----
bwd <- fread(file.path(tdir, "tab_backward_elimination.csv"))
bwd_long <- melt(bwd, id.vars = c("step", "removed", "k"),
  measure.vars = c("cv_AUC_highrisk", "cv_AUC_fusion",
                   "cv_AUC_conj", "cv_AUC_bm", "cv_AUC_mean"),
  variable.name = "outcome", value.name = "AUC")
bwd_long[, outcome := gsub("cv_AUC_", "", outcome)]
bwd_long[, outcome := factor(outcome,
  levels = c("highrisk", "fusion", "conj", "bm", "mean"),
  labels = c("High-risk ARG", "MDR-VF fusion", "Conjugative", "Biocide/metal", "Mean"))]

pF <- ggplot(bwd_long, aes(k, AUC, color = outcome)) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "#BCAAA4", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  scale_color_manual(values = outcome_colors) +
  scale_x_reverse(breaks = 1:10) +
  coord_cartesian(ylim = c(0.55, 1.0)) +
  labs(x = "Dimensions remaining (k)", y = "CV AUC", color = NULL, tag = "F") +
  theme_pub +
  guides(color = guide_legend(nrow = 2))

# ---- assemble ----
fig <- (pA | pB) / (pC | pD) / (pE | pF)

ggsave(file.path(fdir, "Figure4_dimensionality.pdf"), fig,
       width = 7.2, height = 9, units = "in")
ggsave(file.path(fdir, "Figure4_dimensionality.png"), fig,
       width = 7.2, height = 9, units = "in", dpi = 300)
cat("Saved Figure4_dimensionality.pdf/png to", fdir, "\n")
