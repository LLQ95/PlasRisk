#!/usr/bin/env Rscript
# =============================================================================
# Figure 4: Dimensionality analysis
# (A) All-subsets mean AUC vs number of dimensions
# (B) Forward/backward selection paths
# (C) Full vs lite per-outcome AUC comparison
# (D) S_BM ablation: AUC drop for biocide/metal outcome
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
tdir   <- file.path(resdir, "tables")
fdir   <- file.path(resdir, "figures")
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)

warm <- c("#E76F51", "#F4A261", "#E9C46A", "#2A9D8F", "#6D4C41")
outcome_colors <- c("High-risk ARG" = "#E76F51", "MDR-VF fusion" = "#F4A261",
                    "Conjugative" = "#2A9D8F", "Biocide/metal" = "#6D4C41")

theme_pub <- theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 12, color = "#3E2723"),
        axis.text = element_text(color = "#3E2723"),
        plot.margin = margin(4, 6, 4, 4))

# ---- (A) all-subsets ----
sub_file <- file.path(tdir, "tab_allsubsets.csv")
pA <- ggplot() + theme_void() + labs(tag = "A")
if (file.exists(sub_file)) {
  sub <- fread(sub_file)
  kcol <- intersect(c("k", "n_dim", "ndim"), names(sub))[1]
  acol <- intersect(c("mean_auc", "auc", "mean_AUC"), names(sub))[1]
  if (!is.na(kcol) && !is.na(acol)) {
    setnames(sub, c(kcol, acol), c("k", "auc"))
    best <- sub[, .(best = max(auc)), by = k]
    pA <- ggplot(sub, aes(k, auc)) +
      geom_jitter(width = 0.15, alpha = 0.25, color = "#F4A261", size = 0.8) +
      geom_line(data = best, aes(k, best), color = "#9C2C2C", linewidth = 0.7) +
      geom_point(data = best, aes(k, best), color = "#9C2C2C", size = 1.6) +
      geom_vline(xintercept = 5, linetype = "dashed", color = "#2A9D8F") +
      annotate("text", 5.2, min(sub$auc), label = "5-dim lite",
               hjust = 0, color = "#2A9D8F", size = 3) +
      scale_x_continuous(breaks = 1:10) +
      labs(x = "Number of dimensions", y = "Mean CV AUC", tag = "A") +
      theme_pub
  }
}

# ---- (B) selection paths ----
sel_file <- file.path(tdir, "tab_selection_paths.csv")
pB <- ggplot() + theme_void() + labs(tag = "B")
if (file.exists(sel_file)) {
  sel <- fread(sel_file)
  if (all(c("step", "auc", "method") %in% names(sel))) {
    pB <- ggplot(sel, aes(step, auc, color = method)) +
      geom_line(linewidth = 0.7) + geom_point(size = 1.5) +
      scale_color_manual(values = c("forward" = "#E76F51", "backward" = "#2A9D8F")) +
      scale_x_continuous(breaks = 1:10) +
      labs(x = "Step", y = "Mean CV AUC", color = NULL, tag = "B") +
      theme_pub
  }
}

# ---- (C) full vs lite per outcome ----
cmp_file <- file.path(tdir, "tab_full_vs_lite.csv")
pC <- ggplot() + theme_void() + labs(tag = "C")
if (file.exists(cmp_file)) {
  cmp <- fread(cmp_file)
  if (all(c("outcome", "full_auc", "lite_auc") %in% names(cmp))) {
    cl <- melt(cmp, id.vars = "outcome", measure.vars = c("full_auc", "lite_auc"),
               variable.name = "model", value.name = "auc")
    cl[, model := factor(model, levels = c("full_auc", "lite_auc"),
                         labels = c("10-dim full", "5-dim lite"))]
    pC <- ggplot(cl, aes(outcome, auc, fill = model)) +
      geom_col(position = position_dodge(0.8), width = 0.7) +
      scale_fill_manual(values = c("10-dim full" = "#9C2C2C", "5-dim lite" = "#F4A261")) +
      coord_cartesian(ylim = c(0.8, 1.0)) +
      labs(x = NULL, y = "AUC", fill = NULL, tag = "C") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8))
  }
}

# ---- (D) S_BM ablation ----
abl_file <- file.path(tdir, "tab_SBM_weight_sensitivity.csv")
pD <- ggplot() + theme_void() + labs(tag = "D")
if (file.exists(abl_file)) {
  abl <- fread(abl_file)
  if (all(c("S_BM_weight", "AUC_bm") %in% names(abl))) {
    pD <- ggplot(abl, aes(S_BM_weight, AUC_bm)) +
      geom_line(color = "#6D4C41", linewidth = 0.7) +
      geom_point(color = "#E76F51", size = 1.5) +
      geom_vline(xintercept = 0.211, linetype = "dashed", color = "#2A9D8F") +
      annotate("text", 0.211, min(abl$AUC_bm), label = "final 0.211",
               hjust = -0.1, size = 3, color = "#2A9D8F") +
      labs(x = expression("S"[BM]*" weight"), y = "AUC (biocide/metal outcome)", tag = "D") +
      theme_pub
  }
}

fig <- (pA | pB) / (pC | pD)
ggsave(file.path(fdir, "Figure4_dimensionality.pdf"), fig, width = 7.2, height = 7,
       units = "in")
ggsave(file.path(fdir, "Figure4_dimensionality.png"), fig, width = 7.2, height = 7,
       units = "in", dpi = 300)
cat("Saved Figure4_dimensionality.pdf/png to", fdir, "\n")
