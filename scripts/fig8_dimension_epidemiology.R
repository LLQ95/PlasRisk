#!/usr/bin/env Rscript
# =============================================================================
# Figure 8: Genomic epidemiology of the four core risk dimensions
# (A) Heatmap of mean dimension scores across major replicons
# (B) Temporal trends of ARG, BMG, and co-occurrence
# (C) Dimension correlation across replicons
# (D) Replicon risk landscape (S_ARG vs S_BM, labeled with ggrepel)
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

rep_risk <- fread(file.path(tdir, "tab_replicon_risk_10dimensions.csv"))
bm_trend <- fread(file.path(tdir, "tab_temporal_bm_arg_trend.csv"))

# ---- (A) Heatmap of dimension scores across replicons ----
rep_top <- rep_risk[n_PSC >= 200][order(-S10_norm)][1:20]
heat_dims <- c("mean_S_ARG", "mean_S_VF", "mean_S_MOB", "mean_S_BM")
dim_labels <- c(S_ARG = "ARG", S_VF = "VF", S_MOB = "MOB", S_BM = "BM")

heat_long <- melt(rep_top, id.vars = "replicon_primary",
  measure.vars = heat_dims, variable.name = "Dimension", value.name = "Score")
heat_long[, Dimension := factor(Dimension, levels = heat_dims, labels = dim_labels)]
heat_long[, replicon_primary := factor(replicon_primary,
  levels = rev(rep_top$replicon_primary))]

pA <- ggplot(heat_long, aes(Dimension, replicon_primary, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#FFF8E7", high = "#9C2C2C",
                      limits = c(0, 1), name = "Mean\nscore") +
  labs(x = NULL, y = NULL, fill = "Mean\nscore", tag = "A") +
  theme_pub +
  theme(axis.text.y = element_text(size = 7),
        legend.key.height = unit(0.7, "cm"))

# ---- (B) Temporal trends ----
bm_trend[, period_start := as.numeric(sub("-.*", "", period))]
trend_long <- melt(bm_trend, id.vars = c("period", "period_start", "n"),
  measure.vars = c("pct_arg", "pct_bm", "pct_both"),
  variable.name = "Metric", value.name = "Percent")
trend_long[, Metric := factor(Metric,
  levels = c("pct_arg", "pct_bm", "pct_both"),
  labels = c("ARG+", "BMG+", "ARG+ & BMG+"))]

pB <- ggplot(trend_long, aes(period_start, Percent, color = Metric, group = Metric)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("ARG+" = "#E76F51", "BMG+" = "#2A9D8F",
                                "ARG+ & BMG+" = "#6D4C41")) +
  scale_x_continuous(breaks = seq(1970, 2020, 10)) +
  labs(x = "Year period", y = "Percent of PSCs (%)", color = NULL, tag = "B") +
  theme_pub +
  guides(color = guide_legend(nrow = 2))

# ---- (C) Dimension correlation ----
rep_corr <- rep_risk[n_PSC >= 50, .(S_ARG = mean_S_ARG, S_VF = mean_S_VF,
                                     S_MOB = mean_S_MOB, S_BM = mean_S_BM,
                                     S_SIZE = mean_S_SIZE, S_HOST = mean_S_HOST)]
cor_mat <- cor(rep_corr[, .(S_ARG, S_VF, S_MOB, S_BM, S_SIZE, S_HOST)],
               use = "pairwise.complete.obs")
cor_dt <- as.data.table(as.table(cor_mat))
setnames(cor_dt, c("V1", "V2", "r"))
cor_dt[, V1 := factor(V1, levels = rev(colnames(cor_mat)))]
cor_dt[, V2 := factor(V2, levels = colnames(cor_mat))]

pC <- ggplot(cor_dt, aes(V2, V1, fill = r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.8, color = "#3E2723") +
  scale_fill_gradient2(low = "#2A9D8F", mid = "white", high = "#E76F51",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "Pearson r", tag = "C") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        legend.key.height = unit(0.7, "cm"))

# ---- (D) Replicon risk landscape with ggrepel ----
rep_land <- rep_risk[n_PSC >= 100]
rep_land[, label := ifelse(n_PSC >= 1000 | mean_S_ARG > 0.5, replicon_primary, "")]

pD <- ggplot(rep_land, aes(mean_S_ARG, mean_S_BM)) +
  geom_point(aes(size = n_PSC, color = mean_S_MOB), alpha = 0.75) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 30,
                  color = "#3E2723", box.padding = 0.35, min.segment.length = 0,
                  segment.color = "#BCAAA4", segment.size = 0.3) +
  scale_size_continuous(name = "n (PSC)", range = c(1.5, 6)) +
  scale_color_gradient(low = "#FFE8D6", high = "#E76F51", name = "Mean\nS_MOB") +
  labs(x = "Mean S_ARG", y = "Mean S_BM", tag = "D") +
  theme_pub +
  guides(size = guide_legend(nrow = 2), color = guide_colorbar(barwidth = 4))

# ---- assemble ----
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1.2, 1))

ggsave(file.path(fdir, "Figure8_dimension_epidemiology.pdf"), fig,
       width = 7.5, height = 8.5, units = "in")
ggsave(file.path(fdir, "Figure8_dimension_epidemiology.png"), fig,
       width = 7.5, height = 8.5, units = "in", dpi = 300)
cat("Saved Figure8_dimension_epidemiology.pdf/png to", fdir, "\n")
