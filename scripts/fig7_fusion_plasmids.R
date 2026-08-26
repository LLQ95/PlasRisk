#!/usr/bin/env Rscript
# =============================================================================
# Figure 7: MDR-VF fusion plasmid characteristics
# (A) Bubble plot: fusion rate by replicon (labeled with ggrepel)
# (B) Fusion plasmid feature profile
# (C) Replicon distribution of fusion plasmids
# (D) ARG vs VF burden across replicons (fusion rate as color)
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

fus_rep <- fread(file.path(tdir, "tab_fusion_by_replicon.csv"))
fus_sum <- fread(file.path(tdir, "tab_fusion_plasmid_summary.csv"))
rep_vr  <- fread(file.path(tdir, "tab_replicon_virulence_resistance.csv"))

# ---- (A) Bubble plot: fusion rate by replicon ----
fus_plot <- fus_rep[order(-N)][1:15]
fus_plot[, label := replicon_primary]
fus_plot[, replicon_primary := factor(replicon_primary,
  levels = rev(fus_plot$replicon_primary))]

pA <- ggplot(fus_plot, aes(fus_rate, replicon_primary)) +
  geom_point(aes(size = N, color = pct), alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", fus_rate)),
            hjust = -0.6, size = 2.5, color = "#3E2723") +
  scale_size_continuous(name = "Fusion\nplasmids (N)", range = c(1.5, 6)) +
  scale_color_gradient(low = "#FFE8D6", high = "#E76F51",
                       name = "% of all\nfusion plasmids") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(x = "Fusion rate within replicon (%)", y = NULL, tag = "A") +
  theme_pub +
  theme(axis.text.y = element_text(size = 7.5),
        legend.box = "vertical",
        legend.key.height = unit(0.4, "cm"))

# ---- (B) Fusion plasmid feature profile ----
feat <- data.table(
  Feature = c("Mean length (kb)", "Mean ARG count", "Mean VF count",
              "Mean IS count", "Conjugative (%)", "Integron (%)", "Human host (%)"),
  Value   = c(fus_sum$mean_len, fus_sum$mean_arg, fus_sum$mean_vf,
              fus_sum$mean_is, fus_sum$pct_conj, fus_sum$pct_integron, fus_sum$pct_human)
)
feat[, Feature := factor(Feature, levels = rev(Feature))]
feat_group <- ifelse(feat$Feature %in% c("Conjugative (%)", "Integron (%)", "Human host (%)"),
                     "Percentage", "Count/Mean")

pB <- ggplot(feat, aes(Value, Feature)) +
  geom_col(aes(fill = feat_group), width = 0.65, color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", Value)), hjust = -0.1, size = 2.5) +
  scale_fill_manual(values = c("Count/Mean" = "#F4A261", "Percentage" = "#2A9D8F"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Value", y = NULL, tag = "B") +
  theme_pub +
  theme(axis.text.y = element_text(size = 8),
        legend.key.height = unit(0.4, "cm"))

# ---- (C) Replicon distribution ----
fus_dist <- fus_rep[order(-N)][1:12]
fus_dist[, replicon_primary := factor(replicon_primary,
  levels = rev(fus_dist$replicon_primary))]

pC <- ggplot(fus_dist, aes(N, replicon_primary)) +
  geom_col(fill = "#E76F51", width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(aes(label = N), hjust = -0.3, size = 2.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of fusion plasmids", y = NULL, tag = "C") +
  theme_pub +
  theme(axis.text.y = element_text(size = 7.5))

# ---- (D) ARG vs VF burden with ggrepel ----
rep_vr_f <- rep_vr[n >= 100]
rep_vr_f[, label := ifelse(pct_mdr_vf > 5 | n >= 2000, replicon_primary, "")]

pD <- ggplot(rep_vr_f, aes(mean_arg, mean_vf)) +
  geom_point(aes(size = n, color = pct_mdr_vf), alpha = 0.75) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15,
                  color = "#3E2723", box.padding = 0.35, min.segment.length = 0,
                  segment.color = "#BCAAA4", segment.size = 0.3) +
  scale_size_continuous(name = "n (PSC)", range = c(1.5, 6)) +
  scale_color_gradient(low = "#FFE8D6", high = "#9C2C2C",
                       name = "MDR-VF\nfusion (%)") +
  labs(x = "Mean ARG count per PSC", y = "Mean VF count per PSC", tag = "D") +
  theme_pub +
  guides(size = guide_legend(nrow = 2), color = guide_colorbar(barwidth = 4))

# ---- assemble ----
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1.1, 1))

ggsave(file.path(fdir, "Figure7_fusion_plasmids.pdf"), fig,
       width = 7.5, height = 8, units = "in")
ggsave(file.path(fdir, "Figure7_fusion_plasmids.png"), fig,
       width = 7.5, height = 8, units = "in", dpi = 300)
cat("Saved Figure7_fusion_plasmids.pdf/png to", fdir, "\n")
