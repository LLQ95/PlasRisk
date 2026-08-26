#!/usr/bin/env Rscript
# =============================================================================
# Figure 3: Biocide/metal resistance-ARG-replicon co-occurrence network
# (A) Proportion of ARG-carrying plasmids co-carrying each BMG category
# (B) Tripartite ARG-BMG-Replicon co-occurrence network (top 12 replicons only)
# (C) Odds ratio heatmap (QacEdelta excluded)
# (D) BMG-ARG co-carriage vs high-risk ARG rate by replicon (bubble plot)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(igraph)
  library(ggraph)
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

bm_colors <- c(
  "QAC/Disinfectant" = "#E76F51",
  "Zn/Cd/Ni"         = "#F4A261",
  "Mercury (mer)"    = "#9C6C4A",
  "Cu/Ag"            = "#2A9D8F",
  "Arsenic (ars)"    = "#6D4C41"
)

# ---- (A) Proportion bar chart ----
bm_rates <- data.table(
  bm_category = factor(c("QAC/Disinfectant", "Zn/Cd/Ni", "Mercury (mer)", "Cu/Ag", "Arsenic (ars)"),
    levels = c("QAC/Disinfectant", "Zn/Cd/Ni", "Mercury (mer)", "Cu/Ag", "Arsenic (ars)")),
  pct = c(9.3, 7.3, 3.2, 2.5, 2.0)
)

pA <- ggplot(bm_rates, aes(bm_category, pct, fill = bm_category)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 2.8) +
  scale_fill_manual(values = bm_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "% of ARG+ plasmids", tag = "A") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5))

# ---- (B) Tripartite ARG-BMG-Replicon network ----
cooc <- fread(file.path(tdir, "tab_bm_arg_cooccurrence.csv"))
cooc_filt <- cooc[arg_family != "QacEdelta"]
sig <- cooc_filt[p_adj < 0.001 & odds_ratio > 2 & n_both >= 20]

rep_bm <- fread(file.path(tdir, "tab_bm_arg_by_replicon.csv"))
# Keep only major replicons with substantial BMG co-occurrence to avoid clutter
rep_bm_filt <- rep_bm[n >= 1000 & (n_qac >= 50 | n_mer >= 50)][
  order(-(n_qac + n_mer))][1:12]
rep_bm_edges <- rbindlist(list(
  rep_bm_filt[n_qac >= 10, .(from = replicon_primary, to = "QAC/Disinfectant",
                              weight = n_qac, edge_type = "Replicon-BMG")],
  rep_bm_filt[n_mer >= 10, .(from = replicon_primary, to = "Mercury (mer)",
                              weight = n_mer, edge_type = "Replicon-BMG")]
))
arg_bm_edges <- sig[, .(from = arg_family, to = bm_category,
                         weight = n_both, edge_type = "ARG-BMG")]
all_edges <- rbindlist(list(arg_bm_edges, rep_bm_edges), fill = TRUE)

nodes <- unique(rbindlist(list(
  data.table(name = unique(arg_bm_edges$from), type = "ARG"),
  data.table(name = unique(all_edges$to[all_edges$to %in% names(bm_colors)]), type = "BMG"),
  data.table(name = unique(rep_bm_edges$from), type = "Replicon")
)))

g <- graph_from_data_frame(all_edges[, .(from, to, weight, edge_type)],
                           directed = FALSE, vertices = nodes)
set.seed(42)

pB <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight, alpha = edge_type, color = edge_type)) +
  geom_node_point(aes(color = type, size = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2.5, max.overlaps = 30,
                 color = "#3E2723", bg.color = "white", bg.r = 0.08,
                 segment.color = "#BCAAA4", segment.size = 0.2) +
  scale_edge_width(range = c(0.2, 2.5), guide = "none") +
  scale_edge_alpha_manual(values = c("ARG-BMG" = 0.5, "Replicon-BMG" = 0.2), guide = "none") +
  scale_edge_color_manual(values = c("ARG-BMG" = "#E76F51", "Replicon-BMG" = "#2A9D8F"),
                          guide = "none") +
  scale_color_manual(values = c("ARG" = "#E76F51", "BMG" = "#F4A261", "Replicon" = "#2A9D8F"),
                     labels = c("ARG family", "BMG category", "Replicon"), name = NULL) +
  scale_size_manual(values = c("ARG" = 3.5, "BMG" = 6, "Replicon" = 4.5), guide = "none") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(tag = "B") +
  theme_void(base_size = 9) +
  theme(plot.tag = element_text(face = "bold", size = 12, color = "#3E2723"),
        legend.position = "bottom", legend.text = element_text(size = 8),
        plot.margin = margin(2, 2, 2, 2))

# ---- (C) OR heatmap ----
heat <- cooc_filt[p_adj < 0.001]
heat[, OR_show := pmin(odds_ratio, 15)]
heat[, arg_family := factor(arg_family, levels = rev(sort(unique(arg_family))))]
heat[, bm_category := factor(bm_category,
  levels = c("QAC/Disinfectant", "Zn/Cd/Ni", "Mercury (mer)", "Cu/Ag", "Arsenic (ars)"))]

pC <- ggplot(heat, aes(bm_category, arg_family, fill = OR_show)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(OR_show > 1, sprintf("%.1f", OR_show), "")),
            size = 2.2, color = "#3E2723") +
  scale_fill_gradient(low = "#FFF8E7", high = "#9C2C2C",
                      limits = c(0, 15), na.value = "#FAF0E6",
                      name = "OR (cap 15)") +
  labs(x = NULL, y = NULL, tag = "C") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5),
        axis.text.y = element_text(size = 7.5),
        legend.key.height = unit(0.7, "cm"))

# ---- (D) BMG-ARG co-carriage vs high-risk ARG rate by replicon ----
psc_all <- fread(file.path(tdir, "tab_psc_final_scores.csv"),
                 select = c("replicon_primary", "y_highrisk"))
rep_hr <- psc_all[!is.na(replicon_primary), .(
  n_hr = .N, hr_rate = mean(y_highrisk) * 100
), by = replicon_primary]
bubble <- merge(rep_bm[n >= 50], rep_hr, by = "replicon_primary")
top_hr <- bubble[order(-hr_rate)][1:8, replicon_primary]
bubble[, label := ifelse(replicon_primary %in% top_hr, replicon_primary, "")]

pD <- ggplot(bubble, aes(pct_arg_bm, hr_rate)) +
  geom_point(aes(size = n, color = pct_bm), alpha = 0.75) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  color = "#3E2723", box.padding = 0.3, min.segment.length = 0,
                  segment.color = "#BCAAA4", segment.size = 0.3) +
  scale_size_continuous(name = "n (PSC)", range = c(1.5, 6)) +
  scale_color_gradient(low = "#FFE8D6", high = "#E76F51", name = "% BMG+") +
  labs(x = "BMG-ARG co-carriage rate among ARG+ plasmids (%)",
       y = "High-risk ARG carrier rate (%)", tag = "D") +
  theme_pub +
  guides(size = guide_legend(nrow = 2), color = guide_colorbar(barwidth = 4))

# ---- assemble ----
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.1))

ggsave(file.path(fdir, "Figure3_bm_arg_network.pdf"), fig,
       width = 7.5, height = 8.5, units = "in")
ggsave(file.path(fdir, "Figure3_bm_arg_network.png"), fig,
       width = 7.5, height = 8.5, units = "in", dpi = 300)
cat("Saved Figure3_bm_arg_network.pdf/png to", fdir, "\n")
