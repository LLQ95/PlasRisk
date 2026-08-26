#!/usr/bin/env Rscript
# =============================================================================
# Figure 6: High-risk ARG distribution across plasmid replicons
# (A) Top replicons by high-risk ARG carrier rate
# (B) Bubble plot: replicon size vs high-risk ARG rate (labeled with ggrepel)
# (C) Bipartite network: high-risk ARG families x replicons
# (D) Heatmap: high-risk ARG family prevalence across replicons
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

hr_args <- c("NDM", "KPC", "OXA", "VIM", "IMP", "mcr", "CTX-M", "SHV",
             "CMY", "Qnr", "Tet(X4)", "ArmA", "RmtB", "RmtF", "RmtG")

psc <- fread(file.path(tdir, "tab_psc_final_scores.csv"),
             select = c("replicon_primary", "S_final", "y_highrisk", "y_fusion"))
rep_sum <- fread(file.path(tdir, "tab_replicon_summary.csv"))
edges <- fread(file.path(tdir, "tab_arg_replicon_edges.csv"))

# ---- (A) Top replicons by high-risk ARG rate ----
rep_hr <- psc[!is.na(replicon_primary), .(
  n = .N, n_hr = sum(y_highrisk),
  hr_rate = mean(y_highrisk) * 100,
  mean_risk = mean(S_final)
), by = replicon_primary]
rep_hr <- merge(rep_hr, rep_sum[, .(replicon_primary, mean_len, pct_conj)],
                by = "replicon_primary", all.x = TRUE)
rep_hr <- rep_hr[n >= 100][order(-hr_rate)]
rep_top <- rep_hr[1:15]
rep_top[, replicon_primary := factor(replicon_primary,
  levels = rev(rep_top$replicon_primary))]

pA <- ggplot(rep_top, aes(hr_rate, replicon_primary)) +
  geom_col(aes(fill = mean_risk), width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", hr_rate)), hjust = -0.15, size = 2.5) +
  scale_fill_gradient(low = "#FFE8D6", high = "#9C2C2C", name = "Mean\nPlasRisk") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "High-risk ARG carrier rate (%)", y = NULL, tag = "A") +
  theme_pub +
  theme(axis.text.y = element_text(size = 7.5),
        legend.key.height = unit(0.5, "cm"))

# ---- (B) Bubble plot with ggrepel labels ----
rep_bubble <- rep_hr[n >= 100 & !is.na(mean_len) & !is.na(pct_conj)]
rep_bubble[, label := ifelse(n >= 500 | hr_rate > 10, replicon_primary, "")]

pB <- ggplot(rep_bubble, aes(mean_len, hr_rate)) +
  geom_point(aes(size = n, color = pct_conj), alpha = 0.75) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 30,
                  color = "#3E2723", box.padding = 0.4, min.segment.length = 0,
                  segment.color = "#BCAAA4", segment.size = 0.3,
                  max.time = 2, max.iter = 10000) +
  scale_size_continuous(name = "n (PSC)", range = c(1.5, 6)) +
  scale_color_gradient(low = "#FFE8D6", high = "#E76F51", name = "% Conjugative") +
  scale_x_log10() +
  labs(x = "Mean plasmid length (kb, log scale)",
       y = "High-risk ARG rate (%)", tag = "B") +
  theme_pub +
  guides(size = guide_legend(nrow = 2), color = guide_colorbar(barwidth = 4))

# ---- (C) Bipartite network ----
hr_edges <- edges[to %in% hr_args]
top_reps_net <- rep_hr[n >= 200][order(-n)][1:15, replicon_primary]
hr_edges <- hr_edges[from %in% top_reps_net & weight >= 10]

nodes <- unique(rbindlist(list(
  data.table(name = unique(hr_edges$from), type = "Replicon"),
  data.table(name = unique(hr_edges$to), type = "ARG")
)))

g <- graph_from_data_frame(hr_edges[, .(from, to, weight)],
                           directed = FALSE, vertices = nodes)
set.seed(42)

pC <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight, alpha = weight), color = "#BCAAA4") +
  geom_node_point(aes(color = type, size = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2.8, max.overlaps = 30,
                 color = "#3E2723", bg.color = "white", bg.r = 0.08,
                 segment.color = "#BCAAA4", segment.size = 0.2) +
  scale_edge_width(range = c(0.3, 3), guide = "none") +
  scale_edge_alpha(range = c(0.3, 0.8), guide = "none") +
  scale_color_manual(values = c("Replicon" = "#2A9D8F", "ARG" = "#E76F51"),
                     labels = c("ARG family", "Replicon"), name = NULL) +
  scale_size_manual(values = c("Replicon" = 5, "ARG" = 4), guide = "none") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(tag = "C") +
  theme_void(base_size = 9) +
  theme(plot.tag = element_text(face = "bold", size = 12, color = "#3E2723"),
        legend.position = "bottom", legend.text = element_text(size = 8),
        plot.margin = margin(2, 2, 2, 2))

# ---- (D) Heatmap ----
rep_totals <- hr_edges[, .(total = sum(weight)), by = from]
hr_heat <- merge(hr_edges, rep_totals, by = "from")
hr_heat[, pct := weight / total * 100]
hr_heat[, from := factor(from, levels = top_reps_net)]
hr_heat[, to := factor(to, levels = sort(unique(to)))]

pD <- ggplot(hr_heat, aes(to, from, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pct > 5, sprintf("%.0f", pct), "")),
            size = 2.2, color = "#3E2723") +
  scale_fill_gradient(low = "#FFF8E7", high = "#9C2C2C", name = "% of ARG\ncarriage") +
  labs(x = NULL, y = NULL, tag = "D") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5),
        axis.text.y = element_text(size = 7.5),
        legend.key.height = unit(0.7, "cm"))

# ---- assemble ----
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.1))

ggsave(file.path(fdir, "Figure6_highrisk_arg_replicon.pdf"), fig,
       width = 7.5, height = 8.5, units = "in")
ggsave(file.path(fdir, "Figure6_highrisk_arg_replicon.png"), fig,
       width = 7.5, height = 8.5, units = "in", dpi = 300)
cat("Saved Figure6_highrisk_arg_replicon.pdf/png to", fdir, "\n")
