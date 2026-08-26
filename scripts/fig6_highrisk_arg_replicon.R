#!/usr/bin/env Rscript
# =============================================================================
# Figure 6: High-risk ARG carriage across plasmid replicons
# (A) Top replicons by high-risk ARG carrier rate (n >= 50)
# (B) High-risk rate vs mean plasmid length (bubble size = PSC count)
# (C) ARG family - replicon bipartite network (top associations)
# (D) ARG family carriage heatmap across top replicons
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

has_igraph <- requireNamespace("igraph", quietly = TRUE)
has_ggraph <- requireNamespace("ggraph", quietly = TRUE)
if (has_igraph) library(igraph)
if (has_ggraph) library(ggraph)

args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
tdir   <- file.path(resdir, "tables")
fdir   <- file.path(resdir, "figures")
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)

theme_pub <- theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 12, color = "#3E2723"),
        axis.text = element_text(color = "#3E2723"),
        plot.margin = margin(4, 6, 4, 4))

rep <- fread(file.path(tdir, "tab_replicon_summary.csv"))
hrcol  <- intersect(c("highrisk_rate", "hr_rate", "high_risk_rate"), names(rep))[1]
ncol_  <- intersect(c("n", "N", "n_psc"), names(rep))[1]
namecol<- intersect(c("replicon_primary", "replicon", "Replicon"), names(rep))[1]
lencol <- intersect(c("mean_length", "avg_length", "mean_len"), names(rep))[1]
if (is.na(hrcol) || is.na(ncol_) || is.na(namecol)) stop("required columns missing in tab_replicon_summary")
setnames(rep, c(hrcol, ncol_, namecol), c("hr_rate", "n", "replicon"))
if (is.na(lencol)) rep[, mean_length := NA_real_] else setnames(rep, lencol, "mean_length")

# ---- (A) top replicons by high-risk rate ----
top <- rep[n >= 50][order(-hr_rate)][1:18]
top[, replicon := factor(replicon, levels = rev(replicon))]
pA <- ggplot(top, aes(hr_rate * 100, replicon)) +
  geom_col(fill = "#E76F51", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", hr_rate * 100)), hjust = -0.1, size = 2.6) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "High-risk ARG carrier rate (%)", y = NULL, tag = "A") +
  theme_pub + theme(axis.text.y = element_text(size = 7.5))

# ---- (B) rate vs length bubble ----
bub <- rep[n >= 50 & !is.na(mean_length)]
bub[, lab := ifelse(hr_rate > 0.25 | (mean_length < 15000 & hr_rate > 0.3),
                    as.character(replicon), "")]
pB <- ggplot(bub, aes(mean_length / 1000, hr_rate * 100, size = n)) +
  geom_point(alpha = 0.55, color = "#E76F51") +
  geom_text_repel(aes(label = lab), size = 2.4, max.overlaps = 15,
                  color = "#3E2723") +
  scale_size_continuous(range = c(1, 7), name = "PSC count") +
  scale_x_log10() +
  labs(x = expression("Mean plasmid length (kb, "*log[10]*")"),
       y = "High-risk ARG carrier rate (%)", tag = "B") +
  theme_pub

# ---- (C) ARG family - replicon network ----
net_file <- file.path(tdir, "tab_arg_replicon_assoc.csv")
pC <- ggplot() + theme_void() + labs(tag = "C")
if (file.exists(net_file) && has_igraph && has_ggraph) {
  net <- fread(net_file)
  ac <- intersect(c("ARG_family", "arg_family", "ARG"), names(net))[1]
  rc <- intersect(c("replicon", "replicon_primary"), names(net))[1]
  if (!is.na(ac) && !is.na(rc)) {
    setnames(net, c(ac, rc), c("ARG", "REP"))
    d <- net[1:60]
    g <- graph_from_data_frame(d[, .(from = ARG, to = REP)], directed = FALSE)
    V(g)$kind <- ifelse(V(g)$name %in% d$REP, "Replicon", "ARG family")
    set.seed(42)
    pC <- ggraph(g, layout = "fr") +
      geom_edge_link(alpha = 0.3, color = "#F4A261") +
      geom_node_point(aes(color = kind), size = 3) +
      geom_node_text(aes(label = name), repel = TRUE, size = 2.3, max.overlaps = 30) +
      scale_color_manual(values = c("ARG family" = "#E76F51", "Replicon" = "#2A9D8F")) +
      labs(color = NULL, tag = "C") +
      theme_void() + theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 12))
  }
}

# ---- (D) ARG family carriage heatmap ----
heat_file <- file.path(tdir, "tab_replicon_arg_family.csv")
pD <- ggplot() + theme_void() + labs(tag = "D")
if (file.exists(heat_file)) {
  h <- fread(heat_file)
  if (all(c("replicon", "ARG_family", "rate") %in% names(h))) {
    reps <- rep[n >= 50][order(-hr_rate)][1:10]$replicon
    hd <- h[replicon %in% reps]
    hd[, replicon := factor(replicon, levels = rev(reps))]
    pD <- ggplot(hd, aes(ARG_family, replicon, fill = rate * 100)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#FFF8E7", mid = "#F4A261", high = "#9C2C2C",
                           midpoint = 50, name = "Carriage\nrate (%)") +
      labs(x = NULL, y = NULL, tag = "D") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5),
            axis.text.y = element_text(size = 7.5))
  }
}

fig <- (pA | pB) / (pC | pD)
ggsave(file.path(fdir, "Figure6_highrisk_arg_replicon.pdf"), fig, width = 9, height = 9,
       units = "in")
ggsave(file.path(fdir, "Figure6_highrisk_arg_replicon.png"), fig, width = 9, height = 9,
       units = "in", dpi = 300)
cat("Saved Figure6_highrisk_arg_replicon.pdf/png to", fdir, "\n")
