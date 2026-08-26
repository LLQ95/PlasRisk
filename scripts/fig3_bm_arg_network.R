#!/usr/bin/env Rscript
# =============================================================================
# Figure 3: BMG-ARG co-occurrence network
# (A) Global ARG family x BMG category bipartite network (OR > 2, BH p_adj < 0.001)
# (B) High-risk ARG-focused subnetwork (top associations)
# (C) OR bar chart for strongest non-overlapping associations
# (D) Replicon-level BMG-ARG co-carriage (moved from Fig 2D):
#     high-risk ARG carrier rate by replicon, colored by BMG co-carriage
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
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

warm <- c("#E76F51", "#F4A261", "#E9C46A", "#2A9D8F", "#6D4C41",
          "#9C2C2C", "#457B9D", "#BC6C25", "#606C38", "#8D6E63")

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

# ---- load association table (corrected, non-overlapping genes) ----
assoc_file <- file.path(tdir, "tab_bm_arg_associations_corrected.csv")
if (!file.exists(assoc_file)) assoc_file <- file.path(tdir, "tab_bm_arg_associations.csv")
assoc <- fread(assoc_file)

# column flexibility
col_or <- intersect(c("OR", "odds_ratio", "or"), names(assoc))[1]
col_p  <- intersect(c("p_adj", "padj", "fdr", "p.value"), names(assoc))[1]
col_arg <- intersect(c("ARG_family", "arg_family", "ARG"), names(assoc))[1]
col_bmg <- intersect(c("BMG_category", "bmg_category", "BMG"), names(assoc))[1]
if (is.na(col_or)) { assoc[, OR := 1]; col_or <- "OR" }
if (is.na(col_p))  { assoc[, p_adj := 1]; col_p <- "p_adj" }
if (is.na(col_arg)) stop("Cannot find ARG family column")
if (is.na(col_bmg)) stop("Cannot find BMG category column")
setnames(assoc, c(col_or, col_p, col_arg, col_bmg),
         c("OR", "p_adj", "ARG_family", "BMG_category"))

sig <- assoc[OR > 2 & p_adj < 0.001][order(-OR)]

# ---- (A) full bipartite network ----
plot_network <- function(d, tag, title_lab = NULL) {
  if (!has_igraph || !has_ggraph || nrow(d) == 0) {
    return(ggplot() + theme_void() +
      annotate("text", 0.5, 0.5, label = "igraph/ggraph required", size = 3) +
      labs(tag = tag))
  }
  edges <- d[, .(from = ARG_family, to = BMG_category,
                 width = pmin(log10(OR), 2.5))]
  g <- graph_from_data_frame(edges, directed = FALSE)
  V(g)$type <- V(g)$name %in% d$BMG_category
  V(g)$kind <- ifelse(V(g)$type, "BMG", "ARG")
  set.seed(42)
  ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = width), alpha = 0.35, color = "#E76F51") +
    scale_edge_width(range = c(0.2, 2)) +
    geom_node_point(aes(color = kind, size = ifelse(kind == "BMG", 6, 4)),
                    show.legend = FALSE) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.6, max.overlaps = 30,
                   color = "#3E2723") +
    scale_color_manual(values = c("ARG" = "#E76F51", "BMG" = "#2A9D8F")) +
    labs(tag = tag, color = NULL, title = title_lab) +
    theme_void() +
    theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 12),
          plot.title = element_text(size = 9, hjust = 0.5))
}

pA <- plot_network(sig, "A")

# ---- (B) high-risk ARG subnetwork ----
highrisk_args <- c("NDM", "KPC", "IMP", "VIM", "OXA", "MCR", "CTX-M", "SHV",
                   "TetX", "Carbapenemase", "mcr")
hr_pat <- paste(highrisk_args, collapse = "|")
sub_hr <- sig[grepl(hr_pat, ARG_family, ignore.case = TRUE)]
if (nrow(sub_hr) > 24) sub_hr <- sub_hr[1:24]
pB <- plot_network(sub_hr, "B")

# ---- (C) strongest distinct associations bar chart ----
top_bar <- sig[1:min(15, nrow(sig))]
top_bar[, pair := factor(paste(ARG_family, BMG_category, sep = " - "),
                         levels = rev(paste(ARG_family, BMG_category, sep = " - ")))]
pC <- ggplot(top_bar, aes(OR, pair)) +
  geom_col(fill = "#E76F51", alpha = 0.85, width = 0.7) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Odds ratio (log10-scaled axis)", y = NULL, tag = "C") +
  scale_x_log10() +
  theme_pub +
  theme(axis.text.y = element_text(size = 7))

# ---- (D) replicon high-risk ARG rate colored by BMG co-carriage ----
rep_file <- file.path(tdir, "tab_replicon_summary.csv")
pD <- ggplot() + theme_void() + labs(tag = "D")
if (file.exists(rep_file)) {
  rep <- fread(rep_file)
  hrcol  <- intersect(c("highrisk_rate", "hr_rate", "high_risk_rate", "arg_rate"), names(rep))[1]
  bmcol  <- intersect(c("bm_rate", "bmg_rate", "bm_carriage"), names(rep))[1]
  ncol_  <- intersect(c("n", "N", "n_psc"), names(rep))[1]
  namecol<- intersect(c("replicon_primary", "replicon", "Replicon"), names(rep))[1]
  if (!is.na(hrcol) && !is.na(ncol_) && !is.na(namecol)) {
    setnames(rep, c(hrcol, ncol_, namecol), c("hr_rate", "n", "replicon"))
    if (is.na(bmcol)) rep[, bm_rate := 0] else setnames(rep, bmcol, "bm_rate")
    d <- rep[n >= 100][order(-hr_rate)][1:min(20, .SD$n |> length())]
    d <- d[1:min(20, nrow(d))]
    d[, replicon := factor(replicon, levels = rev(replicon))]
    d[, bm_frac := ifelse(bm_rate > 0, pmin(bm_rate, 1), 0)]
    pD <- ggplot(d, aes(hr_rate, replicon)) +
      geom_col(aes(fill = bm_frac), width = 0.7) +
      scale_fill_gradient2(low = "#F4A261", mid = "#E76F51", high = "#9C2C2C",
                           midpoint = 0.5, limits = c(0, 1),
                           name = "BMG\nco-carriage") +
      labs(x = "High-risk ARG carrier rate", y = NULL, tag = "D") +
      theme_pub +
      theme(axis.text.y = element_text(size = 7))
  }
}

fig <- (pA | pB) / (pC | pD)

ggsave(file.path(fdir, "Figure3_bm_arg_network.pdf"), fig,
       width = 9, height = 9, units = "in")
ggsave(file.path(fdir, "Figure3_bm_arg_network.png"), fig,
       width = 9, height = 9, units = "in", dpi = 300)
cat("Saved Figure3_bm_arg_network.pdf/png to", fdir, "\n")
cat(sprintf("Significant associations plotted: %d (of %d total)\n", nrow(sig), nrow(assoc)))
