#!/usr/bin/env Rscript
# =============================================================================
# pipdb_22_bm_arg_coselection.R
# Biocide/metal resistance gene (BMG) - ARG co-selection network analysis
#
# Outputs:
#   tab_bm_arg_cooccurrence.csv    - Fisher exact tests for ARG family x BM category
#   tab_bm_arg_summary.csv         - overall co-occurrence rates
#   tab_bm_arg_by_replicon.csv     - BM-ARG co-occurrence by major replicon
#   tab_bm_arg_highrisk.csv        - BMG co-occurrence with high-risk ARGs
#   fig32_bm_arg_network.pdf/png   - bipartite ARG-BM co-occurrence network
#   fig33_bm_arg_coselection_bar   - co-occurrence rate bar chart
#   fig34_bm_arg_heatmap.pdf/png   - ARG family x BM category OR heatmap
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
})

needed <- c("igraph", "ggraph")
for (pkg in needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}
suppressPackageStartupMessages({ library(igraph); library(ggraph) })

# ---- paths ----
args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) >= 1) args[1] else "results"
fig_dir <- file.path(res_dir, "figures", "descriptive")
tab_dir <- file.path(res_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, name, w = 10, h = 8) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot, width = w, height = h, dpi = 300, limitsize = FALSE)
  cat("  saved", name, "\n")
}

theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA),
        legend.position = "right")

cat("=== Loading data ===\n")
master_cols <- c("id", "plasmid_acc", "replicon_primary", "length_avg",
                 "aro_name", "n_arg", "n_who_arg",
                 "drugclass", "vf_name", "n_vf",
                 "gene_bacmet", "n_metal", "mobility_class",
                 "has_integron", "host_human", "host_animal",
                 "habitat_human", "habitat_animal", "habitat_env")
risk_cols <- c("id", "n_high_risk_arg")

d <- fread(file.path(res_dir, "psc_master.tsv"), select = master_cols,
           na.strings = c("\\N", "", "-", "NA"))
r <- fread(file.path(res_dir, "psc_risk_scores.tsv"), select = risk_cols,
           na.strings = c("\\N", "", "-", "NA"))
d <- merge(d, r, by = "id", all.x = TRUE)
cat(sprintf("  %d PSCs loaded\n", nrow(d)))

# =============================================================================
# 1. Parse BM genes into functional categories
# =============================================================================
cat("=== Parsing BM gene categories ===\n")

bm_categories <- list(
  Mercury = c("mer[A-Z]?", "merR\\d?", "merE", "merD"),
  QAC_Disinf = c("qac", "emr", "mdfA", "norA", "acr", "mex", "sugE", "bcr", "blt", "mtrR", "adeI", "adeL", "vcaM", "mdt", "farA", "wtpC", "ttgV", "actP", "amvA", "oprN", "vexE", "evg", "bae", "gadX", "kdpE"),
  Arsenic = c("ars[A-Z]?", "aioE", "aio[AB]"),
  Copper = c("cop[A-Z]?", "pco[A-Z]?", "csoR", "copZ", "mmco"),
  Silver = c("sil[A-Z]?", "silR", "silS"),
  Zinc_Cadmium = c("czc[A-Z]?", "znt", "zitB", "zra", "znu", "cad", "cnr", "nik", "mnt", "fec", "mod", "smf", "cor")
)

classify_bm <- function(gene_str) {
  if (is.na(gene_str) || gene_str == "" || gene_str == "\\N") return(character(0))
  genes <- trimws(unlist(strsplit(gene_str, ",")))
  genes <- genes[genes != "" & !is.na(genes)]
  cats <- character(0)
  for (g in genes) {
    gl <- tolower(g)
    if (any(vapply(bm_categories$Mercury, function(p) grepl(paste0("^", p), gl), logical(1)))) {
      cats <- c(cats, "Mercury (mer)")
    } else if (any(vapply(bm_categories$QAC_Disinf, function(p) grepl(paste0("^", p), gl), logical(1)))) {
      cats <- c(cats, "QAC/Disinfectant")
    } else if (any(vapply(bm_categories$Arsenic, function(p) grepl(paste0("^", p), gl), logical(1)))) {
      cats <- c(cats, "Arsenic (ars)")
    } else if (any(vapply(bm_categories$Copper, function(p) grepl(paste0("^", p), gl), logical(1)))) {
      cats <- c(cats, "Copper (cop/pco)")
    } else if (any(vapply(bm_categories$Silver, function(p) grepl(paste0("^", p), gl), logical(1)))) {
      cats <- c(cats, "Silver (sil)")
    } else if (any(vapply(bm_categories$Zinc_Cadmium, function(p) grepl(paste0("^", p), gl), logical(1)))) {
      cats <- c(cats, "Zn/Cd/Ni")
    } else {
      cats <- c(cats, "Other/Unknown")
    }
  }
  unique(cats)
}

cat("  Classifying BM genes ...\n")
d[, bm_cats := lapply(gene_bacmet, classify_bm)]
d[, has_bm := n_metal > 0 & !is.na(n_metal)]
d[, has_mer := vapply(bm_cats, function(x) "Mercury (mer)" %in% x, logical(1))]
d[, has_qac := vapply(bm_cats, function(x) "QAC/Disinfectant" %in% x, logical(1))]
d[, has_ars := vapply(bm_cats, function(x) "Arsenic (ars)" %in% x, logical(1))]
d[, has_cop := vapply(bm_cats, function(x) any(c("Copper (cop/pco)", "Silver (sil)") %in% x), logical(1))]
d[, has_zn := vapply(bm_cats, function(x) "Zn/Cd/Ni" %in% x, logical(1))]

cat(sprintf("  BMG+: %d (%.1f%%)\n", sum(d$has_bm), 100*mean(d$has_bm)))
cat(sprintf("  Mercury: %d (%.1f%%)\n", sum(d$has_mer), 100*mean(d$has_mer)))
cat(sprintf("  QAC/Disinfectant: %d (%.1f%%)\n", sum(d$has_qac), 100*mean(d$has_qac)))
cat(sprintf("  Arsenic: %d (%.1f%%)\n", sum(d$has_ars), 100*mean(d$has_ars)))
cat(sprintf("  Cu/Ag: %d (%.1f%%)\n", sum(d$has_cop), 100*mean(d$has_cop)))
cat(sprintf("  Zn/Cd/Ni: %d (%.1f%%)\n", sum(d$has_zn), 100*mean(d$has_zn)))

# =============================================================================
# 2. Parse ARG families
# =============================================================================
cat("=== Parsing ARG families ===\n")

arg_family_map <- list(
  TEM = c("TEM"), CTX_M = c("CTX-M"), SHV = c("SHV"), CMY = c("CMY"),
  KPC = c("KPC"), NDM = c("NDM"), VIM = c("VIM"), IMP = c("IMP"),
  OXA = c("OXA"), AAC_APH = c("AAC", "APH", "ANT", "Aad"),
  Tet = c("Tet", "TetR", "tet"), Sul = c("Sul", "sul"),
  Dfr = c("Dfr", "dfr"), Qnr = c("Qnr", "qnr"),
  MCR = c("MCR", "mcr"), Van = c("van", "Van"),
  Mph_Msr = c("Mph", "Msr", "mph", "msr"),
  Cat_Flo = c("Cat", "cat", "Flo", "flo", "Cml", "cml"),
  Erm = c("Erm", "erm"), QacEdelta = c("QacEdelta", "qacEdelta")
)

classify_arg <- function(gene_str) {
  if (is.na(gene_str) || gene_str == "" || gene_str == "\\N") return(character(0))
  genes <- trimws(unlist(strsplit(gene_str, ",")))
  genes <- genes[genes != "" & !is.na(genes)]
  fams <- character(0)
  for (g in genes) {
    for (fam in names(arg_family_map)) {
      patterns <- arg_family_map[[fam]]
      if (any(vapply(patterns, function(p) grepl(paste0("^", p), g), logical(1)))) {
        fams <- c(fams, fam)
        break
      }
    }
  }
  unique(fams)
}

d[, arg_fams := lapply(aro_name, classify_arg)]
d[, has_arg := n_arg > 0 & !is.na(n_arg)]
d[, has_hr_arg := n_high_risk_arg > 0 & !is.na(n_high_risk_arg)]

# =============================================================================
# 3. Overall co-occurrence summary
# =============================================================================
cat("=== Overall co-occurrence summary ===\n")

n_total <- nrow(d)
n_arg <- sum(d$has_arg)
n_bm <- sum(d$has_bm)
n_both <- sum(d$has_arg & d$has_bm)
n_arg_only <- sum(d$has_arg & !d$has_bm)
n_bm_only <- sum(!d$has_arg & d$has_bm)
n_neither <- sum(!d$has_arg & !d$has_bm)

summary_dt <- data.table(
  category = c("Total PSCs", "ARG+ only", "BMG+ only", "ARG+ & BMG+ (co-occurrence)", "Neither"),
  n = c(n_total, n_arg_only, n_bm_only, n_both, n_neither),
  pct = c(100, 100*n_arg_only/n_total, 100*n_bm_only/n_total, 100*n_both/n_total, 100*n_neither/n_total)
)

pct_bm_among_arg <- 100 * n_both / n_arg
pct_arg_among_bm <- 100 * n_both / n_bm

cat(sprintf("  ARG+ plasmids: %d (%.1f%%)\n", n_arg, 100*n_arg/n_total))
cat(sprintf("  BMG+ plasmids: %d (%.1f%%)\n", n_bm, 100*n_bm/n_total))
cat(sprintf("  Co-occurrence: %d (%.1f%% of all)\n", n_both, 100*n_both/n_total))
cat(sprintf("  Among ARG+: %.1f%% also carry BMGs\n", pct_bm_among_arg))
cat(sprintf("  Among BMG+: %.1f%% also carry ARGs\n", pct_arg_among_bm))

ct <- matrix(c(n_both, n_arg_only, n_bm_only, n_neither), nrow = 2,
             dimnames = list(c("ARG+", "ARG-"), c("BMG+", "BMG-")))
ft <- fisher.test(ct)
cat(sprintf("  Overall Fisher OR = %.2f (95%% CI %.2f-%.2f), p = %.2e\n",
            ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))

summary_dt <- rbind(summary_dt, data.table(
  category = c("BMG% among ARG+", "ARG% among BMG+", "Fisher OR (overall)"),
  n = c(NA, NA, NA),
  pct = c(pct_bm_among_arg, pct_arg_among_bm, ft$estimate)
))
fwrite(summary_dt, file.path(tab_dir, "tab_bm_arg_summary.csv"))

# =============================================================================
# 4. ARG family x BM category Fisher exact tests
# =============================================================================
cat("=== ARG family x BM category co-occurrence tests ===\n")

bm_cols <- c("Mercury (mer)" = "has_mer",
             "QAC/Disinfectant" = "has_qac",
             "Arsenic (ars)" = "has_ars",
             "Cu/Ag" = "has_cop",
             "Zn/Cd/Ni" = "has_zn")

all_arg_fams <- names(arg_family_map)
for (fam in all_arg_fams) {
  d[, (fam) := vapply(arg_fams, function(x) fam %in% x, logical(1))]
}

arg_counts <- sapply(all_arg_fams, function(f) sum(d[[f]]))
test_fams <- names(arg_counts[arg_counts >= 100])
cat(sprintf("  Testing %d ARG families (>=100 occurrences) x %d BM categories\n",
            length(test_fams), length(bm_cols)))

results <- list()
for (fam in test_fams) {
  for (bm_name in names(bm_cols)) {
    bm_col <- bm_cols[[bm_name]]
    a <- sum(d[[fam]] & d[[bm_col]])
    b <- sum(d[[fam]] & !d[[bm_col]])
    c_ <- sum(!d[[fam]] & d[[bm_col]])
    dd <- sum(!d[[fam]] & !d[[bm_col]])
    if (a < 10) next
    ct2 <- matrix(c(a, b, c_, dd), nrow = 2)
    ft2 <- tryCatch(fisher.test(ct2), error = function(e) NULL)
    if (is.null(ft2)) next
    results[[length(results) + 1]] <- data.table(
      arg_family = fam, bm_category = bm_name,
      n_both = a, n_arg_only = b, n_bm_only = c_, n_neither = dd,
      odds_ratio = as.numeric(ft2$estimate),
      or_low = ft2$conf.int[1], or_high = ft2$conf.int[2],
      p_value = ft2$p.value,
      enrichment = (a / (a + b)) / ((a + c_) / n_total)
    )
  }
}

res <- rbindlist(results)
res[, p_adj := p.adjust(p_value, method = "BH")]
setorder(res, -odds_ratio)
fwrite(res, file.path(tab_dir, "tab_bm_arg_cooccurrence.csv"))
cat(sprintf("  %d significant pairs (BH p_adj < 0.001, OR > 2)\n",
            sum(res$p_adj < 0.001 & res$odds_ratio > 2)))
print(head(res[order(-odds_ratio)], 20))

# =============================================================================
# 5. By replicon analysis
# =============================================================================
cat("=== BM-ARG co-occurrence by major replicon ===\n")

rep_stats <- d[!is.na(replicon_primary), .(
  n = .N, n_arg = sum(has_arg), n_bm = sum(has_bm),
  n_both = sum(has_arg & has_bm), n_hr_arg = sum(has_hr_arg),
  n_mer = sum(has_mer), n_qac = sum(has_qac)
), by = replicon_primary][n >= 50][order(-n_both)]

rep_stats[, pct_arg_bm := 100 * n_both / n_arg]
rep_stats[, pct_bm := 100 * n_bm / n]
setorder(rep_stats, -pct_arg_bm)
fwrite(rep_stats, file.path(tab_dir, "tab_bm_arg_by_replicon.csv"))
cat("  Top 15 replicons by BM-ARG co-occurrence rate:\n")
print(head(rep_stats, 15))

# =============================================================================
# 6. Visualization
# =============================================================================
cat("=== Visualization ===\n")

bm_long <- data.table(
  BM_category = c("Mercury (mer)", "QAC/Disinfectant", "Arsenic (ars)", "Cu/Ag", "Zn/Cd/Ni"),
  among_ARG = c(
    100 * sum(d$has_arg & d$has_mer) / n_arg,
    100 * sum(d$has_arg & d$has_qac) / n_arg,
    100 * sum(d$has_arg & d$has_ars) / n_arg,
    100 * sum(d$has_arg & d$has_cop) / n_arg,
    100 * sum(d$has_arg & d$has_zn) / n_arg
  )
)

p_bar <- ggplot(bm_long, aes(x = reorder(BM_category, -among_ARG), y = among_ARG)) +
  geom_col(aes(fill = BM_category), width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", among_ARG)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Mercury (mer)" = "#762a83",
                                "QAC/Disinfectant" = "#1b7837",
                                "Arsenic (ars)" = "#b35806",
                                "Cu/Ag" = "#c51b7d",
                                "Zn/Cd/Ni" = "#053061")) +
  labs(x = NULL, y = "% of ARG-carrying plasmids\nalso carrying BMG category",
       title = "Biocide/metal gene co-occurrence with ARGs") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_bar, "fig33_bm_arg_coselection_bar", w = 7, h = 5)

res_sig <- res[p_adj < 0.001 & odds_ratio > 1.5]
if (nrow(res_sig) > 0) {
  res_sig[, or_plot := pmin(odds_ratio, 15)]
  p_hm <- ggplot(res_sig, aes(x = bm_category, y = arg_family, fill = or_plot)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", odds_ratio)), size = 3, color = "white") +
    scale_fill_gradientn(colours = c("#f7f7f7", "#fdd49e", "#fc8d59", "#d7301f", "#7f0000"),
                         name = "Odds Ratio", limits = c(1, 15), na.value = "#f7f7f7") +
    labs(x = "Biocide/metal resistance category", y = "ARG family",
         title = "ARG-BMG co-occurrence odds ratios",
         subtitle = "Fisher's exact test, BH-adjusted p < 0.001") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p_hm, "fig34_bm_arg_heatmap", w = 8, h = 7)
}

net_edges <- res[p_adj < 0.001 & odds_ratio > 2 & n_both >= 20]
if (nrow(net_edges) > 0) {
  nodes_arg <- data.table(name = unique(net_edges$arg_family), type = "ARG")
  nodes_bm <- data.table(name = unique(net_edges$bm_category), type = "BMG")
  nodes <- rbind(nodes_arg, nodes_bm)
  g <- graph_from_data_frame(net_edges[, .(arg_family, bm_category, odds_ratio, n_both)],
                              vertices = nodes, directed = FALSE)
  arg_tot <- net_edges[, .(tot = sum(n_both)), by = arg_family]
  V(g)$size <- ifelse(V(g)$type == "ARG",
                       arg_tot[match(V(g)$name, arg_family), sqrt(tot)] * 1.5, 12)
  V(g)$size[is.na(V(g)$size)] <- 12
  set.seed(42)
  p_net <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = odds_ratio, alpha = after_stat(index)),
                   edge_colour = "grey60", show.legend = TRUE) +
    geom_node_point(aes(color = type, size = size), show.legend = FALSE) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.5, max.overlaps = 20) +
    scale_edge_width(range = c(0.5, 4), name = "Odds Ratio") +
    scale_color_manual(values = c("ARG" = "#d73027", "BMG" = "#1b7837")) +
    labs(title = "ARG-biocide/metal resistance gene co-occurrence network",
         subtitle = "Edges: Fisher's exact BH p < 0.001, OR > 2, n >= 20") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"))
  save_plot(p_net, "fig32_bm_arg_network", w = 10, h = 8)
}

# =============================================================================
# 7. High-risk ARG co-occurrence with BM categories
# =============================================================================
cat("=== High-risk ARG co-occurrence with BMGs ===\n")

hr_bm <- data.table(
  BM_category = c("Mercury (mer)", "QAC/Disinfectant", "Arsenic (ars)", "Cu/Ag", "Zn/Cd/Ni", "Any BMG"),
  pct_among_hr = c(
    100 * sum(d$has_hr_arg & d$has_mer) / max(sum(d$has_hr_arg), 1),
    100 * sum(d$has_hr_arg & d$has_qac) / max(sum(d$has_hr_arg), 1),
    100 * sum(d$has_hr_arg & d$has_ars) / max(sum(d$has_hr_arg), 1),
    100 * sum(d$has_hr_arg & d$has_cop) / max(sum(d$has_hr_arg), 1),
    100 * sum(d$has_hr_arg & d$has_zn) / max(sum(d$has_hr_arg), 1),
    100 * sum(d$has_hr_arg & d$has_bm) / max(sum(d$has_hr_arg), 1)
  ),
  pct_among_nonhr_arg = c(
    100 * sum(d$has_arg & !d$has_hr_arg & d$has_mer) / max(sum(d$has_arg & !d$has_hr_arg), 1),
    100 * sum(d$has_arg & !d$has_hr_arg & d$has_qac) / max(sum(d$has_arg & !d$has_hr_arg), 1),
    100 * sum(d$has_arg & !d$has_hr_arg & d$has_ars) / max(sum(d$has_arg & !d$has_hr_arg), 1),
    100 * sum(d$has_arg & !d$has_hr_arg & d$has_cop) / max(sum(d$has_arg & !d$has_hr_arg), 1),
    100 * sum(d$has_arg & !d$has_hr_arg & d$has_zn) / max(sum(d$has_arg & !d$has_hr_arg), 1),
    100 * sum(d$has_arg & !d$has_hr_arg & d$has_bm) / max(sum(d$has_arg & !d$has_hr_arg), 1)
  )
)
fwrite(hr_bm, file.path(tab_dir, "tab_bm_arg_highrisk.csv"))
print(hr_bm)

cat("\n=== Done ===\n")
cat("Tables:\n")
cat("  tab_bm_arg_summary.csv       - overall co-occurrence rates\n")
cat("  tab_bm_arg_cooccurrence.csv  - ARG family x BM category Fisher tests\n")
cat("  tab_bm_arg_by_replicon.csv   - co-occurrence by major replicon\n")
cat("  tab_bm_arg_highrisk.csv      - BMG co-occurrence with high-risk ARGs\n")
cat("Figures:\n")
cat("  fig32_bm_arg_network         - bipartite ARG-BMG network\n")
cat("  fig33_bm_arg_coselection_bar - co-occurrence rate bar chart\n")
cat("  fig34_bm_arg_heatmap         - OR heatmap\n")
