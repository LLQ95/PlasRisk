#!/usr/bin/env Rscript
# =============================================================================
# Figure 7: MDR-VF fusion plasmids
# (A) Fusion rate by replicon
# (B) Feature comparison fusion vs non-fusion (length, IS, integron, conjugative)
# (C) Absolute fusion counts by replicon
# (D) Mean VF count vs mean ARG count per replicon
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
  theme(strip.background = element_blank(), legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 12, color = "#3E2723"),
        axis.text = element_text(color = "#3E2723"),
        plot.margin = margin(4, 6, 4, 4))

rep <- fread(file.path(tdir, "tab_replicon_summary.csv"))
namecol <- intersect(c("replicon_primary", "replicon", "Replicon"), names(rep))[1]
ncol_   <- intersect(c("n", "N", "n_psc"), names(rep))[1]
fuscol  <- intersect(c("fusion_rate", "vf_arg_rate", "fusion"), names(rep))[1]
argcol  <- intersect(c("mean_arg", "arg_count", "mean_narg"), names(rep))[1]
vfcol   <- intersect(c("mean_vf", "vf_count", "mean_nvf"), names(rep))[1]
setnames(rep, c(namecol, ncol_), c("replicon", "n"))
if (is.na(fuscol)) rep[, fusion_rate := NA_real_] else setnames(rep, fuscol, "fusion_rate")
if (is.na(argcol)) rep[, mean_arg := NA_real_] else setnames(rep, argcol, "mean_arg")
if (is.na(vfcol))  rep[, mean_vf := NA_real_] else setnames(rep, vfcol, "mean_vf")

# ---- (A) fusion rate by replicon ----
top <- rep[n >= 100 & !is.na(fusion_rate)][order(-fusion_rate)][1:15]
top[, replicon := factor(replicon, levels = rev(replicon))]
pA <- ggplot(top, aes(fusion_rate * 100, replicon)) +
  geom_col(fill = "#F4A261", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", fusion_rate * 100)), hjust = -0.1, size = 2.6) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "MDR-VF fusion rate (%)", y = NULL, tag = "A") +
  theme_pub + theme(axis.text.y = element_text(size = 7.5))

# ---- (B) feature comparison ----
feat_file <- file.path(tdir, "tab_fusion_features.csv")
pB <- ggplot() + theme_void() + labs(tag = "B")
if (file.exists(feat_file)) {
  feat <- fread(feat_file)
  if (all(c("feature", "fusion", "nonfusion") %in% names(feat))) {
    fl <- melt(feat, id.vars = "feature", variable.name = "group", value.name = "value")
    fl[, group := factor(group, levels = c("nonfusion", "fusion"),
                         labels = c("Non-fusion", "Fusion"))]
    pB <- ggplot(fl, aes(feature, value, fill = group)) +
      geom_col(position = position_dodge(0.8), width = 0.7) +
      scale_fill_manual(values = c("Non-fusion" = "#A8DADC", "Fusion" = "#E76F51")) +
      labs(x = NULL, y = "Mean / proportion", fill = NULL, tag = "B") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7.5))
  }
}

# ---- (C) absolute counts ----
rep[, n_fusion := round(fusion_rate * n)]
cn <- rep[!is.na(n_fusion)][order(-n_fusion)][1:15]
cn[, replicon := factor(replicon, levels = rev(replicon))]
pC <- ggplot(cn, aes(n_fusion, replicon)) +
  geom_col(fill = "#2A9D8F", width = 0.7) +
  geom_text(aes(label = n_fusion), hjust = -0.1, size = 2.6) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of fusion plasmids", y = NULL, tag = "C") +
  theme_pub + theme(axis.text.y = element_text(size = 7.5))

# ---- (D) VF vs ARG scatter ----
sc <- rep[n >= 100 & !is.na(mean_arg) & !is.na(mean_vf)]
sc[, lab := ifelse(mean_vf > quantile(mean_vf, 0.9, na.rm = TRUE) |
                     mean_arg > quantile(mean_arg, 0.9, na.rm = TRUE),
                   as.character(replicon), "")]
pD <- ggplot(sc, aes(mean_arg, mean_vf, size = n)) +
  geom_point(alpha = 0.55, color = "#F4A261") +
  geom_text_repel(aes(label = lab), size = 2.4, max.overlaps = 15) +
  scale_size_continuous(range = c(1, 7), name = "PSC count") +
  labs(x = "Mean ARG count", y = "Mean VF count", tag = "D") +
  theme_pub

fig <- (pA | pB) / (pC | pD)
ggsave(file.path(fdir, "Figure7_fusion_plasmids.pdf"), fig, width = 9, height = 8,
       units = "in")
ggsave(file.path(fdir, "Figure7_fusion_plasmids.png"), fig, width = 9, height = 8,
       units = "in", dpi = 300)
cat("Saved Figure7_fusion_plasmids.pdf/png to", fdir, "\n")
