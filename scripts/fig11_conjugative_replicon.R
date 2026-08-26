#!/usr/bin/env Rscript
# =============================================================================
# Figure 11: Conjugative capacity vs resistance cargo across replicons
# (A) Top conjugative replicons (n >= 100)
# (B) Conjugative rate vs high-risk ARG rate (bubble = n)
# (C) Conjugative rate vs BMG carriage (bubble = n)
# (D) Conjugative rate vs MDR-VF fusion rate (bubble = n)
# Spearman correlations are annotated with plotmath (ASCII-safe on R 4.1.x).
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
conjcol <- intersect(c("conj_rate", "conjugative_rate", "mob_rate"), names(rep))[1]
hrcol   <- intersect(c("highrisk_rate", "hr_rate", "high_risk_rate"), names(rep))[1]
bmcol   <- intersect(c("bm_rate", "bmg_rate"), names(rep))[1]
fuscol  <- intersect(c("fusion_rate", "vf_arg_rate"), names(rep))[1]
setnames(rep, c(namecol, ncol_), c("replicon", "n"))
if (is.na(conjcol)) stop("conjugative rate column not found")
setnames(rep, conjcol, "conj_rate")
if (is.na(hrcol))  rep[, hr_rate := NA_real_] else setnames(rep, hrcol, "hr_rate")
if (is.na(bmcol))  rep[, bm_rate := NA_real_] else setnames(rep, bmcol, "bm_rate")
if (is.na(fuscol)) rep[, fus_rate := NA_real_] else setnames(rep, fuscol, "fus_rate")

d <- rep[n >= 100 & !is.na(conj_rate)]

# ---- (A) top conjugative replicons ----
top <- d[order(-conj_rate)][1:15]
top[, replicon := factor(replicon, levels = rev(replicon))]
pA <- ggplot(top, aes(conj_rate * 100, replicon)) +
  geom_col(fill = "#2A9D8F", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", conj_rate * 100)), hjust = -0.1, size = 2.6) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Conjugative rate (%)", y = NULL, tag = "A") +
  theme_pub + theme(axis.text.y = element_text(size = 7.5))

# ---- scatter helper: label only high-axis points to avoid clutter ----
scatter_panel <- function(yvar, ylab, color, tag, topn = 8) {
  dd <- d[!is.na(get(yvar))]
  rho <- suppressWarnings(cor(dd$conj_rate, dd[[yvar]], method = "spearman"))
  dd[, score := conj_rate + get(yvar)]
  dd[, lab := ""]
  idx <- order(-dd$score)[1:min(topn, nrow(dd))]
  dd$lab[idx] <- dd$replicon[idx]
  rho_lab <- sprintf("rho == %.2f", rho)
  ggplot(dd, aes(conj_rate * 100, get(yvar) * 100, size = n)) +
    geom_point(alpha = 0.55, color = color) +
    geom_text_repel(aes(label = lab), size = 2.4, max.overlaps = 20, color = "#3E2723") +
    annotate("text", -Inf, Inf, hjust = -0.1, vjust = 1.5,
             label = rho_lab, parse = TRUE, size = 3, color = "#3E2723") +
    scale_size_continuous(range = c(1, 7), name = "PSC count") +
    labs(x = "Conjugative rate (%)", y = ylab, tag = tag) +
    theme_pub
}

pB <- scatter_panel("hr_rate",  "High-risk ARG rate (%)", "#E76F51", "B")
pC <- scatter_panel("bm_rate",  "BMG carriage (%)",       "#6D4C41", "C")
pD <- scatter_panel("fus_rate", "MDR-VF fusion rate (%)", "#F4A261", "D")

fig <- (pA | pB) / (pC | pD)
ggsave(file.path(fdir, "Figure11_conjugative_replicon.pdf"), fig,
       width = 10, height = 8.5, units = "in")
ggsave(file.path(fdir, "Figure11_conjugative_replicon.png"), fig,
       width = 10, height = 8.5, units = "in", dpi = 300)
cat("Saved Figure11_conjugative_replicon.pdf/png to", fdir, "\n")

cat(sprintf("Spearman rho: conj vs HR=%.3f, conj vs BMG=%.3f, conj vs fusion=%.3f\n",
            suppressWarnings(cor(d$conj_rate, d$hr_rate, method = "spearman")),
            suppressWarnings(cor(d$conj_rate, d$bm_rate, method = "spearman")),
            suppressWarnings(cor(d$conj_rate, d$fus_rate, method = "spearman"))))
