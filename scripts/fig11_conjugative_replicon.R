#!/usr/bin/env Rscript
# =============================================================================
# Figure 11: Conjugative potential and resistance burden across replicons
# (A) Top 15 replicons by conjugative rate (colored by mean PlasRisk score)
# (B) Bubble: conjugative rate vs high-risk ARG carrier rate
# (C) Bubble: conjugative rate vs BMG carrier rate
# (D) Heatmap: four outcome rates across top conjugative replicons
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

# ---- Load and aggregate per-replicon outcomes ----
psc <- fread(file.path(tdir, "tab_psc_final_scores.csv"),
             select = c("replicon_primary", "S_final", "S_MOB",
                        "y_highrisk", "y_fusion", "y_conj", "y_bm"))
rep_sum <- fread(file.path(tdir, "tab_replicon_summary.csv"),
                 select = c("replicon_primary", "mean_len"))

rep <- psc[!is.na(replicon_primary), .(
  n          = .N,
  conj_rate  = mean(y_conj) * 100,
  hr_rate    = mean(y_highrisk) * 100,
  fus_rate   = mean(y_fusion) * 100,
  bm_rate    = mean(y_bm) * 100,
  mean_risk  = mean(S_final),
  mean_mob   = mean(S_MOB)
), by = replicon_primary]

rep <- merge(rep, rep_sum, by = "replicon_primary", all.x = TRUE)
rep <- rep[n >= 100 & !is.na(mean_len)]

# Spearman correlations for annotation
rho_hr  <- cor(rep$conj_rate, rep$hr_rate,  method = "spearman")
rho_bm  <- cor(rep$conj_rate, rep$bm_rate,  method = "spearman")
rho_fus <- cor(rep$conj_rate, rep$fus_rate, method = "spearman")
med_conj <- median(rep$conj_rate)
med_hr   <- median(rep$hr_rate)
med_bm   <- median(rep$bm_rate)

# ---- (A) Top 15 replicons by conjugative rate ----
rep_conj <- rep[order(-conj_rate)][1:15]
rep_conj[, replicon_primary := factor(replicon_primary,
  levels = rev(rep_conj$replicon_primary))]

pA <- ggplot(rep_conj, aes(conj_rate, replicon_primary)) +
  geom_col(aes(fill = mean_risk), width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", conj_rate)), hjust = -0.15, size = 2.5) +
  scale_fill_gradient(low = "#FFE8D6", high = "#9C2C2C",
                      name = "Mean\nPlasRisk") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = "Conjugative rate (%)", y = NULL, tag = "A") +
  theme_pub +
  theme(axis.text.y = element_text(size = 7.5),
        legend.key.height = unit(0.5, "cm"))

# ---- (B) Bubble: conjugative vs high-risk ARG ----
# Label only top replicons on either axis
top_conj_b <- rep[order(-conj_rate)][1:6, replicon_primary]
top_hr     <- rep[order(-hr_rate)][1:6, replicon_primary]
rep[, label_b := ifelse(
  replicon_primary %in% union(top_conj_b, top_hr), replicon_primary, "")]

pB <- ggplot(rep, aes(conj_rate, hr_rate)) +
  geom_hline(yintercept = med_hr, linetype = "dashed",
             color = "#BCAAA4", linewidth = 0.3) +
  geom_vline(xintercept = med_conj, linetype = "dashed",
             color = "#BCAAA4", linewidth = 0.3) +
  geom_point(aes(size = n, color = mean_len), alpha = 0.75) +
  geom_text_repel(aes(label = label_b), size = 2.5, max.overlaps = 40,
                  color = "#3E2723", box.padding = 0.35, min.segment.length = 0,
                  segment.color = "#BCAAA4", segment.size = 0.3,
                  max.time = 2, max.iter = 10000) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("rho == %.2f", rho_hr),
           parse = TRUE, size = 2.8, color = "#6D4C41") +
  scale_size_continuous(name = "n (PSC)", range = c(1.5, 6)) +
  scale_color_gradient(low = "#FFE8D6", high = "#2A9D8F",
                       name = "Mean length\n(kb)") +
  labs(x = "Conjugative rate (%)",
       y = "High-risk ARG carrier rate (%)", tag = "B") +
  theme_pub +
  guides(size = guide_legend(nrow = 2), color = guide_colorbar(barwidth = 4))

# ---- (C) Bubble: conjugative vs BMG carrier ----
# Label only top replicons on either axis; require n >= 200 for BMG leaders
# to avoid piling up labels among numerous small ~100% BMG replicons
top_conj_c <- rep[order(-conj_rate)][1:6, replicon_primary]
top_bm     <- rep[n >= 200][order(-bm_rate)][1:6, replicon_primary]
rep[, label_c := ifelse(
  replicon_primary %in% union(top_conj_c, top_bm), replicon_primary, "")]

pC <- ggplot(rep, aes(conj_rate, bm_rate)) +
  geom_hline(yintercept = med_bm, linetype = "dashed",
             color = "#BCAAA4", linewidth = 0.3) +
  geom_vline(xintercept = med_conj, linetype = "dashed",
             color = "#BCAAA4", linewidth = 0.3) +
  geom_point(aes(size = n, color = mean_risk), alpha = 0.75) +
  geom_text_repel(aes(label = label_c), size = 2.5, max.overlaps = 40,
                  color = "#3E2723", box.padding = 0.35, min.segment.length = 0,
                  segment.color = "#BCAAA4", segment.size = 0.3,
                  max.time = 2, max.iter = 10000) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("rho == %.2f", rho_bm),
           parse = TRUE, size = 2.8, color = "#6D4C41") +
  scale_size_continuous(name = "n (PSC)", range = c(1.5, 6)) +
  scale_color_gradient(low = "#FFE8D6", high = "#9C2C2C",
                       name = "Mean\nPlasRisk") +
  labs(x = "Conjugative rate (%)",
       y = "BMG carrier rate (%)", tag = "C") +
  theme_pub +
  guides(size = guide_legend(nrow = 2), color = guide_colorbar(barwidth = 4))

# ---- (D) Heatmap: four outcomes across top conjugative replicons ----
heat_reps <- rep[order(-conj_rate)][1:15, replicon_primary]
heat_dt <- melt(
  rep[replicon_primary %in% heat_reps,
      .(replicon_primary, Conjugative = conj_rate,
        `High-risk ARG` = hr_rate, `MDR-VF fusion` = fus_rate,
        `Biocide/metal` = bm_rate)],
  id.vars = "replicon_primary", variable.name = "Outcome", value.name = "rate"
)
# Column-scaled fill for visual contrast; raw values annotated
heat_dt[, rate_scaled := rate / max(rate), by = Outcome]
heat_dt[, replicon_primary := factor(replicon_primary,
  levels = rev(heat_reps))]

pD <- ggplot(heat_dt, aes(Outcome, replicon_primary, fill = rate_scaled)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", rate)), size = 2.3, color = "#3E2723") +
  scale_fill_gradient(low = "#FFF8E7", high = "#E76F51",
                      name = "Column-\nscaled", guide = "none") +
  labs(x = NULL, y = NULL, tag = "D") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5),
        axis.text.y = element_text(size = 7.5))

# ---- assemble ----
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1))

ggsave(file.path(fdir, "Figure11_conjugative_replicon.pdf"), fig,
       width = 7.5, height = 8.5, units = "in")
ggsave(file.path(fdir, "Figure11_conjugative_replicon.png"), fig,
       width = 7.5, height = 8.5, units = "in", dpi = 300)
cat("Saved Figure11_conjugative_replicon.pdf/png to", fdir, "\n")
cat(sprintf("Spearman rho: conj vs HR=%.3f, conj vs BMG=%.3f, conj vs fusion=%.3f\n",
            rho_hr, rho_bm, rho_fus))
