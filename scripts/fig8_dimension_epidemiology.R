#!/usr/bin/env Rscript
# =============================================================================
# Figure 8: Dimension epidemiology
# (A) Dimension-score heatmap across top replicons
# (B) Temporal trends 1970-2020 (ARG, BMG, co-carriage)
# (C) Dimension correlation matrix
# (D) Risk landscape bubble (S_ARG vs S_BM)
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

dims <- c("S_ARG", "S_VF", "S_MOB", "S_HOST", "S_REP", "S_SIZE",
          "S_BM", "S_GEO", "S_HAB", "S_GROW")

# ---- (A) heatmap ----
rep <- fread(file.path(tdir, "tab_replicon_risk_10dimensions.csv"))
namecol <- intersect(c("replicon_primary", "replicon", "Replicon"), names(rep))[1]
setnames(rep, namecol, "replicon")
have_dims <- intersect(dims, names(rep))
scorecol <- intersect(c("S_final", "S", "mean_score", "risk_score"), names(rep))[1]
if (!is.na(scorecol)) setorder(rep, -get(scorecol))
top_reps <- rep$replicon[1:20]
hm <- melt(rep[replicon %in% top_reps], id.vars = "replicon",
           measure.vars = have_dims, variable.name = "dimension", value.name = "score")
hm[, replicon := factor(replicon, levels = rev(top_reps))]
pA <- ggplot(hm, aes(dimension, replicon, fill = score)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#FFF8E7", mid = "#F4A261", high = "#9C2C2C",
                       midpoint = 0.4, name = "Mean\nscore") +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
        axis.text.y = element_text(size = 7))

# ---- (B) temporal trends ----
temp_file <- file.path(tdir, "tab_temporal_outcomes.csv")
pB <- ggplot() + theme_void() + labs(tag = "B")
if (file.exists(temp_file)) {
  tt <- fread(temp_file)
  if (all(c("year", "arg_rate", "bm_rate") %in% names(tt))) {
    if (!"both_rate" %in% names(tt)) tt[, both_rate := arg_rate * bm_rate]
    tl <- melt(tt, id.vars = "year",
               measure.vars = intersect(c("arg_rate", "bm_rate", "both_rate"), names(tt)),
               variable.name = "series", value.name = "rate")
    tl[, series := factor(series, levels = c("arg_rate", "bm_rate", "both_rate"),
                          labels = c("ARG carriage", "BMG carriage", "Both"))]
    pB <- ggplot(tl, aes(year, rate * 100, color = series)) +
      geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
      scale_color_manual(values = c("ARG carriage" = "#E76F51",
                                    "BMG carriage" = "#2A9D8F",
                                    "Both" = "#6D4C41")) +
      labs(x = "Year", y = "Carriage (%)", color = NULL, tag = "B") +
      theme_pub
  }
}

# ---- (C) correlation matrix ----
psc <- fread(file.path(tdir, "tab_psc_final_scores.csv"), select = have_dims)
cm <- cor(psc, use = "pairwise.complete.obs", method = "spearman")
ml <- as.data.table(as.table(cm))
setnames(ml, c("V1", "V2", "N"), c("dim1", "dim2", "r"))
ml[, dim1 := factor(dim1, levels = dims)]
ml[, dim2 := factor(dim2, levels = rev(dims))]
pC <- ggplot(ml, aes(dim1, dim2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.2) +
  scale_fill_gradient2(low = "#2A9D8F", mid = "white", high = "#9C2C2C",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman r") +
  labs(x = NULL, y = NULL, tag = "C") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7))

# ---- (D) risk landscape ----
land <- rep
land[, lab := ifelse(get(have_dims[1]) > quantile(get(have_dims[1]), 0.9, na.rm = TRUE) |
                       get("S_BM") > quantile(get("S_BM"), 0.9, na.rm = TRUE),
                     as.character(replicon), "")]
pD <- ggplot(land, aes(S_ARG, S_BM)) +
  geom_point(aes(size = n), alpha = 0.5, color = "#E76F51") +
  geom_text_repel(aes(label = lab), size = 2.3, max.overlaps = 12) +
  scale_size_continuous(range = c(1, 6), name = "PSC count") +
  labs(x = expression("S"[ARG]), y = expression("S"[BM]), tag = "D") +
  theme_pub

fig <- (pA | pB) / (pC | pD)
ggsave(file.path(fdir, "Figure8_dimension_epidemiology.pdf"), fig, width = 10, height = 9,
       units = "in")
ggsave(file.path(fdir, "Figure8_dimension_epidemiology.png"), fig, width = 10, height = 9,
       units = "in", dpi = 300)
cat("Saved Figure8_dimension_epidemiology.pdf/png to", fdir, "\n")
