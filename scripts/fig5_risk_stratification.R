#!/usr/bin/env Rscript
# =============================================================================
# Figure 5: Risk stratification and clinical utility
# (A) Grade distribution
# (B) Observed high-risk rate by grade
# (C) Calibration curves (full vs lite)
# (D) Decision-curve analysis across four outcomes
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
tdir   <- file.path(resdir, "tables")
fdir   <- file.path(resdir, "figures")
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)

grade_colors <- c("A" = "#9C2C2C", "B" = "#E76F51", "C" = "#F4A261",
                  "D" = "#A8DADC", "E" = "#2A9D8F")
outcome_colors <- c("High-risk ARG" = "#E76F51", "MDR-VF fusion" = "#F4A261",
                    "Conjugative" = "#2A9D8F", "Biocide/metal" = "#6D4C41")

theme_pub <- theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 12, color = "#3E2723"),
        axis.text = element_text(color = "#3E2723"),
        plot.margin = margin(4, 6, 4, 4))

psc <- fread(file.path(tdir, "tab_psc_final_scores.csv"))
gcol <- intersect(c("grade", "Grade", "final_grade"), names(psc))[1]
if (is.na(gcol)) stop("grade column not found")
setnames(psc, gcol, "grade")
psc[, grade := factor(grade, levels = c("A", "B", "C", "D", "E"))]

# ---- (A) grade distribution ----
gdist <- psc[, .N, by = grade][order(grade)]
gdist[, pct := N / sum(N) * 100]
pA <- ggplot(gdist, aes(grade, pct, fill = grade)) +
  geom_col(width = 0.75, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = grade_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Risk grade", y = "Proportion of PSCs (%)", tag = "A") +
  theme_pub

# ---- (B) observed high-risk rate by grade ----
grate <- psc[, .(obs = mean(y_highrisk) * 100, n = .N), by = grade][order(grade)]
pB <- ggplot(grate, aes(grade, obs, fill = grade)) +
  geom_col(width = 0.75, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", obs)), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = grade_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Risk grade", y = "Observed high-risk ARG rate (%)", tag = "B") +
  theme_pub

# ---- (C) calibration full vs lite ----
cal_file <- file.path(tdir, "tab_calibration.csv")
pC <- ggplot() + theme_void() + labs(tag = "C")
if (file.exists(cal_file)) {
  cal <- fread(cal_file)
  if (all(c("pred", "obs", "model") %in% names(cal))) {
    pC <- ggplot(cal, aes(pred, obs, color = model)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#BCAAA4") +
      geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
      scale_color_manual(values = c("full" = "#9C2C2C", "lite" = "#2A9D8F")) +
      labs(x = "Mean predicted risk", y = "Observed high-risk rate",
           color = NULL, tag = "C") +
      theme_pub
  }
}

# ---- (D) decision curves ----
dca <- fread(file.path(tdir, "tab_dca_all_outcomes.csv"))
pD <- ggplot(dca, aes(threshold)) +
  geom_line(aes(y = NB_model, color = outcome), linewidth = 0.7) +
  geom_line(aes(y = NB_treat_all), linetype = "dashed", color = "#BCAAA4", linewidth = 0.5) +
  scale_color_manual(values = outcome_colors) +
  scale_x_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "Risk threshold", y = "Net benefit", color = NULL, tag = "D") +
  theme_pub

fig <- (pA | pB) / (pC | pD)
ggsave(file.path(fdir, "Figure5_risk_stratification.pdf"), fig, width = 7.2, height = 7,
       units = "in")
ggsave(file.path(fdir, "Figure5_risk_stratification.png"), fig, width = 7.2, height = 7,
       units = "in", dpi = 300)
cat("Saved Figure5_risk_stratification.pdf/png to", fdir, "\n")
