#!/usr/bin/env Rscript
# =============================================================================
# Figure 9: Global country-level prevalence maps for the four outcomes
# Uses the 'maps' package (geom_polygon) so that sf/scatterpie are not required
# on older R installations. Input: tab_country_outcomes.csv with columns
# country, highrisk_rate, fusion_rate, conj_rate, bm_rate (rates in [0,1]).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(maps)
})

args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
tdir   <- file.path(resdir, "tables")
fdir   <- file.path(resdir, "figures")
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)

world <- as.data.table(map_data("world"))

# ISO-name fallback: map PIPdb country strings to map::world region names
name_map <- c(
  "USA" = "USA", "United States" = "USA", "UK" = "UK",
  "United Kingdom" = "UK", "South Korea" = "South Korea",
  "Republic of Korea" = "South Korea", "Czech Republic" = "Czech Republic",
  "Democratic Republic of the Congo" = "Democratic Republic of the Congo",
  "Republic of the Congo" = "Republic of Congo",
  "Ivory Coast" = "Ivory Coast", "Myanmar" = "Myanmar"
)

cfile <- file.path(tdir, "tab_country_outcomes.csv")
if (!file.exists(cfile)) {
  # graceful fallback: empty panel so the script does not fail
  cdat <- data.table(country = character(), rate = numeric(), outcome = character())
} else {
  cin <- fread(cfile)
  rate_cols <- c(highrisk_rate = "highrisk_rate", fusion_rate = "fusion_rate",
                 conj_rate = "conj_rate", bm_rate = "bm_rate")
  rate_cols <- rate_cols[rate_cols %in% names(cin)]
  cdat <- melt(cin, id.vars = "country", measure.vars = rate_cols,
               variable.name = "outcome", value.name = "rate")
  cdat[, region := ifelse(country %in% names(name_map), name_map[country], country)]
}

outcome_labs <- c(highrisk_rate = "High-risk ARG",
                  fusion_rate   = "MDR-VF fusion",
                  conj_rate     = "Conjugative",
                  bm_rate       = "Biocide/metal")
outcome_cols <- c("#E76F51", "#F4A261", "#2A9D8F", "#6D4C41")

theme_map <- theme_void(base_size = 9) +
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 11),
        plot.title = element_text(size = 9, hjust = 0.5),
        legend.key.width = unit(0.8, "cm"))

make_map <- function(var, col, lab, tag) {
  d <- cdat[outcome == var]
  m <- merge(world, d, by = "region", all.x = TRUE)
  setorder(m, order)
  # crop spurious horizontal lines by clipping longitudes to [-180,180]
  ggplot(m, aes(long, lat, group = group)) +
    geom_polygon(aes(fill = rate * 100), color = "white", linewidth = 0.1) +
    scale_fill_gradient2(low = "#FFF8E7", mid = col, high = "#9C2C2C",
                         midpoint = if (var == "bm_rate") 35 else 12,
                         na.value = "#EDE7E0", name = "%") +
    coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-58, 85)) +
    labs(title = lab, tag = tag) + theme_map
}

panels <- list(
  make_map("highrisk_rate", "#E76F51", "High-risk ARG", "A"),
  make_map("fusion_rate",   "#F4A261", "MDR-VF fusion", "B"),
  make_map("conj_rate",     "#2A9D8F", "Conjugative", "C"),
  make_map("bm_rate",       "#6D4C41", "Biocide/metal", "D")
)
fig <- wrap_plots(panels, ncol = 2)

ggsave(file.path(fdir, "Figure9_global_outcomes_maps.pdf"), fig,
       width = 9, height = 6.5, units = "in")
ggsave(file.path(fdir, "Figure9_global_outcomes_maps.png"), fig,
       width = 9, height = 6.5, units = "in", dpi = 300)
cat("Saved Figure9_global_outcomes_maps.pdf/png to", fdir, "\n")
