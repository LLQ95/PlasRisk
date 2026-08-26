#!/usr/bin/env Rscript
# =============================================================================
# Figure 10: World map with grade-distribution pies per country
# Uses 'maps' polygons plus hand-computed trigonometric pie grobs so that the
# figure renders on R 4.1.x without scatterpie/sf. Input:
# tab_country_grades.csv with columns country, A,B,C,D,E (counts or fractions).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(maps)
})

args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
tdir   <- file.path(resdir, "tables")
fdir   <- file.path(resdir, "figures")
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)

grade_cols <- c(A = "#9C2C2C", B = "#E76F51", C = "#F4A261",
                D = "#A8DADC", E = "#2A9D8F")

world <- as.data.table(map_data("world"))
# country centroids (long, lat) from maps::map where possible, else manual
cent <- as.data.table(map("world", plot = FALSE, fill = TRUE)[c("names", "x", "y")])

cfile <- file.path(tdir, "tab_country_grades.csv")
if (!file.exists(cfile)) stop("tab_country_grades.csv not found")
cg <- fread(cfile)

# normalize grade columns to proportions
gcols <- c("A", "B", "C", "D", "E")
gcols <- gcols[gcols %in% names(cg)]
cg[, total := rowSums(.SD), .SDcols = gcols]
cg <- cg[total > 0]
cg[, (gcols) := lapply(.SD, function(x) x / total), .SDcols = gcols]

# match country names to centroids
name_map <- c("USA" = "USA", "United States" = "USA", "UK" = "UK",
              "United Kingdom" = "UK", "South Korea" = "South Korea",
              "Republic of Korea" = "South Korea")
cg[, nm := ifelse(country %in% names(name_map), name_map[country], country)]
cent[, country := gsub(":.*", "", names)]
cg <- merge(cg, cent[, .(country, x, y)], by.x = "nm", by.y = "country", all.x = TRUE)
cg <- cg[!is.na(x) & !is.na(y)]

# build pie wedges as polygons
radius <- 3.5  # degrees
make_pies <- function(d) {
  wedges <- list()
  for (i in seq_len(nrow(d))) {
    cum <- 0
    for (g in gcols) {
      a0 <- cum * 2 * pi - pi / 2
      cum <- cum + d[[g]][i]
      a1 <- cum * 2 * pi - pi / 2
      theta <- seq(a0, a1, length.out = pmax(3, ceiling((a1 - a0) * 12)))
      wedges[[length(wedges) + 1]] <- data.table(
        x = c(d$x[i], d$x[i] + radius * cos(theta) / 1.3, d$x[i]),
        y = c(d$y[i], d$y[i] + radius * sin(theta), d$y[i]),
        grade = g, idx = i)
    }
  }
  rbindlist(wedges)
}

# only label countries with enough plasmids to keep the map readable
ncol_ <- intersect(c("total", "n", "N"), names(cg))[1]
if (!is.na(ncol_)) cg <- cg[total >= quantile(total, 0.5)]
pies <- make_pies(cg)

p <- ggplot(world, aes(long, lat, group = group)) +
  geom_polygon(fill = "#EDE7E0", color = "white", linewidth = 0.1) +
  geom_polygon(data = pies, aes(x, y, group = interaction(idx, grade), fill = grade),
               inherit.aes = FALSE, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = grade_cols, name = "Grade") +
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-58, 85)) +
  theme_void(base_size = 9) +
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 11)) +
  labs(tag = "A")

ggsave(file.path(fdir, "Figure10_grade_pie_map.pdf"), p,
       width = 10, height = 6.2, units = "in")
ggsave(file.path(fdir, "Figure10_grade_pie_map.png"), p,
       width = 10, height = 6.2, units = "in", dpi = 300)
cat("Saved Figure10_grade_pie_map.pdf/png to", fdir, "\n")
cat(sprintf("Pies drawn for %d countries\n", nrow(cg)))
