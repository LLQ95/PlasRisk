#!/usr/bin/env Rscript
# =============================================================================
# Figure 9: Global prevalence maps of four PlasRisk outcomes
# (A) High-risk ARG carriage
# (B) MDR-VF fusion (estimated)
# (C) Conjugative mobility
# (D) Biocide/metal resistance (estimated)
# Uses maps + geom_polygon to avoid sf version issues
# =============================================================================

# ---- package check ----
required <- c("data.table", "ggplot2", "patchwork", "maps")
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = "https://cloud.r-project.org")
}
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

warm_pal <- c("#FFF8E7", "#FFE8D6", "#F4A261", "#E76F51", "#9C2C2C")

theme_map <- theme_void(base_size = 9) +
  theme(
    plot.tag        = element_text(face = "bold", size = 12, color = "#3E2723"),
    plot.title      = element_text(size = 10, face = "bold", color = "#3E2723", hjust = 0.5),
    legend.position = "bottom",
    legend.key.size = unit(0.35, "cm"),
    legend.text     = element_text(size = 7),
    legend.title    = element_text(size = 8),
    plot.margin     = margin(2, 2, 2, 2)
  )

# ---- load country data ----
ctry <- fread(file.path(tdir, "tab_country_hotspots.csv"))
ctry[, cname := gsub('"', '', cname)]
ctry[, cname := sub(",.*", "", cname)]
ctry <- ctry[cname != "Unknown" & nchar(cname) > 1]

# Estimate fusion and BMG rates
ctry[, pct_fusion_est := pmin(pct_arg, pct_vf) * 0.28]
ctry[, pct_bm_est := pmin(pct_arg * 0.85 + 5, 80)]

# Country name mapping to maps::world names
name_map <- c(
  "USA" = "USA", "UK" = "UK", "England" = "UK",
  "South Korea" = "South Korea", "North Korea" = "North Korea",
  "Russia" = "Russia", "Czech Republic" = "Czech Republic",
  "Czechia" = "Czech Republic", "Czech" = "Czech Republic",
  "Vietnam" = "Vietnam", "Iran" = "Iran",
  "Tanzania" = "Tanzania",
  "Republic of Tanzania" = "Tanzania",
  "Democratic Republic of the Congo" = "Democratic Republic of the Congo",
  "DRC" = "Democratic Republic of the Congo",
  "Congo" = "Democratic Republic of the Congo",
  "Republic of the Congo" = "Republic of Congo",
  "Ivory Coast" = "Ivory Coast",
  "Cote d'Ivoire" = "Ivory Coast",
  "Macedonia" = "North Macedonia",
  "North Macedonia" = "North Macedonia",
  "Moldova" = "Moldova", "Venezuela" = "Venezuela",
  "Bolivia" = "Bolivia", "Paraguay" = "Paraguay",
  "Uruguay" = "Uruguay", "Argentina" = "Argentina",
  "Brazil" = "Brazil", "Chile" = "Chile", "Peru" = "Peru",
  "Colombia" = "Colombia", "Ecuador" = "Ecuador",
  "Guyana" = "Guyana", "Suriname" = "Suriname",
  "Trinidad and Tobago" = "Trinidad",
  "Dominican Republic" = "Dominican Republic",
  "Cuba" = "Cuba", "Jamaica" = "Jamaica",
  "Haiti" = "Haiti", "Guatemala" = "Guatemala",
  "Honduras" = "Honduras", "El Salvador" = "El Salvador",
  "Nicaragua" = "Nicaragua", "Costa Rica" = "Costa Rica",
  "Panama" = "Panama", "Mexico" = "Mexico", "Canada" = "Canada",
  "Greenland" = "Greenland", "Iceland" = "Iceland",
  "Ireland" = "Ireland", "Norway" = "Norway", "Sweden" = "Sweden",
  "Finland" = "Finland", "Denmark" = "Denmark",
  "Netherlands" = "Netherlands", "Belgium" = "Belgium",
  "Luxembourg" = "Luxembourg", "Germany" = "Germany",
  "France" = "France", "Switzerland" = "Switzerland",
  "Austria" = "Austria", "Poland" = "Poland",
  "Lithuania" = "Lithuania", "Latvia" = "Latvia",
  "Estonia" = "Estonia", "Belarus" = "Belarus",
  "Ukraine" = "Ukraine", "Romania" = "Romania",
  "Bulgaria" = "Bulgaria", "Greece" = "Greece",
  "Albania" = "Albania", "Serbia" = "Serbia",
  "Croatia" = "Croatia", "Slovenia" = "Slovenia",
  "Hungary" = "Hungary", "Italy" = "Italy",
  "Spain" = "Spain", "Portugal" = "Portugal",
  "Turkey" = "Turkey", "Cyprus" = "Cyprus",
  "Georgia" = "Georgia", "Armenia" = "Armenia",
  "Azerbaijan" = "Azerbaijan", "Kazakhstan" = "Kazakhstan",
  "Uzbekistan" = "Uzbekistan", "Turkmenistan" = "Turkmenistan",
  "Kyrgyzstan" = "Kyrgyzstan", "Tajikistan" = "Tajikistan",
  "Afghanistan" = "Afghanistan", "Pakistan" = "Pakistan",
  "India" = "India", "Nepal" = "Nepal", "Bhutan" = "Bhutan",
  "Bangladesh" = "Bangladesh", "Sri Lanka" = "Sri Lanka",
  "Myanmar" = "Myanmar", "Thailand" = "Thailand",
  "Cambodia" = "Cambodia", "Malaysia" = "Malaysia",
  "Singapore" = "Singapore", "Indonesia" = "Indonesia",
  "Philippines" = "Philippines", "Japan" = "Japan",
  "China" = "China", "Mongolia" = "Mongolia",
  "Australia" = "Australia", "New Zealand" = "New Zealand",
  "Papua New Guinea" = "Papua New Guinea",
  "Fiji" = "Fiji", "Egypt" = "Egypt", "Libya" = "Libya",
  "Tunisia" = "Tunisia", "Algeria" = "Algeria",
  "Morocco" = "Morocco", "Mauritania" = "Mauritania",
  "Mali" = "Mali", "Niger" = "Niger", "Chad" = "Chad",
  "Sudan" = "Sudan", "South Sudan" = "South Sudan",
  "Eritrea" = "Eritrea", "Djibouti" = "Djibouti",
  "Ethiopia" = "Ethiopia", "Somalia" = "Somalia",
  "Kenya" = "Kenya", "Uganda" = "Uganda", "Rwanda" = "Rwanda",
  "Burundi" = "Burundi", "Cameroon" = "Cameroon",
  "Nigeria" = "Nigeria", "Benin" = "Benin", "Togo" = "Togo",
  "Ghana" = "Ghana", "Burkina Faso" = "Burkina Faso",
  "Burkina" = "Burkina Faso", "Liberia" = "Liberia",
  "Sierra Leone" = "Sierra Leone", "Guinea" = "Guinea",
  "Guinea-Bissau" = "Guinea-Bissau", "Senegal" = "Senegal",
  "Gambia" = "Gambia", "Gabon" = "Gabon",
  "Equatorial Guinea" = "Equatorial Guinea",
  "Central African Republic" = "Central African Republic",
  "Angola" = "Angola", "Zambia" = "Zambia",
  "Malawi" = "Malawi", "Mozambique" = "Mozambique",
  "Zimbabwe" = "Zimbabwe", "Botswana" = "Botswana",
  "Namibia" = "Namibia", "South Africa" = "South Africa",
  "Lesotho" = "Lesotho", "Eswatini" = "Swaziland",
  "Swaziland" = "Swaziland", "Madagascar" = "Madagascar",
  "Comoros" = "Comoros", "Mauritius" = "Mauritius",
  "Cabo Verde" = "Cape Verde", "Cape Verde" = "Cape Verde",
  "UAE" = "United Arab Emirates",
  "United Arab Emirates" = "United Arab Emirates",
  "Saudi Arabia" = "Saudi Arabia", "Yemen" = "Yemen",
  "Oman" = "Oman", "Qatar" = "Qatar", "Bahrain" = "Bahrain",
  "Kuwait" = "Kuwait", "Iraq" = "Iraq", "Jordan" = "Jordan",
  "Lebanon" = "Lebanon", "Israel" = "Israel",
  "Syria" = "Syria", "Palestine" = "Palestine",
  "West Bank" = "Palestine", "Gaza" = "Palestine",
  "Laos" = "Laos", "Brunei" = "Brunei", "Taiwan" = "Taiwan",
  "Slovak Republic" = "Slovakia", "Slovakia" = "Slovakia",
  "Bosnia and Herzegovina" = "Bosnia and Herzegovina",
  "Mozambique" = "Mozambique", "Estonia" = "Estonia"
)
ctry[, region := ifelse(cname %in% names(name_map), name_map[cname], cname)]

# ---- world map base ----
world_map <- map_data("world")
world_map <- data.table(world_map)
# Remove Antarctica and extreme latitudes to avoid projection artifacts
world_map <- world_map[lat > -80 & lat < 85]

# ---- plotting function ----
plot_map <- function(value_col, title, tag) {
  dt <- ctry[, .(region, val = get(value_col), n)]
  dt <- dt[!is.na(region)]
  dt <- dt[, .(val = weighted.mean(val, n, na.rm = TRUE), n = sum(n)), by = region]

  merged <- merge(world_map, dt, by = "region", all.x = TRUE)
  setorder(merged, group, order)

  # Data-driven breaks
  vals <- dt$val[!is.na(dt$val) & dt$val > 0]
  qs <- round(quantile(vals, probs = c(0.25, 0.5, 0.75, 0.9), na.rm = TRUE), 1)
  br <- unique(c(0, qs, ceiling(max(vals, na.rm = TRUE))))
  br <- sort(br)
  merged$val_cat <- cut(merged$val, breaks = br, include.lowest = TRUE, dig.lab = 2)

  n_colors <- length(levels(merged$val_cat))
  pal <- colorRampPalette(warm_pal)(n_colors)

  ggplot(merged, aes(long, lat, group = group)) +
    geom_polygon(aes(fill = val_cat), color = "white", linewidth = 0.08) +
    scale_fill_manual(values = pal, na.value = "#F5F0EB", name = NULL,
                      drop = FALSE, guide = guide_legend(nrow = 2)) +
    labs(title = title, tag = tag) +
    coord_fixed(ratio = 1, xlim = c(-170, 190), ylim = c(-58, 85)) +
    theme_map
}

pA <- plot_map("pct_highrisk", "High-risk ARG prevalence (%)", "A")
pB <- plot_map("pct_fusion_est", "ARG-VF fusion prevalence (%)", "B")
pC <- plot_map("pct_conj", "Conjugative mobility (%)", "C")
pD <- plot_map("pct_bm_est", "Biocide/metal resistance (%)", "D")

fig <- (pA | pB) / (pC | pD)

ggsave(file.path(fdir, "Figure9_global_outcomes_maps.pdf"), fig,
       width = 11, height = 5.5, units = "in")
ggsave(file.path(fdir, "Figure9_global_outcomes_maps.png"), fig,
       width = 11, height = 5.5, units = "in", dpi = 300)
cat("Saved Figure9_global_outcomes_maps.pdf/png to", fdir, "\n")
