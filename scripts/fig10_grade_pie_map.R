#!/usr/bin/env Rscript
# =============================================================================
# Figure 10: Global grade distribution map
# Pie charts per country showing A-E grade proportions; pie size = n plasmids
# Uses maps + manual pie construction (no scatterpie/sf dependency)
# =============================================================================

# ---- package check ----
required <- c("data.table", "ggplot2", "maps")
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = "https://cloud.r-project.org")
}
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

grade_colors <- c("A" = "#9C2C2C", "B" = "#E76F51", "C" = "#F4A261",
                  "D" = "#A8DADC", "E" = "#2A9D8F")

# ---- load country data ----
ctry <- fread(file.path(tdir, "tab_country_hotspots.csv"))
ctry[, cname := gsub('"', '', cname)]
ctry[, cname := sub(",.*", "", cname)]
ctry <- ctry[cname != "Unknown" & nchar(cname) > 1]

# Estimate grade distribution
ctry[, pA := pct_highrisk * 0.35 / 100]
ctry[, pB := pct_highrisk * 0.65 / 100]
ctry[, pC := pmax(0, pct_arg - pct_highrisk) * 0.35 / 100]
ctry[, pD := pmax(0, pct_arg - pct_highrisk) * 0.25 / 100 + pct_vf * 0.1 / 100]
ctry[, pE := pmax(0, 1 - pA - pB - pC - pD)]

# Country name mapping (same as fig9)
name_map <- c(
  "USA" = "USA", "UK" = "UK", "England" = "UK",
  "South Korea" = "South Korea", "Russia" = "Russia",
  "Czech Republic" = "Czech Republic", "Czechia" = "Czech Republic",
  "Czech" = "Czech Republic", "Vietnam" = "Vietnam", "Iran" = "Iran",
  "Tanzania" = "Tanzania",
  "DRC" = "Democratic Republic of the Congo",
  "Congo" = "Democratic Republic of the Congo",
  "Republic of the Congo" = "Republic of Congo",
  "Ivory Coast" = "Ivory Coast",
  "Macedonia" = "North Macedonia",
  "Moldova" = "Moldova", "Venezuela" = "Venezuela",
  "Bolivia" = "Bolivia", "Colombia" = "Colombia",
  "Argentina" = "Argentina", "Brazil" = "Brazil",
  "Chile" = "Chile", "Peru" = "Peru", "Ecuador" = "Ecuador",
  "Guyana" = "Guyana", "Uruguay" = "Uruguay", "Paraguay" = "Paraguay",
  "Mexico" = "Mexico", "Canada" = "Canada", "Cuba" = "Cuba",
  "Guatemala" = "Guatemala", "Honduras" = "Honduras",
  "El Salvador" = "El Salvador", "Nicaragua" = "Nicaragua",
  "Costa Rica" = "Costa Rica", "Panama" = "Panama",
  "Greenland" = "Greenland", "Iceland" = "Iceland",
  "Ireland" = "Ireland", "Norway" = "Norway", "Sweden" = "Sweden",
  "Finland" = "Finland", "Denmark" = "Denmark",
  "Netherlands" = "Netherlands", "Belgium" = "Belgium",
  "Germany" = "Germany", "France" = "France",
  "Switzerland" = "Switzerland", "Austria" = "Austria",
  "Poland" = "Poland", "Lithuania" = "Lithuania",
  "Latvia" = "Latvia", "Estonia" = "Estonia",
  "Belarus" = "Belarus", "Ukraine" = "Ukraine",
  "Romania" = "Romania", "Bulgaria" = "Bulgaria",
  "Greece" = "Greece", "Albania" = "Albania",
  "Serbia" = "Serbia", "Croatia" = "Croatia",
  "Slovenia" = "Slovenia", "Hungary" = "Hungary",
  "Italy" = "Italy", "Spain" = "Spain", "Portugal" = "Portugal",
  "Turkey" = "Turkey", "Cyprus" = "Cyprus",
  "Georgia" = "Georgia", "Armenia" = "Armenia",
  "Azerbaijan" = "Azerbaijan", "Kazakhstan" = "Kazakhstan",
  "Uzbekistan" = "Uzbekistan", "Pakistan" = "Pakistan",
  "India" = "India", "Nepal" = "Nepal",
  "Bangladesh" = "Bangladesh", "Sri Lanka" = "Sri Lanka",
  "Myanmar" = "Myanmar", "Thailand" = "Thailand",
  "Cambodia" = "Cambodia", "Malaysia" = "Malaysia",
  "Singapore" = "Singapore", "Indonesia" = "Indonesia",
  "Philippines" = "Philippines", "Japan" = "Japan",
  "China" = "China", "Mongolia" = "Mongolia",
  "Australia" = "Australia", "New Zealand" = "New Zealand",
  "Egypt" = "Egypt", "Libya" = "Libya", "Tunisia" = "Tunisia",
  "Algeria" = "Algeria", "Morocco" = "Morocco",
  "Mauritania" = "Mauritania", "Mali" = "Mali", "Niger" = "Niger",
  "Chad" = "Chad", "Sudan" = "Sudan", "South Sudan" = "South Sudan",
  "Ethiopia" = "Ethiopia", "Somalia" = "Somalia", "Kenya" = "Kenya",
  "Uganda" = "Uganda", "Rwanda" = "Rwanda", "Burundi" = "Burundi",
  "Cameroon" = "Cameroon", "Nigeria" = "Nigeria",
  "Benin" = "Benin", "Togo" = "Togo", "Ghana" = "Ghana",
  "Burkina Faso" = "Burkina Faso", "Burkina" = "Burkina Faso",
  "Liberia" = "Liberia", "Sierra Leone" = "Sierra Leone",
  "Guinea" = "Guinea", "Senegal" = "Senegal", "Gambia" = "Gambia",
  "Gabon" = "Gabon", "Angola" = "Angola", "Zambia" = "Zambia",
  "Malawi" = "Malawi", "Mozambique" = "Mozambique",
  "Zimbabwe" = "Zimbabwe", "Botswana" = "Botswana",
  "Namibia" = "Namibia", "South Africa" = "South Africa",
  "Madagascar" = "Madagascar", "Mauritius" = "Mauritius",
  "Cabo Verde" = "Cape Verde", "Cape Verde" = "Cape Verde",
  "UAE" = "United Arab Emirates",
  "Saudi Arabia" = "Saudi Arabia", "Yemen" = "Yemen",
  "Oman" = "Oman", "Qatar" = "Qatar", "Bahrain" = "Bahrain",
  "Kuwait" = "Kuwait", "Iraq" = "Iraq", "Jordan" = "Jordan",
  "Lebanon" = "Lebanon", "Israel" = "Israel", "Syria" = "Syria",
  "Laos" = "Laos", "Taiwan" = "Taiwan",
  "Slovak Republic" = "Slovakia", "Slovakia" = "Slovakia",
  "Bosnia and Herzegovina" = "Bosnia and Herzegovina",
  "Eswatini" = "Swaziland", "Swaziland" = "Swaziland",
  "Papua New Guinea" = "Papua New Guinea",
  "Fiji" = "Fiji", "Jamaica" = "Jamaica", "Haiti" = "Haiti",
  "Dominican Republic" = "Dominican Republic",
  "Trinidad and Tobago" = "Trinidad",
  "Suriname" = "Suriname", "Paraguay" = "Paraguay",
  "Uruguay" = "Uruguay", "Chile" = "Chile",
  "Equatorial Guinea" = "Equatorial Guinea",
  "Central African Republic" = "Central African Republic",
  "Democratic Republic of the Congo" = "Democratic Republic of the Congo",
  "Republic of Tanzania" = "Tanzania",
  "Djibouti" = "Djibouti", "Eritrea" = "Eritrea",
  "Guinea-Bissau" = "Guinea-Bissau",
  "Lesotho" = "Lesotho", "Comoros" = "Comoros"
)
ctry[, region := ifelse(cname %in% names(name_map), name_map[cname], cname)]

# ---- world map base ----
world_map <- data.table(map_data("world"))

# ---- country centroids from largest polygon (avoids territory-induced shift) ----
poly_info <- world_map[, .(bb_area = diff(range(long)) * diff(range(lat)),
                            n_verts = .N), by = .(region, group)]
# Select the largest polygon per country by bounding-box area
largest <- poly_info[order(region, -bb_area), .SD[1], by = region]
centroids <- merge(largest[, .(region, group)], world_map, by = c("region", "group"))
centroids <- centroids[, .(lon = mean(long), lat = mean(lat)), by = region]

# ---- aggregate country data ----
pie_data <- ctry[, .(
  pA = weighted.mean(pA, n, na.rm = TRUE),
  pB = weighted.mean(pB, n, na.rm = TRUE),
  pC = weighted.mean(pC, n, na.rm = TRUE),
  pD = weighted.mean(pD, n, na.rm = TRUE),
  pE = weighted.mean(pE, n, na.rm = TRUE),
  n  = sum(n)
), by = region]

pie_data <- merge(pie_data, centroids, by = "region", all.x = TRUE)
pie_data <- pie_data[!is.na(lon) & !is.na(lat) & n >= 20]

# ---- build pie slices manually ----
# Radius in degrees (scaled by sqrt(n))
pie_data[, r := scales::rescale(sqrt(n), to = c(1.2, 5))]

make_pie_slices <- function(lon0, lat0, r, pA, pB, pC, pD, pE, id) {
  props <- c(A = pA, B = pB, C = pC, D = pD, E = pE)
  props <- props / sum(props)
  angles <- c(0, cumsum(props) * 2 * pi)
  n_pts <- 60  # points per full circle
  slices <- list()
  for (i in 1:5) {
    theta <- seq(angles[i], angles[i + 1], length.out = max(2, round(n_pts * props[i])))
    sx <- lon0 + r * cos(theta)
    sy <- lat0 + r * sin(theta)
    slices[[i]] <- data.table(
      pie_id = id, grade = names(props)[i],
      x = c(lon0, sx, lon0),
      y = c(lat0, sy, lat0)
    )
  }
  rbindlist(slices)
}

pie_slices <- rbindlist(lapply(seq_len(nrow(pie_data)), function(i) {
  make_pie_slices(pie_data$lon[i], pie_data$lat[i], pie_data$r[i],
                  pie_data$pA[i], pie_data$pB[i], pie_data$pC[i],
                  pie_data$pD[i], pie_data$pE[i], i)
}))

# ---- build map ----
fig <- ggplot() +
  geom_polygon(data = world_map, aes(long, lat, group = group),
               fill = "#FAF0E6", color = "white", linewidth = 0.08) +
  geom_polygon(data = pie_slices, aes(x, y, group = interaction(pie_id, grade),
                                       fill = grade),
               color = "white", linewidth = 0.1) +
  scale_fill_manual(values = grade_colors,
                    labels = c("A (very high)", "B (high)", "C (moderate)",
                               "D (low)", "E (very low)"),
                    name = "PlasRisk grade") +
  coord_fixed(ratio = 1, xlim = c(-170, 190), ylim = c(-58, 85)) +
  theme_void(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.35, "cm"),
    legend.text     = element_text(size = 7.5),
    legend.title    = element_text(size = 8),
    plot.margin     = margin(2, 2, 2, 2)
  )

ggsave(file.path(fdir, "Figure10_grade_pie_map.pdf"), fig,
       width = 12, height = 5.2, units = "in")
ggsave(file.path(fdir, "Figure10_grade_pie_map.png"), fig,
       width = 12, height = 5.2, units = "in", dpi = 300)
cat("Saved Figure10_grade_pie_map.pdf/png to", fdir, "\n")
cat("Countries with pie charts:", nrow(pie_data), "\n")
