#!/usr/bin/env Rscript
# =============================================================================
# pipdb_13_descriptive_epidemiology.R
# =============================================================================
# Comprehensive descriptive epidemiology & visualization of PIPdb plasmids.
# Produces publication-quality figures + summary statistics tables.
#
# Usage:
#   Rscript pipdb_13_descriptive_epidemiology.R [results_dir]
#
# Output: results/figures/descriptive/  (PDF + PNG, 300 dpi)
#         results/tables/               (CSV summary tables)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(ggrepel)
  library(scales)
  library(dplyr)
  library(tidyr)
  library(maps)
})

# ---- args ----
res_dir <- if (length(commandArgs(trailingOnly=TRUE))>=1) commandArgs(trailingOnly=TRUE)[1] else "results"
fig_dir <- file.path(res_dir, "figures", "descriptive")
tab_dir <- file.path(res_dir, "tables")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)

cat("Loading data ...\n")

# ---- load & merge ----
master_cols <- c("plasmid_acc","replicon_primary","pmlst","length_avg","length_min","length_max",
  "isolate_mark","host_rank2","country","n_country","year_min","year_max",
  "year_mid","single_year","species_name","genus_name","phylum_name","gram_stain",
  "aro_name","n_arg","n_who_arg","n_drugclass","drugclass",
  "vf_name","vf_category","n_vf","mobility_class","n_is","n_is_family","is_density_per_kb",
  "n_metal","annual_growth_rate","combined_risk_index","has_integron",
  "host_class","host_human","habitat_human","habitat_animal","habitat_env")
risk_cols <- c("plasmid_acc","n_high_risk_arg","S_ARG","S_MOB","S_HOST","S_REP","S_SIZE",
               "S_GEO","S_HAB","S_GROW","risk_score","Q_grade")

m <- fread(file.path(res_dir,"psc_master.tsv"), select=master_cols, na.strings=c("\\N","-","","NA"))
r <- fread(file.path(res_dir,"psc_risk_scores.tsv"), select=risk_cols)
d <- merge(m, r, by="plasmid_acc", all.x=TRUE)
rm(m,r); gc()

# ---- clean ----
d[, length_avg := as.numeric(length_avg)]
d[, year_mid := as.numeric(year_mid)]
d[, n_arg := as.numeric(n_arg)]
d[, n_vf := as.numeric(n_vf)]
d[, n_country := as.numeric(n_country)]
d[, annual_growth_rate := as.numeric(annual_growth_rate)]
d[, n_is := as.numeric(n_is)]
d[, n_metal := as.numeric(n_metal)]
d[, is_density_per_kb := as.numeric(is_density_per_kb)]
d <- d[!is.na(year_mid) & year_mid>=1920 & year_mid<=2025]
d <- d[!is.na(length_avg) & length_avg>=500]
bool_cols <- c("has_integron","host_human","habitat_human","habitat_animal","habitat_env")
for (bc in bool_cols) {
  if (bc %in% names(d) && is.character(d[[bc]]))
    d[, (bc) := toupper(substr(get(bc),1,1))=="T"]
}
d[, era := cut(year_mid, breaks=c(1920,1970,1990,2000,2010,2015,2020,2026),
               labels=c("<=1970","1971-1990","1991-2000","2001-2010","2011-2015","2016-2020","2021-2025"),
               right=TRUE, include.lowest=TRUE)]
rep_tab <- d[!is.na(replicon_primary), .N, by=replicon_primary][order(-N)]
top_reps <- rep_tab$replicon_primary[1:min(20,nrow(rep_tab))]
top20 <- top_reps
d[, rep_top := ifelse(replicon_primary %in% top_reps, replicon_primary, "Other")]
d[, rep_top := factor(rep_top, levels=c(top_reps,"Other"))]

cat(sprintf("  %s PSCs after filtering\n", format(nrow(d), big.mark=",")))

# ---- theme ----
theme_pub <- theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="grey90"),
        legend.key.size=unit(0.4,"cm"),
        plot.title=element_text(face="bold", size=12))
palette_rep <- c(brewer.pal(12,"Set3"), brewer.pal(8,"Set2"), brewer.pal(8,"Set1"), "grey80")
names(palette_rep) <- c(top_reps,"Other")

save_plot <- function(p, name, w=10, h=7) {
  ggsave(file.path(fig_dir,paste0(name,".pdf")), p, width=w, height=h, units="in")
  ggsave(file.path(fig_dir,paste0(name,".png")), p, width=w, height=h, units="in", dpi=300)
  cat(sprintf("  saved %s\n", name))
}

# =============================================================================
# FIG 1: Temporal trends
# =============================================================================
cat("\n=== Fig 1: Temporal trends ===\n")
era_tot <- d[, .(total=.N), by=era]
era_rep <- d[replicon_primary %in% top_reps, .N, by=.(era,replicon_primary)]
era_rep <- merge(era_rep, era_tot, by="era")
era_rep[, prop := N/total]
sig_tests <- rbindlist(lapply(top_reps, function(rp) {
  tab <- table(d$era[d$replicon_primary %in% top_reps],
               d$replicon_primary[d$replicon_primary %in% top_reps] == rp)
  if(ncol(tab)==2 && nrow(tab)>=2) {
    ft <- fisher.test(tab, simulate.p.value=TRUE, B=10000)
    chisq <- suppressWarnings(chisq.test(tab))
    data.table(replicon_primary=rp, p_fisher=ft$p.value, p_chisq=chisq$p.value)
  } else data.table(replicon_primary=rp, p_fisher=NA_real_, p_chisq=NA_real_)
}))
sig_labels <- c("***","**","*","")
sig_tests[, sig := ifelse(is.na(p_fisher), "",
  sig_labels[pmin(4, findInterval(p_fisher, c(0.001,0.01,0.05,1))+1)])]
fwrite(sig_tests, file.path(tab_dir,"tab_temporal_significance.csv"))
trends <- rbindlist(lapply(top_reps, function(rp) {
  dat <- d[replicon_primary %in% top_reps, .(year_mid, is_rep = replicon_primary==rp)]
  fit <- glm(is_rep ~ year_mid, data=dat, family=binomial)
  data.table(replicon_primary=rp, slope=coef(fit)[2], p=summary(fit)$coefficients[2,4], OR=exp(coef(fit)[2]))
}))
fwrite(trends, file.path(tab_dir,"tab_temporal_trends.csv"))
p1a <- ggplot(era_rep, aes(x=era, y=prop, fill=replicon_primary, group=replicon_primary)) +
  geom_area(position="stack", alpha=0.85, colour="white", linewidth=0.2) +
  scale_fill_manual(values=palette_rep, name="Replicon") +
  scale_y_continuous(labels=percent_format(), expand=c(0,0)) +
  labs(x="Collection era", y="Proportion of PSCs", title="Temporal dynamics of major plasmid replicons") +
  theme_pub + theme(legend.position="right")
top12 <- rep_tab$replicon_primary[1:12]
sig_lab <- sig_tests[match(top12, replicon_primary)]
rep12_labels <- setNames(paste0(top12, sig_lab$sig), top12)
p1b <- ggplot(era_rep[replicon_primary %in% top12],
              aes(x=era, y=prop, colour=replicon_primary, group=replicon_primary)) +
  geom_line(linewidth=0.8) + geom_point(size=1.8) +
  scale_colour_manual(values=palette_rep[top12], name="Replicon (Fisher p)", labels=rep12_labels) +
  scale_y_continuous(labels=percent_format()) +
  labs(x="Collection era", y="Proportion", title="Per-replicon trends") +
  theme_pub + theme(legend.position="right")
save_plot(p1a/p1b, "fig01_temporal_trends", 10, 10)

# =============================================================================
# FIG 2: World map with replicon pie charts by country
# =============================================================================
cat("\n=== Fig 2: Geographic distribution ===\n")
cty <- d[!is.na(country) & country!="" & !is.na(replicon_primary), .N, by=.(country,replicon_primary)]
cty_tot <- cty[, .(total=sum(N)), by=country]
cty <- merge(cty, cty_tot, by="country")
major_countries <- cty_tot[total>=100][order(-total)]$country
cat(sprintf("  %d countries with >=100 PSCs\n", length(major_countries)))
cty[, country := sub("^U\\.K\\.$","UK",country)]
cty[, country := sub("^United States$","USA",country)]
centroids <- data.table(
  country=c("China","USA","UK","Germany","France","Japan","India","Brazil","Australia","Canada","Italy","Spain","Netherlands",
    "South Korea","Denmark","Switzerland","Belgium","Sweden","Thailand","Vietnam","Egypt","Nigeria","Kenya","South Africa","Mexico","Argentina",
    "Russia","Poland","Portugal","Norway","Finland","Ireland","Austria","Greece","Turkey","Iran","Pakistan","Bangladesh","Indonesia","Malaysia",
    "Philippines","Singapore","New Zealand","Colombia","Chile","Peru","Czech Republic","Hungary","Israel","Saudi Arabia","Taiwan","Hong Kong",
    "Iceland","Croatia","Slovenia","Lithuania","Estonia","Latvia","Ukraine","Romania","Bulgaria","Serbia","Slovakia","Luxembourg"),
  lon=c(104,-95,-2,10,2,138,78,-51,134,-106,12,-3,5,127,10,8,4,16,101,108,30,8,37,25,-102,-64,
    90,19,-8,10,25,-8,14,22,35,53,69,90,113,102,122,104,174,-74,-71,-76,15,19,35,45,121,114,
    -19,16,15,24,26,24,32,25,25,21,19,6),
  lat=c(35,38,54,51,46,36,22,-10,-25,56,42,40,52,36,56,47,50,62,15,16,27,9,-1,-29,23,-34,
    61,52,39,62,64,53,47,39,39,33,30,24,-2,4,13,1,-41,4,-35,-10,50,47,31,24,24,22,
    65,45,46,55,59,57,49,46,43,44,48,50))
centroids <- unique(centroids, by="country")
name_map <- c("USA"="USA","United States"="USA","U.K."="UK","United Kingdom"="UK",
  "South Korea"="South Korea","Korea"="South Korea","Czech Republic"="Czech Republic","Czechia"="Czech Republic",
  "Russia"="Russia","Russian Federation"="Russia","Vietnam"="Vietnam","Viet Nam"="Vietnam")
cty[, cname := ifelse(country %in% names(name_map), name_map[country], country)]
if (requireNamespace("sf", quietly=TRUE)) {
  rdata_pre <- "world_maps.RData"
  if (!file.exists(rdata_pre)) {
    for (p in c("../world_maps.RData","/db/student/metagenome/pipdb/contigs/world_maps.RData")) {
      if (file.exists(p)) { rdata_pre <- p; break }
    }
  }
  if (file.exists(rdata_pre)) {
    e <- new.env(); load(rdata_pre, envir=e)
    if ("world_GBD" %in% ls(e)) {
      wsf <- e$world_GBD
      missing_cty <- setdiff(unique(cty[country %in% major_countries]$cname), centroids$country)
      if (length(missing_cty)>0) {
        cat(sprintf("  Auto-computing centroids for %d countries...\n", length(missing_cty)))
        for (cn in missing_cty) {
          idx <- which(tolower(as.character(wsf[[2]]))==tolower(cn) | tolower(as.character(wsf[[1]]))==tolower(cn))
          if (length(idx)>0) {
            cent <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_union(wsf$geometry[idx]))))
            centroids <- rbind(centroids, data.table(country=cn, lon=as.numeric(cent[1,1]), lat=as.numeric(cent[1,2])))
          }
        }
      }
    }
  }
}
cty_map <- merge(cty[country %in% major_countries], centroids, by.x="cname", by.y="country", allow.cartesian=TRUE)
cty_map <- cty_map[!is.na(lon)]
top5_reps <- rep_tab$replicon_primary[1:min(5,nrow(rep_tab))]
cty_map[, rep5 := ifelse(replicon_primary %in% top5_reps, replicon_primary, "Others")]
cty_pie <- cty_map[, .(N=sum(N)), by=.(country,lon,lat,rep5)]
cty_pie[, total:=sum(N), by=country]
cty_pie[, frac:=N/total]
cty_pie[, rep5 := factor(rep5, levels=c(top5_reps,"Others"))]
world <- NULL; world_sf <- NULL
rdata_path <- "world_maps.RData"
if (!file.exists(rdata_path)) {
  for (p in c("../world_maps.RData","/db/student/metagenome/pipdb/contigs/world_maps.RData")) {
    if (file.exists(p)) { rdata_path <- p; break }
  }
}
if (file.exists(rdata_path)) {
  cat(sprintf("  Loading world map from: %s\n", rdata_path))
  objs <- load(rdata_path)
  cat(sprintf("  RData contains: %s\n", paste(objs, collapse=", ")))
  for (obj in objs) {
    o <- get(obj)
    if (inherits(o, "sf") || (is.data.frame(o) && "geometry" %in% names(o))) { world_sf <- o; break }
    if (is.data.frame(o)) {
      cn <- tolower(names(o))
      if (any(cn %in% c("long","longitude","lon","x")) && any(cn %in% c("lat","latitude","y"))) {
        world <- as.data.frame(o)
        for (i in seq_along(names(world))) {
          cl <- tolower(names(world)[i])
          if (cl %in% c("long","longitude","lon","x")) names(world)[i] <- "long"
          if (cl %in% c("lat","latitude","y")) names(world)[i] <- "lat"
          if (cl %in% c("group","piece")) names(world)[i] <- "group"
        }
        break
      }
    }
  }
  if (!is.null(world_sf) && requireNamespace("sf", quietly=TRUE)) {
    cat(sprintf("  Converting sf map to polygon data.frame...\n"))
    world_sf <- suppressWarnings(sf::st_make_valid(world_sf))
    geoms <- sf::st_geometry(world_sf)
    pieces <- vector("list", length(geoms)); k <- 0
    for (i in seq_along(geoms)) {
      g <- geoms[i]
      gtype <- as.character(sf::st_geometry_type(g))
      if (gtype == "GEOMETRYCOLLECTION") {
        g <- tryCatch(sf::st_collection_extract(g, "POLYGON"), error=function(e) NULL)
        if (is.null(g) || length(g)==0) next
        gtype <- as.character(sf::st_geometry_type(g))
      }
      if (!gtype %in% c("POLYGON","MULTIPOLYGON")) next
      if (sf::st_is_empty(g)) next
      coords <- sf::st_coordinates(g)
      if (nrow(coords)==0) next
      k <- k+1
      l2 <- if ("L2" %in% colnames(coords)) coords[,"L2"] else 0
      pieces[[k]] <- data.frame(long=as.numeric(coords[,"X"]), lat=as.numeric(coords[,"Y"]),
                                group=paste(i, coords[,"L1"], l2, sep="_"), stringsAsFactors=FALSE)
    }
    world <- do.call(rbind, pieces[seq_len(k)])
    world_sf <- NULL
    cat(sprintf("  Converted: %d polygons, %d coordinates\n", k, nrow(world)))
  }
}
if (is.null(world)) { cat("  Using ggplot2::map_data('world')\n"); world <- ggplot2::map_data("world") }

make_pie_polygons <- function(cty_pie, r_scale=35, n_arc=60) {
  cty_pie <- cty_pie[order(country, rep5)]
  cty_pie[, cum_start := (shift(cumsum(frac), fill=0)), by=country]
  cty_pie[, cum_end := cumsum(frac), by=country]
  pies <- vector("list", nrow(cty_pie))
  for (i in seq_len(nrow(cty_pie))) {
    cx <- cty_pie$lon[i]; cy <- cty_pie$lat[i]
    r  <- sqrt(cty_pie$total[i]) / r_scale
    a0 <- pi/2 - cty_pie$cum_end[i]   * 2*pi
    a1 <- pi/2 - cty_pie$cum_start[i] * 2*pi
    theta <- seq(a1, a0, length.out=n_arc)
    pies[[i]] <- data.frame(
      x = c(cx, cx + r*cos(theta), cx),
      y = c(cy, cy + r*sin(theta), cy),
      group = paste0(cty_pie$country[i],"_",cty_pie$rep5[i]),
      category = as.character(cty_pie$rep5[i]),
      country = cty_pie$country[i], stringsAsFactors=FALSE)
  }
  do.call(rbind, pies)
}
pie_polys <- make_pie_polygons(cty_pie, r_scale=35)
cat(sprintf("  Pie polygons: %d slices for %d countries\n", nrow(cty_pie), uniqueN(cty_pie$country)))
pie5_cols <- palette_rep[top5_reps]; pie5_cols["Others"] <- "grey60"
size_vals <- c(1000, 5000, 10000, 20000)
size_key <- do.call(rbind, lapply(size_vals, function(n) {
  r <- sqrt(n)/35; theta <- seq(0, 2*pi, length.out=80)
  data.frame(x=-155 + r*cos(theta), y=-48 + r*sin(theta), group=paste0("key_",n), n=as.character(n))
}))
p2 <- ggplot() +
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill="white", colour="grey40", linewidth=0.15) +
  geom_polygon(data=pie_polys, aes(x=x, y=y, group=group, fill=category), colour="white", linewidth=0.15) +
  geom_path(data=size_key, aes(x=x, y=y, group=group), colour="grey30", linewidth=0.3) +
  annotate("text", x=-155, y=-48-sqrt(20000)/35-2, label="20k", size=2.5, colour="grey30") +
  annotate("text", x=-155, y=-48-sqrt(10000)/35-2, label="10k", size=2.5, colour="grey30") +
  annotate("text", x=-155, y=-48-sqrt(5000)/35-2, label="5k", size=2.5, colour="grey30") +
  annotate("text", x=-155, y=-48-sqrt(1000)/35-2, label="1k", size=2.5, colour="grey30") +
  annotate("text", x=-155, y=-48+sqrt(20000)/35+3, label="n PSCs", size=2.8, fontface="bold") +
  scale_fill_manual(values=pie5_cols, name="Replicon type", breaks=c(top5_reps,"Others")) +
  coord_quickmap(xlim=c(-170,180), ylim=c(-55,75), expand=FALSE) +
  labs(title="Geographic distribution of major plasmid replicons",
       subtitle="Pie charts show top 5 replicons + Others; radius proportional to sqrt(n PSCs)") +
  theme_void(base_size=10) +
  theme(plot.title=element_text(face="bold", size=12),
        legend.position.inside=c(0.98,0.35), legend.justification=c(1,0.5),
        legend.background=element_rect(fill=alpha("white",0.7), colour=NA))
cty_top <- cty_tot[total>=200][order(-total)][1:min(25,.N)]
cty_bar <- cty_pie[country %in% cty_top$country]
cty_bar[, country := factor(country, levels=cty_top$country)]
p2b <- ggplot(cty_bar, aes(x=country, y=frac, fill=rep5)) +
  geom_col(position="stack", width=0.8) +
  scale_fill_manual(values=pie5_cols, name="Replicon") +
  scale_y_continuous(labels=percent_format(), expand=c(0,0)) +
  coord_flip() + labs(x=NULL, y="Proportion of PSCs", title="Replicon composition by country (top 25)") +
  theme_pub + theme(legend.position="bottom")
save_plot(p2 / p2b, "fig02_world_map_pies", 12, 10)

# =============================================================================
# FIG 3: ARG count & plasmid length per replicon
# =============================================================================
cat("\n=== Fig 3: ARG burden & plasmid size ===\n")
rep_summary <- d[!is.na(replicon_primary), .(
  n=.N, mean_arg=mean(n_arg,na.rm=TRUE), median_arg=median(n_arg,na.rm=TRUE),
  pct_arg_carrier=100*mean(n_arg>0,na.rm=TRUE), mean_len=mean(length_avg,na.rm=TRUE)/1000,
  median_len=median(length_avg,na.rm=TRUE)/1000, mean_vf=mean(n_vf,na.rm=TRUE),
  pct_conj=100*mean(mobility_class %in% c("conjugative_complete","conjugative_likely"),na.rm=TRUE),
  mean_risk=mean(risk_score,na.rm=TRUE)), by=replicon_primary][order(-n)]
fwrite(rep_summary, file.path(tab_dir,"tab_replicon_summary.csv"))
top30 <- rep_summary$replicon_primary[1:30]
ds <- d[replicon_primary %in% top30]
arg_order <- rep_summary[match(top30,replicon_primary)][order(mean_arg)]$replicon_primary
len_order <- rep_summary[match(top30,replicon_primary)][order(mean_len)]$replicon_primary
p3a <- ggplot(ds, aes(x=factor(replicon_primary, levels=arg_order), y=n_arg, fill=replicon_primary)) +
  geom_violin(scale="width", alpha=0.7, colour=NA) + geom_boxplot(width=0.15, outlier.size=0.3, alpha=0.9) +
  scale_fill_manual(values=palette_rep[top30], guide="none") +
  stat_summary(fun=mean, geom="point", shape=23, size=1.5, fill="white") +
  coord_flip() + scale_y_sqrt(breaks=c(0,1,3,10,30,100)) +
  labs(x=NULL, y="Number of ARGs (sqrt scale)", title="ARG burden per replicon (ordered by mean ARG)") + theme_pub
p3b <- ggplot(ds, aes(x=factor(replicon_primary, levels=len_order), y=length_avg/1000, fill=replicon_primary)) +
  geom_violin(scale="width", alpha=0.7, colour=NA) + geom_boxplot(width=0.15, outlier.size=0.3, alpha=0.9) +
  scale_fill_manual(values=palette_rep[top30], guide="none") +
  stat_summary(fun=mean, geom="point", shape=23, size=1.5, fill="white") +
  coord_flip() + scale_y_log10(breaks=c(1,3,10,30,100,300)) +
  labs(x=NULL, y="Plasmid length kb (log scale)", title="Plasmid size per replicon (ordered by mean length)") + theme_pub
p3c <- ggplot(rep_summary[n>=50], aes(x=mean_len, y=mean_arg, size=n, colour=mean_risk)) +
  geom_point(alpha=0.7) +
  geom_text_repel(data=rep_summary[n>=200][1:20], aes(label=replicon_primary), size=2.8, max.overlaps=25, colour="black") +
  scale_size_continuous(range=c(1,10), labels=comma, name="n PSCs") +
  scale_colour_gradient2(low="steelblue", mid="white", high="firebrick", midpoint=0.3, name="Mean risk\nscore") +
  scale_x_log10(breaks=c(1,3,10,30,100,300)) +
  labs(x="Mean plasmid length (kb, log)", y="Mean ARG count", title="ARG burden vs plasmid size across replicons") + theme_pub
save_plot((p3a|p3b)/p3c + plot_layout(heights=c(2,1)), "fig03_arg_length_per_replicon", 14, 14)

# =============================================================================
# FIG 4: Size evolution
# =============================================================================
cat("\n=== Fig 4: Size evolution ===\n")
size_trends <- rbindlist(lapply(top_reps, function(rp) {
  dat <- d[replicon_primary==rp]
  fit <- lm(length_avg/1000 ~ year_mid, data=dat)
  data.table(replicon_primary=rp, slope_kb_per_year=coef(fit)[2], p=summary(fit)$coefficients[2,4],
             r2=summary(fit)$r.squared, mean_len=mean(dat$length_avg/1000), n=nrow(dat))
}))[order(-slope_kb_per_year)]
size_trends[, sig := p<0.05]
fwrite(size_trends, file.path(tab_dir,"tab_size_trends.csv"))
print(size_trends)
p4a <- ggplot(size_trends, aes(x=reorder(replicon_primary,slope_kb_per_year), y=slope_kb_per_year, fill=sig)) +
  geom_col() + coord_flip() + scale_fill_manual(values=c("grey70","firebrick"), name="p<0.05") +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x=NULL, y="Change in length (kb/year)", title="Plasmid size change over time per replicon") + theme_pub
p4b <- ggplot(d[replicon_primary %in% top12], aes(x=era, y=length_avg/1000, fill=replicon_primary)) +
  geom_boxplot(outlier.size=0.2, linewidth=0.3) + scale_fill_manual(values=palette_rep[top12], name="Replicon") +
  scale_y_log10(breaks=c(1,3,10,30,100,300)) + labs(x="Collection era", y="Plasmid length kb (log)") + theme_pub +
  theme(axis.text.x=element_text(angle=30, hjust=1))
p4c <- ggplot(d[replicon_primary %in% top12], aes(x=year_mid, y=length_avg/1000, colour=replicon_primary)) +
  geom_smooth(method="loess", se=TRUE, linewidth=0.8, alpha=0.1) + scale_colour_manual(values=palette_rep[top12], name="Replicon") +
  scale_y_log10(breaks=c(1,3,10,30,100,300)) + labs(x="Year", y="Plasmid length kb (log)", title="Smoothed size trajectories") + theme_pub
save_plot((p4a|p4b)/p4c + plot_layout(heights=c(2,1)), "fig04_size_evolution", 14, 14)

# =============================================================================
# FIG 5: High-risk ARG heatmap
# =============================================================================
cat("\n=== Fig 5: High-risk ARGs ===\n")
arg_patterns <- list("mcr"="mcr","blaNDM"="NDM","blaKPC"="KPC","blaOXA-48-like"="OXA-48|OXA-181|OXA-232|OXA-244",
  "blaCTX-M"="CTX-M","blaCMY"="CMY","blaSHV"="SHV","blaTEM"="TEM","blaIMP"="IMP","blaVIM"="VIM",
  "tet(X)"="tet\\(X\\)|tetX","van"="van[A-Z]","cfr"="cfr","optrA"="optrA","poxtA"="poxtA",
  "qnr"="qnr","aac(6')-Ib-cr"="aac.6.-Ib-cr","rmt"="rmt[A-H]","fos"="fos[A-Z]")
for (nm in names(arg_patterns)) d[, (nm) := grepl(arg_patterns[[nm]], aro_name, ignore.case=TRUE, perl=TRUE)]
arg_cols <- names(arg_patterns)
arg_rep <- d[!is.na(replicon_primary) & replicon_primary %in% top30,
             c(list(n=.N), lapply(.SD, mean, na.rm=TRUE)), by=replicon_primary, .SDcols=arg_cols]
arg_long <- melt(arg_rep, id.vars=c("replicon_primary","n"), variable.name="ARG", value.name="prevalence")
arg_long[, ARG := factor(ARG, levels=names(arg_patterns))]
arg_rep[, total_hr := rowSums(.SD), .SDcols=arg_cols]
arg_ord <- arg_rep[order(total_hr), replicon_primary]
arg_long[, replicon_primary := factor(replicon_primary, levels=arg_ord)]
p5 <- ggplot(arg_long, aes(x=ARG, y=replicon_primary, fill=prevalence*100)) +
  geom_tile(colour="white", linewidth=0.3) +
  geom_text(aes(label=ifelse(prevalence>=0.02, sprintf("%.0f%%",prevalence*100),"")), size=2.5, colour="white") +
  scale_fill_gradientn(colours=c("grey95","#fff7bc","#fec44f","#d95f0e","#b30000"),
                       values=scales::rescale(c(0,1,5,20,50)), name="% PSCs\ncarrying", na.value="grey95") +
  labs(x="High-risk ARG family", y=NULL, title="Prevalence of high-risk ARGs across major replicons") +
  theme_pub + theme(axis.text.x=element_text(angle=45, hjust=1))
arg_era <- d[, c(list(total=.N), lapply(.SD, sum, na.rm=TRUE)), by=era, .SDcols=arg_cols]
arg_era_long <- melt(arg_era, id.vars=c("era","total"), variable.name="ARG", value.name="count")
arg_era_long[, rate := count/total*100]
p5b <- ggplot(arg_era_long[ARG %in% c("mcr","blaNDM","blaCTX-M","blaKPC","blaOXA-48-like","tet(X)","cfr","van","qnr","optrA")],
              aes(x=era, y=rate, colour=ARG, group=ARG)) +
  geom_line(linewidth=0.9) + geom_point(size=2) + scale_colour_brewer(palette="Set1") +
  labs(x="Collection era", y="% PSCs carrying", title="Temporal trends of high-risk ARG families") +
  theme_pub + theme(axis.text.x=element_text(angle=30,hjust=1))
save_plot(p5/p5b + plot_layout(heights=c(2,1)), "fig05_highrisk_arg", 12, 13)

# =============================================================================
# FIG 6: Virulence vs resistance
# =============================================================================
cat("\n=== Fig 6: Virulence vs resistance ===\n")
d[, plasmid_type := case_when(n_arg>0 & n_vf>0 ~ "Resistance+Virulence (MDR-VF)",
  n_arg>0 & n_vf==0 ~ "Resistance only", n_arg==0 & n_vf>0 ~ "Virulence only", TRUE ~ "Neither")]
d[, plasmid_type := factor(plasmid_type, levels=c("Resistance+Virulence (MDR-VF)","Resistance only","Virulence only","Neither"))]
type_tab <- d[, .N, by=plasmid_type][order(-N)]; type_tab[, pct := N/sum(N)*100]
fwrite(type_tab, file.path(tab_dir,"tab_plasmid_types.csv")); print(type_tab)
rep_vr <- d[!is.na(replicon_primary) & replicon_primary %in% top30, .(
  n=.N, pct_mdr_vf=100*mean(plasmid_type=="Resistance+Virulence (MDR-VF)"),
  pct_res=100*mean(plasmid_type=="Resistance only"), pct_vf=100*mean(plasmid_type=="Virulence only"),
  mean_arg=mean(n_arg), mean_vf=mean(n_vf), mean_len=mean(length_avg)/1000,
  pct_conj=100*mean(mobility_class %in% c("conjugative_complete","conjugative_likely")),
  growth=mean(annual_growth_rate,na.rm=TRUE)), by=replicon_primary][order(-pct_mdr_vf)]
fwrite(rep_vr, file.path(tab_dir,"tab_replicon_virulence_resistance.csv"))
p6a <- ggplot(rep_vr, aes(x=reorder(replicon_primary,pct_mdr_vf), y=pct_mdr_vf, fill=pct_conj)) +
  geom_col() + coord_flip() + scale_fill_gradient(low="grey80", high="firebrick", name="%\nconjugative") +
  labs(x=NULL, y="% PSCs carrying both ARGs and VFs", title="MDR-virulence fusion plasmids per replicon") + theme_pub
p6b <- ggplot(rep_vr, aes(x=mean_arg, y=mean_vf, size=n, colour=pct_conj)) +
  geom_point(alpha=0.7) + geom_text_repel(aes(label=replicon_primary), size=2.5, max.overlaps=20, colour="black") +
  geom_vline(xintercept=mean(rep_vr$mean_arg), linetype="dashed", colour="grey50") +
  geom_hline(yintercept=mean(rep_vr$mean_vf), linetype="dashed", colour="grey50") +
  scale_size_continuous(range=c(1,8), name="n PSCs") + scale_colour_gradient(low="steelblue", high="firebrick", name="%\nconjugative") +
  labs(x="Mean ARG count", y="Mean VF count", title="Virulence-resistance spectrum") + theme_pub
evo_comp <- d[plasmid_type %in% c("Resistance only","Virulence only","Resistance+Virulence (MDR-VF)"),
              .(plasmid_type, year_mid, length_avg, mobility_class, n_country, annual_growth_rate, risk_score, is_density_per_kb)]
kw_year <- kruskal.test(year_mid~plasmid_type, data=evo_comp)
kw_len <- kruskal.test(length_avg~plasmid_type, data=evo_comp)
kw_grow <- kruskal.test(annual_growth_rate~plasmid_type, data=evo_comp)
cat(sprintf("\nKruskal-Wallis: year p=%.2e, length p=%.2e, growth p=%.2e\n", kw_year$p.value, kw_len$p.value, kw_grow$p.value))
p6c <- ggplot(evo_comp, aes(x=plasmid_type, y=year_mid, fill=plasmid_type)) + geom_violin(alpha=0.7) + geom_boxplot(width=0.15, outlier.size=0.3) +
  scale_fill_brewer(palette="Set1", guide="none") + labs(x=NULL, y="Collection year") + theme_pub + theme(axis.text.x=element_text(angle=20,hjust=1))
p6d <- ggplot(evo_comp, aes(x=plasmid_type, y=length_avg/1000, fill=plasmid_type)) + geom_violin(alpha=0.7) + geom_boxplot(width=0.15, outlier.size=0.3) +
  scale_fill_brewer(palette="Set1", guide="none") + scale_y_log10() + labs(x=NULL, y="Length kb (log)") + theme_pub + theme(axis.text.x=element_text(angle=20,hjust=1))
p6e <- ggplot(evo_comp, aes(x=plasmid_type, y=annual_growth_rate, fill=plasmid_type)) + geom_violin(alpha=0.7) + geom_boxplot(width=0.15, outlier.size=0.3) +
  scale_fill_brewer(palette="Set1", guide="none") + coord_cartesian(ylim=c(-0.5,1)) + labs(x=NULL, y="Annual growth rate") + theme_pub + theme(axis.text.x=element_text(angle=20,hjust=1))
save_plot((p6a|p6b)/(p6c|p6d|p6e) + plot_layout(heights=c(1.2,1)), "fig06_virulence_resistance", 14, 14)

# =============================================================================
# FIG 7: Virulence gene composition
# =============================================================================
cat("\n=== Fig 7: Virulence gene composition ===\n")
vf_dt <- d[n_vf>0 & !is.na(vf_name), .(vf=trimws(unlist(strsplit(vf_name,",")))), by=.(plasmid_acc,replicon_primary,plasmid_type)]
vf_dt <- vf_dt[vf!="" & !is.na(vf)]
vf_top <- vf_dt[, .N, by=vf][order(-N)][1:30]
fwrite(vf_top, file.path(tab_dir,"tab_top_virulence_genes.csv"))
p7a <- ggplot(vf_top, aes(x=reorder(vf,N), y=N)) + geom_col(fill="steelblue", alpha=0.8) + coord_flip() +
  scale_y_continuous(labels=comma) + labs(x=NULL, y="Number of PSCs", title="Top 30 virulence genes") + theme_pub
vf_cat <- d[n_vf>0 & !is.na(vf_category), .(cat=trimws(unlist(strsplit(vf_category,",")))),by=plasmid_acc]
vf_cat <- vf_cat[cat!=""]
vf_cat_tab <- vf_cat[, .N, by=cat][order(-N)]; vf_cat_tab[, pct := N/sum(N)*100]
fwrite(vf_cat_tab, file.path(tab_dir,"tab_vf_categories.csv"))
p7b <- ggplot(vf_cat_tab[1:15], aes(x=reorder(cat,N), y=N)) + geom_col(fill="coral", alpha=0.8) + coord_flip() +
  scale_y_continuous(labels=comma) + labs(x=NULL, y="Number of PSCs", title="VF categories (top 15)") + theme_pub
vf_rep <- vf_dt[replicon_primary %in% top_reps, .N, by=.(replicon_primary,vf)]
vf_rep_top <- vf_rep[vf %in% vf_top$vf[1:15]]
p7c <- ggplot(vf_rep_top, aes(x=vf, y=replicon_primary, fill=N)) + geom_tile(colour="white") +
  scale_fill_gradientn(colours=c("grey95","#deebf7","#3182bd","#08519c"), name="Count") +
  labs(x="Virulence gene", y=NULL) + theme_pub + theme(axis.text.x=element_text(angle=45,hjust=1))
save_plot(p7a/p7b/p7c, "fig07_virulence_genes", 11, 14)

# =============================================================================
# FIG 8: Mobility
# =============================================================================
cat("\n=== Fig 8: Mobility ===\n")
mob_rep <- d[!is.na(replicon_primary) & replicon_primary %in% top30, .N, by=.(replicon_primary,mobility_class)]
mob_rep[, total:=sum(N), by=replicon_primary]; mob_rep[, pct:=N/total*100]
mob_levels <- c("conjugative_complete","conjugative_likely","mobilizable","non-mobilizable")
mob_rep[, mobility_class:=factor(mobility_class, levels=mob_levels)]
conj_order <- mob_rep[mobility_class %in% c("conjugative_complete","conjugative_likely"), .(conj_pct=sum(pct)), by=replicon_primary][order(conj_pct)]$replicon_primary
mob_rep[, replicon_primary := factor(replicon_primary, levels=conj_order)]
p8 <- ggplot(mob_rep, aes(x=replicon_primary, y=pct, fill=mobility_class)) + geom_col() + coord_flip() +
  scale_fill_manual(values=c("conjugative_complete"="#b30000","conjugative_likely"="#fc8d59","mobilizable"="#fee08b","non-mobilizable"="#91cf60"), name="Mobility") +
  scale_y_continuous(labels=percent_format(scale=1)) + labs(x=NULL, y="% PSCs", title="Mobility class per replicon") + theme_pub
mob_era <- d[, .N, by=.(era,mobility_class)]; mob_era[, total:=sum(N), by=era]; mob_era[, pct:=N/total*100]
mob_era[, mobility_class:=factor(mobility_class, levels=mob_levels)]
p8b <- ggplot(mob_era, aes(x=era, y=pct, fill=mobility_class)) + geom_area(alpha=0.85, colour="white", linewidth=0.2) +
  scale_fill_manual(values=c("#b30000","#fc8d59","#fee08b","#91cf60"), name="Mobility") +
  scale_y_continuous(labels=percent_format(scale=1)) + labs(x="Collection era", y="% PSCs", title="Mobility trends over time") + theme_pub
save_plot(p8/p8b + plot_layout(heights=c(3,1)), "fig08_mobility", 10, 12)

# =============================================================================
# FIG 9: Host range
# =============================================================================
cat("\n=== Fig 9: Host range ===\n")
sp <- d[!is.na(species_name), .N, by=species_name][order(-N)][1:20]
fwrite(sp, file.path(tab_dir,"tab_top_species.csv"))
p9a <- ggplot(sp, aes(x=reorder(species_name,N), y=N)) + geom_col(fill="darkolivegreen4", alpha=0.8) + coord_flip() +
  scale_y_continuous(labels=comma) + labs(x=NULL, y="Number of PSCs", title="Top 20 host species") + theme_pub
d[, species_clean := sub(",.*","",species_name)]
host_range <- d[!is.na(replicon_primary) & replicon_primary %in% top30,
  .(n_species=uniqueN(species_clean[!is.na(species_clean)]), n_genus=uniqueN(genus_name[!is.na(genus_name)]), n=.N),
  by=replicon_primary][order(-n_species)]
p9b <- ggplot(host_range, aes(x=reorder(replicon_primary,n_species), y=n_species)) + geom_col(fill="darkolivegreen3", alpha=0.8) + coord_flip() +
  labs(x=NULL, y="Number of host species", title="Host range breadth per replicon") + theme_pub
gram <- d[!is.na(gram_stain) & replicon_primary %in% top20, .N, by=.(replicon_primary,gram_stain)]
gram[, total:=sum(N), by=replicon_primary]
p9c <- ggplot(gram, aes(x=reorder(replicon_primary,-total), y=N, fill=gram_stain)) + geom_col(position="fill") + coord_flip() +
  scale_fill_manual(values=c("positive"="#4daf4a","negative"="#377eb8"), name="Gram stain") +
  scale_y_continuous(labels=percent_format()) + labs(x=NULL, y="Proportion", title="Gram-stain distribution") + theme_pub
save_plot(p9a/p9b/p9c, "fig09_host_range", 10, 14)

# =============================================================================
# FIG 10: IS, integrons, metal resistance
# =============================================================================
cat("\n=== Fig 10: IS, integrons, metal resistance ===\n")
mobil <- d[!is.na(replicon_primary) & replicon_primary %in% top30, .(
  n=.N, pct_integron=100*mean(has_integron==TRUE,na.rm=TRUE), mean_is=mean(n_is,na.rm=TRUE),
  is_density=mean(is_density_per_kb,na.rm=TRUE), pct_metal=100*mean(n_metal>0,na.rm=TRUE),
  mean_metal=mean(n_metal,na.rm=TRUE)), by=replicon_primary][order(-pct_integron)]
fwrite(mobil, file.path(tab_dir,"tab_mobilization_potential.csv"))
p10a <- ggplot(mobil, aes(x=reorder(replicon_primary,pct_integron), y=pct_integron)) + geom_col(fill="#756bb1", alpha=0.8) + coord_flip() +
  labs(x=NULL, y="% PSCs with integron") + theme_pub
p10b <- ggplot(mobil, aes(x=reorder(replicon_primary,is_density), y=is_density)) + geom_col(fill="#9e9ac8", alpha=0.8) + coord_flip() +
  labs(x=NULL, y="IS density (per kb)") + theme_pub
p10c <- ggplot(mobil, aes(x=reorder(replicon_primary,pct_metal), y=pct_metal)) + geom_col(fill="#bcbddc", alpha=0.8) + coord_flip() +
  labs(x=NULL, y="% PSCs with metal resistance") + theme_pub
int_era <- d[, .(pct=100*mean(has_integron==TRUE,na.rm=TRUE)), by=era]
p10d <- ggplot(int_era, aes(x=era, y=pct, group=1)) + geom_line(linewidth=1, colour="#756bb1") + geom_point(size=2.5) +
  labs(x="Collection era", y="% PSCs with integron") + theme_pub
save_plot(p10a/p10b/p10c/p10d + plot_layout(heights=c(1,1,1,0.7)), "fig10_is_integron_metal", 10, 18)

# =============================================================================
# FIG 11: Risk scores
# =============================================================================
cat("\n=== Fig 11: Risk scores ===\n")
p11a <- ggplot(d[!is.na(risk_score)], aes(x=risk_score, fill=Q_grade)) +
  geom_histogram(bins=80, alpha=0.9, colour="white", linewidth=0.1) +
  scale_fill_manual(values=c(Q1="#b30000",Q2="#fc8d59",Q3="#fee08b",Q4="#91cf60"), name="Risk grade") +
  labs(x="Composite risk score", y="Number of PSCs", title="Distribution of risk scores") + theme_pub
q_era <- d[!is.na(Q_grade), .N, by=.(era,Q_grade)]; q_era[, total:=sum(N), by=era]; q_era[, pct:=N/total*100]
q_era[, Q_grade:=factor(Q_grade, levels=c("Q1","Q2","Q3","Q4"))]
p11b <- ggplot(q_era, aes(x=era, y=pct, fill=Q_grade)) + geom_col(position="dodge") +
  scale_fill_manual(values=c(Q1="#b30000",Q2="#fc8d59",Q3="#fee08b",Q4="#91cf60"), name="Risk grade") +
  labs(x="Collection era", y="% PSCs", title="Risk grade over time") + theme_pub + theme(axis.text.x=element_text(angle=30,hjust=1))
q1_rep <- d[!is.na(replicon_primary) & replicon_primary %in% top30, .(pct_q1=100*mean(Q_grade=="Q1"), n=.N), by=replicon_primary][order(-pct_q1)]
p11c <- ggplot(q1_rep, aes(x=reorder(replicon_primary,pct_q1), y=pct_q1)) + geom_col(fill="#b30000", alpha=0.8) + coord_flip() +
  labs(x=NULL, y="% Q1 (highest risk)") + theme_pub
save_plot(p11a/p11b/p11c, "fig11_risk_scores", 10, 14)

# =============================================================================
# Summary
# =============================================================================
overall <- data.table(
  metric=c("Total PSCs","With replicon typed","With ARGs","With VFs","With integron",
           "Conjugative (complete)","Conjugative (likely)","Mobilizable","Non-mobilizable",
           "Q1","Q2","Q3","Q4","Median length (kb)","Mean ARG count","Mean VF count",
           "Distinct replicons","Distinct species","Distinct countries","Year range"),
  value=c(format(nrow(d),big.mark=","), format(nrow(d[!is.na(replicon_primary)]),big.mark=","),
    format(nrow(d[n_arg>0]),big.mark=","), format(nrow(d[n_vf>0]),big.mark=","),
    format(nrow(d[has_integron==TRUE]),big.mark=","),
    format(nrow(d[mobility_class=="conjugative_complete"]),big.mark=","),
    format(nrow(d[mobility_class=="conjugative_likely"]),big.mark=","),
    format(nrow(d[mobility_class=="mobilizable"]),big.mark=","),
    format(nrow(d[mobility_class=="non-mobilizable"]),big.mark=","),
    format(nrow(d[Q_grade=="Q1"]),big.mark=","), format(nrow(d[Q_grade=="Q2"]),big.mark=","),
    format(nrow(d[Q_grade=="Q3"]),big.mark=","), format(nrow(d[Q_grade=="Q4"]),big.mark=","),
    round(median(d$length_avg,na.rm=TRUE)/1000,1), round(mean(d$n_arg,na.rm=TRUE),2),
    round(mean(d$n_vf,na.rm=TRUE),2), uniqueN(d$replicon_primary[!is.na(d$replicon_primary)]),
    uniqueN(d$species_clean[!is.na(d$species_clean)]), uniqueN(d$country[!is.na(d$country)&d$country!=""]),
    paste0(min(d$year_mid),"-",max(d$year_mid))))
fwrite(overall, file.path(tab_dir,"tab_overall_summary.csv"))
cat(sprintf("\nDONE. Figures in %s\nTables in %s\n", fig_dir, tab_dir))
