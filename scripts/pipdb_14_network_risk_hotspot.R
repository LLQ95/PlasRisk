#!/usr/bin/env Rscript
# =============================================================================
# pipdb_14_network_risk_hotspot.R
# Priority analyses:
#   Fig 16: ARG-replicon bipartite network
#   Fig 17: ARG-VF co-occurrence network
#   Fig 18: Risk score validation (Q-grade stratification, ROC, sensitivity)
#   Fig 19: Spatial hotspot analysis (Getis-Ord Gi*)
#   Fig 20: MDR-VF fusion plasmid comparative genomics
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(RColorBrewer); library(ggrepel); library(scales)
})
needed <- c("igraph","ggraph","spdep","sp","pROC")
for (pkg in needed) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg, repos="https://cloud.r-project.org", quiet=TRUE)
}
suppressPackageStartupMessages({ library(igraph); library(ggraph); library(pROC) })

args <- commandArgs(trailingOnly=TRUE)
res_dir <- if (length(args)>=1) args[1] else "results"
fig_dir <- file.path(res_dir, "figures", "descriptive")
tab_dir <- file.path(res_dir, "tables")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)

save_plot <- function(plot, name, w=10, h=8) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot, width=w, height=h, limitsize=FALSE)
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot, width=w, height=h, dpi=300, limitsize=FALSE)
  cat("  saved", name, "\n")
}
theme_pub <- theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey90", colour=NA), legend.position="right")

cat("Loading data ...\n")
master_cols <- c("id","plasmid_acc","replicon_primary","length_avg","species_name",
  "genus_name","country","year_min","year_max","year_mid","single_year",
  "aro_name","n_arg","n_who_arg","drugclass","n_drugclass",
  "vf_name","vf_category","n_vf","mobility_class",
  "has_integron","n_is","n_is_family","n_metal","annual_growth_rate",
  "host_human","host_animal","habitat_human","habitat_animal","habitat_env")
risk_cols <- c("id","n_high_risk_arg","risk_score","Q_grade")
d <- fread(file.path(res_dir, "psc_master.tsv"), select=master_cols, na.strings=c("\\N","","-","NA"))
r <- fread(file.path(res_dir, "psc_risk_scores.tsv"), select=risk_cols, na.strings=c("\\N","","-","NA"))
d <- merge(d, r, by="id", all.x=TRUE)

# Override with data-driven final weights if available
final_scores_file <- file.path(tab_dir, "tab_psc_final_scores.csv")
if (file.exists(final_scores_file)) {
  cat("  Loading data-driven final scores ...\n")
  fs <- fread(final_scores_file, select=c("id","S_final"))
  d <- merge(d, fs, by="id", all.x=TRUE)
  d[!is.na(S_final), risk_score := S_final]
  d[, S_final := NULL]
  d[!is.na(risk_score), Q_grade := {
    qs <- quantile(risk_score, probs=c(0.25,0.50,0.75), na.rm=TRUE)
    ifelse(risk_score >= qs[3], "Q1", ifelse(risk_score >= qs[2], "Q2", ifelse(risk_score >= qs[1], "Q3", "Q4")))
  }]
}
cat(sprintf("  %s PSCs loaded\n", format(nrow(d), big.mark=",")))
for (col in c("has_integron","host_human","host_animal","habitat_human","habitat_animal","habitat_env","single_year")) {
  if (col %in% names(d)) d[, (col):=as.logical(get(col))]
}
dt <- d[!is.na(replicon_primary) & replicon_primary!=""]
dt[, era := cut(year_mid, breaks=c(0,1990,2000,2010,2015,2020,2030),
  labels=c("<=1990","1991-2000","2001-2010","2011-2015","2016-2020","2021-2023"))]
rep_tab <- dt[, .N, by=replicon_primary][order(-N)]
top20 <- rep_tab$replicon_primary[1:min(20,nrow(rep_tab))]
top30 <- rep_tab$replicon_primary[1:min(30,nrow(rep_tab))]

# =============================================================================
# FIG 16: ARG-replicon bipartite network
# =============================================================================
cat("\n=== Fig 16: ARG-replicon bipartite network ===\n")
arg_dt <- dt[n_arg>0 & !is.na(aro_name), .(arg=trimws(unlist(strsplit(aro_name, ",")))), by=.(plasmid_acc, replicon_primary)]
arg_dt <- arg_dt[arg!="" & !is.na(arg)]
arg_family_map <- list("TEM"="TEM","SHV"="SHV","CTX-M"="CTX-M","CMY"="CMY","NDM"="NDM","KPC"="KPC",
  "OXA"="OXA","VIM"="VIM","IMP"="IMP","mcr"="mcr","tet(A)"="tet(A)","tet(B)"="tet(B)","tet(M)"="tet(M)",
  "tet(X)"="tet(X)","Qnr"="Qnr","AAC"="AAC/APH","APH"="AAC/APH","ANT"="AAC/APH","Aad"="Aad",
  "Sul"="Sul","Dfr"="Dfr","Erm"="Erm","Mph"="Mph","Qac"="Qac","cfr"="cfr","optrA"="optrA","van"="van")
arg_dt[, arg_family := arg]
for (pat in names(arg_family_map)) arg_dt[grepl(pat, arg, ignore.case=TRUE), arg_family := arg_family_map[[pat]]]
edges_ar <- arg_dt[replicon_primary %in% top20, .(weight=uniqueN(plasmid_acc)), by=.(replicon_primary, arg_family)]
setnames(edges_ar, c("replicon_primary","arg_family"), c("from","to"))
edges_ar_f <- edges_ar[weight>=50]
g_ar <- graph_from_data_frame(edges_ar_f, directed=FALSE)
V(g_ar)$type <- V(g_ar)$name %in% edges_ar_f$to
V(g_ar)$node_type <- ifelse(V(g_ar)$type, "ARG family", "Replicon")
deg <- degree(g_ar)
V(g_ar)$size <- sqrt(deg) * 1.5
V(g_ar)$label <- ifelse(deg >= quantile(deg, 0.6), V(g_ar)$name, "")
E(g_ar)$width <- sqrt(edges_ar_f$weight) / 8
rep_proj <- bipartite_projection(g_ar, which="false")
comm <- cluster_louvain(rep_proj)
mod_val <- modularity(comm)
cat(sprintf("  Replicon projection modularity: %.3f (%d communities)\n", mod_val, length(comm)))
btw <- betweenness(g_ar)
fwrite(edges_ar[order(-weight)], file.path(tab_dir,"tab_arg_replicon_edges.csv"))
fwrite(data.table(node=names(deg), type=V(g_ar)$node_type, degree=as.numeric(deg), betweenness=as.numeric(btw)),
       file.path(tab_dir,"tab_arg_replicon_network_metrics.csv"))
set.seed(42)
p16 <- ggraph(g_ar, layout="fr") +
  geom_edge_link(aes(width=width, alpha=weight), colour="grey60", show.legend=FALSE) +
  geom_node_point(aes(size=size, fill=node_type), shape=21, colour="grey20", stroke=0.3) +
  geom_node_text(aes(label=label), size=2.8, repel=TRUE, max.overlaps=30) +
  scale_fill_manual(values=c("ARG family"="#e41a1c", "Replicon"="#377eb8"), name="Node type") +
  scale_size_continuous(range=c(2,10), guide="none") +
  labs(title="ARG family-replicon bipartite network", subtitle=sprintf("modularity=%.2f; %d nodes, %d edges", mod_val, vcount(g_ar), ecount(g_ar))) +
  theme_void(base_size=11) + theme(plot.title=element_text(face="bold"), legend.position="bottom")
deg_proj <- degree(rep_proj)
V(rep_proj)$size <- sqrt(deg_proj) * 2
p16b <- ggraph(rep_proj, layout="fr") +
  geom_edge_link(alpha=0.3, colour="grey70") +
  geom_node_point(aes(size=size, fill=factor(membership(comm))), shape=21, colour="grey20", stroke=0.3) +
  geom_node_text(aes(label=name), size=3, repel=TRUE, max.overlaps=25) +
  scale_fill_brewer(palette="Set2", name="Community") + scale_size_continuous(range=c(3,12), guide="none") +
  labs(title="Replicon co-occurrence projection") + theme_void(base_size=11) + theme(plot.title=element_text(face="bold"), legend.position="bottom")
save_plot(p16/p16b, "fig16_arg_replicon_network", 14, 14)

# =============================================================================
# FIG 17: ARG-VF co-occurrence network
# =============================================================================
cat("\n=== Fig 17: ARG-VF co-occurrence network ===\n")
vf_dt <- dt[n_vf>0 & !is.na(vf_name), .(vf=trimws(unlist(strsplit(vf_name, ",")))), by=.(plasmid_acc, replicon_primary)]
vf_dt <- vf_dt[vf!="" & !is.na(vf)]
vf_top <- vf_dt[, .N, by=vf][order(-N)][1:40]$vf
arg_top <- arg_dt[, .N, by=arg_family][order(-N)][1:25]$arg_family
pl_arg <- unique(arg_dt[arg_family %in% arg_top, .(plasmid_acc, arg_family)])
pl_vf  <- unique(vf_dt[vf %in% vf_top, .(plasmid_acc, vf)])
pl_both <- merge(pl_arg, pl_vf, by="plasmid_acc", allow.cartesian=TRUE)
cooc <- pl_both[, .N, by=.(arg_family, vf)]; setnames(cooc, "N", "both")
n_arg_pl <- pl_arg[, .(n_arg=uniqueN(plasmid_acc)), by=arg_family]
n_vf_pl  <- pl_vf[, .(n_vf=uniqueN(plasmid_acc)), by=vf]
total_pl <- uniqueN(dt$plasmid_acc)
cooc <- merge(merge(cooc, n_arg_pl, by="arg_family"), n_vf_pl, by="vf")
cat("  Running Fisher's exact tests...\n")
fisher_results <- rbindlist(lapply(seq_len(nrow(cooc)), function(i) {
  a <- cooc$both[i]; b <- cooc$n_arg[i]-a; c <- cooc$n_vf[i]-a; d <- total_pl-a-b-c
  ft <- fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")
  data.table(arg_family=cooc$arg_family[i], vf=cooc$vf[i], both=a, n_arg=cooc$n_arg[i], n_vf=cooc$n_vf[i],
             odds_ratio=unname(ft$estimate), p=ft$p.value)
}))
fisher_results[, p_adj := p.adjust(p, method="BH")]
fisher_results[, enrichment := log2(odds_ratio)]
fisher_results[is.infinite(enrichment), enrichment := 10]
setorder(fisher_results, p_adj)
fwrite(fisher_results, file.path(tab_dir,"tab_arg_vf_cooccurrence.csv"))
sig_pairs <- fisher_results[p_adj<0.001 & both>=20]
cat(sprintf("  %d significant ARG-VF associations (BH p<0.001, n>=20)\n", nrow(sig_pairs)))
if (nrow(sig_pairs) > 0) {
  edges_av <- sig_pairs[, .(from=arg_family, to=vf, weight=both, enrichment=pmax(pmin(enrichment,10),-10))]
  g_av <- graph_from_data_frame(edges_av, directed=FALSE)
  V(g_av)$node_type <- ifelse(V(g_av)$name %in% sig_pairs$arg_family, "ARG family", "Virulence factor")
  V(g_av)$size <- degree(g_av)
  E(g_av)$width <- sqrt(edges_av$weight)/6
  set.seed(123)
  p17 <- ggraph(g_av, layout="fr") +
    geom_edge_link(aes(width=width, alpha=weight), edge_colour="grey60", show.legend=FALSE) +
    geom_node_point(aes(size=size, fill=node_type), shape=21, colour="grey20", stroke=0.3) +
    geom_node_text(aes(label=name), size=2.8, repel=TRUE, max.overlaps=35) +
    scale_fill_manual(values=c("ARG family"="#e41a1c", "Virulence factor"="#4daf4a"), name="Node type") +
    scale_size_continuous(range=c(2,10), guide="none") +
    labs(title="ARG-virulence factor co-occurrence network", subtitle=sprintf("BH p<0.001; %d nodes, %d edges", vcount(g_av), ecount(g_av))) +
    theme_void(base_size=11) + theme(plot.title=element_text(face="bold"), legend.position="bottom")
} else { p17 <- ggplot() + annotate("text", x=0, y=0, label="No significant pairs") + theme_void() }
top_pairs <- fisher_results[order(-enrichment)][1:25]
top_pairs[, pair := paste(arg_family, vf, sep=" - ")]
p17b <- ggplot(top_pairs, aes(x=reorder(pair, enrichment), y=enrichment, fill=log10(both+1))) +
  geom_col() + coord_flip() + scale_fill_gradient(low="grey80", high="#d73027", name="log10(n)") +
  labs(x=NULL, y="log2(Odds Ratio)", title="Top 25 enriched ARG-VF pairs") + theme_pub
save_plot(p17/p17b + plot_layout(heights=c(1.5,1)), "fig17_arg_vf_network", 14, 16)

# =============================================================================
# FIG 18: Risk score validation
# =============================================================================
cat("\n=== Fig 18: Risk score validation ===\n")
dv <- d[!is.na(risk_score) & !is.na(Q_grade)]
q_metrics <- dv[, .(n=.N, mean_arg=mean(n_arg,na.rm=TRUE), mean_highrisk=mean(n_high_risk_arg,na.rm=TRUE),
  pct_arg=100*mean(n_arg>0,na.rm=TRUE), pct_highrisk=100*mean(n_high_risk_arg>0,na.rm=TRUE),
  pct_conj=100*mean(mobility_class %in% c("conjugative_complete","conjugative_likely"),na.rm=TRUE),
  pct_integron=100*mean(has_integron==TRUE,na.rm=TRUE), pct_human=100*mean(host_human==TRUE,na.rm=TRUE),
  mean_vf=mean(n_vf,na.rm=TRUE), mean_growth=mean(annual_growth_rate,na.rm=TRUE), mean_risk=mean(risk_score,na.rm=TRUE)),
  by=Q_grade][order(Q_grade)]
fwrite(q_metrics, file.path(tab_dir,"tab_risk_qgrade_metrics.csv")); print(q_metrics)
q_long <- melt(q_metrics, id.vars="Q_grade", measure.vars=c("mean_arg","mean_highrisk","pct_conj","pct_integron","pct_human","mean_vf","mean_growth"))
p18a <- ggplot(q_long, aes(x=Q_grade, y=value, fill=Q_grade)) + geom_col(alpha=0.85) +
  facet_wrap(~variable, scales="free_y", ncol=4) +
  scale_fill_manual(values=c(Q1="#b30000",Q2="#fc8d59",Q3="#fee08b",Q4="#91cf60"), guide="none") +
  labs(x="Risk quartile", y="Value", title="Risk quartile stratification") + theme_pub
dv[, has_highrisk := n_high_risk_arg > 0]; dv[, has_arg := n_arg > 0]
dv[, has_conj := mobility_class %in% c("conjugative_complete","conjugative_likely")]
roc_hr <- roc(dv$has_highrisk, dv$risk_score, quiet=TRUE)
roc_arg <- roc(dv$has_arg, dv$risk_score, quiet=TRUE)
roc_conj <- roc(dv$has_conj, dv$risk_score, quiet=TRUE)
roc_df <- rbind(
  data.table(FPR=1-roc_hr$specificities, TPR=roc_hr$sensitivities, Outcome=sprintf("High-risk ARG (AUC=%.3f)", as.numeric(auc(roc_hr)))),
  data.table(FPR=1-roc_arg$specificities, TPR=roc_arg$sensitivities, Outcome=sprintf("Any ARG (AUC=%.3f)", as.numeric(auc(roc_arg)))),
  data.table(FPR=1-roc_conj$specificities, TPR=roc_conj$sensitivities, Outcome=sprintf("Conjugative (AUC=%.3f)", as.numeric(auc(roc_conj)))))
p18b <- ggplot(roc_df, aes(x=FPR, y=TPR, colour=Outcome)) + geom_line(linewidth=0.8) +
  geom_abline(intercept=0, slope=1, linetype="dashed", colour="grey50") + scale_colour_brewer(palette="Set1") + coord_equal() +
  labs(x="False positive rate", y="True positive rate", title="ROC validation") + theme_pub +
  theme(legend.position.inside=c(0.75,0.25), legend.background=element_rect(fill=alpha("white",0.8)))
p18d <- ggplot(dv[sample(.N, min(50000,.N))], aes(x=risk_score, y=n_arg+0.5, colour=Q_grade)) +
  geom_point(alpha=0.15, size=0.5) + geom_smooth(method="lm", colour="black", linewidth=0.6) +
  scale_colour_manual(values=c(Q1="#b30000",Q2="#fc8d59",Q3="#fee08b",Q4="#91cf60")) +
  scale_y_log10(breaks=c(0.5,1,2,5,10,30), labels=c(0,1,2,5,10,30)) +
  labs(x="Composite risk score", y="ARG count (log)", title="Risk score vs observed ARG burden") + theme_pub
save_plot((p18a)/(p18b)/(p18d), "fig18_risk_validation", 14, 16)

# =============================================================================
# FIG 19: Spatial hotspot analysis
# =============================================================================
cat("\n=== Fig 19: Spatial hotspot analysis ===\n")
tryCatch({
name_map <- c("U.K."="UK","United States"="USA","USA"="USA","UK"="UK","South Korea"="South Korea","Korea"="South Korea",
  "Czech Republic"="Czech Republic","Czechia"="Czech Republic","Russia"="Russia","Vietnam"="Vietnam","Tanzania"="Tanzania",
  "Myanmar"="Myanmar","Burma"="Myanmar","Laos"="Laos","Iran"="Iran","Taiwan"="Taiwan","Hong Kong"="Hong Kong","Singapore"="Singapore")
dt[, cname := ifelse(country %in% names(name_map), name_map[country], country)]
cty <- dt[!is.na(cname) & cname!="", .(n=.N, n_arg=sum(n_arg>0,na.rm=TRUE), n_highrisk=sum(n_high_risk_arg>0,na.rm=TRUE),
  n_mcr=sum(grepl("mcr",aro_name,ignore.case=TRUE),na.rm=TRUE), n_ndm=sum(grepl("NDM",aro_name),na.rm=TRUE),
  n_ctxm=sum(grepl("CTX-M",aro_name),na.rm=TRUE), n_kpc=sum(grepl("KPC",aro_name),na.rm=TRUE),
  n_vf=sum(n_vf>0,na.rm=TRUE), n_conj=sum(mobility_class %in% c("conjugative_complete","conjugative_likely"),na.rm=TRUE),
  mean_risk=mean(risk_score,na.rm=TRUE)), by=cname]
cty[, `:=`(pct_arg=100*n_arg/n, pct_highrisk=100*n_highrisk/n, pct_mcr=100*n_mcr/n, pct_ndm=100*n_ndm/n,
  pct_ctxm=100*n_ctxm/n, pct_kpc=100*n_kpc/n, pct_vf=100*n_vf/n, pct_conj=100*n_conj/n)]
cty <- cty[n>=50]
fwrite(cty[order(-pct_highrisk)], file.path(tab_dir,"tab_country_hotspots.csv"))
world <- NULL
rdata_path <- "world_maps.RData"
if (!file.exists(rdata_path)) for (p in c("../world_maps.RData","/db/student/metagenome/pipdb/contigs/world_maps.RData"))
  if (file.exists(p)) { rdata_path <- p; break }
if (file.exists(rdata_path)) {
  e <- new.env(); load(rdata_path, envir=e)
  if ("world_GBD" %in% ls(e)) {
    wsf <- e$world_GBD; geoms <- sf::st_geometry(wsf)
    pieces <- vector("list", length(geoms)); k <- 0
    for (i in seq_along(geoms)) {
      g <- geoms[i]; gtype <- as.character(sf::st_geometry_type(g))
      if (gtype=="GEOMETRYCOLLECTION") { g <- tryCatch(sf::st_collection_extract(g,"POLYGON"), error=function(e) NULL); if (is.null(g)||length(g)==0) next }
      gt <- as.character(sf::st_geometry_type(g)); if (!gt %in% c("POLYGON","MULTIPOLYGON")) next
      if (sf::st_is_empty(g)) next; coords <- sf::st_coordinates(g); if (nrow(coords)==0) next
      k <- k+1; l2 <- if ("L2" %in% colnames(coords)) coords[,"L2"] else 0
      pieces[[k]] <- data.frame(long=as.numeric(coords[,"X"]), lat=as.numeric(coords[,"Y"]),
        group=paste(i, coords[,"L1"], l2, sep="_"), region=as.character(wsf[[2]][i]), stringsAsFactors=FALSE)
    }
    world <- do.call(rbind, pieces[seq_len(k)])
  }
}
if (is.null(world)) world <- ggplot2::map_data("world")
map_data <- merge(world, cty, by.x="region", by.y="cname", all.x=TRUE)
map_data <- map_data[order(map_data$group, map_data$long), ]
if (requireNamespace("spdep", quietly=TRUE)) {
  cat("  Computing Getis-Ord Gi* hotspots...\n")
  cty_centroids <- aggregate(cbind(long,lat) ~ region, data=world, FUN=mean)
  cty_sp <- merge(cty_centroids, cty, by.x="region", by.y="cname"); setDT(cty_sp)
  cty_sp <- cty_sp[!is.na(pct_highrisk)]
  if (nrow(cty_sp) > 10) {
    coords <- as.matrix(cty_sp[, c("long","lat")])
    nb <- spdep::knearneigh(coords, k=5); nb <- spdep::knn2nb(nb)
    lw <- spdep::nb2listw(nb, style="W", zero.policy=TRUE)
    gi <- spdep::localG(cty_sp$pct_highrisk, lw, zero.policy=TRUE)
    cty_sp$gi <- as.numeric(gi); cty_sp$gi_p <- 2*pnorm(-abs(cty_sp$gi))
    cty_sp$hotspot <- ifelse(cty_sp$gi>0 & cty_sp$gi_p<0.05, "Hotspot (95%)", ifelse(cty_sp$gi<0 & cty_sp$gi_p<0.05, "Coldspot (95%)", "Not significant"))
    fwrite(cty_sp[, .(region, n, pct_highrisk, gi, gi_p, hotspot)], file.path(tab_dir,"tab_getis_ord_hotspots.csv"))
  }
} else cty_sp <- NULL
p19a <- ggplot() +
  geom_polygon(data=map_data, aes(x=long, y=lat, group=group, fill=pct_highrisk), colour="white", linewidth=0.1) +
  scale_fill_gradientn(colours=c("#ffffcc","#fed976","#fd8d3c","#e31a1c","#800026"), na.value="grey90", name="% high-risk ARGs") +
  coord_quickmap(xlim=c(-170,180), ylim=c(-55,75), expand=FALSE) +
  labs(title="Global prevalence of high-risk ARG carriage") + theme_void(base_size=10) + theme(plot.title=element_text(face="bold"))
if (!is.null(cty_sp) && "hotspot" %in% names(cty_sp)) {
  p19b <- ggplot() +
    geom_polygon(data=world, aes(x=long, y=lat, group=group), fill="grey92", colour="white", linewidth=0.1) +
    geom_point(data=cty_sp, aes(x=long, y=lat, size=n, colour=hotspot), alpha=0.8) +
    scale_colour_manual(values=c("Hotspot (95%)"="#d73027","Coldspot (95%)"="#4575b4","Not significant"="grey60"), name="Gi* hotspot") +
    scale_size_continuous(range=c(1,6), name="n PSCs", labels=comma) +
    coord_quickmap(xlim=c(-170,180), ylim=c(-55,75), expand=FALSE) +
    labs(title="Getis-Ord Gi* spatial hotspots") + theme_void(base_size=10) + theme(plot.title=element_text(face="bold"))
} else { p19b <- ggplot() + geom_polygon(data=world, aes(x=long,y=lat,group=group), fill="grey92", colour="white") + theme_void() }
top_cty <- cty[order(-pct_highrisk)][1:min(20,.N)]
p19c <- ggplot(top_cty, aes(x=reorder(cname,pct_highrisk), y=pct_highrisk, fill=mean_risk)) + geom_col() + coord_flip() +
  scale_fill_gradient(low="#fee08b", high="#b30000", name="Mean risk") + labs(x=NULL, y="% high-risk ARGs", title="Top 20 countries") + theme_pub
save_plot((p19a/p19b)|p19c + plot_layout(widths=c(1.2,1)), "fig19_spatial_hotspots", 16, 12)
}, error=function(e) cat("  [SKIPPED] Fig 19 failed:", conditionMessage(e), "\n"))

# =============================================================================
# FIG 20: MDR-VF fusion plasmids
# =============================================================================
cat("\n=== Fig 20: MDR-VF fusion plasmids ===\n")
dt[, plasmid_type := fifelse(n_arg>0 & n_vf>0, "MDR-VF", fifelse(n_arg>0 & n_vf==0, "Resistance only", fifelse(n_arg==0 & n_vf>0, "Virulence only", "Neither")))]
fusion <- dt[plasmid_type=="MDR-VF"]
cat(sprintf("  %s MDR-VF fusion plasmids\n", format(nrow(fusion), big.mark=",")))
fus_rep <- fusion[, .N, by=replicon_primary][order(-N)][1:20]
all_rep_pct <- dt[, .(total=.N, n_fus=sum(plasmid_type=="MDR-VF")), by=replicon_primary]
all_rep_pct[, fus_rate := 100*n_fus/total]
fus_rep <- merge(fus_rep, all_rep_pct[, .(replicon_primary, fus_rate)], by="replicon_primary")
p20a <- ggplot(fus_rep, aes(x=reorder(replicon_primary,N), y=N, fill=fus_rate)) + geom_col() + coord_flip() +
  scale_fill_gradient(low="#fc8d59", high="#b30000", name="% fusion") + labs(x=NULL, y="Count", title="MDR-VF by replicon") + theme_pub
fus_arg <- fusion[n_arg>0 & !is.na(aro_name), .(arg=trimws(unlist(strsplit(aro_name,",")))), by=.(plasmid_acc, replicon_primary)]
fus_arg <- fus_arg[arg!="" & !is.na(arg)]
fus_arg_top <- fus_arg[, .N, by=arg][order(-N)][1:20]
p20b <- ggplot(fus_arg_top, aes(x=reorder(arg,N), y=N)) + geom_col(fill="#e41a1c", alpha=0.8) + coord_flip() +
  labs(x=NULL, y="Count", title="Top ARGs on MDR-VF plasmids") + theme_pub
fus_vf <- fusion[n_vf>0 & !is.na(vf_name), .(vf=trimws(unlist(strsplit(vf_name,",")))), by=replicon_primary]
fus_vf <- fus_vf[vf!="" & !is.na(vf)]
fus_vf_top <- fus_vf[, .N, by=vf][order(-N)][1:20]
p20c <- ggplot(fus_vf_top, aes(x=reorder(vf,N), y=N)) + geom_col(fill="#4daf4a", alpha=0.8) + coord_flip() +
  labs(x=NULL, y="Count", title="Top VFs on MDR-VF plasmids") + theme_pub
mob_comp <- dt[plasmid_type %in% c("MDR-VF","Resistance only","Virulence only"), .N, by=.(plasmid_type, mobility_class)]
mob_comp[, total:=sum(N), by=plasmid_type]; mob_comp[, pct:=100*N/total]
mob_comp[, mobility_class:=factor(mobility_class, levels=c("conjugative_complete","conjugative_likely","mobilizable","non-mobilizable"))]
p20e <- ggplot(mob_comp, aes(x=plasmid_type, y=pct, fill=mobility_class)) + geom_col() +
  scale_fill_manual(values=c("conjugative_complete"="#b30000","conjugative_likely"="#fc8d59","mobilizable"="#fee08b","non-mobilizable"="#91cf60"), name="Mobility") +
  labs(x=NULL, y="% of plasmids", title="Mobility by plasmid type") + theme_pub
fus_era <- dt[, .(total=.N, fus=sum(plasmid_type=="MDR-VF")), by=era]
fus_era[, pct_fus := 100*fus/total]
p20g <- ggplot(fus_era[!is.na(era)], aes(x=era, y=pct_fus, group=1)) +
  geom_line(linewidth=1, colour="#b30000") + geom_point(size=3, colour="#b30000") +
  labs(x="Collection era", y="% MDR-VF", title="Temporal trend of MDR-VF proportion") + theme_pub
fwrite(fusion[, .(n=.N, mean_len=mean(length_avg/1000,na.rm=TRUE), mean_arg=mean(n_arg,na.rm=TRUE),
  mean_vf=mean(n_vf,na.rm=TRUE), pct_conj=100*mean(mobility_class %in% c("conjugative_complete","conjugative_likely"),na.rm=TRUE),
  pct_integron=100*mean(has_integron==TRUE,na.rm=TRUE), pct_human=100*mean(host_human==TRUE,na.rm=TRUE),
  n_species=uniqueN(species_name), n_countries=uniqueN(cname), n_replicons=uniqueN(replicon_primary))],
  file.path(tab_dir,"tab_fusion_plasmid_summary.csv"))
save_plot((p20a|p20b|p20c)/(p20e|p20g), "fig20_fusion_plasmids", 18, 14)

cat("\n=== All analyses complete ===\n")
