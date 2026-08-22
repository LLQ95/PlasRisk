#!/usr/bin/env Rscript
# =============================================================================
# pipdb_16_risk_10dimensions.R
# 10-dimension risk assessment: adds S_BM (biocides/metals) to the 9-dim model.
# Compares 9-dim vs 10-dim rankings, RF importance, sensitivity, visualization.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(RColorBrewer); library(scales)
})
for (pkg in c("randomForest","pROC")) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg, repos="https://cloud.r-project.org", quiet=TRUE)
}
suppressPackageStartupMessages({ library(randomForest); library(pROC) })

args <- commandArgs(trailingOnly=TRUE)
res_dir <- if (length(args)>=1) args[1] else "results"
fig_dir <- file.path(res_dir, "figures", "descriptive")
tab_dir <- file.path(res_dir, "tables")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE); dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)
save_plot <- function(plot, name, w=10, h=8) {
  ggsave(file.path(fig_dir, paste0(name,".pdf")), plot, width=w, height=h, limitsize=FALSE)
  ggsave(file.path(fig_dir, paste0(name,".png")), plot, width=w, height=h, dpi=300, limitsize=FALSE)
  cat("  saved", name, "\n")
}
theme_pub <- theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey90", colour=NA))

# 9-dim weights (original, normalized)
W9 <- c(S_ARG=0.30, S_VF=0.30, S_MOB=0.20, S_HOST=0.15, S_REP=0.10, S_SIZE=0.10, S_GEO=0.05, S_HAB=0.05, S_GROW=0.05)
W9 <- W9 / sum(W9)
# 10-dim weights: add S_BM=0.10, redistribute proportionally
W10 <- c(W9 * 0.9, S_BM=0.10)
cat(sprintf("9-dim weight sum: %.1f\n", sum(W9))); cat(sprintf("10-dim weight sum: %.1f\n", sum(W10)))

cat("Loading data ...\n")
master_cols <- c("id","plasmid_acc","replicon_primary","length_avg","species_name","country","n_country",
  "year_min","year_max","year_mid","aro_name","n_arg","n_who_arg","vf_name","vf_category","n_vf",
  "mobility_class","has_t4cp","has_relaxase","has_oriT","has_auxiliary","has_integron","n_is","n_metal",
  "gene_bacmet","annual_growth_rate","host_human","host_animal","habitat_human","habitat_animal","habitat_env")
risk_cols <- c("id","S_ARG","S_MOB","S_HOST","S_REP","S_SIZE","S_GEO","S_HAB","S_GROW","n_high_risk_arg","risk_score","Q_grade")
d <- fread(file.path(res_dir,"psc_master.tsv"), select=master_cols, na.strings=c("\\N","","-","NA"))
r <- fread(file.path(res_dir,"psc_risk_scores.tsv"), select=risk_cols, na.strings=c("\\N","","-","NA"))
d <- merge(d, r, by="id", all.x=TRUE)
cat(sprintf("  %s PSCs loaded\n", format(nrow(d), big.mark=",")))
for (col in c("has_integron","has_t4cp","has_relaxase","has_oriT","has_auxiliary","host_human","host_animal","habitat_human","habitat_animal","habitat_env"))
  if (col %in% names(d)) d[, (col):=as.logical(get(col))]

# S_VF
cat("Computing S_VF ...\n")
d[, vf_base := 0.0]; d[n_vf > 0 & !is.na(n_vf), vf_base := 0.30]
d[, vf_count_bonus := pmin(n_vf * 0.03, 0.40)]
d[, has_exotoxin := grepl("Exotoxin", vf_category, ignore.case=TRUE)]
d[, has_effector := grepl("Effector delivery", vf_category, ignore.case=TRUE)]
d[is.na(has_exotoxin), has_exotoxin := FALSE]; d[is.na(has_effector), has_effector := FALSE]
d[, vf_category_bonus := 0.15*has_exotoxin + 0.15*has_effector]
d[, S_VF := pmin(vf_base + vf_count_bonus + vf_category_bonus, 1.0)]
d[n_vf==0 | is.na(n_vf), S_VF := 0.0]

# S_BM
cat("Computing S_BM ...\n")
d[, n_metal := fifelse(is.na(n_metal), 0L, n_metal)]
d[, has_mer := as.integer(grepl("mer|mercury", gene_bacmet, ignore.case=TRUE))]
d[, has_qac := as.integer(grepl("qac|quaternary|disinfectant", gene_bacmet, ignore.case=TRUE))]
d[, has_ars := as.integer(grepl("ars|cop|sil|czc|cad|pco", gene_bacmet, ignore.case=TRUE))]
d[is.na(has_mer), has_mer := 0L]; d[is.na(has_qac), has_qac := 0L]; d[is.na(has_ars), has_ars := 0L]
d[, S_BM := 0.0]
d[n_metal > 0, S_BM := pmin(0.25 + pmin(n_metal*0.04, 0.35) + 0.15*has_mer + 0.15*has_qac + 0.10*has_ars, 1.0)]
cat(sprintf("  S_BM: mean=%.3f, %d PSCs with S_BM>0 (%.1f%%)\n", mean(d$S_BM,na.rm=TRUE), sum(d$S_BM>0,na.rm=TRUE), 100*mean(d$S_BM>0,na.rm=TRUE)))
cat(sprintf("  Mercury: %.1f%%, Qac: %.1f%%, As/Cu/Ag: %.1f%%\n", 100*mean(d$has_mer,na.rm=TRUE), 100*mean(d$has_qac,na.rm=TRUE), 100*mean(d$has_ars,na.rm=TRUE)))

# 9-dim and 10-dim composite scores
comp9 <- c("S_ARG","S_VF","S_MOB","S_HOST","S_REP","S_SIZE","S_GEO","S_HAB","S_GROW")
comp10 <- c(comp9, "S_BM")
for (col in comp10) d[is.na(get(col)), (col):=0]
d[, S9 := as.matrix(.SD) %*% W9[comp9], .SDcols=comp9]
d[, S10 := as.matrix(.SD) %*% W10[comp10], .SDcols=comp10]

# Per-replicon tables
cat("Building per-replicon tables ...\n")
dt <- d[!is.na(replicon_primary) & replicon_primary!=""]
rep9 <- dt[, .(S9=mean(S9,na.rm=TRUE)), by=replicon_primary][order(-S9)]; rep9[, rank9 := .I]
rep10 <- dt[, .(S10=mean(S10,na.rm=TRUE), n_PSC=.N,
  mean_S_ARG=mean(S_ARG,na.rm=TRUE), mean_S_VF=mean(S_VF,na.rm=TRUE), mean_S_BM=mean(S_BM,na.rm=TRUE),
  mean_S_MOB=mean(S_MOB,na.rm=TRUE), mean_S_HOST=mean(S_HOST,na.rm=TRUE), mean_S_REP=mean(S_REP,na.rm=TRUE),
  mean_S_SIZE=mean(S_SIZE,na.rm=TRUE), mean_S_GEO=mean(S_GEO,na.rm=TRUE), mean_S_HAB=mean(S_HAB,na.rm=TRUE), mean_S_GROW=mean(S_GROW,na.rm=TRUE)),
  by=replicon_primary][order(-S10)]; rep10[, rank10 := .I]
rep_compare <- merge(rep9, rep10, by="replicon_primary")
rep_compare[, rank_change := rank9 - rank10]
setorder(rep_compare, rank10)
num_cols <- names(rep_compare)[sapply(rep_compare, is.numeric)]
rep_compare[, (num_cols) := lapply(.SD, function(x) round(x,4)), .SDcols=num_cols]
fwrite(rep_compare, file.path(tab_dir,"tab_replicon_risk_10dimensions.csv"))
cat(sprintf("  Table saved: %d replicons (%d major)\n", nrow(rep_compare), sum(rep_compare$n_PSC>=50)))

# S_BM analysis
cat("\n=== S_BM analysis ===\n")
bm_by_rep <- dt[, .(n_PSC=.N, mean_S_BM=mean(S_BM,na.rm=TRUE), pct_BM=100*mean(S_BM>0,na.rm=TRUE),
  pct_mer=100*mean(has_mer,na.rm=TRUE), pct_qac=100*mean(has_qac,na.rm=TRUE), pct_ars=100*mean(has_ars,na.rm=TRUE)),
  by=replicon_primary][n_PSC>=50][order(-mean_S_BM)]
cat("Top 10 replicons by S_BM:\n"); print(head(bm_by_rep, 10))

# S_BM temporal trend
d[, era := cut(year_mid, breaks=c(0,1970,1990,2000,2010,2015,2020,2030),
  labels=c("<=1970","1971-1990","1991-2000","2001-2010","2011-2015","2016-2020",">2020"), include.lowest=TRUE)]
bm_temporal <- d[!is.na(era), .(n_PSC=.N, mean_S_BM=mean(S_BM,na.rm=TRUE), pct_BM=100*mean(S_BM>0,na.rm=TRUE)), by=era]
p25a <- ggplot(bm_temporal, aes(x=era, y=mean_S_BM, group=1)) +
  geom_line(linewidth=0.8, colour="#756bb1") + geom_point(size=3, colour="#756bb1") +
  labs(x="Isolation era", y="Mean S_BM", title="Temporal trend of biocide/metal resistance") + theme_pub
p25b <- ggplot(bm_temporal, aes(x=era, y=pct_BM, group=1)) +
  geom_col(fill="#bcbddc", colour="#756bb1", linewidth=0.3) +
  labs(x="Isolation era", y="% PSCs with S_BM > 0", title="Prevalence of BM resistance genes") + theme_pub
p25c <- ggplot(d[S_BM>0], aes(x=S_BM)) +
  geom_histogram(bins=50, fill="#756bb1", colour="white", linewidth=0.1) +
  facet_wrap(~replicon_primary[replicon_primary %in% head(bm_by_rep$replicon_primary,6)]) +
  labs(x="S_BM", y="Count", title="S_BM distribution by top replicons") + theme_pub
p25d <- ggplot(bm_by_rep[1:20], aes(x=reorder(replicon_primary, mean_S_BM), y=mean_S_BM)) +
  geom_col(fill="#756bb1", colour="white") + coord_flip() +
  labs(x=NULL, y="Mean S_BM", title="Top 20 replicons by S_BM") + theme_pub
save_plot(p25a/p25b, "fig25_SBM_distribution", 10, 8)

# 9-dim vs 10-dim comparison
cat("\n=== 9-dim vs 10-dim comparison ===\n")
major <- rep_compare[n_PSC>=50]
rho_all <- cor(rep_compare$rank9, rep_compare$rank10, method="spearman")
rho_major <- cor(major$rank9, major$rank10, method="spearman")
cat(sprintf("  Rank correlation (all): rho=%.4f\n", rho_all))
cat(sprintf("  Rank correlation (major): rho=%.4f\n", rho_major))
p26a <- ggplot(rep_compare, aes(x=rank9, y=rank10, size=n_PSC)) +
  geom_point(alpha=0.4, colour="grey50") +
  geom_point(data=major, aes(x=rank9, y=rank10, size=n_PSC), colour="#756bb1", alpha=0.7) +
  geom_abline(intercept=0, slope=1, linetype="dashed", colour="red") +
  scale_size_continuous(range=c(1,8), name="n PSC") +
  labs(x="Rank (9-dim)", y="Rank (10-dim)", title="Rank comparison: 9-dim vs 10-dim") + theme_pub
p26b <- ggplot(major, aes(x=rank_change)) +
  geom_histogram(binwidth=1, fill="#756bb1", colour="white") +
  geom_vline(xintercept=0, linetype="dashed", colour="red") +
  labs(x="Rank change (9-dim -> 10-dim)", y="Number of replicons", title="Rank shift distribution") + theme_pub
save_plot(p26a|p26b, "fig26_rank_comparison", 12, 6)

# RF importance for 10 dimensions
cat("\n=== RF importance (10-dim) ===\n")
ml_dt <- dt[!is.na(S_ARG), .(S_ARG, S_VF, S_MOB, S_HOST, S_REP, S_SIZE, S_BM, S_GEO, S_HAB, S_GROW,
  n_high_risk_arg, n_arg, n_vf, mobility_class, n_metal,
  highrisk=factor(ifelse(n_high_risk_arg>0,"Yes","No")), mdr_vf=factor(ifelse(n_arg>0 & n_vf>0,"Yes","No")))]
ml_dt <- ml_dt[complete.cases(ml_dt)]
set.seed(42); ml_sample <- ml_dt[sample(.N, min(100000,.N))]
ml_sample[, conj := factor(ifelse(mobility_class %in% c("conjugative_complete","conjugative_likely"),"Yes","No"))]
ml_sample[, bm := factor(ifelse(n_metal>0,"Yes","No"))]
rf1 <- randomForest(highrisk ~ S_ARG+S_VF+S_MOB+S_HOST+S_REP+S_SIZE+S_BM+S_GEO+S_HAB+S_GROW, data=ml_sample, ntree=300, importance=TRUE)
rf2 <- randomForest(mdr_vf ~ S_ARG+S_VF+S_MOB+S_HOST+S_REP+S_SIZE+S_BM+S_GEO+S_HAB+S_GROW, data=ml_sample, ntree=300, importance=TRUE)
rf3 <- randomForest(conj ~ S_ARG+S_VF+S_MOB+S_HOST+S_REP+S_SIZE+S_BM+S_GEO+S_HAB+S_GROW, data=ml_sample, ntree=300, importance=TRUE)
rf4 <- randomForest(bm ~ S_ARG+S_VF+S_MOB+S_HOST+S_REP+S_SIZE+S_BM+S_GEO+S_HAB+S_GROW, data=ml_sample, ntree=300, importance=TRUE)
cat(sprintf("  AUC highrisk=%.3f, mdrvf=%.3f, conj=%.3f, bm=%.3f\n",
  as.numeric(auc(roc(ml_sample$highrisk, predict(rf1,type="prob")[,2], quiet=TRUE))),
  as.numeric(auc(roc(ml_sample$mdr_vf, predict(rf2,type="prob")[,2], quiet=TRUE))),
  as.numeric(auc(roc(ml_sample$conj, predict(rf3,type="prob")[,2], quiet=TRUE))),
  as.numeric(auc(roc(ml_sample$bm, predict(rf4,type="prob")[,2], quiet=TRUE)))))
imp1 <- as.data.table(importance(rf1), keep.rownames="feature"); setnames(imp1,"MeanDecreaseGini","MDG_highrisk")
imp2 <- as.data.table(importance(rf2), keep.rownames="feature"); setnames(imp2,"MeanDecreaseGini","MDG_mdrvf")
imp3 <- as.data.table(importance(rf3), keep.rownames="feature"); setnames(imp3,"MeanDecreaseGini","MDG_conj")
imp4 <- as.data.table(importance(rf4), keep.rownames="feature"); setnames(imp4,"MeanDecreaseGini","MDG_bm")
imp <- Reduce(function(x,y) merge(x,y,by="feature"), list(imp1[,.(feature,MDG_highrisk)], imp2[,.(feature,MDG_mdrvf)], imp3[,.(feature,MDG_conj)], imp4[,.(feature,MDG_bm)]))
imp[, `:=`(MDG_highrisk_norm=MDG_highrisk/sum(MDG_highrisk), MDG_mdrvf_norm=MDG_mdrvf/sum(MDG_mdrvf), MDG_conj_norm=MDG_conj/sum(MDG_conj), MDG_bm_norm=MDG_bm/sum(MDG_bm))]
imp[, mean_MDG_norm := (MDG_highrisk_norm+MDG_mdrvf_norm+MDG_conj_norm+MDG_bm_norm)/4]
imp[, weight_10dim := W10[feature]]; imp[, weight_10dim_norm := weight_10dim/sum(weight_10dim)]
setorder(imp, -mean_MDG_norm)
fwrite(imp, file.path(tab_dir,"tab_rf_importance_10dim.csv"))
imp_long <- melt(imp, id.vars="feature", measure.vars=c("weight_10dim_norm","mean_MDG_norm"), variable.name="source", value.name="importance")
imp_long[, source := factor(source, levels=c("weight_10dim_norm","mean_MDG_norm"), labels=c("Expert weight","RF MDG"))]
imp_long[, feature := factor(feature, levels=imp$feature)]
p27 <- ggplot(imp_long, aes(x=feature, y=importance, fill=source)) +
  geom_col(position="dodge", colour="white", linewidth=0.3) + coord_flip() +
  scale_fill_manual(values=c("#756bb1","#fc8d59"), name=NULL) +
  labs(x=NULL, y="Normalized importance", title="Expert weights vs RF feature importance (10-dim)") + theme_pub
save_plot(p27, "fig27_rf_importance_10dim", 10, 7)

# S_BM weight sensitivity
cat("\n=== S_BM weight sensitivity ===\n")
set.seed(123); bm_weights <- seq(0, 0.30, by=0.02); sens_results <- list()
for (bw in bm_weights) {
  w_adj <- c(W9 * (1-bw), S_BM=bw)
  s_adj <- as.matrix(dt[, .SD, .SDcols=comp10]) %*% w_adj
  tmp <- data.table(replicon_primary=dt$replicon_primary, s=as.numeric(s_adj))[, .(s=mean(s)), by=replicon_primary]
  tmp <- merge(tmp, dt[, .N, by=replicon_primary], by="replicon_primary")[N>=50][order(-s)]
  tmp[, rank_adj := .I]
  tmp <- merge(tmp, rep10[, .(replicon_primary, rank10)], by="replicon_primary")
  rho <- cor(tmp$rank10, tmp$rank_adj, method="spearman")
  top10 <- length(intersect(rep10[rank10<=10]$replicon_primary, tmp[rank_adj<=10]$replicon_primary))
  sens_results[[length(sens_results)+1]] <- data.table(S_BM_weight=bw, spearman_rho=rho, top10_overlap=top10)
}
sens_dt <- rbindlist(sens_results)
fwrite(sens_dt, file.path(tab_dir,"tab_SBM_weight_sensitivity.csv"))
p28 <- ggplot(sens_dt, aes(x=S_BM_weight, y=spearman_rho)) +
  geom_line(linewidth=0.8, colour="#756bb1") + geom_point(size=2, colour="#756bb1") +
  geom_vline(xintercept=0.10, linetype="dashed", colour="red") +
  annotate("text", x=0.12, y=min(sens_dt$spearman_rho), label="chosen: 0.10", hjust=0, colour="red", size=3) +
  labs(x="S_BM weight", y="Spearman rho vs 10-dim ranking", title="S_BM weight sensitivity") + theme_pub
save_plot(p28, "fig28_SBM_sensitivity", 8, 5)

# 10-dim heatmap for top replicons
top30 <- rep10[n_PSC>=50]$replicon_primary[1:min(30, sum(rep10$n_PSC>=50))]
heat_dt <- melt(rep10[replicon_primary %in% top30], id.vars="replicon_primary",
  measure.vars=paste0("mean_", comp10), variable.name="component", value.name="score")
heat_dt[, component := gsub("mean_", "", component)]
heat_dt[, replicon_primary := factor(replicon_primary, levels=rev(top30))]
heat_dt[, component := factor(component, levels=comp10)]
p29 <- ggplot(heat_dt, aes(x=component, y=replicon_primary, fill=score)) +
  geom_tile(colour="white", linewidth=0.3) +
  scale_fill_gradientn(colours=c("#f7fcf5","#c7e9c0","#74c476","#238b45","#00441b"), limits=c(0,1), name="Score") +
  scale_x_discrete(position="top") + labs(x=NULL, y=NULL, title="10-dimension risk component scores (top 30 replicons)") +
  theme_pub + theme(axis.text.x=element_text(angle=45, hjust=0, face="bold"), axis.text.y=element_text(size=8))
save_plot(p29, "fig29_10dim_heatmap", 12, 10)
cat("\n=== All done ===\n")
