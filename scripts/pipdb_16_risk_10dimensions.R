#!/usr/bin/env Rscript
# =============================================================================
# pipdb_16_risk_10dimensions.R
# Evaluate adding S_BM (biocide & metal resistance) as a 10th risk dimension.
# Compares 9-dim expert vs 10-dim data-driven weights, RF validation, sensitivity.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(RColorBrewer); library(scales); library(ggrepel)
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

W9 <- c(S_ARG=0.30, S_VF=0.30, S_MOB=0.20, S_HOST=0.15, S_REP=0.10, S_SIZE=0.10, S_GEO=0.05, S_HAB=0.05, S_GROW=0.05)
wf <- file.path(tab_dir, "tab_weight_comparison.csv")
if (file.exists(wf)) {
  wc <- fread(wf); W10 <- wc$Final; names(W10) <- wc$Component
  cat("Loaded data-driven final weights\n")
} else {
  W10 <- c(S_ARG=0.2448, S_VF=0.1096, S_MOB=0.2041, S_HOST=0.0282, S_REP=0.0030,
           S_SIZE=0.1808, S_BM=0.2112, S_GEO=0.0015, S_HAB=0.0022, S_GROW=0.0147)
}
W9_dd <- W10[names(W10) != "S_BM"]; W9_dd <- W9_dd / sum(W9_dd)

cat("Loading data ...\n")
master_cols <- c("id","plasmid_acc","replicon_primary","length_avg","species_name","country","year_min","year_max","year_mid",
  "aro_name","n_arg","n_who_arg","vf_name","vf_category","n_vf","mobility_class","has_integron","n_is",
  "gene_bacmet","n_metal","annual_growth_rate","host_human","host_animal")
risk_cols <- c("id","S_ARG","S_MOB","S_HOST","S_REP","S_SIZE","S_GEO","S_HAB","S_GROW","n_high_risk_arg","risk_score","Q_grade")
d <- fread(file.path(res_dir,"psc_master.tsv"), select=master_cols, na.strings=c("\\N","","-","NA"))
r <- fread(file.path(res_dir,"psc_risk_scores.tsv"), select=risk_cols, na.strings=c("\\N","","-","NA"))
d <- merge(d, r, by="id", all.x=TRUE)
cat(sprintf("  %s PSCs loaded\n", format(nrow(d), big.mark=",")))
for (col in c("has_integron","host_human","host_animal")) if (col %in% names(d)) d[, (col):=as.logical(get(col))]

# S_VF
cat("Computing S_VF ...\n")
d[, vf_base := 0.0]; d[n_vf > 0 & !is.na(n_vf), vf_base := 0.30]
d[, vf_count_bonus := pmin(n_vf * 0.03, 0.40)]
d[, has_exotoxin := grepl("Exotoxin", vf_category, ignore.case=TRUE)]
d[, has_effector := grepl("Effector delivery", vf_category, ignore.case=TRUE)]
d[is.na(has_exotoxin), has_exotoxin := FALSE]; d[is.na(has_effector), has_effector := FALSE]
d[, S_VF := pmin(vf_base + vf_count_bonus + 0.15*has_exotoxin + 0.15*has_effector, 1.0)]
d[n_vf==0 | is.na(n_vf), S_VF := 0.0]

# S_BM
cat("Computing S_BM ...\n")
d[, has_bm := n_metal > 0 & !is.na(n_metal)]
d[, bm_base := ifelse(has_bm, 0.25, 0.0)]
d[, bm_count_bonus := pmin(pmax(n_metal, 0) * 0.04, 0.35)]
d[, has_mer := grepl("mer|mercury", gene_bacmet, ignore.case=TRUE)]
d[, has_qac := grepl("qac|quaternary|disinfectant", gene_bacmet, ignore.case=TRUE)]
d[, has_ars_cop := grepl("ars|cop|sil|czc|cad|pco", gene_bacmet, ignore.case=TRUE)]
d[is.na(has_mer), has_mer := FALSE]; d[is.na(has_qac), has_qac := FALSE]; d[is.na(has_ars_cop), has_ars_cop := FALSE]
d[, S_BM := pmin(bm_base + bm_count_bonus + 0.15*has_mer + 0.15*has_qac + 0.10*has_ars_cop, 1.0)]
d[!has_bm | is.na(n_metal), S_BM := 0.0]
cat(sprintf("  S_BM: mean=%.3f, %.1f%% PSCs with S_BM>0\n", mean(d$S_BM,na.rm=TRUE), 100*mean(d$S_BM>0,na.rm=TRUE)))

s9_cols <- c("S_ARG","S_VF","S_MOB","S_HOST","S_REP","S_SIZE","S_GEO","S_HAB","S_GROW")
s10_cols <- c(s9_cols, "S_BM")
for (col in s10_cols) d[is.na(get(col)), (col):=0]
d[, S9 := S_ARG*W9["S_ARG"] + S_VF*W9["S_VF"] + S_MOB*W9["S_MOB"] + S_HOST*W9["S_HOST"] +
  S_REP*W9["S_REP"] + S_SIZE*W9["S_SIZE"] + S_GEO*W9["S_GEO"] + S_HAB*W9["S_HAB"] + S_GROW*W9["S_GROW"]]
d[, S10 := S_ARG*W10["S_ARG"] + S_VF*W10["S_VF"] + S_MOB*W10["S_MOB"] + S_HOST*W10["S_HOST"] +
  S_REP*W10["S_REP"] + S_SIZE*W10["S_SIZE"] + S_BM*W10["S_BM"] + S_GEO*W10["S_GEO"] + S_HAB*W10["S_HAB"] + S_GROW*W10["S_GROW"]]
d[, S9_norm := S9/sum(W9)]; d[, S10_norm := S10/sum(W10)]

# Per-replicon tables
cat("Building per-replicon tables ...\n")
dt <- d[!is.na(replicon_primary) & replicon_primary!=""]
rep10_basic <- dt[, .(n_PSC=.N, n_plasmid=uniqueN(plasmid_acc), mean_length_kb=mean(length_avg/1000,na.rm=TRUE),
  mean_ARG=mean(n_arg,na.rm=TRUE), mean_VF=mean(n_vf,na.rm=TRUE), mean_metal=mean(n_metal,na.rm=TRUE),
  pct_ARG=100*mean(n_arg>0,na.rm=TRUE), pct_VF=100*mean(n_vf>0,na.rm=TRUE), pct_BM=100*mean(has_bm,na.rm=TRUE),
  pct_mer=100*mean(has_mer,na.rm=TRUE), pct_qac=100*mean(has_qac,na.rm=TRUE),
  pct_human=100*mean(host_human==TRUE,na.rm=TRUE),
  pct_conj=100*mean(mobility_class %in% c("conjugative_complete","conjugative_likely"),na.rm=TRUE)), by=replicon_primary]
rep10_s <- dt[, lapply(.SD, mean, na.rm=TRUE), by=replicon_primary, .SDcols=s10_cols]
setnames(rep10_s, s10_cols, paste0("mean_", s10_cols))
rep10 <- merge(rep10_basic, rep10_s, by="replicon_primary")
rep10[, S9_total := mean_S_ARG*W9["S_ARG"] + mean_S_VF*W9["S_VF"] + mean_S_MOB*W9["S_MOB"] + mean_S_HOST*W9["S_HOST"] +
  mean_S_REP*W9["S_REP"] + mean_S_SIZE*W9["S_SIZE"] + mean_S_GEO*W9["S_GEO"] + mean_S_HAB*W9["S_HAB"] + mean_S_GROW*W9["S_GROW"]]
rep10[, S10_total := mean_S_ARG*W10["S_ARG"] + mean_S_VF*W10["S_VF"] + mean_S_MOB*W10["S_MOB"] + mean_S_HOST*W10["S_HOST"] +
  mean_S_REP*W10["S_REP"] + mean_S_SIZE*W10["S_SIZE"] + mean_S_BM*W10["S_BM"] + mean_S_GEO*W10["S_GEO"] + mean_S_HAB*W10["S_HAB"] + mean_S_GROW*W10["S_GROW"]]
rep10[, S9_norm := S9_total/sum(W9)]; rep10[, S10_norm := S10_total/sum(W10)]
setorder(rep10, -S10_norm); rep10[, rank10_all := .I]; rep10[, major := n_PSC >= 50]
rep9 <- rep10[major==TRUE][order(-S9_norm)]; rep9[, rank9 := .I]
rep10 <- merge(rep10, rep9[, .(replicon_primary, rank9)], by="replicon_primary", all.x=TRUE)
rep10[major==TRUE, rank10 := frank(-S10_norm, ties.method="first")]
rep10[, rank_change := rank9 - rank10]
setorder(rep10, rank10_all)
num_cols <- names(rep10)[sapply(rep10, is.numeric)]
rep10[, (num_cols) := lapply(.SD, function(x) round(x, 4)), .SDcols=num_cols]
fwrite(rep10, file.path(tab_dir,"tab_replicon_risk_10dimensions.csv"))

# Fig 25
cat("\n=== Fig 25: S_BM analysis ===\n")
top20 <- rep10[major==TRUE][order(rank10)][1:20]$replicon_primary
bm_dt <- dt[replicon_primary %in% top20]; bm_dt[, replicon_primary := factor(replicon_primary, levels=rev(top20))]
p25a <- ggplot(bm_dt, aes(x=replicon_primary, y=S_BM)) +
  geom_violin(fill="#fdae61", alpha=0.6, scale="width") + geom_boxplot(width=0.12, outlier.size=0.1, alpha=0.7) + coord_flip() +
  labs(x=NULL, y="S_BM score", title="Biocide/metal resistance score by replicon") + theme_pub
set.seed(42); samp <- d[sample(.N, min(100000,.N))]
p25b <- ggplot(samp, aes(x=S_ARG, y=S_BM)) + geom_point(alpha=0.08, size=0.5, colour="#d73027") +
  geom_smooth(method="lm", colour="black", linewidth=0.6, se=FALSE) +
  annotate("text", x=0.05, y=0.9, label=sprintf("Spearman rho = %.3f", cor(samp$S_ARG, samp$S_BM, method="spearman", use="complete.obs")), hjust=0, size=3.5) +
  labs(x="S_ARG", y="S_BM", title="S_ARG vs S_BM correlation") + theme_pub
save_plot(p25a|p25b, "fig25_SBM_distribution", 14, 8)

# Fig 26
cat("\n=== Fig 26: 9-dim vs 10-dim comparison ===\n")
major <- rep10[major==TRUE]
rank_cor_major <- cor(major$rank9, major$rank10, method="spearman", use="complete.obs")
cat(sprintf("  Rank correlation (major): rho=%.4f\n", rank_cor_major))
p26a <- ggplot(major, aes(x=rank9, y=rank10, colour=S10_norm)) + geom_point(alpha=0.7, size=2) +
  geom_abline(intercept=0, slope=1, linetype="dashed", colour="grey50") +
  geom_text_repel(data=major[abs(rank_change)>=5 | rank10<=15], aes(label=replicon_primary), size=2.8, max.overlaps=25) +
  scale_colour_gradient(low="#fee08b", high="#b30000", name="S10") +
  labs(x="Rank (9-dim)", y="Rank (10-dim)", title=sprintf("9-dim vs 10-dim ranking (rho=%.4f)", rank_cor_major)) + theme_pub
chg <- major[abs(rank_change)>=3][order(rank_change)]
if (nrow(chg)>0) {
  chg[, direction := ifelse(rank_change>0, "Up (higher risk with S_BM)", "Down")]
  chg[, replicon_primary := factor(replicon_primary, levels=chg$replicon_primary)]
  p26b <- ggplot(chg, aes(x=replicon_primary, y=rank_change, fill=direction)) + geom_col() + coord_flip() +
    scale_fill_manual(values=c("Up (higher risk with S_BM)"="#d73027","Down"="#4575b4"), name="Direction") +
    labs(x=NULL, y="Rank change", title="Replicons with rank change >=3") + theme_pub
} else p26b <- ggplot() + annotate("text",x=0,y=0,label="No rank changes >=3") + theme_void()
save_plot(p26a/p26b, "fig26_rank_comparison", 12, 12)

# Fig 27: RF
cat("\n=== Fig 27: RF validation with S_BM ===\n")
ml_dt <- dt[, .(S_ARG, S_VF, S_MOB, S_HOST, S_REP, S_SIZE, S_GEO, S_HAB, S_GROW, S_BM,
  n_high_risk_arg, n_arg, n_vf, n_metal, mobility_class,
  highrisk=factor(ifelse(n_high_risk_arg>0, "Yes","No")), mdr_vf=factor(ifelse(n_arg>0 & n_vf>0, "Yes","No")),
  bm_carrier=factor(ifelse(n_metal>0, "Yes","No")))]
ml_dt <- ml_dt[complete.cases(ml_dt)]
set.seed(42); ml_sample <- ml_dt[sample(.N, min(100000,.N))]
rf1 <- randomForest(highrisk ~ S_ARG+S_VF+S_MOB+S_HOST+S_REP+S_SIZE+S_GEO+S_HAB+S_GROW+S_BM, data=ml_sample, ntree=300, importance=TRUE)
rf1_auc <- as.numeric(auc(roc(ml_sample$highrisk, predict(rf1,type="prob")[,", quiet=TRUE)))
rf2 <- randomForest(mdr_vf ~ S_ARG+S_VF+S_MOB+S_HOST+S_REP+S_SIZE+S_GEO+S_HAB+S_GROW+S_BM, data=ml_sample, ntree=300, importance=TRUE)
rf2_auc <- as.numeric(auc(roc(ml_sample$mdr_vf, predict(rf2,type="prob")[,", quiet=TRUE)))
ml_sample[, conj := factor(ifelse(mobility_class %in% c("conjugative_complete","conjugative_likely"),"Yes","No"))]
rf3 <- randomForest(conj ~ S_ARG+S_VF+S_MOB+S_HOST+S_REP+S_SIZE+S_GEO+S_HAB+S_GROW+S_BM, data=ml_sample, ntree=300, importance=TRUE)
rf3_auc <- as.numeric(auc(roc(ml_sample$conj, predict(rf3,type="prob")[,", quiet=TRUE)))
rf4 <- randomForest(bm_carrier ~ S_ARG+S_VF+S_MOB+S_HOST+S_REP+S_SIZE+S_GEO+S_HAB+S_GROW, data=ml_sample, ntree=300, importance=TRUE)
rf4_auc <- as.numeric(auc(roc(ml_sample$bm_carrier, predict(rf4,type="prob")[,", quiet=TRUE)))
cat(sprintf("  RF AUCs: highrisk=%.3f, mdrvf=%.3f, conj=%.3f, bm=%.3f\n", rf1_auc, rf2_auc, rf3_auc, rf4_auc))
imp_list <- lapply(list(rf1,rf2,rf3,rf4), function(rf) { imp <- as.data.table(importance(rf), keep.rownames="feature"); setnames(imp, "MeanDecreaseGini", "MDG"); imp[, MDG_norm := MDG/sum(MDG)]; imp[, .(feature, MDG_norm)] })
imp <- Reduce(function(x,y) merge(x,y,by="feature",all=TRUE), imp_list)
setnames(imp, c("feature","MDG_highrisk","MDG_mdrvf","MDG_conj","MDG_bm"))
imp[, mean_MDG := rowMeans(.SD, na.rm=TRUE), .SDcols=c("MDG_highrisk","MDG_mdrvf","MDG_conj")]
imp[, user_weight_10 := W10[feature]]; imp[is.na(user_weight_10), user_weight_10 := 0]
imp[, user_weight_norm := user_weight_10/sum(user_weight_10)]
setorder(imp, -mean_MDG)
fwrite(imp, file.path(tab_dir,"tab_rf_importance_10dim.csv"))
imp_long <- melt(imp, id.vars="feature", measure.vars=c("user_weight_norm","MDG_highrisk","MDG_mdrvf","MDG_conj"), variable.name="method", value.name="importance")
imp_long[, method := factor(method, levels=c("user_weight_norm","MDG_highrisk","MDG_mdrvf","MDG_conj"),
  labels=c("Expert weights (10-dim)", sprintf("RF: high-risk ARG (AUC=%.3f)",rf1_auc), sprintf("RF: MDR-VF (AUC=%.3f)",rf2_auc), sprintf("RF: conjugative (AUC=%.3f)",rf3_auc)))]
imp_long[, feature := factor(feature, levels=imp$feature)]
p27 <- ggplot(imp_long, aes(x=feature, y=importance, fill=method)) + geom_col(position="dodge", alpha=0.85) +
  scale_fill_brewer(palette="Set2", name="Method") + labs(x=NULL, y="Relative importance", title="Expert weights vs RF feature importance (10 dim)") +
  theme_pub + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.text=element_text(size=8))
save_plot(p27, "fig27_rf_importance_10dim", 12, 6)

# Fig 28: S_BM weight sensitivity
cat("\n=== Fig 28: S_BM weight sensitivity ===\n")
set.seed(123); bm_weights <- seq(0, 0.30, by=0.02)
sens_res <- rbindlist(lapply(bm_weights, function(w) {
  w_try <- W10; w_try["S_BM"] <- w
  other_sum <- sum(w_try[names(w_try) != "S_BM"])
  w_try[names(w_try) != "S_BM"] <- w_try[names(w_try) != "S_BM"] * (1 - w) / other_sum
  s_try <- as.numeric(as.matrix(dt[, .SD, .SDcols=s10_cols]) %*% w_try[s10_cols])
  tmp <- data.table(replicon_primary=dt$replicon_primary, s=s_try)[, .(s=mean(s)), by=replicon_primary]
  tmp <- merge(tmp, rep10[, .(replicon_primary, major, rank10)], by="replicon_primary")[major==TRUE][order(-s)]
  tmp[, rank := .I]
  data.table(bm_weight=w, spearman=cor(tmp$rank, tmp$rank10, method="spearman"),
    top10_overlap=length(intersect(tmp[1:10]$replicon_primary, rep10[major==TRUE][rank10<=10]$replicon_primary)),
    top20_overlap=length(intersect(tmp[1:20]$replicon_primary, rep10[major==TRUE][rank10<=20]$replicon_primary)))
}))
fwrite(sens_res, file.path(tab_dir,"tab_SBM_weight_sensitivity.csv"))
bm_final <- as.numeric(W10["S_BM"])
p28a <- ggplot(sens_res, aes(x=bm_weight, y=spearman)) + geom_line(linewidth=0.8, colour="#0570b0") + geom_point(size=2) +
  geom_vline(xintercept=bm_final, linetype="dashed", colour="red") +
  annotate("text", x=bm_final, y=Inf, vjust=2, hjust=-0.2, label=sprintf("final = %.3f", bm_final), size=3, colour="red") +
  labs(x="S_BM weight", y="Spearman rho", title="Ranking stability across S_BM weights") + theme_pub
p28b <- ggplot(sens_res, aes(x=bm_weight)) +
  geom_line(aes(y=top10_overlap/10, colour="Top-10"), linewidth=0.8) + geom_line(aes(y=top20_overlap/20, colour="Top-20"), linewidth=0.8) +
  geom_point(aes(y=top10_overlap/10, colour="Top-10"), size=2) + geom_point(aes(y=top20_overlap/20, colour="Top-20"), size=2) +
  scale_colour_brewer(palette="Set1", name="Overlap") + scale_y_continuous(labels=percent, limits=c(0,1)) +
  labs(x="S_BM weight", y="Overlap", title="Top-N overlap") + theme_pub
save_plot(p28a|p28b, "fig28_SBM_sensitivity", 12, 4.5)

# Fig 29: 10-dim heatmap
cat("\n=== Fig 29: 10-dim visualization ===\n")
top30 <- rep10[major==TRUE][order(rank10)][1:30]$replicon_primary
heat10 <- melt(rep10[replicon_primary %in% top30], id.vars="replicon_primary", measure.vars=paste0("mean_",s10_cols), variable.name="component", value.name="score")
heat10[, component := gsub("mean_","",component)]; heat10[, component := factor(component, levels=s10_cols)]
heat10[, replicon_primary := factor(replicon_primary, levels=rev(top30))]
heat10[, weight := W10[component]]; heat10[, contribution := score*weight]
p29a <- ggplot(heat10, aes(x=component, y=replicon_primary, fill=score)) +
  geom_tile(colour="white", linewidth=0.3) + geom_text(aes(label=sprintf("%.2f",score)), size=1.8) +
  scale_fill_gradientn(colours=c("#f7fcf5","#c7e9c0","#74c476","#238b45","#00441b"), limits=c(0,1), name="Score") +
  scale_x_discrete(position="top") + labs(x=NULL, y=NULL, title="10-dimension risk component scores (top 30)") +
  theme_pub + theme(axis.text.x=element_text(angle=45, hjust=0, face="bold"), axis.text.y=element_text(size=8))
p29b <- ggplot(heat10, aes(x=reorder(replicon_primary, contribution, sum), y=contribution, fill=component)) + geom_col() + coord_flip() +
  scale_fill_brewer(palette="Set3", name="Component") + labs(x=NULL, y="Weighted contribution", title="Stacked 10-dim contributions") +
  theme_pub + theme(legend.position="bottom")
save_plot(p29a/p29b, "fig29_10dim_heatmap", 14, 16)
cat("\n=== All done ===\n")
