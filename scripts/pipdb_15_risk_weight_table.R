#!/usr/bin/env Rscript
# =============================================================================
# pipdb_15_risk_weight_table.R
# Per-replicon composite risk table with data-driven final weights (10 dimensions).
# Includes: S_VF, S_BM computation, RF feature-importance validation,
#           weight sensitivity analysis, logistic regression, visualization.
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
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)
save_plot <- function(plot, name, w=10, h=8) {
  ggsave(file.path(fig_dir, paste0(name,".pdf")), plot, width=w, height=h, limitsize=FALSE)
  ggsave(file.path(fig_dir, paste0(name,".png")), plot, width=w, height=h, dpi=300, limitsize=FALSE)
  cat("  saved", name, "\n")
}
theme_pub <- theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey90", colour=NA))

comp_cols <- c("S_ARG","S_VF","S_MOB","S_HOST","S_REP","S_SIZE","S_BM","S_GEO","S_HAB","S_GROW")
wf <- file.path(tab_dir, "tab_weight_comparison.csv")
if (file.exists(wf)) {
  wc <- fread(wf); W <- wc$Final; names(W) <- wc$Component
  cat("Loaded data-driven final weights from tab_weight_comparison.csv\n")
} else {
  W <- c(S_ARG=0.2448, S_VF=0.1096, S_MOB=0.2041, S_HOST=0.0282, S_REP=0.0030,
         S_SIZE=0.1808, S_BM=0.2112, S_GEO=0.0015, S_HAB=0.0022, S_GROW=0.0147)
  cat("Using built-in final consensus weights\n")
}
W_norm <- W / sum(W)
cat(sprintf("Weight sum = %.4f\n", sum(W))); print(round(W,4))

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

# Composite S
for (col in comp_cols) d[is.na(get(col)), (col):=0]
d[, S_total_raw := S_ARG*W["S_ARG"] + S_VF*W["S_VF"] + S_MOB*W["S_MOB"] + S_HOST*W["S_HOST"] +
  S_REP*W["S_REP"] + S_SIZE*W["S_SIZE"] + S_BM*W["S_BM"] + S_GEO*W["S_GEO"] + S_HAB*W["S_HAB"] + S_GROW*W["S_GROW"]]
d[, S_total_norm := S_total_raw / sum(W)]

# Per-replicon table
cat("Building per-replicon table ...\n")
dt <- d[!is.na(replicon_primary) & replicon_primary!=""]
rep_table <- dt[, .(n_PSC=.N, n_plasmid=uniqueN(plasmid_acc), n_species=uniqueN(species_name), n_country=uniqueN(country),
  mean_length_kb=mean(length_avg/1000,na.rm=TRUE), mean_ARG=mean(n_arg,na.rm=TRUE), mean_VF=mean(n_vf,na.rm=TRUE),
  pct_ARG=100*mean(n_arg>0,na.rm=TRUE), pct_VF=100*mean(n_vf>0,na.rm=TRUE),
  pct_highrisk_ARG=100*mean(n_high_risk_arg>0,na.rm=TRUE),
  pct_conjugative=100*mean(mobility_class %in% c("conjugative_complete","conjugative_likely"),na.rm=TRUE),
  pct_integron=100*mean(has_integron==TRUE,na.rm=TRUE), pct_human=100*mean(host_human==TRUE,na.rm=TRUE),
  mean_growth=mean(annual_growth_rate,na.rm=TRUE),
  S_ARG=mean(S_ARG,na.rm=TRUE), S_VF=mean(S_VF,na.rm=TRUE), S_MOB=mean(S_MOB,na.rm=TRUE),
  S_HOST=mean(S_HOST,na.rm=TRUE), S_REP=mean(S_REP,na.rm=TRUE), S_SIZE=mean(S_SIZE,na.rm=TRUE),
  S_BM=mean(S_BM,na.rm=TRUE), S_GEO=mean(S_GEO,na.rm=TRUE), S_HAB=mean(S_HAB,na.rm=TRUE), S_GROW=mean(S_GROW,na.rm=TRUE),
  S_total_norm=mean(S_total_norm,na.rm=TRUE)), by=replicon_primary]
setorder(rep_table, -S_total_norm)
rep_table[, rank := NA_integer_]; rep_table[n_PSC >= 50, rank := .I]; rep_table[, major := n_PSC >= 50]
num_cols <- names(rep_table)[sapply(rep_table, is.numeric)]
rep_table[, (num_cols) := lapply(.SD, function(x) round(x, 4)), .SDcols=num_cols]
fwrite(rep_table, file.path(tab_dir,"tab_replicon_risk_weight_table.csv"))
cat(sprintf("  Table saved: %d replicons (%d major)\n", nrow(rep_table), sum(rep_table$major)))

# Save per-PSC final scores (for use by other scripts)
fwrite(d[, .(id, S_total_raw, S_total_norm)], file.path(tab_dir,"tab_psc_final_scores.csv"))

# Visualization
top30 <- rep_table[major==TRUE]$replicon_primary[1:min(30,sum(rep_table$major))]
heat_dt <- melt(rep_table[replicon_primary %in% top30], id.vars="replicon_primary", measure.vars=comp_cols, variable.name="component", value.name="score")
heat_dt[, replicon_primary := factor(replicon_primary, levels=rev(top30))]
heat_dt[, component := factor(component, levels=comp_cols)]
heat_dt[, weight := W[as.character(component)]]
heat_dt[, contribution := score * weight]
p21a <- ggplot(heat_dt, aes(x=component, y=replicon_primary, fill=score)) +
  geom_tile(colour="white", linewidth=0.3) + geom_text(aes(label=sprintf("%.2f",score)), size=2.2) +
  scale_fill_gradientn(colours=c("#f7fcf5","#c7e9c0","#74c476","#238b45","#00441b"), limits=c(0,1), name="Score") +
  scale_x_discrete(position="top") + labs(x=NULL, y=NULL, title="Per-replicon risk component scores") +
  theme_pub + theme(axis.text.x=element_text(angle=45, hjust=0, face="bold"), axis.text.y=element_text(size=8))
p21c <- ggplot(heat_dt, aes(x=reorder(replicon_primary, contribution, sum), y=contribution, fill=component)) +
  geom_col() + coord_flip() + scale_fill_brewer(palette="Set1", name="Component") +
  labs(x=NULL, y="Weighted contribution", title="Stacked weighted risk contributions") + theme_pub + theme(legend.position="bottom")
save_plot(p21a, "fig21a_risk_component_heatmap", 12, 10)
save_plot(p21c, "fig21b_risk_stacked_bar", 10, 10)

# ML validation
cat("\n=== ML validation ===\n")
ml_dt <- dt[!is.na(S_ARG), .(S_ARG, S_VF, S_MOB, S_HOST, S_REP, S_SIZE, S_BM, S_GEO, S_HAB, S_GROW,
  n_high_risk_arg, n_arg, n_vf, mobility_class,
  highrisk=factor(ifelse(n_high_risk_arg>0, "Yes", "No")), mdr_vf=factor(ifelse(n_arg>0 & n_vf>0, "Yes", "No")))]
ml_dt <- ml_dt[complete.cases(ml_dt)]
set.seed(42)
ml_sample <- ml_dt[sample(.N, min(100000, .N))]
cat("  Training RF (high-risk ARG)...\n")
rf1 <- randomForest(highrisk ~ S_ARG + S_VF + S_MOB + S_HOST + S_REP + S_SIZE + S_BM + S_GEO + S_HAB + S_GROW,
                    data=ml_sample, ntree=300, importance=TRUE)
rf1_auc <- as.numeric(auc(roc(ml_sample$highrisk, predict(rf1, type="prob")[,", quiet=TRUE)))
cat(sprintf("  RF high-risk ARG AUC = %.3f\n", rf1_auc))
cat("  Training RF (MDR-VF)...\n")
rf2 <- randomForest(mdr_vf ~ S_ARG + S_VF + S_MOB + S_HOST + S_REP + S_SIZE + S_BM + S_GEO + S_HAB + S_GROW,
                    data=ml_sample, ntree=300, importance=TRUE)
rf2_auc <- as.numeric(auc(roc(ml_sample$mdr_vf, predict(rf2, type="prob")[,", quiet=TRUE)))
cat(sprintf("  RF MDR-VF AUC = %.3f\n", rf2_auc))
ml_sample[, conj := factor(ifelse(mobility_class %in% c("conjugative_complete","conjugative_likely"),"Yes","No"))]
cat("  Training RF (conjugative)...\n")
rf3 <- randomForest(conj ~ S_ARG + S_VF + S_MOB + S_HOST + S_REP + S_SIZE + S_BM + S_GEO + S_HAB + S_GROW,
                    data=ml_sample, ntree=300, importance=TRUE)
rf3_auc <- as.numeric(auc(roc(ml_sample$conj, predict(rf3, type="prob")[,", quiet=TRUE)))
cat(sprintf("  RF conjugative AUC = %.3f\n", rf3_auc))
imp1 <- as.data.table(importance(rf1), keep.rownames="feature"); setnames(imp1, "MeanDecreaseGini", "MDG_highrisk")
imp2 <- as.data.table(importance(rf2), keep.rownames="feature"); setnames(imp2, "MeanDecreaseGini", "MDG_mdrvf")
imp3 <- as.data.table(importance(rf3), keep.rownames="feature"); setnames(imp3, "MeanDecreaseGini", "MDG_conj")
imp <- merge(merge(imp1[, .(feature, MDG_highrisk)], imp2[, .(feature, MDG_mdrvf)], by="feature"), imp3[, .(feature, MDG_conj)], by="feature")
imp[, `:=`(MDG_highrisk_norm=MDG_highrisk/sum(MDG_highrisk), MDG_mdrvf_norm=MDG_mdrvf/sum(MDG_mdrvf), MDG_conj_norm=MDG_conj/sum(MDG_conj))]
imp[, final_weight := W[feature]]; imp[, final_weight_norm := final_weight/sum(final_weight)]
imp[, mean_MDG_norm := (MDG_highrisk_norm+MDG_mdrvf_norm+MDG_conj_norm)/3]
setorder(imp, -mean_MDG_norm)
fwrite(imp, file.path(tab_dir,"tab_rf_feature_importance.csv"))

# Sensitivity analysis
cat("\n=== Sensitivity analysis ===\n")
set.seed(123); n_pert <- 100; rank_cor <- numeric(n_pert); top10_overlap <- numeric(n_pert)
original_ranks <- rep_table[major==TRUE, .(replicon_primary, rank)]
rep_n <- dt[, .N, by=replicon_primary]
for (i in seq_len(n_pert)) {
  w_pert <- W * runif(length(W), 0.7, 1.3); w_pert <- w_pert / sum(w_pert)
  s_pert <- as.matrix(dt[, .SD, .SDcols=comp_cols]) %*% w_pert
  tmp <- data.table(replicon_primary=dt$replicon_primary, s=as.numeric(s_pert))[, .(s=mean(s)), by=replicon_primary]
  tmp <- merge(tmp, rep_n, by="replicon_primary")[N>=50][order(-s)]; tmp[, rank_pert := .I]
  tmp <- merge(tmp, original_ranks, by="replicon_primary")
  rank_cor[i] <- cor(tmp$rank, tmp$rank_pert, method="spearman")
  top10_overlap[i] <- length(intersect(original_ranks[rank<=10]$replicon_primary, tmp[rank_pert<=10]$replicon_primary))
}
cat(sprintf("  Spearman rho: mean=%.4f; Top-10 overlap: mean=%.1f/10\n", mean(rank_cor), mean(top10_overlap)))
sens_dt <- data.table(iteration=seq_len(n_pert), spearman_rho=rank_cor, top10_overlap=top10_overlap)
fwrite(sens_dt, file.path(tab_dir,"tab_weight_sensitivity.csv"))

# Logistic regression
cat("\n=== Logistic regression ===\n")
glm1 <- glm(highrisk ~ S_ARG + S_VF + S_MOB + S_HOST + S_REP + S_SIZE + S_BM + S_GEO + S_HAB + S_GROW,
            data=ml_sample, family=binomial)
glm_coef <- as.data.table(summary(glm1)$coefficients, keep.rownames="feature")
setnames(glm_coef, c("feature","Estimate","Std.Error","z","p"))
glm_coef[, OR := exp(Estimate)]; glm_coef[, OR_low := exp(Estimate-1.96*Std.Error)]; glm_coef[, OR_high := exp(Estimate+1.96*Std.Error)]
glm_coef[, sig := ifelse(p<0.001,"***", ifelse(p<0.01,"**", ifelse(p<0.05,"*","ns")))]
fwrite(glm_coef, file.path(tab_dir,"tab_logistic_coefficients.csv"))
cat("\n=== All done ===\n")
