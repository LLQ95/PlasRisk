#!/usr/bin/env Rscript
# =============================================================================
# pipdb_19_case_study.R
# PlasRisk application case study: scoring real plasmids + cross-validation.
# Cases: pNDM-1 (NDM carbapenemase), pHNSHP45 (mcr-1), pColE1-like (cryptic).
# =============================================================================
required_pkgs <- c("data.table", "ggplot2", "patchwork")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, repos="https://cloud.r-project.org", dependencies=TRUE, quiet=TRUE)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) >= 1) args[1] else "results"
fig_dir <- file.path(res_dir, "figures", "descriptive")
tab_dir <- file.path(res_dir, "tables")
case_dir <- file.path(res_dir, "case_study")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE); dir.create(tab_dir, showWarnings=FALSE, recursive=TRUE); dir.create(case_dir, showWarnings=FALSE, recursive=TRUE)

comp_cols <- c("S_ARG","S_VF","S_MOB","S_HOST","S_REP","S_SIZE","S_BM","S_GEO","S_HAB","S_GROW")
wf <- file.path(tab_dir, "tab_weight_comparison.csv")
if (file.exists(wf)) { wc <- fread(wf); final_w <- wc$Final; names(final_w) <- wc$Component
} else { final_w <- c(S_ARG=0.2448, S_VF=0.1096, S_MOB=0.2041, S_HOST=0.0282, S_REP=0.0030, S_SIZE=0.1808, S_BM=0.2112, S_GEO=0.0015, S_HAB=0.0022, S_GROW=0.0147) }

cat("=== PlasRisk case study ===\n")
case1 <- data.table(seq_id="pNDM-1 (NC_019050)", length_bp=50000, n_arg=5, n_vf=0, n_bm=3, replicon="IncN",
  has_t4cp=TRUE, has_relaxase=TRUE, has_oriT=TRUE, has_aux=TRUE, high_risk="blaNDM",
  description="NDM-1 carbapenemase plasmid, IncN, conjugative")
case2 <- data.table(seq_id="pHNSHP45 (NC_022535)", length_bp=33000, n_arg=2, n_vf=0, n_bm=1, replicon="IncX4",
  has_t4cp=FALSE, has_relaxase=TRUE, has_oriT=TRUE, has_aux=FALSE, high_risk="mcr-1",
  description="MCR-1 colistin resistance plasmid, IncX4, mobilizable")
case3 <- data.table(seq_id="pColE1-like (cryptic)", length_bp=6500, n_arg=0, n_vf=1, n_bm=0, replicon="Col(MG828)",
  has_t4cp=FALSE, has_relaxase=FALSE, has_oriT=FALSE, has_aux=FALSE, high_risk="",
  description="Cryptic colicin plasmid, no ARGs, non-mobilizable")
cases <- rbindlist(list(case1, case2, case3))

rep_lookup <- fread(text="replicon,S_REP,S_GEO,S_HAB,S_GROW,S_HOST
IncX3,0.85,0.80,0.80,0.85,0.80
IncN,0.80,0.85,0.80,0.75,0.85
IncFII,0.65,0.90,0.85,0.70,0.85
IncFII(K),0.70,0.80,0.75,0.80,0.80
IncI1,0.60,0.85,0.80,0.65,0.80
IncX4,0.75,0.80,0.75,0.90,0.75
IncR,0.70,0.70,0.70,0.65,0.70
IncQ1,0.65,0.80,0.75,0.50,0.90
ColRNAI,0.35,0.75,0.70,0.55,0.65
Col440I,0.30,0.70,0.65,0.45,0.60
Col(MG828),0.10,0.40,0.40,0.30,0.40
ColpVC,0.05,0.30,0.30,0.25,0.30
Unknown,0.30,0.30,0.30,0.30,0.50")

score_case <- function(case, rep_lookup) {
  rl <- rep_lookup[replicon == case$replicon]; if (nrow(rl)==0) rl <- rep_lookup[replicon=="Unknown"]
  S_ARG <- min(0.25 + min(case$n_arg*0.05, 0.35) + 0.20*(case$n_arg>0) + 0.20*(nchar(case$high_risk)>0), 1.0)
  S_VF <- if (case$n_vf==0) 0 else min(0.30 + min(case$n_vf*0.03, 0.40), 1.0)
  n_mob <- sum(c(case$has_t4cp, case$has_relaxase, case$has_oriT, case$has_aux))
  S_MOB <- c(0.05, 0.30, 0.55, 0.85, 1.0)[n_mob + 1]
  S_SIZE <- 1/(1+exp(-(case$length_bp/1000-30)/15))
  S_BM <- if (case$n_bm==0) 0 else min(0.25 + min(case$n_bm*0.04, 0.35) + 0.15*(case$n_bm>=2), 1.0)
  comps <- c(S_ARG=S_ARG, S_VF=S_VF, S_MOB=S_MOB, S_HOST=rl$S_HOST, S_REP=rl$S_REP, S_SIZE=S_SIZE,
    S_BM=S_BM, S_GEO=rl$S_GEO, S_HAB=rl$S_HAB, S_GROW=rl$S_GROW)
  S_norm <- sum(comps * final_w[names(comps)])
  grade <- if (S_norm>=0.60) c("A","Very High") else if (S_norm>=0.45) c("B","High") else if (S_norm>=0.30) c("C","Moderate") else if (S_norm>=0.15) c("D","Low") else c("E","Minimal")
  data.table(seq_id=case$seq_id, t(comps), S_norm=round(S_norm,4), grade=grade[1], grade_label=grade[2])
}
results <- rbindlist(lapply(seq_len(nrow(cases)), function(i) score_case(as.list(cases[i]), rep_lookup)))
cat("\n  Case study results:\n"); print(results[, .(seq_id, S_ARG, S_VF, S_MOB, S_BM, S_REP, S_SIZE, S_norm, grade)])
fwrite(results, file.path(tab_dir, "tab_case_study_scores.csv"))

# Population distribution
cat("\n  Loading PIPdb scores for distribution ...\n")
master <- fread(file.path(res_dir, "psc_master.tsv"), na.strings=c("\\N","","NA"))
risk <- fread(file.path(res_dir, "psc_risk_scores.tsv"), na.strings=c("\\N","","NA"))
d <- merge(master, risk, by="id", suffixes=c("",".r"))
d[,n_vf:=fifelse(is.na(n_vf),0L,n_vf)]; d[,n_metal:=fifelse(is.na(n_metal),0L,n_metal)]
d[,S_VF:=0.0]; d[n_vf>0,S_VF:=pmin(0.30+pmin(n_vf*0.03,0.40)+0.15*grepl("Exotoxin",vf_category)+0.15*grepl("Effector delivery",vf_category),1.0)]
d[,has_mer:=as.integer(grepl("mer",gene_bacmet,ignore.case=TRUE))]; d[,has_qac:=as.integer(grepl("qac",gene_bacmet,ignore.case=TRUE))]; d[,has_ars:=as.integer(grepl("ars|cop|sil",gene_bacmet,ignore.case=TRUE))]
d[is.na(has_mer),has_mer:=0L]; d[is.na(has_qac),has_qac:=0L]; d[is.na(has_ars),has_ars:=0L]
d[,S_BM:=0.0]; d[n_metal>0,S_BM:=pmin(0.25+pmin(n_metal*0.04,0.35)+0.15*has_mer+0.15*has_qac+0.10*has_ars,1.0)]
for (cc in comp_cols) d[[cc]] <- fifelse(is.na(d[[cc]]),0.0,as.numeric(d[[cc]]))
d[,S_norm:=as.matrix(.SD)%*%final_w[comp_cols],.SDcols=comp_cols]
cat(sprintf("  Population: %d PSCs, median S_norm=%.3f\n", nrow(d), median(d$S_norm,na.rm=TRUE)))
for (i in seq_len(nrow(results))) { pct <- mean(d$S_norm <= results$S_norm[i], na.rm=TRUE)*100; cat(sprintf("  %s: percentile=%.1f%%\n", results$seq_id[i], pct)); results[i, percentile := round(pct,1)] }

# Figure
cat("\n  Generating case study figure ...\n")
p_dist <- ggplot(d[!is.na(S_norm)], aes(x=S_norm)) +
  geom_histogram(bins=80, fill="grey70", colour="white", linewidth=0.1) +
  geom_vline(data=results, aes(xintercept=S_norm, colour=seq_id), linetype="dashed", linewidth=0.8) +
  scale_colour_manual(values=c("#d94801","#2171b5","#2ca02c")) +
  labs(title="PlasRisk score distribution (PIPdb)", x="PlasRisk S_norm", y="Number of PSCs") +
  theme_bw(base_size=11) + theme(plot.title=element_text(face="bold"), legend.position="none")
res_long <- melt(results, id.vars=c("seq_id","grade","grade_label","percentile"), measure.vars=comp_cols, variable.name="Component", value.name="Score")
res_long[, seq_id := factor(seq_id, levels=results$seq_id)]; res_long[, Component := factor(Component, levels=rev(comp_cols))]
p_heat <- ggplot(res_long, aes(x=seq_id, y=Component, fill=Score)) +
  geom_tile(colour="white", linewidth=1) + geom_text(aes(label=sprintf("%.2f",Score)), size=3, fontface="bold") +
  scale_fill_gradient2(low="#f7f7f7", mid="#fdae6b", high="#b30000", midpoint=0.5, limits=c(0,1)) +
  scale_x_discrete(labels=function(x) gsub(" \\(.*","",x)) +
  labs(title="Component score breakdown", x=NULL, y=NULL, fill="Score") +
  theme_bw(base_size=11) + theme(plot.title=element_text(face="bold"), axis.text.x=element_text(angle=25, hjust=1, face="bold"), panel.grid=element_blank())
grade_tbl <- results[, .(seq_id=gsub(" \\(.*","",seq_id), S_norm, grade, percentile)]
grade_tbl[, label := sprintf("S=%.3f  Grade %s  Pctl=%.0f%%", S_norm, grade, percentile)]
p_grade <- ggplot(grade_tbl, aes(x=reorder(seq_id,-S_norm), y=S_norm, fill=grade)) +
  geom_col(width=0.6, colour="white") + geom_text(aes(label=label), hjust=1.1, size=3.5, fontface="bold") +
  scale_fill_manual(values=c(A="#b30000",B="#fc8d59",C="#fdcc8a",D="#e0f3db",E="#f7f7f7")) +
  coord_flip() + ylim(0,1) + labs(title="Risk grade assignment", x=NULL, y="S_norm") +
  theme_bw(base_size=11) + theme(plot.title=element_text(face="bold"), legend.position="none")
comb <- (p_dist / (p_heat | p_grade)) + plot_annotation(title="PlasRisk case study: scoring clinically relevant plasmids",
  subtitle="pNDM-1 (NDM), pHNSHP45 (mcr-1), pColE1-like (cryptic)",
  theme=theme(plot.title=element_text(face="bold",size=14), plot.subtitle=element_text(size=10))) + plot_layout(heights=c(1,1.2))
ggsave(file.path(fig_dir, "fig26_case_study.png"), comb, width=14, height=10, dpi=300)
ggsave(file.path(fig_dir, "fig26_case_study.pdf"), comb, width=14, height=10)
cat("  saved fig26_case_study\n")

# Example FASTA files for CLI testing
for (i in seq_len(nrow(cases))) {
  cs <- cases[i]; fa_file <- file.path(case_dir, paste0(gsub(" \\(.*","",cs$seq_id), ".fasta"))
  con <- file(fa_file, "w"); writeLines(sprintf(">%s %s", gsub(" \\(.*","",cs$seq_id), cs$description), con)
  set.seed(i); seq_str <- paste(sample(c("A","T","G","C"), cs$length_bp, replace=TRUE), collapse="")
  for (j in seq(1, nchar(seq_str), by=80)) writeLines(substr(seq_str, j, min(j+79, nchar(seq_str))), con)
  close(con); cat(sprintf("  Wrote %s\n", fa_file))
}
fwrite(grade_tbl, file.path(tab_dir, "tab_case_study_summary.csv"))
cat("\nDone. Case study complete.\n")
