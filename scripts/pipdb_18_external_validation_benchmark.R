#!/usr/bin/env Rscript
# pipdb_18_external_validation_benchmark.R
# External independent validation (NCBI plasmids) + benchmarking vs PIPdb ordinal, 9-dim, single dims.
required_pkgs <- c("data.table", "pROC", "ggplot2", "patchwork")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, repos="https://cloud.r-project.org", dependencies=TRUE, quiet=TRUE)
suppressPackageStartupMessages({ library(data.table); library(pROC); library(ggplot2); library(patchwork) })
fasta_lengths <- function(path) { lines <- readLines(path, warn=FALSE); h <- grep("^>", lines); ids <- sub("^>(\\S+).*","\\1",lines[h]); s <- h+1; e <- c(h[-1]-1, length(lines)); setNames(mapply(function(a,b) sum(nchar(lines[a:b])), s, e), ids) }
args <- commandArgs(trailingOnly=TRUE); res_dir <- if (length(args)>=1) args[1] else "results"
fig_dir <- file.path(res_dir,"figures","descriptive"); tab_dir <- file.path(res_dir,"tables"); ext_dir <- file.path(res_dir,"external_validation")
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE); dir.create(tab_dir, showWarnings=FALSE, recursive=TRUE); dir.create(ext_dir, showWarnings=FALSE, recursive=TRUE)

cat("=== A. External validation ===\n")
external_plasmids <- data.table(
  accession=c("NC_019050","NC_024795","NC_019148","NC_022535","NC_017837","NC_020065","NC_019288","NC_023091","NC_019051","NC_021502","NC_025192","NC_031453","NC_019289","NC_022082","NC_017710","NC_019376","NC_024817","NC_025193","NC_032349","NC_040315","NC_001372","NC_002122","NC_000964","NC_001371","NC_000913","NC_002483","NC_005207","NC_001370","NC_007800","NC_013956","NC_010411","NC_004463","NC_007708","NC_010557","NC_011746","NC_013505","NC_016034","NC_017726","NC_019014","NC_020085"),
  label=c(rep("high_risk",20),rep("low_risk",20)), y=c(rep(1L,20),rep(0L,20)))
cat(sprintf("  External test set: %d plasmids\n", nrow(external_plasmids)))
fasta_path <- file.path(ext_dir, "external_plasmids.fasta")
if (!file.exists(fasta_path)) {
  cat("  Downloading from NCBI ...\n")
  for (acc in external_plasmids$accession) {
    out_fa <- file.path(ext_dir, paste0(acc,".fasta"))
    if (!file.exists(out_fa)) { system(sprintf('curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=fasta&retmode=text" -o "%s"', acc, out_fa)); Sys.sleep(0.4) }
  }
  fa_files <- list.files(ext_dir, pattern="\\.fasta$", full.names=TRUE); fa_files <- fa_files[!grepl("external_plasmids", fa_files)]
  if (length(fa_files)>0) system(sprintf('cat %s > "%s"', paste(shQuote(fa_files),collapse=" "), fasta_path))
}
annotate_ext <- function(fa, ed) {
  res <- list()
  for (db in c("card","vfdb","plasmidfinder")) {
    out <- file.path(ed, sprintf("abricate_%s.tsv", db))
    if (!file.exists(out)) system2("abricate", args=c("--db",db,"--minid","75","--mincov","50",fa), stdout=out, stderr=NULL)
    if (file.exists(out) && file.info(out)$size>100) res[[db]] <- fread(out)
  }
  res
}
score_ext <- function(annot, ep, cc, w) {
  lens <- fasta_lengths(fasta_path); sc <- data.table(accession=ep$accession); sc[, length_bp := lens[accession]]
  if (!is.null(annot$card)) { ac <- annot$card[,.N,by=.(SEQUENCE=sub(" .*","",SEQUENCE))]; sc <- merge(sc,ac,by.x="accession",by.y="SEQUENCE",all.x=TRUE); setnames(sc,"N","n_arg"); sc[is.na(n_arg),n_arg:=0L]; sc[,n_highrisk:=0L]; hr <- annot$card[grepl("NDM|KPC|mcr|OXA-48|CTX-M|tet.X|VIM|IMP",GENE)]; if (nrow(hr)>0) { hc <- hr[,uniqueN(GENE),by=.(SEQUENCE=sub(" .*","",SEQUENCE))]; sc[hc,n_highrisk:=i.V1,on=.(accession=SEQUENCE)] } } else sc[,`:=`(n_arg=0L,n_highrisk=0L)]
  if (!is.null(annot$vfdb)) { vc <- annot$vfdb[,.N,by=.(SEQUENCE=sub(" .*","",SEQUENCE))]; sc <- merge(sc,vc,by.x="accession",by.y="SEQUENCE",all.x=TRUE); setnames(sc,"N","n_vf"); sc[is.na(n_vf),n_vf:=0L] } else sc[,n_vf:=0L]
  if (!is.null(annot$plasmidfinder)) { rp <- annot$plasmidfinder[,.(replicon=GENE[1]),by=.(SEQUENCE=sub(" .*","",SEQUENCE))]; sc <- merge(sc,rp,by.x="accession",by.y="SEQUENCE",all.x=TRUE); sc[is.na(replicon),replicon:="Unknown"] } else sc[,replicon:="Unknown"]
  sc[,S_ARG:=pmin(0.25+pmin(n_arg*0.05,0.35)+0.20*(n_arg>0)+0.20*(n_highrisk>0),1.0)]
  sc[,S_VF:=0.0]; sc[n_vf>0,S_VF:=pmin(0.30+pmin(n_vf*0.03,0.40),1.0)]
  sc[,`:=`(S_MOB=0.30,S_HOST=0.50,S_REP=0.30,S_SIZE=1/(1+exp(-(length_bp/1000-30)/15)),S_BM=0.0,S_GEO=0.30,S_HAB=0.30,S_GROW=0.30)]
  sc[,S_final:=as.matrix(.SD)%*%w,.SDcols=cc]; merge(sc, ep[,.(accession,label,y)], by="accession")
}
wf <- file.path(tab_dir,"tab_weight_comparison.csv")
if (file.exists(wf)) { wc <- fread(wf); final_w <- wc$Final; names(final_w) <- wc$Component } else final_w <- c(S_ARG=0.2448,S_VF=0.1096,S_MOB=0.2041,S_HOST=0.0282,S_REP=0.0030,S_SIZE=0.1808,S_BM=0.2112,S_GEO=0.0015,S_HAB=0.0022,S_GROW=0.0147)
comp_cols <- names(final_w)
if (file.exists(fasta_path) && Sys.which("abricate")!="") {
  annot <- annotate_ext(fasta_path, ext_dir); ext_scores <- score_ext(annot, external_plasmids, comp_cols, final_w)
  fwrite(ext_scores, file.path(tab_dir,"tab_external_validation_scores.csv"))
  ext_roc <- roc(ext_scores$y, ext_scores$S_final, quiet=TRUE); ext_auc <- as.numeric(auc(ext_roc))
  cat(sprintf("  External AUC = %.3f\n", ext_auc))
} else { cat("  abricate/FASTA unavailable; skipping external scoring.\n"); ext_auc <- NA }

cat("\n=== B. ROC/PR benchmarking ===\n")
master <- fread(file.path(res_dir,"psc_master.tsv"), na.strings=c("\\N","","NA"))
risk <- fread(file.path(res_dir,"psc_risk_scores.tsv"), na.strings=c("\\N","","NA"))
d <- merge(master, risk, by="id", suffixes=c("",".r"))
d[,n_vf:=fifelse(is.na(n_vf),0L,n_vf)]; d[,n_metal:=fifelse(is.na(n_metal),0L,n_metal)]; d[,n_arg:=fifelse(is.na(n_arg),0L,n_arg)]
d[,S_VF:=0.0]; d[n_vf>0,S_VF:=pmin(0.30+pmin(n_vf*0.03,0.40)+0.15*grepl("Exotoxin",vf_category)+0.15*grepl("Effector delivery",vf_category),1.0)]
d[,has_mer:=as.integer(grepl("mer",gene_bacmet,ignore.case=TRUE))]; d[,has_qac:=as.integer(grepl("qac",gene_bacmet,ignore.case=TRUE))]; d[,has_ars:=as.integer(grepl("ars|cop|sil",gene_bacmet,ignore.case=TRUE))]
d[is.na(has_mer),has_mer:=0L]; d[is.na(has_qac),has_qac:=0L]; d[is.na(has_ars),has_ars:=0L]
d[,S_BM:=0.0]; d[n_metal>0,S_BM:=pmin(0.25+pmin(n_metal*0.04,0.35)+0.15*has_mer+0.15*has_qac+0.10*has_ars,1.0)]
for (cc in comp_cols) d[[cc]] <- fifelse(is.na(d[[cc]]),0.0,as.numeric(d[[cc]]))
d[,y_highrisk:=as.integer(n_high_risk_arg>0)]; d[,y_fusion:=as.integer(n_arg>0 & n_vf>0)]
dd <- d[!is.na(S_ARG)]; set.seed(42); ti <- sample(seq_len(nrow(dd)),min(100000,nrow(dd))); dt <- dd[ti]
Xc <- as.matrix(dt[,.SD,.SDcols=comp_cols])
models <- list("PlasRisk (10-dim)"=Xc%*%final_w, "PlasRisk (9-dim)"={w9<-final_w; w9["S_BM"]<-0; w9<-w9/sum(w9); Xc%*%w9}, "PIPdb ordinal"=dt$risk_score/max(dt$risk_score,na.rm=TRUE), "S_ARG only"=dt$S_ARG, "S_MOB only"=dt$S_MOB, "ARG count"=dt$n_arg/max(dt$n_arg,na.rm=TRUE))
rl <- list(); pl <- list(); ba <- data.table()
for (nm in names(models)) {
  sc <- as.numeric(models[[nm]]); sc[is.na(sc)] <- 0
  ro <- roc(dt$y_highrisk, sc, quiet=TRUE); a <- as.numeric(auc(ro))
  pr <- data.table(recall=ro$sensitivities, precision=ro$sensitivities/(ro$sensitivities+(1-ro$specificities)*sum(dt$y_highrisk==0)/sum(dt$y_highrisk==1))); pr[,Model:=nm]; pl[[nm]] <- pr
  rl[[nm]] <- data.table(FPR=1-ro$specificities, TPR=ro$sensitivities, Model=nm, AUC=a)
  ba <- rbind(ba, data.table(Model=nm, AUC_highrisk=round(a,4), AUC_fusion=round(as.numeric(auc(roc(dt$y_fusion,sc,quiet=TRUE))),4)))
}
print(ba); fwrite(ba, file.path(tab_dir,"tab_benchmark_auc.csv"))
roc_df <- rbindlist(rl); roc_df[,Label:=sprintf("%s (AUC=%.3f)",Model,AUC)]
p_roc <- ggplot(roc_df,aes(x=FPR,y=TPR,colour=Label))+geom_line(linewidth=0.7)+geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey60")+scale_colour_brewer(palette="Set1")+labs(title="ROC: high-risk ARG prediction",x="FPR",y="TPR")+theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"),legend.position=c(0.72,0.22),legend.title=element_blank(),legend.text=element_text(size=8))
pr_df <- rbindlist(pl)
p_pr <- ggplot(pr_df,aes(x=recall,y=precision,colour=Model))+geom_line(linewidth=0.7)+scale_colour_brewer(palette="Set1")+ylim(0,1)+labs(title="Precision-Recall",x="Recall",y="Precision")+theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"),legend.position="none")

cat("\n=== C. Tool comparison table ===\n")
tool_compare <- data.table(
  Feature=c("Input","Output","ARG","VF","Replicon","Mobility","Biocide/metal","Risk score","Dimensions","Continuous","Validated","Conda","CLI","Language","Scale","License"),
  PlasmidFinder=c("FASTA","Replicon","No","No","Yes","No","No","No","0","N/A","N/A","Yes","Yes","Perl/Python","PF DB","GPL"),
  MOB_suite=c("FASTA","Mobility+relaxase","No","No","Yes","Yes","No","No","0","N/A","N/A","Yes","Yes","Python","Internal","Apache"),
  RFPlasmid=c("FASTA","Plasmid/chrom","No","No","No","No","No","No","0","N/A","RF","Yes","Yes","Python/R","Train","GPL"),
  Plasmer=c("FASTA/reads","Plasmid prob","No","No","No","No","No","No","0","N/A","ML","Yes","Yes","Python","PLSDB","MIT"),
  PIPdb=c("Web","Ordinal 1-5","Yes","Yes","Yes","Yes","Yes","Yes(ordinal)","6(equal)","No","No","No","No","Web/PHP","792k PSCs","Web"),
  PlasRisk=c("FASTA","Score A-E","Yes","Yes","Yes","Yes","Yes","Yes(continuous)","10(weighted)","Yes[0,1]","RF+LASSO+LORO","Yes","Yes","Python","792k PSCs","MIT"))
fwrite(tool_compare, file.path(tab_dir,"tab_tool_comparison.csv"))

cat("\n=== E. Figures ===\n")
if (exists("ext_auc") && !is.na(ext_auc)) {
  p_ext <- ggplot(data.table(FPR=1-ext_roc$specificities,TPR=ext_roc$sensitivities),aes(x=FPR,y=TPR))+geom_line(linewidth=1,colour="#d94801")+geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey60")+annotate("text",x=0.6,y=0.2,label=sprintf("External AUC=%.3f",ext_auc),size=5,fontface="bold")+labs(title="External validation (NCBI)",x="FPR",y="TPR")+theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
} else p_ext <- ggplot()+annotate("text",x=0.5,y=0.5,label="External validation\n(requires abricate+internet)",size=5)+theme_void()+labs(title="External validation")+theme(plot.title=element_text(face="bold"))
comb <- (p_roc|p_pr)/p_ext+plot_annotation(title="PlasRisk benchmarking and external validation",theme=theme(plot.title=element_text(face="bold",size=14)))
ggsave(file.path(fig_dir,"fig25_benchmark_validation.png"),comb,width=14,height=10,dpi=300)
ggsave(file.path(fig_dir,"fig25_benchmark_validation.pdf"),comb,width=14,height=10)
cat("  saved fig25_benchmark_validation\nDone.\n")
