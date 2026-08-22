#!/usr/bin/env Rscript
# pipdb_12_temporal_signal.R — root-to-tip regression temporal signal check
suppressPackageStartupMessages({ library(ape); library(data.table); library(ggplot2) })
evo_dir <- commandArgs(trailingOnly=TRUE)[1]
if (is.na(evo_dir)) evo_dir <- "results/evolution"

check_one <- function(rep_dir) {
  rep <- basename(rep_dir)
  td <- file.path(rep_dir, "tip_dates.tsv")
  tree_candidates <- c(
    file.path(rep_dir,"tree","snippy",paste0(rep,"_ml.treefile")),
    file.path(rep_dir,"tree",paste0(rep,"_core.treefile")),
    file.path(rep_dir,"tree",paste0(rep,"_snp.treefile")))
  tree_file <- tree_candidates[file.exists(tree_candidates)][1]
  if (is.na(tree_file) || !file.exists(td)) return(NULL)
  tr <- read.tree(tree_file)
  dates <- fread(td)
  id_col <- if ("plasmid_seq_id" %in% names(dates)) "plasmid_seq_id" else "seq_id"
  dates <- dates[, .(id=get(id_col), year)]
  dates[, id:=gsub("\\.","_",id)]
  dates <- dates[!duplicated(id)]
  common <- intersect(tr$tip.label, dates$id)
  if (length(common) < 10) return(NULL)
  tr <- keep.tip(tr, common)
  d <- setNames(dates[id %in% common]$year, dates[id %in% common]$id)
  d <- d[tr$tip.label]
  rtt_tr <- rtt(tr, d, objective="rms")
  rtt_dist <- node.depth.edgelength(rtt_tr)[1:length(rtt_tr$tip.label)]
  fit <- lm(rtt_dist ~ d)
  sm <- summary(fit)
  data.table(replicon=rep, n_tips=length(common), year_min=min(d), year_max=max(d),
    R2=round(sm$r.squared,4), p_value=signif(pf(sm$fstatistic[1],sm$fstatistic[2],sm$fstatistic[3],lower.tail=FALSE),3),
    slope=signif(coef(fit)[2],4), tmcra_year=round(-coef(fit)[1]/coef(fit)[2],1),
    has_temporal_signal=sm$r.squared>0.1 & pf(sm$fstatistic[1],sm$fstatistic[2],sm$fstatistic[3],lower.tail=FALSE)<0.05)
}

dirs <- list.dirs(evo_dir, recursive=FALSE)
res <- rbindlist(lapply(dirs, check_one), fill=TRUE)
if (nrow(res)==0) { cat("No trees found.\n") } else {
  fwrite(res, file.path(evo_dir,"temporal_signal_summary.tsv"), sep="\t")
  print(res)
}
