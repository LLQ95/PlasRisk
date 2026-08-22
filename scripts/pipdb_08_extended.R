#!/usr/bin/env Rscript
# pipdb_08_extended.R — Extended analyses: MDR-VF fusion, host range, mobility trends
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
res_dir <- commandArgs(trailingOnly=TRUE)[1]
if (is.na(res_dir)) res_dir <- "results"
fig_dir <- file.path(res_dir,"figures"); dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
tab_dir <- file.path(res_dir,"tables"); dir.create(tab_dir, showWarnings=FALSE, recursive=TRUE)

d <- fread(file.path(res_dir,"psc_risk_scores.tsv"))

# MDR-VF fusion
d[, has_ARG := n_arg > 0]
d[, has_VF := n_vf > 0]
d[, plasmid_type := ifelse(has_ARG & has_VF, "Resistance+Virulence (MDR-VF)",
                    ifelse(has_ARG, "Resistance only",
                    ifelse(has_VF, "Virulence only", "Neither")))]
fusion <- d[, .N, by=plasmid_type][, pct := N/sum(N)*100]
fwrite(fusion, file.path(tab_dir,"tab_mdr_vf_fusion.csv"))

# Mobility temporal trend
mob_trend <- d[!is.na(year_mid) & single_year==TRUE, .N, by=.(year_mid, mobility_class)]
mob_trend[, total := sum(N), by=year_mid]
mob_trend[, frac := N/total]
fwrite(mob_trend, file.path(tab_dir,"tab_mobility_temporal.csv"))

# Host range by replicon
host_range <- d[!is.na(replicon_primary), .(
  n = .N,
  n_species = uniqueN(species_name),
  n_country = mean(n_country, na.rm=TRUE),
  pct_human = mean(host_human, na.rm=TRUE)*100,
  pct_conj = mean(mobility_class %in% c("conjugative_complete","conjugative_likely"))*100
), by=replicon_primary][n>=50][order(-n)]
fwrite(host_range, file.path(tab_dir,"tab_host_range.csv"))

cat("Extended analysis tables saved.\n")
