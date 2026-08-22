#!/usr/bin/env Rscript
# pipdb_04_figures.R — Publication figures for PIPdb plasmid risk & evolution
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(RColorBrewer); library(ggrepel)
})
res_dir <- if (length(commandArgs(trailingOnly=TRUE))>=1) commandArgs(trailingOnly=TRUE)[1] else "results"
fig_dir <- file.path(res_dir,"figures"); dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
theme_set(theme_bw(base_size=11))

d <- fread(file.path(res_dir,"psc_risk_scores.tsv"), select=c(
  "plasmid_acc","replicon_primary","pmlst","length_avg","species_name","genus_name",
  "phylum_name","gram_stain","host_class","host_human","n_country","country",
  "year_mid","single_year","n_arg","n_who_arg","n_high_risk_arg","n_vf",
  "mobility_class","has_integron","n_is","n_metal","is_density_per_kb",
  "annual_growth_rate","combined_risk_index",
  "S_ARG","S_MOB","S_HOST","S_REP","S_SIZE","S_GEO","S_HAB","S_GROW",
  "risk_score","Q_grade","Q_tree"))

# Fig 1: risk overview
p1 <- ggplot(d, aes(risk_score, fill=Q_grade)) +
  geom_histogram(bins=120, color=NA) +
  scale_fill_brewer(palette="RdYlBu_r", direction=-1) +
  labs(x="Plasmid risk score", y="Number of PSCs", fill="Grade") +
  theme(legend.position=c(0.99,0.99), legend.justification=c(1,1))
p2 <- d %>% count(Q_grade, mobility_class) %>% group_by(Q_grade) %>% mutate(frac=n/sum(n)) %>%
  ggplot(aes(Q_grade, frac, fill=mobility_class)) + geom_col() +
  scale_fill_brewer(palette="Set2") + labs(x="Q grade", y="Fraction", fill="Mobility")
p3 <- d %>% group_by(Q_grade) %>% summarise(across(c(S_ARG,S_MOB,S_HOST,S_REP,S_GEO,S_GROW), mean)) %>%
  pivot_longer(-Q_grade) %>%
  ggplot(aes(Q_grade, value, fill=name)) + geom_col(position="dodge") +
  scale_fill_brewer(palette="Dark2") + labs(x="Q grade", y="Mean sub-score", fill="Component")
(p1/p2/p3) + plot_annotation(title="PIPdb plasmid risk assessment")
ggsave(file.path(fig_dir,"Fig1_risk_overview.pdf"), width=10, height=13)

# Fig 2: replicon ranking
rep <- fread(file.path(res_dir,"replicon_risk.tsv"))
top <- rep %>% filter(n>=30) %>% arrange(desc(mean_score)) %>% head(30)
p4 <- ggplot(top, aes(reorder(replicon_primary, mean_score), mean_score, fill=pct_Q1)) +
  geom_col() + coord_flip() + scale_fill_viridis_c() +
  labs(x=NULL, y="Mean risk score", fill="% Q1")
p5 <- ggplot(top, aes(reorder(replicon_primary, mean_score), pct_conjugative, size=n, color=mean_ARG)) +
  geom_point() + coord_flip() + scale_color_viridis_c() +
  labs(x=NULL, y="% conjugative", size="N PSCs", color="Mean S_ARG")
(p4|p5) + plot_annotation(title="Highest-risk replicon types")
ggsave(file.path(fig_dir,"Fig2_replicon_ranking.pdf"), width=14, height=10)

# Fig 4: temporal trends
dt <- d[!is.na(year_mid) & single_year==TRUE]
dt[, decade:=floor(year_mid/10)*10]
trend <- dt[, .(n=.N, pct_Q1=mean(Q_grade=="Q1")*100,
  pct_ARG=mean(n_arg>0)*100,
  pct_conj=mean(mobility_class%in%c("conjugative_complete","conjugative_likely"))*100),
  by=decade][order(decade)]
p7 <- ggplot(trend[decade>=1980], aes(factor(decade), pct_ARG, group=1)) +
  geom_line(linewidth=1) + geom_point(size=2) + labs(x="Decade", y="% PSCs carrying ARGs")
p8 <- ggplot(trend[decade>=1980], aes(factor(decade), pct_Q1, group=1)) +
  geom_line(linewidth=1, color="firebrick") + geom_point(size=2) + labs(x="Decade", y="% Q1")
p9 <- ggplot(trend[decade>=1980], aes(factor(decade), n, group=1)) +
  geom_col(fill="steelblue") + scale_y_log10() + labs(x="Decade", y="N PSCs (log)")
(p7/p8/p9) + plot_annotation(title="Temporal trends")
ggsave(file.path(fig_dir,"Fig4_temporal_trends.pdf"), width=8, height=12)

cat("\nAll figures in", fig_dir, "\n")
