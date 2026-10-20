library(tidyverse)
library(data.table)


setwd("comp_bio/capture_c_adipocytes/main_capture_analysis")


gene_bodies <- fread("comp_bio/capture_c_adipocytes/main_capture_analysis/M1_sig_mpcc_gene_bodies.bed")


unique_genes <- distinct(gene_bodies, V12, .keep_all = TRUE)

write_tsv(as.data.frame(unique_genes$V10), "enrichr.txt")

