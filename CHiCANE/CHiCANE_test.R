### ------------------------ ###
# Procedure: Test CHiCANE
# Project: Capture C
# User: CZJS
# Date: APR-2025
### ------------------------ ###


### ------------------------ ###
# Description
# Takes chinputs made in making_chinputs...R and creates chicanes and saves them
# before replicate chicanes are compared and tested, then saved
### ------------------------ ###

# Setup -------------------------------------------------------------------

#r has to be 4.4

library(chicane)
library(tidyverse)
library(data.table)
library(purrr)

setwd("comp_bio/capture_c_adipocytes/CHiCANE")

# test with package data --------------------------------------------------------------------

# start from pre-processed data

data(bre80)
results <- chicane(interactions = bre80)



hmmmmm <- get.trans.counts(bre80[, .(bait.chr, target.chr, bait.id, target.id, count)])
whatisthis <- bre80[, .(bait.chr, target.chr, bait.id, target.id, count)]



# test with own data --------------------------------------------------------------------


restriction_fragments <- fread("../main_capture_analysis/genome.digest.bed") %>% # waittt....lol
    dplyr::rename(chr = V1, start = V2, stop = V3)


M1 <- fread("../main_capture_analysis/DAAP6mADs_M1_full.chinput")


#restriction_fragments <- fread("genome.digest.bed") %>% # waittt....lol
#    dplyr::rename(chr = V1, start = V2, stop = V3)

M1_chicane <- left_join(M1, restriction_fragments, by = c("baitID" = "V4")) %>% 
    dplyr::rename(bait.chr = chr, bait.start = start, bait.end = stop) %>% 
    left_join(., restriction_fragments, by = c("otherEndID" = "V4")) %>% 
    dplyr::rename(target.chr = chr, target.start = start, target.end = stop) %>% 
    dplyr::mutate(bait.to.bait = ifelse(baitID == otherEndID,TRUE, FALSE)) %>% 
    dplyr::select(target.id = otherEndID, bait.id = baitID, target.chr, target.start, target.end, bait.chr, bait.start, bait.end, bait.to.bait, count = N) %>% 
    dplyr::filter(!(target.chr == "chr2" & target.start == 32915071 & target.end == 32918563))
# filter out chr2 32915071 32918563


M1_chicane <- add.covariates(M1_chicane) %>% 
    filter(distance < 1e6) # first run: 5e6, second run 1.25e5

M1_results <- chicane(interactions = M1_chicane, cores = 2, verbose = T)

# what is a good q/p value to filter by? 
# paper mentions q value < 0.05
table(M1_results$q.value < 0.05)

M1_results_sig <- M1_results[M1_results$q.value < 0.05]

write_delim(M1_results_sig, "M1_results_005_peaky_distances.txt")

# reorder for bedtools intersect

M1_results_sig_bedtools <- M1_results_sig %>% 
    dplyr::select(-target.id, -bait.id) %>% 
    filter(distance < 1e6) # filter by distance because there are some truly rogue ones here

write_tsv(M1_results_sig_bedtools, "M1_results_005_bedtools_peaky_distances.txt", col_names = F)


# Check distributions -------------------------------------------------------

m1_chicane <- fread("M1_results_005.txt")

hist(m1_chicane$distance)

# let's have a look at a few of the interactions which are massively far away

far_away <- m1_chicane[which(m1_chicane$distance > 1e6)]

# chr6     27777842   27777935 
# chr6     27902513   27903678
# chr6     27592841   27593503
# chr17     64741712   64742187

m1_chicane[which(m1_chicane$distance > 1e6)]

peakys <- fread("../peaky/M1_all/peaky_output.csv")

hist(abs(peakys$dist[peakys$rjmcmc_pos>=0.1]))



# replicates --------------------------------------------------------------



# manually changed the following to M1, M2, M£ then P1, P2, P3

{replicate <- fread("../main_capture_analysis/DAAP6pADs_P3_processed_full.chinput")


#restriction_fragments <- fread("genome.digest.bed") %>% # waittt....lol
#    dplyr::rename(chr = V1, start = V2, stop = V3)

rep_chicane <- left_join(replicate, restriction_fragments, by = c("baitID" = "V4")) %>% 
    dplyr::rename(bait.chr = chr, bait.start = start, bait.end = stop) %>% 
    left_join(., restriction_fragments, by = c("otherEndID" = "V4")) %>% 
    dplyr::rename(target.chr = chr, target.start = start, target.end = stop) %>% 
    dplyr::mutate(bait.to.bait = ifelse(baitID == otherEndID,TRUE, FALSE)) %>% 
    dplyr::select(target.id = otherEndID, bait.id = baitID, target.chr, target.start, target.end, bait.chr, bait.start, bait.end, bait.to.bait, count = N)# %>% 
#    dplyr::filter(!(target.chr == "chr2" & target.start == 32915071 & target.end == 32918563))
# filter out chr2 32915071 32918563


rep_chicane <- add.covariates(rep_chicane) #%>% 
#    filter(distance < 1e6) # first run: 5e6, second run 1.25e5, 3rd run 1e6
# actually, filter later in case I want to look at the ultra long range regions

rm(replicate)

write_tsv(rep_chicane, "P3_chicane.tsv")}

### compare replicates

#compare.replicates(c("M1_chicane.tsv", "M2_chicane.tsv", "M3_chicane.tsv"), output.directory = "comparisons")

M1_chicane <- fread("P1_chicane.tsv")# %>%  dplyr::filter(distance < 1e6) 
M2_chicane <- fread("P2_chicane.tsv")# %>%  dplyr::filter(distance < 1e6) 
M3_chicane <- fread("P3_chicane.tsv")# %>%  dplyr::filter(distance < 1e6) 

combined_chicane <- combine.replicates(list(M1_chicane, M2_chicane, M3_chicane), method = "sum")

combined_chicane_covar <- add.covariates(combined_chicane) %>% 
    dplyr::filter(distance < 1e6) # first run: 5e6, second run 1.25e5, 3rd run 1e6

# 
# combined_results_tnb <- chicane(interactions = combined_chicane_covar, 
#         verbose = TRUE,
#         cores = 6, distribution = "truncated-negative-binomial")

combined_results <- chicane(interactions = combined_chicane_covar, 
                                verbose = TRUE,
                                cores = 6)



# remotes::install_github("https://github.com/r-forge/countreg", subdir="pkg")
# install.packages("topmodels", repos = "https://R-Forge.R-project.org")



# quick sanity checks
hist(combined_results$distance)

hist(combined_results$q.value)
table(combined_results$q.value < 0.05) # 5753 true! # 8106 true when distance filtering before chicane anything... seems to eliminate nas? 

# saving results
combined_results_sig <- combined_results[combined_results$q.value < 0.05]

write_delim(combined_results_sig, "combined_results_P_blacklisted_only_final.txt")

#combined_results_sig <- fread("../CHiCANE/combined_results_M_blacklisted_only_final.txt")
# reorder for bedtools intersect

combined_results_bedtools <- combined_results_sig %>% 
    dplyr::select(-target.id, -bait.id) %>% 
    filter(distance < 1e6) # filter by distance because there are some truly rogue ones here

#write_tsv(combined_results_bedtools, "combined_results_005_M_blacklisted_only_final_fulltable.bed", col_names = F)

#write_tsv(combined_results_bedtools %>% dplyr::select(target.chr, target.start, target.end), "combined_results_005_P_blacklisted_only_final.bed", col_names = F)

# blacklisted_regions <- fread("../ref_genomes/hg38-blacklist.v2_slop_5kbp.bed") %>% 
#     dplyr::select(-V4)
# write_tsv(blacklisted_regions, "../ref_genomes/blacklist_ucsc.bed", col_names = F)

# plots --------------------------------------------------------------


# Create data frame with the distances and counts
distances <- c(100, 200, 500, 1000, 2000, 5000, 10000, 20000, 30000, 40000, 50000, 80000, 100000)

# Read counts from the files
counts <- c(3869, 3562, 2984, 2435, 1825, 1149, 787, 548, 442, 386, 348, 285, 258)

# Create data frame
data <- data.frame(distance = distances, count = counts)

# Plot
library(ggplot2)
ggplot(data, aes(x = distance, y = count)) +
    geom_line() +
    geom_point() +
    scale_x_log10(breaks = distances, labels = distances) +
    labs(x = "Merge Distance (bp)", y = "Number of Merged Regions", 
         title = "Effect of Merge Distance on Number of Regions") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("merge_distance_plot.png", width = 10, height = 6)

# wrangle clustered chicane --------------------------------------------------------------

M_10000 <- fread("merged_chicane_M_10000_fulltable.bed")
P_10000 <- fread("merged_chicane_P_10000_fulltable.bed")
