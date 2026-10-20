### ------------------------ ###
# Procedure: Make chinputs out of the processed hdf5s from capcruncher
# Project: Capture C
# User: CZJS
# Date: FEB-2025
### ------------------------ ###


library("data.table")
library("tidyverse")
library("rhdf5")

setwd("comp_bio/capture_c_adipocytes/main_capture_analysis/")

### R versioning and global packages are not compatible SMH
# args <- commandArgs(trailingOnly = TRUE)
# print(paste("passage is",args[[1]]))

#passage <- args[[1]]
passage <- "DAAP6pADs_P1"

hdf5_file_to_read <- paste0("comp_bio/capture_c_adipocytes/main_capture_analysis/capcruncher_output/results/", passage, "/", passage, ".hdf5")

print(paste("reading ", hdf5_file_to_read))

### what's the structure of the cooler? 
#foo <- h5ls(hdf5_file_to_read)


### get the baitIDs chromosome positions
key_fragments <- fread("comp_bio/capture_c_adipocytes/main_capture_analysis/unique_restriction_fragments.bed") %>% 
    mutate(name = paste0("/", V4, "/pixels")) %>% 
    mutate(key = paste0(V1, V3))

corrected_key_fragments <-  key_fragments %>%
    mutate(chr = V1,
           start = V2+4,
           stop = V3-4)


# read in full restriction fragment rmap...
#restriction_fragments <- fread("comp_bio/capture_c_adipocytes/chicago/trial_run_chicago/testing_r_chicago/restriction_fragments.rmap") #this restriction fragment map is wrong!!!
restriction_fragments <- fread("genome.digest.bed") %>% # waittt....lol
    dplyr::rename(chr = V1, start = V2, stop = V3)

#there seems to be a consistent 4 base offset....
# which makes sense because dpnii is a 4 base cutter
# chr start stop
# chr10 113021010 113022756
# chr10 113021014 113022760
# 
# chr10 113021014 113022752
# chr10 113022756 113023534

corrected_key_fragments <-  key_fragments %>%
    mutate(chr = V1,
           start = V2+4,
           stop = V3-4)


master_fragment_key <- merge(corrected_key_fragments, restriction_fragments, by = c("chr", "start", "stop"))

#write_tsv(master_fragment_key, "master_fragment_key.bed")


##################################
### UMIs ----
##################################
M1 <- fread("../main_capture_analysis/DAAP6mADs_M2_full.chinput", nrows = 300000)


#restriction_fragments <- fread("genome.digest.bed") %>% # waittt....lol
#    dplyr::rename(chr = V1, start = V2, stop = V3)

M1_umi <- left_join(M1, restriction_fragments, by = c("baitID" = "V4")) %>% 
    dplyr::rename(chr_bait = chr, start_bait = start, end_bait = stop) %>% 
    left_join(., restriction_fragments, by = c("otherEndID" = "V4")) %>% 
    dplyr::rename(chr_contact = chr, start_contact = start, end_contact = stop) %>% 
    dplyr::select(chr_bait, start_bait, end_bait, chr_contact, start_contact, end_contact, UMIs = N, name = baitID, otherEndID) %>% 
    filter(chr_contact == chr_bait)

baits <- as.list(unique(M1_umi$name))



# for (i in 1:length(baits)) {
#     temp_UMI <- M1_umi[M1_umi$name == baits[i]] %>% dplyr::select(-name)
#     write_tsv(temp_UMI, paste0("umi_counts/", "M2_umi_", baits[i], "_counts.tsv"))
# }

## P1s
P1 <- fread("DAAP6pADs_P2_full.chinput", nrows = 300000)


#restriction_fragments <- fread("genome.digest.bed") %>% # waittt....lol
#    dplyr::rename(chr = V1, start = V2, stop = V3)

P1_umi <- left_join(P1, restriction_fragments, by = c("baitID" = "V4")) %>% 
    dplyr::rename(chr_bait = chr, start_bait = start, end_bait = stop) %>% 
    left_join(., restriction_fragments, by = c("otherEndID" = "V4")) %>% 
    dplyr::rename(chr_contact = chr, start_contact = start, end_contact = stop) %>% 
    dplyr::select(chr_bait, start_bait, end_bait, chr_contact, start_contact, end_contact, UMIs = N, name = baitID, otherEndID) %>% 
    filter(chr_contact == chr_bait)

baits <- as.list(unique(P1_umi$name))



for (i in 1:length(baits)) {
    temp_UMI <- P1_umi[P1_umi$name == baits[i]] %>% dplyr::select(-name)
    write_tsv(temp_UMI, paste0("umi_counts/", "P2_umi_", baits[i], "_counts.tsv"))
}

pADs_P1_43123 <- fread("umi_counts/P1_umi_43123_counts.tsv")
mADs_M1_43123 <- fread("umi_counts/M1_umi_43123_counts.tsv")
pADs_P2_43123 <- fread("umi_counts/P2_umi_43123_counts.tsv")
mADs_M2_43123 <- fread("umi_counts/M2_umi_43123_counts.tsv")
# Fill in missing rows for pADs

master_UMIs <- rbind(mADs_M1_43123, pADs_P1_43123, mADs_M2_43123, pADs_P2_43123) %>% 
    distinct(otherEndID, .keep_all = T)

# M1
missing_mADs_M1 <-setdiff(master_UMIs[,!"UMIs"], mADs_M1_43123[,!"UMIs"])
missing_mADs_M1$UMIs <- 0

M1_filled <- rbind(mADs_M1_43123, missing_mADs_M1)
write_tsv(M1_filled[,!"otherEndID"], "umi_counts/mADs_M1_filled.tsv")

# M2
missing_mADs_M2 <-setdiff(master_UMIs[,!"UMIs"], mADs_M2_43123[,!"UMIs"])
missing_mADs_M2$UMIs <- 0

M2_filled <- rbind(mADs_M2_43123, missing_mADs_M2)
write_tsv(M2_filled[,!"otherEndID"], "umi_counts/mADs_M2_filled.tsv")

# P1
missing_pADs_P1 <-setdiff(master_UMIs[,!"UMIs"], pADs_P1_43123[,!"UMIs"])
missing_pADs_P1$UMIs <- 0

P1_filled <- rbind(pADs_P1_43123, missing_pADs_P1)
write_tsv(P1_filled[,!"otherEndID"], "umi_counts/pADs_P1_filled.tsv")

# P2
missing_pADs_P2 <-setdiff(master_UMIs[,!"UMIs"], pADs_P2_43123[,!"UMIs"])
missing_pADs_P2$UMIs <- 0

P2_filled <- rbind(pADs_P2_43123, missing_pADs_P2)
write_tsv(P2_filled[,!"otherEndID"], "umi_counts/pADs_P2_filled.tsv")

#(head -n 1 mADs_M1_filled.tsv && tail -n +2 mADs_M1_filled.tsv | sort -t$'\t' -k4,4 -k5,5n -k6,6n) > mADs_M1_filled_sorted.tsv