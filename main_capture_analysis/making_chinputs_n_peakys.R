### ------------------------ ###
# Procedure: Make chinputs out of the processed hdf5s from capcruncher
# Project: Capture C
# User: CZJS
# Date: FEB-2025
### ------------------------ ###

### ------------------------ ###
# Description
# This file prcoesses the hdf5 cooler files from capcruncher into chinputs,
# which are used later to make peakys in peaky.R and chicanes in CHiCANE_test.R
# as the files are read and processed, preys within +/-2.5kbp of a blacklisted
# region are removed
# individual and passage level chinputs are made in this file
### ------------------------ ###

library("data.table")
library("tidyverse")
library("rhdf5")

setwd("comp_bio/capture_c_adipocytes/main_capture_analysis/")

### R versioning and global packages are not compatible SMH
# args <- commandArgs(trailingOnly = TRUE)
# print(paste("passage is",args[[1]]))

#passage <- args[[1]]
passage <- "DAAP6mADs_M3"

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
##################################-
### minimum example ----
##################################-

#/chr10_113021939_113021940_rs12258200:T:C_cluster_227_whradjbmi_combined_H3K27ac_ATAC/pixels

#qux <- h5read(file = hdf5_file_to_read, name = "/chr10_113021939_113021940_rs12258200:T:C_cluster_227_whradjbmi_combined_H3K27ac_ATAC/pixels")


#hmm <- h5read(file = hdf5_file_to_read, name = ids_to_read)
#fii <- h5read(file = hdf5_file_to_read, name = "/chr10_122098310_122100128/bins")


# just extracting one bait id! 

#cooler_section <- h5read(hdf5_file_to_read, name = master_fragment_key$name[3])


# cooler ids that match the baitid of interest, either in col1 or col2
# matching_indices_col1 <- which(cooler_section$bin1_id == master_fragment_key$V4.y[3])
# matching_indices_col2 <- which(cooler_section$bin2_id == master_fragment_key$V4.y[3])

# Extract the matching bin1_id values
# clean_cooler_section_bait <- c(cooler_section$bin1_id[matching_indices_col1], cooler_section$bin2_id[matching_indices_col2]) # should all be the same ID
# clean_cooler_section_oe <- c(cooler_section$bin2_id[matching_indices_col1], cooler_section$bin1_id[matching_indices_col2]) # should be the opposite ID
# clean_cooler_section_count <- c(cooler_section$count[matching_indices_col1], cooler_section$count[matching_indices_col2])
# 
# clean_cooler_section <- data.table(baitID = clean_cooler_section_bait, otherEndID = clean_cooler_section_oe, N =  clean_cooler_section_count)

### calculate the distsign
# 
# baitID_rfs_temp <- restriction_fragments[restriction_fragments$V4 %in% clean_cooler_section$baitID]
# otherEndID_rfs_temp <- restriction_fragments[restriction_fragments$V4 %in% clean_cooler_section$otherEndID]
# 
# condition <- baitID_rfs_temp$chr == otherEndID_rfs_temp$chr
# 
# otherEndID_rfs_temp$distSign <- ifelse(
#     condition,
#     ifelse(
#         otherEndID_rfs_temp$stop >= baitID_rfs_temp$start, # is the otherend further along than the bait
#         otherEndID_rfs_temp$start - baitID_rfs_temp$stop, #if yes subtract OE$start from Bait$stop
#         otherEndID_rfs_temp$stop - baitID_rfs_temp$start # if no, subtract OE$stop from Bait$start
#     ),
#     NA
# )

### calculate the OtherEndLength

# otherEndID_rfs_temp$otherEndLen <- otherEndID_rfs_temp$stop - otherEndID_rfs_temp$start
# 
# #baitID	otherEndID	N	otherEndLen	distSign
# check_chinput_form <- merge(clean_cooler_section, otherEndID_rfs_temp, by.x = "otherEndID", by.y = "V4", all = TRUE)
# 
# 
# final_chinput_form <- merge(clean_cooler_section, otherEndID_rfs_temp, by.x = "otherEndID", by.y = "V4", all = TRUE) %>% 
#     select(baitID, otherEndID, N, distSign, otherEndLen)
# 

#write_delim(final_chinput_form, "trial_run_chicago/testing_r_chicago/final.chinput")

##################################-
### extracting multiple bait ids ----
##################################-

chinput_compile <- list()

for (i in seq_along(master_fragment_key$name)) { #seq_along(master_fragment_key$name)
    
    print(paste0("Processing ", master_fragment_key$name[i]))
    
    cooler_section <- h5read(hdf5_file_to_read, name = master_fragment_key$name[i])
    
    
    # cooler ids that match the baitid of interest, either in col1 or col2
    matching_indices_col1 <- which(cooler_section$bin1_id == master_fragment_key$V4.y[i])
    matching_indices_col2 <- which(cooler_section$bin2_id == master_fragment_key$V4.y[i])
    
    # Extract the matching bin1_id values
    clean_cooler_section_bait <- c(cooler_section$bin1_id[matching_indices_col1], cooler_section$bin2_id[matching_indices_col2]) # should all be the same ID
    clean_cooler_section_oe <- c(cooler_section$bin2_id[matching_indices_col1], cooler_section$bin1_id[matching_indices_col2]) # should be the opposite ID
    clean_cooler_section_count <- c(cooler_section$count[matching_indices_col1], cooler_section$count[matching_indices_col2])
    
    clean_cooler_section <- data.table(baitID = clean_cooler_section_bait, otherEndID = clean_cooler_section_oe, N =  clean_cooler_section_count)
    
    ### calculate the distsign
    
    baitID_rfs_temp <- restriction_fragments[restriction_fragments$V4 %in% clean_cooler_section$baitID]
    otherEndID_rfs_temp <- restriction_fragments[restriction_fragments$V4 %in% clean_cooler_section$otherEndID]
    
    condition <- baitID_rfs_temp$chr == otherEndID_rfs_temp$chr
    
    otherEndID_rfs_temp$distSign <- ifelse(
        condition,
        ifelse(
            otherEndID_rfs_temp$stop >= baitID_rfs_temp$start, # is the otherend further along than the bait
            otherEndID_rfs_temp$start - baitID_rfs_temp$stop, #if yes subtract OE$start from Bait$stop
            otherEndID_rfs_temp$stop - baitID_rfs_temp$start # if no, subtract OE$stop from Bait$start
        ),
        NA
    )
    
    ### calculate the OtherEndLength
    
    otherEndID_rfs_temp$otherEndLen <- otherEndID_rfs_temp$stop - otherEndID_rfs_temp$start
    
    #baitID	otherEndID	N	otherEndLen	distSign
    check_chinput_form <- merge(clean_cooler_section, otherEndID_rfs_temp, by.x = "otherEndID", by.y = "V4", all = TRUE)
    
    chinput_compile[[i]] <- check_chinput_form
    
}


# final_chinput_form <- do.call(rbind, chinput_compile) %>%
#     select(baitID, otherEndID, N, distSign, otherEndLen)

final_chinput_form <- do.call(rbind, chinput_compile) %>% 
    relocate(c("chr", "start", "stop"))

# save then blacklist in terminal, then come back

write_tsv(final_chinput_form, paste0(passage, "_unproc.bed"), col_names = F)

#bedtools intersect -b ../ref_genomes/hg38-blacklist.v2_slop_5kbp.bed -a DAAP6pADs_P1_unproc.bed -v > blacklist_5kbp_P1.bed 


blacklisted <- fread("blacklist_5kbp_M3.bed") 

procd <- blacklisted %>% 
    dplyr::select(baitID = V5, otherEndID = V4, N = V6, distSign = V7, otherEndLen = V8)

save_file_name <- paste0(passage, "_processed_full.chinput")

print(paste0("saving as ", save_file_name))

write_delim(procd, save_file_name)

# rinse repeat with the other passages...

##################################-
### groupings ----
##################################-

# group individual samples by grouping baitID-preyID for all samples, then summing N.
# split by passage?

M1_chinput <- fread("DAAP6mADs_M1_processed_full.chinput")
M2_chinput <- fread("DAAP6mADs_M2_processed_full.chinput")
M3_chinput <- fread("DAAP6mADs_M3_processed_full.chinput")

all_mADs_chinput <- bind_rows(M1_chinput, M2_chinput, M3_chinput) %>% 
    group_by(baitID, otherEndID) %>% 
    summarise(N = sum(N, na.rm = TRUE), .groups = 'drop')


#write_delim(all_mADs_chinput, "all_mADs_blacklisted.chinput")

P1_chinput <- fread("DAAP6pADs_P1_processed_full.chinput")
P2_chinput <- fread("DAAP6pADs_P2_processed_full.chinput")
P3_chinput <- fread("DAAP6pADs_P3_processed_full.chinput")

all_pADs_chinput <- bind_rows(P1_chinput, P2_chinput, P3_chinput) %>% 
    group_by(baitID, otherEndID) %>% 
    summarise(N = sum(N, na.rm = TRUE), .groups = 'drop')


#write_delim(all_pADs_chinput, "all_pADs_blacklisted.chinput")
