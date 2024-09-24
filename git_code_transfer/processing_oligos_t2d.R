### ------------------------ ###
# Procedure: Find the guides which are "perfect", reselect snps for regions that aren't, rinse repeat
# Project: Capture C
# User: CZJS
# Date: Aug-2024
### ------------------------ ###

library(tidyverse)
library(data.table)

setwd("/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_t2d/")

t2d_oligos_r1 <- fread("T2D_EUR_hg38_selected_oligos.bed") %>% 
    mutate(associations = str_remove(associations, ",")) %>% 
    dplyr::rename(new_name = unique_name)

t2d_oligos_f1 <- fread("leftright_coords_round1.bed") %>% 
    mutate(associations = str_remove_all(new_name, c("_L|_R|,")))

t2d_oligos_r2 <- fread("t2d_round2_oligo_info.txt") %>% 
    mutate(new_name = str_remove(associations, ","),
           associations = str_remove_all(associations, "_L|_R|,"))

t2d_oligos_f2 <- fread("leftright_coords_round2.bed") %>% 
    mutate(associations = str_remove(new_name, c("_L|_R|,")))



t2d_oligos_reselected_r1 <- fread("t2d_reselected_oligo_info.txt") %>% 
    mutate(associations = str_remove(associations, ",")) 

t2d_oligos_reselected_f1 <- fread("leftright_coords_reselected_round1.bed")%>%
    mutate(associations = str_remove_all(new_name, c("_L|_R|,")))

t2d_oligos_reselected_r2 <- fread("t2d_reselected_oligo_info_round2.txt")%>% 
    mutate(new_name = str_remove(associations, ","),
           associations = str_remove_all(associations, "_L|_R|,"))

t2d_oligos_reselected_f2 <- fread("leftright_coords_reselected_round2.bed") %>% 
    mutate(associations = str_remove(new_name, c("_L|_R|,")))


### !intersections

t2d_perfect_r2 <- t2d_oligos_r2[!t2d_oligos_r2$new_name %in% t2d_oligos_f2$associations,] # perfect in the second round, as they are not in the 3rd round

t2d_perfect_r2_pairs_only <- t2d_oligos_r2[!t2d_oligos_r2$new_name %in% t2d_oligos_f2$associations,] %>% 
    group_by(associations) %>% 
    filter(n() == 2) %>% 
    ungroup() # perfect in the second round, as they are not in the 3rd round, and only has a single pair

# of these ones, if the probes to the left and to the right are perfect, choose the lowest density score, if that's a tie, lowest gc%

t2d_choose <- t2d_perfect_r2 %>% 
    group_by(associations) %>% 
    filter(n() > 2) %>% 
    group_by(associations, new_name) %>% 
    mutate(average_density = mean(density_score),
           average_gc_score = mean(`GC%`)) %>% 
    ungroup() %>% 
    group_by(associations) %>% 
    filter(average_density == min(average_density)) %>% 
    filter(average_gc_score == min(average_gc_score)) %>% 
    ungroup() %>% 
    select(-average_density, -average_gc_score)


t2d_no_doubles_found <- t2d_oligos_r2[t2d_oligos_r2$new_name %in% t2d_oligos_f2$associations,] # or the probes overlap so can't be used....

t2d_perfect_r1 <- t2d_oligos_r1[!t2d_oligos_r1$associations %in% t2d_oligos_f1$associations] # perfect from the first round, because they don't appear in the 2nd round

t2d_redundant_probes <- t2d_perfect_r1[t2d_perfect_r1$associations %in% t2d_perfect_r2$associations] # 0, good

t2d_combined_perfect <- rbind(t2d_perfect_r1, t2d_perfect_r2_pairs_only, t2d_choose)

length(unique(t2d_combined_perfect$associations))

### write separately or...

### export the ones that need new snps 

### same but for the reselected snps


t2d_resel_perfect_r2 <- t2d_oligos_reselected_r2[!t2d_oligos_reselected_r2$new_name %in% t2d_oligos_reselected_f2$associations,] 

t2d_resel_perfect_r2_pairs_only <- t2d_oligos_reselected_r2[!t2d_oligos_reselected_r2$new_name %in% t2d_oligos_reselected_f2$associations,] %>% 
    group_by(associations) %>% 
    filter(n() == 2) %>% 
    ungroup()

# of these ones, if the probes to the left and to the right are perfect, choose the lowest density score, if that's a tie, lowest gc%

t2d_resel_choose <- t2d_resel_perfect_r2 %>% 
    group_by(associations) %>% 
    filter(n() > 2) %>% 
    group_by(associations, new_name) %>% 
    mutate(average_density = mean(density_score),
           average_gc_score = mean(`GC%`)) %>% 
    ungroup() %>% 
    group_by(associations) %>% 
    filter(average_density == min(average_density)) %>% 
    filter(average_gc_score == min(average_gc_score)) %>% 
    ungroup() %>% 
    select(-average_density, -average_gc_score)


t2d_resel_no_doubles_found <- t2d_oligos_reselected_r2[t2d_oligos_reselected_r2$new_name %in% t2d_oligos_reselected_f2$associations,] # or the probes overlap so can't be used....

t2d_resel_perfect_r1 <- t2d_oligos_reselected_r1[!t2d_oligos_reselected_r1$associations %in% t2d_oligos_reselected_f1$associations]

t2d_resel_perfect_r1 <- t2d_resel_perfect_r1 %>% 
    mutate(new_name = associations,
           associations = str_remove_all(associations, "_reselected_snps"))

t2d_resel_perfect_r2 <- t2d_resel_perfect_r2 %>% 
    mutate(associations = str_remove_all(associations, "_reselected_snps"))


t2d_resel_redundant_probes <- t2d_resel_perfect_r1[t2d_resel_perfect_r1$associations %in% t2d_resel_perfect_r2$associations] # 0, good

t2d_resel_combined_perfect <- rbind(t2d_resel_perfect_r1, t2d_resel_perfect_r2_pairs_only, t2d_resel_choose)

length(unique(t2d_resel_combined_perfect$associations))
# checks

hmm <- str_extract(t2d_resel_combined_perfect$associations, "cluster_\\d+")

why <- data.frame(cluster = c(hmm), stringsAsFactors = F)

multiple_times <- why %>%
    group_by(cluster) %>% 
    filter( n() != 2)


#how many overlap with the first set...

overlapping_probes <- t2d_combined_perfect[t2d_combined_perfect$associations %in% t2d_resel_combined_perfect$associations,] #0, good

final_t2d_set <- rbindlist(list("t2d_combined_perfect" = t2d_combined_perfect, "t2d_resel_combined_perfect" = t2d_resel_combined_perfect), idcol = "ID") %>% 
    select(-ID, ID)


### choose oligos if the have just one good one

# which regions are left?

# the full list of regions should be in t2d_oligos_r1

which_left <- t2d_oligos_r1[!t2d_oligos_r1$associations %in% final_t2d_set$associations]

# how many of these have a passable single probe?

which_pass_r1 <- which_left[which_left$density_score < 40 ]

# bet some of these are not unique, but these may be overlapping instead

length(unique.default(which_pass_r1$associations)) #hmm, not this time strangely enough

# how many are left that don't have a good single probe from r1

# not quite, because one probe of a pair could be good enough

which_fail_r1 <- which_left[!which_left$associations %in% which_pass_r1$associations]

# get the uniques

which_fail_r1_unique <- unique(which_fail_r1$associations) #3

# do they exist in round 2, and if yes, do they have a passable single probe?

which_in_r2 <- t2d_oligos_r2[t2d_oligos_r2$associations %in% which_fail_r1_unique]

which_pass_r2 <- which_in_r2[which_in_r2$density_score < 40]

length(unique(which_pass_r2$associations))

# which regions still don't have a probe?

which_fail_r2 <- unique(which_in_r2[!which_in_r2$associations %in% which_pass_r2$associations]$associations) #1 region

# do these failed regions have a single successful probe in the reselected snps? 
t2d_oligos_reselected_r1 <- t2d_oligos_reselected_r1 %>% 
    mutate(new_name = str_remove_all(associations, "_reselected_snps"))

# string extract the region numbers

which_resel_1 <- t2d_oligos_reselected_r1[t2d_oligos_reselected_r1$new_name %in% which_fail_r2] # doesn't work

which_fail_r2_region <-  str_extract(which_fail_r2, "cluster_\\d+")

which_resel_1 <- t2d_oligos_reselected_r1 %>% 
    filter(str_detect(associations, paste(which_fail_r2_region, collapse = "|")))

# do they have a low enough density

which_resel_pass_1 <- which_resel_1[which_resel_1$density_score < 40]

# which ones do not pass

which_resel_fail_1 <-  unique(which_resel_1[!which_resel_1$associations %in% which_resel_pass_1$associations]$associations) 

# do they exist in round 2

which_in_resel2 <- t2d_oligos_reselected_r2[t2d_oligos_reselected_r2$associations %in% which_resel_fail_1]

# do they have a low enough density

which_resel_pass_2 <- which_in_resel2[which_in_resel2$density_score < 40]


# which ones don't have a double or single/ failed both snps and both rounds

which_resel_fail_2 <- unique(which_in_resel2[!which_in_resel2$associations %in% which_resel_pass_2$associations]$associations) 


# compile all singlets and check they are unique

which_single_probes <- rbindlist(list("which_pass_r1" = which_pass_r1, "which_pass_r2" = which_pass_r2, "which_resel_pass_1" = which_resel_pass_1, "which_resel_pass_2" = which_resel_pass_2), idcol = "ID") %>% 
    select(-ID, ID)


final_double_single_list <- rbindlist(list("doubles" = final_t2d_set, "singles" = which_single_probes), idcol = "probe_type") %>% 
    select(-ID, ID)

### check these sizes, somethings wrong here

##### DOUBLES FIRST 

# get the cluster number for the doubles:

doubles_cluster_num <-  str_extract(final_t2d_set$associations, "cluster_\\d+")

singles_cluster_num <-  str_extract(which_single_probes$associations, "cluster_\\d+")

# are there any in the doubles which are not just 2 per region?

doubles <- data.frame(
    cluster = c(doubles_cluster_num),
    stringsAsFactors = FALSE
)

# Use dplyr and stringr to filter rows where col1 contains the target string exactly twice
only_twice <- doubles %>%
    group_by(cluster) %>% 
    filter( n() == 2)

multiple_times <- doubles %>%
    group_by(cluster) %>% 
    filter( n() != 2)

# multiple of 4. Theres a problem with the filtering of LRLR

# choose only the ones from doubles

# ok something weird happened in the reselection but...
# just based on the id key, choose the doubles from the first round


clean_final_round_selection <- final_t2d_set %>% 
    mutate(cluster = str_extract(associations, "cluster_\\d+")) %>% 
    group_by(cluster) %>% 
    filter(n() > 2) %>% 
    filter(ID == "t2d_combined_perfect")


clean_final_round <- final_t2d_set %>% 
    mutate(cluster = str_extract(associations, "cluster_\\d+")) %>% 
    group_by(cluster) %>% 
    filter(n() == 2)

clean_final_doubles <- rbindlist(list("clean_final_round" = clean_final_round, "clean_final_round_selection" = clean_final_round_selection))


##### SINGLES NEXT


# why my god... anyways let's just choose the doubles if they exist, and amongst the singles, the best scoring one...

# cleaning singles

which_single_probes_clean <- which_single_probes %>% 
    filter(`GC%`<= 60) %>% 
    filter(is.na(repeat_class)) %>% 
    mutate(cluster = str_extract(associations, "cluster_\\d+")) %>% 
    group_by(cluster) %>% 
    filter(density_score == min(density_score)) %>% 
    filter(`GC%` == min(`GC%`)) %>% 
    ungroup() 


# get the cluster number for the doubles:

doubles_cluster_num <-  str_extract(clean_final_doubles$associations, "cluster_\\d+")

singles_cluster_num <-  str_extract(which_single_probes_clean$associations, "cluster_\\d+")

# check overlap

overlap_cluster_num <- intersect(doubles_cluster_num, singles_cluster_num)

# there are some non-uniques from the singles:

unique_singles_cluster_num <- unique(singles_cluster_num)

unique_doubles_cluster_num <- unique(doubles_cluster_num)

# choosing doubles...

remaining_singles <- which_single_probes_clean %>%
    filter(!cluster %in% overlap_cluster_num)

##### FINAL SET

final_clean_nondup_oligos <- rbindlist(list("clean_final_doubles" = clean_final_doubles , "remaining_singles" = remaining_singles))

# litttle checks

length(unique(final_clean_nondup_oligos$associations)) # 385
length(unique(t2d_oligos_r1$associations)) # 386

# the one that is missing just can't be made, which is in "which_resel_fail_2"

#write_tsv(final_clean_nondup_oligos, "final_clean_nondup_oligos_t2d.txt")

t2d_final_set <- fread("/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_t2d/final_clean_nondup_oligos_t2d.txt")

whr_final_set <- fread("/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_whr/final_clean_nondup_oligos_whr.txt")

merged_oligos <- rbindlist(list("t2d_final_set" = t2d_final_set , "whr_final_set" = whr_final_set), id = "gwas") %>% 
    select(-gwas, gwas)

#write_tsv(merged_oligos, "/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/final_merged_oligos.txt")

#bedtools intersect -a final_clean_nondup_oligos_whr.bed -b final_clean_nondup_oligos_t2d.bed -wb > whr_t2d_overlap.bed

# read in the overlap file and remove from overall big file

setwd("/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps")


whr_t2d_overlap <- fread("whr_t2d_overlap.bed") %>% 
    mutate(overlap_grp = ceiling(row_number() / 2), 
           t2d_cluster = str_extract(V29,"cluster_\\d+"),
           whr_cluster = str_extract(V13,"cluster_\\d+"))

# extract their respective clusters

whr_overlap_cluster <- str_extract(whr_t2d_overlap$V13,"cluster_\\d+")
t2d_overlap_cluster <- str_extract(whr_t2d_overlap$V29,"cluster_\\d+")

# select their respective rows
whr_t2d_overlap_selected <- whr_t2d_overlap %>% 
    select(overlap_grp, t2d_cluster, whr_cluster)

t2d_overlap_rows <- t2d_final_set[t2d_final_set$cluster %in% t2d_overlap_cluster] %>% 
    left_join(whr_t2d_overlap_selected, by = c("cluster" = "t2d_cluster")) %>% 
    mutate(new_name = paste0(new_name, "_", side_of_fragment)) %>% 
    distinct(new_name, .keep_all = TRUE) %>% 
    select(-whr_cluster)


whr_overlap_rows <-  whr_final_set[whr_final_set$cluster %in% whr_overlap_cluster] %>% 
    left_join(whr_t2d_overlap_selected, by = c("cluster" = "whr_cluster")) %>% 
    mutate(new_name = paste0(new_name, "_", side_of_fragment)) %>% 
    distinct(new_name, .keep_all = TRUE) %>% 
    select(-t2d_cluster)


# compare by lowest density score and or GC%

whr_t2d_overlap_rows <- rbind(whr_overlap_rows, t2d_overlap_rows) %>% 
    group_by(overlap_grp, cluster) %>%
    mutate(average_density = mean(density_score),
           average_gc_score = mean(`GC%`)) %>% 
    ungroup() %>% 
    group_by(overlap_grp) %>% 
    slice_min(order_by = average_density, n = 2, with_ties = FALSE) %>% 
    slice_min(order_by = average_gc_score, n = 2, with_ties = FALSE) %>% 
    ungroup() %>% 
    select(-average_density, -average_gc_score, -overlap_grp)

    
# select rows that don't overlap

t2d_no_overlap_rows <- t2d_final_set[!t2d_final_set$cluster %in% t2d_overlap_cluster]

whr_no_overlap_rows <- whr_final_set[!whr_final_set$cluster %in% whr_overlap_cluster]

no_overlap_merge <- rbind(whr_no_overlap_rows, t2d_no_overlap_rows, whr_t2d_overlap_rows)

#write_tsv(no_overlap_merge, "/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/final_merged_oligos_non_overlapping.txt")
#write_tsv(no_overlap_merge, "/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/final_merged_oligos_non_overlapping.bed")

# final check in bedtools 

t2d_final_combo <- rbind(t2d_no_overlap_rows, whr_t2d_overlap_rows[whr_t2d_overlap_rows$ID ==  "t2d_combined_perfect",])

whr_final_combo <- rbind(whr_no_overlap_rows, whr_t2d_overlap_rows[whr_t2d_overlap_rows$ID ==  "whr_combined_perfect",])

# write_tsv(t2d_final_combo, "/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/t2d_no_overlap_rows.txt" ,col_names = T)
# write_tsv(t2d_final_combo, "/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/t2d_no_overlap_rows.bed" ,col_names = F)
# write_tsv(whr_final_combo, "/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/whr_no_overlap_rows.txt" ,col_names = T)
# write_tsv(whr_final_combo, "/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/whr_no_overlap_rows.bed" ,col_names = F)
