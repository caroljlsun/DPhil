### ------------------------ ###
# Procedure: Find the guides which are "perfect", reselect snps for regions that aren't, rinse repeat
# Project: Capture C
# User: CZJS
# Date: Aug-2024
### ------------------------ ###

library(tidyverse)
library(data.table)

setwd("/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_whr//")

whr_oligos_r1 <- fread("whradjbmi_combined_hg38_selected_oligos.bed") %>% 
    mutate(associations = str_remove(associations, ",")) %>% 
    dplyr::rename(new_name = unique_name)
whr_oligos_f1 <- fread("leftright_coords_round1.bed") %>% 
    mutate(associations = str_remove_all(new_name, c("_L|_R|,")))

whr_oligos_r2 <- fread("whradjbmi_round2_oligo_info.txt") %>% 
    mutate(new_name = str_remove(associations, ","),
           associations = str_remove_all(associations, "_L|_R|,"))

whr_oligos_f2 <- fread("leftright_coords_round2.bed") %>% 
    mutate(associations = str_remove(new_name, c("_L|_R|,")))



whr_oligos_reselected_r1 <- fread("whradjbmi_reselected_oligo_info.txt") %>% 
    mutate(associations = str_remove(associations, ",")) 

whr_oligos_reselected_f1 <- fread("leftright_coords_reselected_round1.bed")%>%
    mutate(associations = str_remove_all(new_name, c("_L|_R|,")))

whr_oligos_reselected_r2 <- fread("whradjbmi_reselected_oligo_info_round2.txt")%>% 
    mutate(new_name = str_remove(associations, ","),
    associations = str_remove_all(associations, "_L|_R|,"))

whr_oligos_reselected_f2 <- fread("leftright_coords_reselected_round2.bed") %>% 
    mutate(associations = str_remove(new_name, c("_L|_R|,")))


### !intersections

whr_perfect_r2 <- whr_oligos_r2[!whr_oligos_r2$new_name %in% whr_oligos_f2$associations,] # perfect in the second round, as they are not in the 3rd round

whr_perfect_r2_pairs_only <- whr_oligos_r2[!whr_oligos_r2$new_name %in% whr_oligos_f2$associations,] %>% 
    group_by(associations) %>% 
    filter(n() == 2) %>% 
    ungroup() # perfect in the second round, as they are not in the 3rd round, and only has a single pair

# of these ones, if the probes to the left and to the right are perfect, choose the lowest density score, if that's a tie, lowest gc%

whr_choose <- whr_perfect_r2 %>% 
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


whr_no_doubles_found <- whr_oligos_r2[whr_oligos_r2$new_name %in% whr_oligos_f2$associations,] # or the probes overlap so can't be used....

whr_perfect_r1 <- whr_oligos_r1[!whr_oligos_r1$associations %in% whr_oligos_f1$associations] # perfect from the first round, because they don't appear in the 2nd round

whr_redundant_probes <- whr_perfect_r1[whr_perfect_r1$associations %in% whr_perfect_r2$associations] # 0, good

whr_combined_perfect <- rbind(whr_perfect_r1, whr_perfect_r2_pairs_only, whr_choose)

length(unique(whr_combined_perfect$associations))

### write separately or...

### export the ones that need new snps 

### same but for the reselected snps


whr_resel_perfect_r2 <- whr_oligos_reselected_r2[!whr_oligos_reselected_r2$new_name %in% whr_oligos_reselected_f2$associations,] 

whr_resel_perfect_r2_pairs_only <- whr_oligos_reselected_r2[!whr_oligos_reselected_r2$new_name %in% whr_oligos_reselected_f2$associations,] %>% 
    group_by(associations) %>% 
    filter(n() == 2) %>% 
    ungroup()

# of these ones, if the probes to the left and to the right are perfect, choose the lowest density score, if that's a tie, lowest gc%

whr_resel_choose <- whr_resel_perfect_r2 %>% 
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


whr_resel_no_doubles_found <- whr_oligos_reselected_r2[whr_oligos_reselected_r2$new_name %in% whr_oligos_reselected_f2$associations,] # or the probes overlap so can't be used....

whr_resel_perfect_r1 <- whr_oligos_reselected_r1[!whr_oligos_reselected_r1$associations %in% whr_oligos_reselected_f1$associations]

whr_resel_perfect_r1 <- whr_resel_perfect_r1 %>% 
    mutate(new_name = associations,
    associations = str_remove_all(associations, "_reselected_snps"))

whr_resel_perfect_r2 <- whr_resel_perfect_r2 %>% 
    mutate(associations = str_remove_all(associations, "_reselected_snps"))


whr_resel_redundant_probes <- whr_resel_perfect_r1[whr_resel_perfect_r1$associations %in% whr_resel_perfect_r2$associations] # 0, good

whr_resel_combined_perfect <- rbind(whr_resel_perfect_r1, whr_resel_perfect_r2_pairs_only, whr_resel_choose)

length(unique(whr_resel_combined_perfect$associations))
# checks

hmm <- str_extract(whr_resel_combined_perfect$associations, "cluster_\\d+")

why <- data.frame(cluster = c(hmm), stringsAsFactors = F)

multiple_times <- why %>%
    group_by(cluster) %>% 
    filter( n() != 2)


#how many overlap with the first set...

overlapping_probes <- whr_combined_perfect[whr_combined_perfect$associations %in% whr_resel_combined_perfect$associations,] #0, good

final_whr_set <- rbindlist(list("whr_combined_perfect" = whr_combined_perfect, "whr_resel_combined_perfect" = whr_resel_combined_perfect), idcol = "ID") %>% 
    select(-ID, ID)


### choose oligos if the have just one good one

# which regions are left?

# the full list of regions should be in whr_oligos_r1

which_left <- whr_oligos_r1[!whr_oligos_r1$associations %in% final_whr_set$associations]

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

which_in_r2 <- whr_oligos_r2[whr_oligos_r2$associations %in% which_fail_r1_unique]

which_pass_r2 <- which_in_r2[which_in_r2$density_score < 40]

length(unique(which_pass_r2$associations))

# which regions still don't have a probe?

which_fail_r2 <- unique(which_in_r2[!which_in_r2$associations %in% which_pass_r2$associations]$associations) #1 region

# do these failed regions have a single successful probe in the reselected snps? 
whr_oligos_reselected_r1 <- whr_oligos_reselected_r1 %>% 
    mutate(new_name = str_remove_all(associations, "_reselected_snps"))

# string extract the region numbers

which_resel_1 <- whr_oligos_reselected_r1[whr_oligos_reselected_r1$new_name %in% which_fail_r2] # doesn't work

which_fail_r2_region <-  str_extract(which_fail_r2, "cluster_\\d+")

which_resel_1 <- whr_oligos_reselected_r1 %>% 
    filter(str_detect(associations, paste(which_fail_r2_region, collapse = "|" )))

# do they have a low enough density

which_resel_pass_1 <- which_resel_1[which_resel_1$density_score < 40]

# which ones do not pass

which_resel_fail_1 <-  unique(which_resel_1[!which_resel_1$associations %in% which_resel_pass_1$associations]$associations) 

# do they exist in round 2

which_in_resel2 <- whr_oligos_reselected_r2[whr_oligos_reselected_r2$associations %in% which_resel_fail_1]

# do they have a low enough density

which_resel_pass_2 <- which_in_resel2[which_in_resel2$density_score < 40]


# which ones don't have a double or single/ failed both snps and both rounds

which_resel_fail_2 <- unique(which_in_resel2[!which_in_resel2$associations %in% which_resel_pass_2$associations]$associations) 


# compile all singlets and check they are unique

which_single_probes <- rbindlist(list("which_pass_r1" = which_pass_r1, "which_pass_r2" = which_pass_r2, "which_resel_pass_1" = which_resel_pass_1, "which_resel_pass_2" = which_resel_pass_2), idcol = "ID") %>% 
    select(-ID, ID)


final_double_single_list <- rbindlist(list("doubles" = final_whr_set, "singles" = which_single_probes), idcol = "probe_type") %>% 
    select(-ID, ID)

### check these sizes, somethings wrong here

##### DOUBLES FIRST 

doubles_cluster_num <-  str_extract(final_whr_set$associations, "cluster_\\d+")

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


clean_final_round_selection <- final_whr_set %>% 
    mutate(cluster = str_extract(associations, "cluster_\\d+")) %>% 
    group_by(cluster) %>% 
    filter(n() > 2) %>% 
    filter(ID == "whr_combined_perfect")


clean_final_round <- final_whr_set %>% 
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

length(unique(final_clean_nondup_oligos$associations)) # 211
length(unique(whr_oligos_r1$associations)) # 198

# 13 regions could not be designed for

# the one that is missing just can't be made, which is in "which_resel_fail_2"

#write_tsv(final_clean_nondup_oligos, "final_clean_nondup_oligos_whr.txt")

