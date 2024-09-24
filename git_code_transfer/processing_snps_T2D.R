library(tidyverse)
library(data.table)

setwd("/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps")

t2d_finemapping <- fread("T2D_EUR_Suzuki_finemapping.bed") 
t2d_snps <- fread("T2D_EUR_hg38_clustered_sorted_final_no_exons_no_TSS.bed")


# also consider that I only need this strategy for multiplet regions?

### isolate a dataset with multiplet regions 

singlets <- t2d_snps %>% 
    group_by(V17) %>%
    filter(n() == 1)


multiplets <- t2d_snps %>% 
    group_by(V17) %>%
    filter(n() != 1)


### then group by region, and select the top pip if it has one 
# so first intersect with finemapping data

all_snp_multiplets <- inner_join(t2d_finemapping, multiplets, by = join_by(V5 == V4))

#hmm <- multiplets[multiplets$V4 %in% t2d_finemapping$V5,]

top_snp_multiplets <- all_snp_multiplets %>%
    group_by(V17, V2.x, V3.x) %>%
    filter(V4 == max(abs(V4))) %>%
    filter(V7.y == min(V7.y)) %>% # there are sometimes multiple rows for the same region because of window/cs reasons
    ungroup() %>% 
    group_by(V17) %>% 
    filter(V4 == max(abs(V4))) %>%
    filter(V7.y == min(V7.y))
    
top_snp_multiplets_unified <- multiplets[multiplets$V4 %in% top_snp_multiplets$V5,]

# 83 top pips selected for 227 multiplet regions
# some of them are still targetting the same region so pick the highest of the region I guess

# which ones don't have a selection yet?

which_empty <- multiplets[!multiplets$V17 %in% top_snp_multiplets$V17,]

### if they don't have one, by the ones with eQTLs

# load in all intersected eqtl datasets

Adipose_Subcutaneous_snps_lifted_intersected_T2D_EUR <- fread("Adipose_Subcutaneous_snps_lifted_intersected_T2D_EUR.bed")
AdipoExpress_snps_lifted_intersected_T2D_EUR <- fread("AdipoExpress_snps_lifted_intersected_T2D_EUR.bed")
FUSION_ge_adipose_naive_snps_lifted_intersected_T2D_EUR <- fread("FUSION_ge_adipose_naive_snps_lifted_intersected_T2D_EUR.bed")
TwinsUK_ge_fat_snps_lifted_intersected_T2D_EUR <-  fread("TwinsUK_ge_fat_snps_lifted_intersected_T2D_EUR.bed")


# merge

eQTLs <- bind_rows(
    Adipose_Subcutaneous_snps_lifted_intersected_T2D_EUR = Adipose_Subcutaneous_snps_lifted_intersected_T2D_EUR,
    AdipoExpress_snps_lifted_intersected_T2D_EUR = AdipoExpress_snps_lifted_intersected_T2D_EUR,
    FUSION_ge_adipose_naive_snps_lifted_intersected_T2D_EUR = FUSION_ge_adipose_naive_snps_lifted_intersected_T2D_EUR,
    TwinsUK_ge_fat_snps_lifted_intersected_T2D_EUR = TwinsUK_ge_fat_snps_lifted_intersected_T2D_EUR,
    .id = "id"
)

eQTL_collapsed <- eQTLs %>% 
    group_by(V1, V2, V3) %>% 
    summarise(genes = paste0(V4, collapse = ","))


# intersect

eQTL_snp_multiplets <- inner_join(eQTL_collapsed, which_empty, by = join_by(V1, V2, V3))


eQTL_unified <- multiplets[multiplets$V4 %in% eQTL_snp_multiplets$V4,]

### if they don't have eQTLs, can't be within 2 kbp of a promoter

# so remove the ones with eqtl data 

which_empty_empty <- which_empty[!which_empty$V17 %in% eQTL_unified$V17,]


T2D_EUR_hg38_2kbp_TSS <- fread("T2D_EUR_hg38_2kbp_TSS.bed")


# get a common key
T2D_EUR_hg38_2kbp_TSS_rsids <- T2D_EUR_hg38_2kbp_TSS %>% 
    mutate(rsids = sub(":.*", "", V4))

# intersect again, but by rsids

which_far_enough <- which_empty_empty[which_empty_empty$V4 %in% T2D_EUR_hg38_2kbp_TSS_rsids$rsids,] # regions = 86

# rank these ones by beta effect

far_enough_beta_ranked <- which_far_enough %>% 
    group_by( V17 ) %>%
    filter(abs(V5) == max(abs(V5))) %>%
    filter(V7 == min(V7)) %>%
    ungroup()

### ok find the regions which don't have any snps that are far enough

which_regions_left <- which_empty_empty[!which_empty_empty$V17 %in% which_far_enough$V17,] # regions = 46

# of these regions, also rank by beta effect...

too_close_beta_ranked <- which_regions_left %>%
    group_by(V17) %>%
    filter(abs(V5) == max(abs(V5))) %>%
    filter(V7 == min(V7)) %>%
    ungroup()

length(unique(top_snp_multiplets$V17))

### think about whether you want an additional few columns for extra information e.g. genes if eqtl derived?

### anyways unify into one data table:

unified <- bind_rows(singlets = singlets,
                     top_snp_multiplets_unified = top_snp_multiplets_unified,
                     eQTL_unified = eQTL_unified,
                     far_enough_beta_ranked = far_enough_beta_ranked,
                     too_close_beta_ranked = too_close_beta_ranked, .id="id") 

unified <- unified %>% select(-id, id)

#write_tsv(unified, "T2D_EUR_hg38_selected_snps.bed")
#write_tsv(unified, "T2D_EUR_hg38_selected_snps_round2.bed")

# maybe write as a bed file...

### which oligos are these snps going to use?

all_oligos <- fread("/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/t2d_no_exons_no_TSS/oligo_info.txt")
# remove that awkward comma...

all_oligos_processed <- all_oligos %>% 
    mutate(unique_name = sub(",", "", associations))
# "intersect"

which_oligos <- all_oligos_processed[all_oligos_processed$unique_name %in% unified$V18,]

#write_tsv(which_oligos, "T2D_EUR_hg38_selected_oligos.bed")


### read snps which failed...

# first the key of their failure:

setwd("/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_t2d/")
t2d_pick_new_snps <- fread("leftright_coords_round2.bed")

# clean key?

t2d_pick_new_snps <- t2d_pick_new_snps %>% 
    mutate(associations = str_remove_all(new_name, "_L_L|_R_R|,"))

# which snp exactly was it?

exact_snp <- unified[unified$V18 %in% t2d_pick_new_snps$associations,]

# ok if it is a singleton, not much can be done for it

# if it is a multiplet, got through the above again 

regions_to_pick_again <- multiplets[multiplets$V17 %in% exact_snp$V17,]

# exclude the snps that have already been through the gauntlet

new_regions_to_pick_again <- regions_to_pick_again[!regions_to_pick_again$V18 %in% t2d_pick_new_snps$associations,]


### ok I'm about to do some real cursed looping

#t2d_snps <- new_regions_to_pick_again


# t2d_snps <- unified %>%
#     ungroup() %>%
#     select(V1, V2, V3, V18) %>%
#     dplyr::rename(chr = V1,
#                   start = V2,
#                   stop = V3,
#                   new_name = V18) %>%
#     mutate(new_name = paste0(new_name, "_reselected_snps"))


#write_tsv(t2d_snps, "t2d_hg38_selected_snps_round2_capready.bed")

