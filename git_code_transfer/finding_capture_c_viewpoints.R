library("tidyverse")
library("readr")
library("data.table")
library("stringr")

path_current <- "~/doctoral_work/2023"
setwd(dir = path_current)

distances <- read_tsv("distances.bed", col_names = F)

# find guide regions that are below 500kbp from each other and shame them

sub_500 <- distances %>% 
    filter(X9 < 500000 & X9>0)
# looks abvout right


# the ones we will keep were manually selected after considering the h3k17ac etc. 
# documentation in powerpoint
# "which regions are >500 kbp from each other?"

regions_to_remove <- c("chr5_77144102_77145371", "chr5_56511757_56513278", 
                       "chr3_12452186_12453928", "chr3_12458357_12459859")


# port intersecting cres and remove regions within 500 kbp of each other

intersections <- read_tsv("cre_intersections.bed", col_names = F) %>% 
    filter(!X4 %in% regions_to_remove)

# separate the -1s

intersections_full <- intersections %>% 
    filter(X6 != -1)

no_cres <- intersections %>% 
    filter(X6 == -1)
# remove whichever ones are duplicates/ choose the most central one/ most overlap

## consider a version which only finds the most central of the overlap only

# those which are duplicates, choose the middle of the guide regions

unique_intersections <- intersections_full %>% 
    distinct(X4, .keep_all = T)

# middle of regions for main dataset

centres <- unique_intersections %>% mutate(reg_centre = ceiling(X6+((X7-X6)/2)) ) # corrected 23.8.23

# middle of regions for regions with no cres

centres_non_cres <- no_cres %>% mutate(reg_centre = ceiling(X2+((X3-X2)/2)) )

#format into perfect bed files

final_centres <- centres[c("X5", "reg_centre", "X4")] %>% 
    rename(chromosome = X5,
           name = X4,
           start = reg_centre) %>% 
    mutate(stop = start+1, .before = name)

final_centres_non_cres <-  centres_non_cres[c("X1", "reg_centre", "X4")] %>% 
    rename(chromosome = X1,
           name = X4,
           start = reg_centre) %>% 
    mutate(stop = start+1, .before = name)

merged_centres <- bind_rows(final_centres, final_centres_non_cres)

# export

#write_tsv(merged_centres, "single_bp_regions.bed")

### after running capsequm, one region could not be made bc of tiny dpnii fragment
### manually replacing chr2:106030681-106030682 with chr2:106030696-106030698

replaced_centres <- merged_centres %>% filter(!start == 106030681) %>% 
    add_row(chromosome = "chr2", start = 106030696 , stop = 106030697, name = "chr2_106029903_106031690
")


write_tsv(replaced_centres, "corrected_replaced_single_bp_regions.bed")

### issues with highly non-specific guides being made
### find probes which are highly repetitive, pull coords
### to their left, and their right
### and redesign



foo <- fread("oligo/successful_attempt/oligo_info.txt")

### extract if no. of alighnments greater than 10

bar <- foo[foo$total_number_of_alignments>10]

### get unique region names e.g. some have both probes repetitive, some just one.

bar_names <- unique(bar$associations)

### get the left and right sequences

buzz <- foo[foo$associations%in%bar_names]

### grab coordinates depending on left or right

buzz_left <- buzz[side_of_fragment=="L"]
buzz_right <- buzz[side_of_fragment=="R"]


### new column with single bp coords
# why 71? 
# because each oligo is 70bp
qux_left <- buzz_left %>% mutate(single_bp_stop = start - 71,
                                 single_bp_start = start - 72,
                                 new_name = paste0(str_remove(associations, ","), "_L"))
qux_right <- buzz_right %>% mutate(single_bp_start = stop + 71,
                                   single_bp_stop = stop + 72,
                                   new_name = paste0(str_remove(associations, ","), "_R"))


### new bed file with new left and right coords


foobar <- qux_left[, list(chr, single_bp_start, single_bp_stop, new_name)] %>% rbind(qux_right[, list(chr, single_bp_start, single_bp_stop, new_name)]) %>% 
    rename(start = single_bp_start,
           stop = single_bp_stop)

### export

write_tsv(foobar, "left_right_single_coords.bed")





### check the new left right guides for quality


foo <- fread("oligo/attempt_051023/19351_2/oligo_info.txt")

### extract if no. of alighnments greater than 10

bar <- foo[foo$density_score<40]

# which have perfect probes?
quack <- foo %>% 
    group_by(associations) %>% 
    filter(density_score<40) %>% 
    filter(n() != 1)

#write_tsv(quack, "perfect_coords_2nd_round.txt")

# which have one good probe?

quam <- foo %>% 
    group_by(associations) %>% 
    filter(density_score<40) %>% 
    filter(n() == 1)

#write_tsv(quam, "semiperfect_coords_2nd_round.txt")

# exclude these from the next list


quam_names <- unique(quam$associations) %>%
    str_remove_all(c("_L|_R")) #remove _L or _R

foo2 <- foo






### get unique region names e.g. some have both probes repetitive, some just one.

bar_names <- unique(bar$associations)

### get the left and right sequences
# these sequences need to be redesigned ahhhhhhhhhhhhhhhhhhhhhh

buzz <- foo[foo$associations%in%bar_names]

#40/56 are not useable, and need to be redone

### grab coordinates depending on left or right
# the left has to go more left, and the right more right

buzz_left <- buzz[side_of_fragment=="L" & str_detect(associations, "_L")]
buzz_right <- buzz[side_of_fragment=="R" & str_detect(associations, "_R")]


### new column with single bp coords
# why 71? 
# because each oligo is 70bp
qux_left <- buzz_left %>% mutate(single_bp_stop = start - 71,
                                 single_bp_start = start - 72,
                                 new_name = paste0(str_remove(associations, ","), "_L"))
qux_right <- buzz_right %>% mutate(single_bp_start = stop + 71,
                                   single_bp_stop = stop + 72,
                                   new_name = paste0(str_remove(associations, ","), "_R"))


### new bed file with new left and right coords


foobar <- qux_left[, list(chr, single_bp_start, single_bp_stop, new_name)] %>% rbind(qux_right[, list(chr, single_bp_start, single_bp_stop, new_name)]) %>% 
    rename(start = single_bp_start,
           stop = single_bp_stop)

### export

#write_tsv(foobar, "left_right_single_coords.bed")

#export the ones which passed only:

write_tsv(bar, "passed_coords_2nd_round.txt")
