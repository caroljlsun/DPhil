### ------------------------ ###
# Procedure: Read 1st round of capsequm guides, rerun if density above 40
# Project: Capture C
# User: CZJS
# Date: Oct-2023
### ------------------------ ###

library("tidyverse")
library("readr")
library("data.table")
library("stringr")

path_current <- "~/doctoral_work/2023/oligo/"
setwd(dir = path_current)


### ------------------------ ###
#
#the L/R oligos for the position you originally selected if they both have density scores <40 ----
#
### ------------------------ ###


foo <- fread("attempt_manual_round_1/oligo_info.txt")


# which have perfect probes?

quack <- foo %>% 
    group_by(associations) %>% 
    filter(density_score<40) %>% 
    filter(n() != 1)%>% 
    mutate(region = str_remove_all(associations, c("_L|_R|,")))

#write_tsv(quack, "final_probes/perfect_1st_round.txt")

# which need to be reran again?

bar <- foo[foo$density_score>40] #the too dense ones

# get unique region names e.g. some have both probes repetitive, some just one

bar_names <- unique(bar$associations)

buzz <- foo[foo$associations%in%bar_names]

# grab coordinates depending on left or right

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

# export


#write_tsv(foobar, paste0( "attempt_manual_round_1/","problematic_region", "1" , ".bed"))


### ------------------------ ###
# the L/R oligos for the +1/-1 restriction fragment if those both have density scores <40----
### ------------------------ ###


foo_2 <- fread("attempt_manual_round_2/oligo_info.txt")


# which have perfect probes?

quack_2 <- foo_2 %>% 
    group_by(associations) %>% 
    filter(density_score<40) %>% 
    filter(n() != 1) %>% 
    mutate(region = str_remove_all(associations, c("_L|_R|,")))


#write_tsv(quack_2, "final_probes/perfect_2nd_round.txt")

### ------------------------ ###
#
# Either the L or the R oligo for the position you originally selected if only one of them has a density score < 40 ----
#
### ------------------------ ###


# only one good probe?

quam <- foo %>% 
    group_by(associations) %>% 
    filter(density_score<40) %>% 
    filter(n() == 1)%>% 
    mutate(region = str_remove_all(associations, c("_L|_R|,")))

#write_tsv(quam, "final_probes/semi_perfect_1st_round.txt")


### ------------------------ ###
#
# Either the L or the R oligo for the +1/-1 restriction fragment if only one of them has a density score <40 ----
#
### ------------------------ ###



# only one good probe?

quam_2 <- foo_2 %>% 
    group_by(associations) %>% 
    filter(density_score<40) %>% 
    filter(n() == 1)%>% 
    mutate(region = str_remove_all(associations, c("_L|_R|,")))

#write_tsv(quam_2, "final_probes/semi_perfect_2nd_round.txt")

### ------------------------ ###
#
# which ones have failed all of this? ----
#
### ------------------------ ###

bar_2 <- foo_2[foo_2$density_score>40] #the too dense ones

# remove the ones with one good probe in the og region


quam_names <- unique(quam$associations) %>%
    str_remove_all(c("_L|_R|,")) #remove _L or _R


bax <- bar_2 %>% filter(!str_detect(associations, paste(quam_names, collapse = "|")))



# remove the ones with one good probe in +/- region

quam_names_2 <- unique(quam_2$associations) %>%
    str_remove_all(c("_L|_R|,")) #remove _L or _R


bax_2 <- bar_2 %>% filter(!str_detect(associations, paste(quam_names_2, collapse = "|")))

# bc bax and bax 2 do not match, we're good?


### ------------------------ ###
#
# merge then export ----
#
### ------------------------ ###


beast <- bind_rows("round_1_both_probes_sub_40" = quack,
                   "round_2_both_probes_sub_40" = quack_2,
                   "round_1_one_probe_sub_40" = quam,
                   "round_2_one_probe_sub_40" = quam_2, .id = "type")


# remove non-uniques.... hopefully with priority to the top of the df?

table(beast$region)

beast_twos <- beast %>% 
    group_by(region) %>% 
    filter(n() == 2)

beast_multiples <- beast %>% 
    group_by(region) %>% 
    filter(n() > 2)

table(beast_multiples$region)

# keep the ones which are completed in round 2

choose_doubles <- beast_multiples %>% 
    filter(type == "round_2_both_probes_sub_40")

choose_doubles_names <- choose_doubles$region

# filter them from the beast

beast_multiples_2 <- beast_multiples %>% filter(!str_detect(region, paste(choose_doubles_names, collapse = "|")))

# remove any more replicates

table(beast_multiples_2$region)

choose_singles <- beast_multiples_2%>% 
    filter(type == "round_1_one_probe_sub_40")

### merge the chosens with the beast_twos


final_beast <- rbind(beast_twos, choose_doubles, choose_singles)

#export


#write_tsv(final_beast, "final_probes/capsequm_oligos.txt")


### ------------------------ ###
#
# rerun for these 2 problematic regions ----
#
### ------------------------ ###

# chr12_124020197_124021765
# chr2_85704617_85706382

#let's see if they exist already...


problematic <- c("chr12_124020197_124021765",
                 "chr2_85704617_85706382")

do_they_exist <- beast %>% filter(str_detect(region, paste(problematic, collapse = "|")))

#yes, but in round 1 only


# make new df
#rerun from the top of the code, not best

buzz <- do_they_exist %>% 
    select(-type, -region) %>% 
    ungroup() %>% 
    as.data.table()

#write into bed file



#write_tsv(foobar, paste0( "attempt_manual_round_1/","leftright_coords_round", "1" , ".bed"))


problematic_foo <- fread("attempt_manual_round_3/oligo_info.txt") %>%
    mutate(region = str_remove_all(associations, c("_L|_R|,")))

#write_tsv(problematic_foo, "final_probes/problematic_regions.txt")



######### fluff

# only one good probe?

quam <- foo %>% 
    group_by(associations) %>% 
    filter(density_score<40) %>% 
    filter(n() == 1)

# exclude these from the next list


quam_names <- unique(quam$associations) %>%
    str_remove_all(c("_L|_R")) #remove _L or _R

foo2 <- foo

foo2 <- foo2[!foo2$associations %in% quam_names]

