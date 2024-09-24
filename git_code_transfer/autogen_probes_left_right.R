library("tidyverse")
library("readr")
library("data.table")
library("stringr")


### ------------------------------------------------- ###
#
# This is a script to automate the capture c probe design search
# 
### ------------------------------------------------- ###


### temp debug params
# args <- list()
# args[[1]] <- "/nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/main_attempt_090824/5890_1"
# args[[2]] <- 2




args <- commandArgs(trailingOnly = TRUE)
print(paste("path is",args[[1]]))
print(paste("round is",args[[2]]))



path_current <- args[[1]]
setwd(dir = path_current)


threshold_for_repeats = 40



### find probes which are highly repetitive, pull coords
### to their left, and their right
### and redesign



foo <- fread("oligo_info.txt")

### extract if density score over 40

bar <- foo[foo$density_score>threshold_for_repeats]

### get unique region names e.g. some have both probes repetitive, some just one.

bar_names <- unique(bar$associations)

### get the left and right sequences

buzz <- foo[foo$associations%in%bar_names]

if (args[[2]] == 1) {
    
    
    #you need this section for the first round
    
    
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
    
    
    
}else{
    
    
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
               stop = single_bp_stop) %>% as.data.frame()
}


### export

number_of_guide_regions_left <- length(bar_names)

write(number_of_guide_regions_left, paste0("leftright_coords_", args[[2]] , ".txt"))

print(foobar)


write_tsv(foobar, paste0("leftright_coords_round", args[[2]] , ".bed"))

######### error debug 
# 
# #you need this section for the first round
# 
# 
# buzz_left <- buzz[side_of_fragment=="L"]
# buzz_right <- buzz[side_of_fragment=="R"]
# 
# 
# ### new column with single bp coords
# # why 71?
# # because each oligo is 70bp
# qux_left <- buzz_left %>% mutate(single_bp_stop = start - 71,
#                                  single_bp_start = start - 72,
#                                  new_name = paste0(str_remove(associations, ","), "_L"))
# qux_right <- buzz_right %>% mutate(single_bp_start = stop + 71,
#                                    single_bp_stop = stop + 72,
#                                    new_name = paste0(str_remove(associations, ","), "_R"))
# 
# 
# ### new bed file with new left and right coords
# 
# 
# foobar <- qux_left[, list(chr, single_bp_start, single_bp_stop, new_name)] %>% rbind(qux_right[, list(chr, single_bp_start, single_bp_stop, new_name)]) %>%
#     rename(start = single_bp_start,
#            stop = single_bp_stop)
