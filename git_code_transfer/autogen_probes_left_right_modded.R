library("tidyverse")
library("readr")
library("data.table")
library("stringr")
library("GenomicRanges")



### ------------------------------------------------- ###
#
# This is a script to automate the capture c probe design search
# 
### ------------------------------------------------- ###

### temp debug params
args <- list()
args[[1]] <- "/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_t2d/"
args[[2]] <- 2

# args <- commandArgs(trailingOnly = TRUE)
# print(paste("path is",args[[1]]))
# print(paste("round is",args[[2]]))

path_current <- args[[1]]
setwd(dir = path_current)

threshold_for_repeats = 40
GC_threshold <- 60


### find probes which are highly repetitive, pull coords
### to their left, and their right
### and redesign
### also remove any that have repeats and have a GC content over 60%

foo <- fread("t2d_reselected_oligo_info_round2.txt")

### extract if density score over 40
bar <- foo[foo$density_score > threshold_for_repeats] %>%
    rbind(foo[foo$`GC%` >= GC_threshold]) %>% 
    rbind(foo[!is.na(repeat_class)])

### get unique region names e.g. some have both probes repetitive, some just one.
bar_names <- unique(bar$associations)

### get the left and right sequences
buzz <- foo[foo$associations %in% bar_names]

### Detect overlapping probes using GenomicRanges

# Convert foobar to a GRanges object
foo_gr <- GRanges(seqnames = foo$chr,
                  ranges = IRanges(start = foo$start, end = foo$stop),
                  new_name = foo$associations)

# Find overlapping regions
overlaps <- findOverlaps(foo_gr, foo_gr) 

# Extract overlapping probes
overlapping_probes <- data.frame(
    chr = seqnames(foo_gr)[queryHits(overlaps)],
    start = start(foo_gr)[queryHits(overlaps)],
    stop = end(foo_gr)[queryHits(overlaps)],
    new_name = mcols(foo_gr)$new_name[queryHits(overlaps)],
    chr2 = seqnames(foo_gr)[subjectHits(overlaps)],
    start2 = start(foo_gr)[subjectHits(overlaps)],
    stop2 = end(foo_gr)[subjectHits(overlaps)],
    new_name2 = mcols(foo_gr)$new_name[subjectHits(overlaps)]
)

# Filter out self-overlaps
#overlapping_probes <- overlapping_probes[overlapping_probes$new_name != overlapping_probes$new_name2, ]

overlapping_probes <- overlapping_probes[
    !(overlapping_probes$chr == overlapping_probes$chr2 &
          overlapping_probes$start == overlapping_probes$start2 &
          overlapping_probes$stop == overlapping_probes$stop2), 
]

# fill info
new_foo <- foo[foo$associations %in% overlapping_probes$new_name,] #selecting rows with overlapping probes

#add only uniques

new_foo_names <- unique(new_foo$associations)


# check for duplicates with density filtered oligos

overlapping_buzz <- buzz[buzz$associations %in% new_foo_names,]

# final df of non duplicate buzz and new foo

non_overlap_foo <- buzz[!buzz$associations %in% overlapping_buzz$associations, ] %>% 
    rbind.data.frame(new_foo)

    

if (args[[2]] == 1) {
    # you need this section for the first round
    buzz_left <- non_overlap_foo[side_of_fragment == "L"]
    buzz_right <- non_overlap_foo[side_of_fragment == "R"]
    
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
        dplyr::rename(start = single_bp_start,
               stop = single_bp_stop)
} else {
    ### grab coordinates depending on left or right
    # the left has to go more left, and the right more right
    buzz_left <- non_overlap_foo[side_of_fragment == "L" & str_detect(associations, "_L")]
    buzz_right <- non_overlap_foo[side_of_fragment == "R" & str_detect(associations, "_R")]
    
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
        dplyr::rename(start = single_bp_start,
               stop = single_bp_stop) %>% as.data.frame()
}

### export
number_of_guide_regions_left <- length(non_overlap_foo$associations)
#write(number_of_guide_regions_left, paste0("leftright_coords_", args[[2]], ".txt"))

#print(foobar)

#write_tsv(foobar, paste0("leftright_coords_round", args[[2]], ".bed"))
#write_tsv(foobar, paste0("leftright_coords_reselected_round", args[[2]], ".bed"))
