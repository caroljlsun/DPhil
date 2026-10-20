# -------------------------------------------------------- #
# Procedure: Documenting processing of files
# Project: Capture C
# User: CZJS
# Date: FEB-2025
# -------------------------------------------------------- #

library(tidyverse)
library(data.table)

setwd("comp_bio/capture_c_adipocytes/main_capture_analysis")

final_oligos <- fread("final_merged_oligos_non_overlapping.bed")

unique(final_oligos$associations)


final_uniques <- final_oligos %>% 
    distinct(associations, .keep_all = T) %>% 
    select(chr, fragment_start, fragment_stop, associations)

write.table(final_uniques, "unique_restriction_fragments.bed", sep = "\t", col.names = F, quote = F, row.names = F)

# may need to remove the colons....