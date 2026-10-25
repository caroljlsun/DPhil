#---------------------------------------------
# Procedure: Script to reiterate peaky
# Project: Peaky
# User: CZJS
# Date: Jan-2025
#---------------------------------------------



# Libraries ---------------------------------------------------------------



#library(Chicago)
library(tidyverse)
#library(devtools)
library(rJavaEnv)
library(peaky)
library(data.table)

java_check_version_rjava()


# Argument definitions ----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
print(paste("results saved under",args[[1]]))
print(paste("max distance for interactions is",args[[2]]))
print(paste("subsampling size is",args[[3]]))
print(paste("number of bins is",args[[4]]))

folder_name <- args[[1]]
max_dist <- args[[2]]
subsample_size <- args[[3]]
n_bins <- args[[4]]


# Main iteration ----------------------------------------------------------


my_base <- "comp_bio/capture_c_adipocytes/peaky/"
save_dirs <- paste0(my_base, folder_name)
my_int_file <- paste0(my_base,"cleaned_interactions.tsv")
my_bin_dir <- paste0(save_dirs,"bins")
my_frag_file <- paste0(my_base,"cleaned_fragments.tsv")

# binning interactions 

bin_interactions_fs(my_int_file, my_frag_file, output_dir = my_bin_dir, max_dist)

# modelling gamlss to different bins

my_fits_dir = paste0(save_dirs,"fits")

for(bin_index in 1:n_bins){
    model_bin_fs(my_bin_dir, bin_index, output_dir = my_fits_dir, subsample_size)
}

# running peaky on each bin per bait

my_baits_dir = paste0(save_dirs,"baits")

split_baits_fs(my_bin_dir,residuals_dir = my_fits_dir, indices=1:5, output_dir = my_baits_dir)
# also plots them?

baitlist = paste0(my_baits_dir,"/baitlist.txt")
my_rjmcmc_dir = paste0(save_dirs,"rjmcmc")
omega_power = -5

# peaky_fs function breaks if the bait is empty (of which there are many)
# best solution would probably be to submit as multiple parallel scripts
# but it seems fast enough to do locally so:
bait_count <- length(list.files(path = my_baits_dir, pattern = "\\.rds$"))

for (i in 1:bait_count) {
    
    skip_to_next <- FALSE
    
    tryCatch(peaky_fs(baitlist,i,output_dir=my_rjmcmc_dir,omega_power=omega_power), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }     
}

## interpret the peaky data

rjmcmc_list(my_rjmcmc_dir)

rjmcmclist = paste0(my_rjmcmc_dir,"/rjmcmclist.txt")
my_baits_rjmcmc_dir = paste0(my_base,"baits_rjmcmc")


for (i in 1:bait_count) {
    
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch(interpret_peaky_fs(rjmcmclist, i, my_baits_dir, my_baits_rjmcmc_dir, omega_power)
             , error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }     
}


## plots

peaky_wrapup(save_dirs)


#### get the chr start stop 

peaky_output <- fread(paste0(save_dirs,"peaky_output.csv"))
fragments_mine <- fread(my_frag_file)

foo <- merge(peaky_output,fragments_mine[,.(ID,start,end)],by.x="baitID",by.y="ID")
fii <- merge(foo,fragments_mine[,.(ID,start,end)],by.x="preyID",by.y="ID", ,suffixes=c(".bait",".prey"))

fii$p.chr<- paste0("chr", fii$p.chr)

qux <- fii %>% 
    dplyr::select(chr = p.chr,
                  start = start.prey,
                  stop = end.prey,
                  MPCC = rjmcmc_pos)


write_tsv(qux, paste0(save_dirs,"plot_all_mpcc.bedGraph"), col_names = F)


### new bigwig but with only "significant" mpccs


bar <- fii %>% 
    filter(rjmcmc_pos >= 0.1) %>% 
    dplyr::select(chr = p.chr,
                  start = start.prey,
                  stop = end.prey,
                  MPCC = rjmcmc_pos)


write_tsv(bar, paste0(save_dirs,"plot_sig_mpcc.bedGraph"), col_names = F)

