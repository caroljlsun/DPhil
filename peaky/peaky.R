
# Setup -------------------------------------------------------------------

#r has to be 4.4

#library(Chicago)
library(tidyverse)
#library(devtools)
library(rJavaEnv)
library(peaky)
library(data.table)
library(purrr)

setwd("comp_bio/capture_c_adipocytes/peaky/")
#options(scipen = 999)

# Library installations ---------------------------------------------------


#install.packages("devtools")
#install.packages("stringi")
#install_github("cqgd/R2BGLiMS")
#install.packages("cowplot")
#install_github("cqgd/pky")


# Java installation -------------------------------------------------------



#https://cran.rstudio.com/web/packages/rJavaEnv/vignettes/rJavaEnv.html
#install.packages("rJavaEnv")
library(rJavaEnv)
#java_quick_install(version = 21)
java_check_version_rjava()


# Import example data -----------------------------------------------------

base = system.file("extdata",package="peaky")

interactions_file = paste0(base,"/counts.tsv")
fragments_file = paste0(base,"/fragments.bed")

interactions = data.table::fread(interactions_file)
fragments = data.table::fread(fragments_file)

print(head(interactions))
print(head(fragments))


# Running example ---------------------------------------------------------


# modified example according to author comments:
# Have you tried limiting the "analysis" to fragments closer by, where there is less noise, and it's easier to model N as a function of distance
BI = bin_interactions(interactions, fragments, bins=5, , max_dist= 500e3)

print(BI$bins)
print(BI$interactions)

models_example = by(BI$interactions, BI$interactions$dist.bin, model_bin, subsample_size=450)

plot(models_example[[3]]$fit)

# BM = model_bin(BI$interactions[dist.bin==2,], subsample_size=9000)
# plot(BM$fit)

residuals = lapply(models_example, "[[", "residuals")
bins = split(BI$interactions, BI$interactions$dist.bin)
BTS = split_baits(bins, residuals)

relevant_bait = BTS[baitID==618421]
zoom = relevant_bait[abs(relevant_bait$dist)<1e6]
plot(x=zoom$dist,
     y=zoom$residual,
     xlab="Distance from bait (bp)",
     ylab="Adjusted readcount",
     main=paste("Bait",unique(zoom$baitID)))


#Isolate the highest peak and adjusted readcounts within 40 kbp
peak_center = zoom$dist[which.max(zoom$residual)]
zoom$distance_from_peak = abs(zoom$dist-peak_center)
single_peak = zoom[distance_from_peak<4e5,]

#Define the distance function and find the best fit for omega
dist_exp = function(offset, strength, omega_power){
    omega = 10^omega_power
    strength * exp(-(abs(offset * omega)))
}

fit_omega = nls(residual~dist_exp(distance_from_peak,strength,omega_power),
                data=single_peak,
                start=list(strength=5, omega_power=-3))

filtered_single_peak <- single_peak[!is.na(single_peak$residual), ]


{ # curlies mean execute the lines together
    plot(single_peak$distance_from_peak, single_peak$residual,
         main=paste0("Decay of signal for an isolated peak\nEstimate of omega = 10^",
                     round(coefficients(fit_omega)["omega_power"],3)),
         xlab="Distance from center of peak (bp)",
         ylab="Adjusted readcount")
    lines(filtered_single_peak$distance_from_peak, fitted.values(fit_omega), col="red", lwd=3)
}

relevant_bait = BTS[baitID==618421]
omega_power = -3


PKS = peaky(relevant_bait, omega_power, iterations=10e3)
P = interpret_peaky(relevant_bait, PKS, omega_power)

par(mfrow=c(3,1))
zoom = P[abs(P$dist)<1e6]
plot(x=zoom$dist, xlab="Distance from bait (bp)",
     y=zoom$residual, ylab="Adjusted readcount")

plot(x=zoom$dist, xlab="Distance from bait (bp)",
     y=zoom$beta_mean, ylab="Mean contact strength",
     col="green")

plot(x=zoom$dist, xlab="Distance from bait (bp)",
     y=zoom$rjmcmc_pos, ylab="MPPC",
     col="blue")


# Import own data ---------------------------------------------------------


interactions_mine = data.table::fread("comp_bio/capture_c_adipocytes/main_capture_analysis/all_pADs_blacklisted.chinput") %>% 
    rename(preyID = otherEndID) %>% 
    select(baitID, preyID, N) %>% 
    na.omit(.)

#write_delim(interactions_mine, "all_P_cleaned_interactions.tsv")

fragments_mine = data.table::fread("restriction_fragments.bed") %>% 
    rename(chrom = V1,
           chromStart = V2, 
           chromEnd = V3, 
           ID = V4) %>% 
    na.omit(.)

#write_delim(fragments_mine, "cleaned_fragments.tsv")

# peaky_iterate(max_dist = 250e3, subsample_size = 15000, n_bins = 6, folder_name = "M1/", interaction_file = "M1_cleaned_interactions.tsv", fragment_file = "cleaned_fragments.tsv")
# Run own data ------------------------------------------------------------
# 
# interactions_mine <- fread("M1_cleaned_interactions.tsv")
# 
#fragments_mine <- fread("../peaky/cleaned_fragments.tsv")

### excluding peaks which are within 250kbp of each other... ### WAIT A MINUTE

# merged_by_250 <- fread("../main_capture_analysis/master_250kbp.bed")
# 
# exclude_by_250 <- merged_by_250 %>% 
#     group_by(V11) %>% 
#     filter(n() == 1)
# 
# interactions_excluded <- semi_join(interactions_mine, exclude_by_250, by = join_by("baitID" == "V10"))
# 
# write_delim(interactions_excluded, "excluded_250_interactions.tsv")
# 
# hmm_by_250 <- merged_by_250 %>% 
#     group_by(V11) %>% 
#     filter(n() != 1)
# 
# clusters <- split(hmm_by_250, hmm_by_250$V11)
# 
# # Find the maximum number of rows in any cluster
# max_rows <- max(sapply(clusters, nrow))
# 
# # Create separate tables, each with one entry from each cluster
# tables <- list()
# for (i in 1:max_rows) {
#     # For each cluster, get the i-th row if it exists
#     tables[[i]] <- rbindlist(lapply(clusters, function(cluster_data) {
#         if (i <= nrow(cluster_data)) {
#             return(cluster_data[i])
#         } else {
#             return(NULL)
#         }
#     }), fill = TRUE)
# }
# 
# # Access individual tables
# table1 <- tables[[1]]  # First entry from each cluster
# table2 <- tables[[2]]  # Second entry from each cluster (if exists)
# table3 <- tables[[3]]  # Third entry from each cluster (if exists)

# try again, but with much more limited sizes
BI_mine = bin_interactions(interactions_mine, fragments_mine, bins=5, max_dist=1.25e5)

#BI_mine = bin_interactions_mine(interactions_mine, fragments_mine, bins=3, max_dist= 25e3)

print(BI_mine$bins)
#print(BI_mine$interactions)

models = by(BI_mine$interactions, BI_mine$interactions$dist.bin, model_bin, subsample_size=7500)

plot(models[[1]]$fit)

table(is.na(models[[1]]$residuals))
table(is.na(models[[4]]$residuals))

termplot(models[[1]]$fit)
plot(resid(models[[1]]))

# check that 'distant' bins contain enough of variance in read counts to make accurate models
#BI_mine$interactions %>% subset (dist.bin == 6) %>% select (N) %>% table %>% t

residuals = lapply(models, "[[", "residuals")
bins = split(BI_mine$interactions, BI_mine$interactions$dist.bin)
BTS = split_baits(bins, residuals)

relevant_bait = BTS[baitID==43123]
zoom = relevant_bait[abs(relevant_bait$dist)<5e5]
plot(x=zoom$dist,
     y=zoom$residual,
     xlab="Distance from bait (bp)",
     ylab="Adjusted readcount",
     main=paste("Bait",unique(zoom$baitID)))


#Isolate the highest peak and adjusted readcounts within 40 kbp
peak_center = zoom$dist[which.max(zoom$residual)]
zoom$distance_from_peak = abs(zoom$dist-peak_center)
single_peak = zoom[distance_from_peak<4e5,]

### HMMMMMMMM ###
#Isolate the highest peak and adjusted readcounts within 125 kbp
peak_center = zoom$dist[which.max(zoom$residual)]
zoom$distance_from_peak = abs(zoom$dist-peak_center)
single_peak = zoom[distance_from_peak<250e3,]



#Define the distance function and find the best fit for omega
dist_exp = function(offset, strength, omega_power){
    omega = 10^omega_power
    strength * exp(-(abs(offset * omega)))
}

fit_omega = nls(residual~dist_exp(distance_from_peak,strength,omega_power),
                data=single_peak,
                start=list(strength=5, omega_power=-3))

filtered_single_peak <- single_peak[!is.na(single_peak$residual), ]


{
    plot(single_peak$distance_from_peak, single_peak$residual,
         main=paste0("Decay of signal for an isolated peak\nEstimate of omega = 10^",
                     round(coefficients(fit_omega)["omega_power"],3)),
         xlab="Distance from center of peak (bp)",
         ylab="Adjusted readcount")
    lines(filtered_single_peak$distance_from_peak, fitted.values(fit_omega), col="red", lwd=3)
}

hmm <- fitted.values(fit_omega)

coef(fit_omega)
relevant_bait = BTS[baitID==43123]
omega_power = -4.2


PKS = peaky(relevant_bait, omega_power, iterations=1e6)
P = interpret_peaky(relevant_bait, PKS, omega_power)

par(mfrow=c(3,1))
zoom = P[abs(P$dist)<1e6]
plot(x=zoom$dist, xlab="Distance from bait (bp)",
     y=zoom$residual, ylab="Adjusted readcount")

plot(x=zoom$dist, xlab="Distance from bait (bp)",
     y=zoom$beta_mean, ylab="Mean contact strength",
     col="green")

plot(x=zoom$dist, xlab="Distance from bait (bp)",
     y=zoom$rjmcmc_pos, ylab="MPPC",
     col="blue")

zoom[which(zoom$rjmcmc_pos == max(zoom$rjmcmc_pos))]

zoom[which(zoom$rjmcmc == max(zoom$rjmcmc))]



# filesystem application --------------------------------------------------


base = system.file("extdata",package="peaky")
interactions_file = paste0(base,"/counts.tsv")

bins_dir = paste0(base,"/bins")
fragments_file = paste0(base,"/fragments.bed")

bin_interactions_fs(interactions_file, fragments_file, output_dir=bins_dir, max_dist= 500e3)

fits_dir = paste0(base,"/fits")

for(bin_index in 1:5){
    model_bin_fs(bins_dir,bin_index,output_dir=fits_dir,subsample_size=450)
}

baits_dir = paste0(base,"/baits")

split_baits_fs(bins_dir,residuals_dir = fits_dir, indices=1:5, output_dir = baits_dir)

baitlist = paste0(baits_dir,"/baitlist.txt")
rjmcmc_dir = paste0(base,"/rjmcmc")
omega_power = -5

for(i in 1:5){
    peaky_fs(baitlist,i,output_dir=rjmcmc_dir,omega_power=omega_power)
}

rjmcmc_list(rjmcmc_dir)

rjmcmclist = paste0(rjmcmc_dir,"/rjmcmclist.txt")
baits_rjmcmc_dir = paste0(base,"/baits_rjmcmc")
interpret_peaky_fs(rjmcmclist,1,baits_dir,baits_rjmcmc_dir,omega_power)

### how do we interpret the resulting rds, and can we modify where it saves so it doesn't land in the personal directory

# my file paths
my_base <- "comp_bio/capture_c_adipocytes/peaky/excluded/"
my_int_file <- paste0("excluded_250_interactions.tsv")
my_bin_dir <- paste0(my_base,"bins")
my_frag_file <- paste0("cleaned_fragments.tsv")

# binning interactions 

bin_interactions_fs(my_int_file, my_frag_file, output_dir=my_bin_dir, max_dist= 250e3)

# modelling gamlss to different bins

my_fits_dir = paste0(my_base,"fits")

for(bin_index in 1:5){
    model_bin_fs(my_bin_dir,bin_index,output_dir=my_fits_dir,subsample_size=8000)
}

# running peaky on each bin per bait

my_baits_dir = paste0(my_base,"baits")

split_baits_fs(my_bin_dir,residuals_dir = my_fits_dir, indices=1:5, output_dir = my_baits_dir)
# also plots them?

baitlist = paste0(my_baits_dir,"/baitlist.txt")
my_rjmcmc_dir = paste0(my_base,"rjmcmc")
omega_power = -4.2

# for(i in 1:10){
#     peaky_fs(baitlist,i,output_dir=my_rjmcmc_dir,omega_power=omega_power)
# }

# function breaks if the bait is empty (of which there are many)
# best solution would probably be to submit as multiple parallel scripts
# but it seems fast enough to do locally so:

for (i in 1:10) {
    
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch(peaky_fs(baitlist,i,output_dir=my_rjmcmc_dir,omega_power=omega_power), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }     
}

## interpret the peaky data

rjmcmc_list(my_rjmcmc_dir)

rjmcmclist = paste0(my_rjmcmc_dir,"/rjmcmclist.txt")
my_baits_rjmcmc_dir = paste0(my_base,"baits_rjmcmc")
#interpret_peaky_fs(rjmcmclist, 1, my_baits_dir, my_baits_rjmcmc_dir, omega_power)


for (i in 1:10) {
    
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch(interpret_peaky_fs(rjmcmclist, i, my_baits_dir, my_baits_rjmcmc_dir, omega_power)
             , error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }     
}


## plots

peaky_wrapup(my_base)


#### get the chr start stop 



peaky_output <- fread("excluded/peaky_output.csv")


foo <- merge(peaky_output,fragments_mine[,.(ID,start,end)],by.x="baitID",by.y="ID")
fii <- merge(foo,fragments_mine[,.(ID,start,end)],by.x="preyID",by.y="ID", ,suffixes=c(".bait",".prey"))

fii$p.chr<- paste0("chr", fii$p.chr)

qux <- fii %>% 
    dplyr::select(chr = p.chr,
                  start = start.prey,
                  stop = end.prey,
                  MPCC = rjmcmc_pos)

qux$MPCC <- as.character(qux$MPCC)

write_tsv(qux, "excluded/plot_all_mpcc.bedGraph", col_names = F)

# LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s plot_all_mpcc.bedGraph > plot_all_mpcc.bedGraph.tmp
# bedtools merge -i plot_all_mpcc.bedGraph.tmp -c 4 -d 0 -o max > plot_all_mpcc.bedGraph.out.bdg
# LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s plot_all_mpcc.bedGraph.out.bdg > plot_all_mpcc.bedGraph.sorted
# bedGraphToBigWig plot_all_mpcc.bedGraph.sorted comp_bio/capture_c_adipocytes/ref_genomes/sizes.hg.38.fai plot_all_mpcc.bigwig


### new bigwig but with only "significant" mpccs


bar <- fii %>% 
    filter(rjmcmc_pos >= 0.1) %>% 
    dplyr::select(chr = p.chr,
                  start = start.prey,
                  stop = end.prey,
                  MPCC = rjmcmc_pos)

bar$MPCC <- as.character(bar$MPCC)


write_tsv(bar, "excluded/plot_sig_mpcc.bedGraph", col_names = F)

LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s plot_sig_mpcc.bedGraph > plot_sig_mpcc.bedGraph.tmp
bedtools merge -i plot_sig_mpcc.bedGraph.tmp -c 4 -d 0 -o max > plot_sig_mpcc.bedGraph.out.bdg
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s plot_sig_mpcc.bedGraph.out.bdg > plot_sig_mpcc.bedGraph.sorted
bedGraphToBigWig plot_sig_mpcc.bedGraph.sorted comp_bio/capture_c_adipocytes/ref_genomes/sizes.hg.38.fai plot_sig_mpcc.bigwig



# function to iterate faster ----------------------------------------------

# peaky_settings <- c(folder_name = "default/", max_dist = 100e3, subsample_size= 450, n_bins = 5)

# peaky_settings <- c(folder_name = "reduce_bins/", max_dist = 250e3, subsample_size= 450, n_bins = 3)


peaky_iterate = function(folder_name = "manually/", max_dist = 250e3, subsample_size= 450, n_bins = 3, interaction_file, fragment_file){
    my_base <- "comp_bio/capture_c_adipocytes/peaky/"
    save_dirs <- paste0(my_base, folder_name)
    my_bin_dir <- paste0(save_dirs,"bins")
    
    my_int_file <- interaction_file
    my_frag_file <- fragment_file
    
    ##### BINS #####
    ### if this chunk has successfully ran before, do not rerun
    if (file.exists(paste0(my_bin_dir,"/bin_",n_bins , ".rds"))) {
        print("bins already exist!")
    }else{
        
        # binning interactions 
        
        bin_interactions_fs(my_int_file, my_frag_file, output_dir = my_bin_dir, max_dist = max_dist, bins = n_bins)
        
    }
    
    ##### FITS for BINS #####
    
    # modelling gamlss to different bins
    
    my_fits_dir = paste0(save_dirs,"fits")
    
    
    ### if this chunk has successfully ran before, do not rerun
    if (file.exists(paste0(my_fits_dir,"/fit_",n_bins , ".rds"))) {
        print("fits for bins already exist!")
        
    }else{
        
        for(bin_index in 1:n_bins){
            model_bin_fs(my_bin_dir, bin_index, output_dir = my_fits_dir, subsample_size = subsample_size)
        }
    }
    
    ##### BAIT SPLITTING #####
    # running peaky on each bin per bait
    
    my_baits_dir = paste0(save_dirs,"baits")
    
    ### if this chunk has successfully ran before, do not rerun
    
    if (file.exists(paste0(my_baits_dir,"/baitlist.txt"))) {
     baitlist <- fread(paste0(my_baits_dir,"/baitlist.txt"))   
    }else{

            split_baits_fs(my_bin_dir,residuals_dir = my_fits_dir, indices=1:5, output_dir = my_baits_dir)
            # also plots them?

        }
    

    
    ##### PEAKY per BAIT #####
    
    
    baitlist = paste0(my_baits_dir,"/baitlist.txt")
    my_rjmcmc_dir = paste0(save_dirs,"rjmcmc")
    omega_power = -4.2
    
    # peaky_fs function breaks if the bait is empty (of which there are many)
    # best solution would probably be to submit as multiple parallel scripts
    # but it seems fast enough to do locally so
    
    bait_count <- length(list.files(path = my_baits_dir, pattern = "\\.rds$"))
    
    if (file.exists(paste0(my_rjmcmc_dir,"/rjmcmc_", word(paste0(last(baitlist)), -1, sep = fixed("_"))))) {
        print("rjmcmc's already exist!")
        
    }else{
        
        for (i in 1:bait_count) {
            
            skip_to_next <- FALSE
            
            tryCatch(peaky_fs(baitlist,i,output_dir=my_rjmcmc_dir,omega_power=omega_power), error = function(e) { skip_to_next <<- TRUE})
            
            if(skip_to_next) { next }     
        }
    }
    
    ##### PEAKY INTERPRET #####
    
    rjmcmc_list(my_rjmcmc_dir)
    
    rjmcmclist = paste0(my_rjmcmc_dir,"/rjmcmclist.txt")
    my_baits_rjmcmc_dir = paste0(save_dirs,"baits_rjmcmc")
    
    
    if (file.exists(paste0(my_baits_rjmcmc_dir,"/bait_rjmcmc_", word(paste0(last(baitlist)), -1, sep = fixed("_"))))) {
        print("interpreted rjmcmc's already exist!")
        
    }else{
        
        for (i in 1:bait_count) {
            
            skip_to_next <- FALSE
            
            # Note that print(b) fails since b doesn't exist
            
            tryCatch(interpret_peaky_fs(rjmcmclist, i, my_baits_dir, my_baits_rjmcmc_dir, omega_power)
                     , error = function(e) { skip_to_next <<- TRUE})
            
            if(skip_to_next) { next }     
        }
        
    }
    ## plots
    
    peaky_wrapup(save_dirs)
    
    
    #### get the chr start stop 
    
    peaky_output <- fread(paste0(save_dirs,"peaky_output.csv"))
    
    foo <- merge(peaky_output,fragments_mine[,.(ID,start,end)],by.x="baitID",by.y="ID")
    fii <- merge(foo,fragments_mine[,.(ID,start,end)],by.x="preyID",by.y="ID", ,suffixes=c(".bait",".prey"))
    
    fii$p.chr<- paste0("chr", fii$p.chr)
    fii$b.chr<- paste0("chr", fii$b.chr)
    
    qux <- fii %>% 
        dplyr::select(chr_prey = p.chr,
                      start_prey = start.prey,
                      stop_prey = end.prey,
                      MPCC = rjmcmc_pos,
                      chr_bait = b.chr,
                      start_bait = start.bait,
                      stop_bait = end.bait)
    
    
    write_tsv(qux, paste0(save_dirs,"plot_all_mpcc.bedGraph"), col_names = F)
    
    # LC_ALL=C sort -k1,1 -k2,2n plot_all_mpcc.bedGraph > sorted_plot_all_mpcc.bedGraph
    # bedGraphToBigWig sorted_plot_all_mpcc.bedGraph ../ref_genomes/hg38.chrom.sizes all_mpcc.bw
    
    
    ### new bigwig but with only "significant" mpccs
    
    
    bar <- fii %>% 
        filter(rjmcmc_pos >= 0.1) %>% 
        dplyr::select(chr_prey = p.chr,
                      start_prey = start.prey,
                      stop_prey = end.prey,
                      MPCC = rjmcmc_pos,
                      chr_bait = b.chr,
                      start_bait = start.bait,
                      stop_bait = end.bait)
    
    
    write_tsv(bar, paste0(save_dirs,"plot_sig_mpcc.bedGraph"), col_names = F)
    
    
}


#peaky_iterate(max_dist = 250e3, subsample_size = 15000, n_bins = 6, folder_name = "M1/", interaction_file = "M1_cleaned_interactions.tsv", fragment_file = "cleaned_fragments.tsv")


peaky_iterate(max_dist = 1.25e5, subsample_size = 7500, n_bins = 5, folder_name = "P_all_blacklisted/", interaction_file = "all_P_cleaned_interactions.tsv", fragment_file = "cleaned_fragments.tsv")


#peaky_iterate(subsample_size = 1000, n_bins = 6)
#rate <- rate_backoff(pause_base = 0.1, pause_min = 0.005, max_times = 4)

rate <- rate_delay(0.1)


# insistent_peaky <- insistently(peaky_iterate(max_dist = 250e3, folder_name = "reduce_bin/", subsample_size = 10000, n_bins = 3), rate)


insistent_peaky <- insistently(peaky_iterate, rate)
set.seed(8)
try(insistent_peaky())


################ rewrite #####


#### get the chr start stop for baits AND preys
### and also get their counts
save_dirs <- "comp_bio/capture_c_adipocytes/peaky/P_all_blacklisted/"
peaky_output <- fread(paste0(save_dirs,"peaky_output.csv"))


foo <- merge(peaky_output,fragments_mine[,.(ID,start,end)],by.x="baitID",by.y="ID")
fii <- merge(foo,fragments_mine[,.(ID,start,end)],by.x="preyID",by.y="ID", ,suffixes=c(".bait",".prey"))

fii$p.chr<- paste0("chr", fii$p.chr)
fii$b.chr<- paste0("chr", fii$b.chr)

qux <- fii %>% 
    dplyr::select(chr_prey = p.chr,
                  start_prey = start.prey,
                  stop_prey = end.prey,
                  MPCC = rjmcmc_pos,
                  chr_bait = b.chr,
                  start_bait = start.bait,
                  stop_bait = end.bait,
                  rf_counts = N)


#write_tsv(qux, paste0(save_dirs,"plot_all_mpcc.bedGraph"), col_names = F)

# LC_ALL=C sort -k1,1 -k2,2n plot_all_mpcc.bedGraph > sorted_plot_all_mpcc.bedGraph
# bedGraphToBigWig sorted_plot_all_mpcc.bedGraph ../ref_genomes/hg38.chrom.sizes all_mpcc.bw


### new bigwig but with only "significant" mpccs


bar <- fii %>% 
    filter(rjmcmc_pos >= 0.1) %>% 
    dplyr::select(chr_prey = p.chr,
                  start_prey = start.prey,
                  stop_prey = end.prey,
                  MPCC = rjmcmc_pos,
                  chr_bait = b.chr,
                  start_bait = start.bait,
                  stop_bait = end.bait,
                  rf_counts = N)


#write_tsv(bar, paste0(save_dirs,"plot_sig_mpcc.bedGraph"), col_names = F)


################ all interactions binning ################ 

bin_interactions_fs("excluded_250_interactions.tsv", my_frag_file, output_dir = "all_int", max_dist = 2e6, bins = 30) # do I want to bin the excluded 250 ones? might be easier to start iwht it but need to revisit later.

bin_1 <- readRDS("all_int/bin_1.rds")

bin_1_mod <- bin_1 %>% 
    mutate(z_score_count = c(scale(N))) %>% 
    mutate(b.end = round((b.length + 2*b.mid)/2), # you have to round here bc floating points
           b.start = b.mid*2 - b.end) %>% 
    mutate(b.end = ifelse(b.length %% 2 == 0, b.end, b.end +1)) %>% 
    mutate(p.end = round((p.length + 2*p.mid)/2), # you have to round here bc floating points
           p.start = p.mid*2 - p.end) %>% 
    mutate(p.start = ifelse(p.length %% 2 == 0, p.start, p.start -1)) 

# b.end b.start b.length
# 2235104 2233828 1277
# 16013428 1601327 152
# 25794502 25793998 505
# master_fragment_key
# 2233828 2235105
# 16013276 16013428
# 25793998 25794503
# 
# p.end p.start
# 2085351 2084923
# 2086005 2085356
# fragments mine
# 2084923 2085351
# 2085355 2086004


# so if b.length is odd, add 1 to b.end...
# threshold so that it picks up preys that peaky picks up

# we only have p/b.mid
# from the source
# 1) mid=round((end+start)/2) #### !!!!! has been rounded!!!!
# 2) length=end-start


# 1)
# end = mid*2-start
# start = mid*2 - end

# 2)
# end = length + start

# 1) + 2)
# end = length + mid*2 - end
# end*2 = length + mid*2
# end = (length + mid*2)/2

# start = mid*2 - end 


sig_mpcc <- fread("excluded/plot_sig_mpcc.bedGraph")

# match by p.chr p.start p.end
# then compare their z-scores
bin_1_mod$p.chr <- paste0("chr", bin_1_mod$p.chr)
bin_1_mod$b.chr <- paste0("chr", bin_1_mod$b.chr)

# add the "chr"s to p.chr
which_overlap_bin_1 <- inner_join(bin_1_mod, sig_mpcc, by = c("p.chr" = "V1", "p.start" = "V2", "p.end" = "V3" ))

corrected_key_fragments <-  key_fragments %>%
    mutate(chr = V1,
           start = V2+4,
           stop = V3-4)

# plots ----------------------------------------------


# Create data frame with the distances and counts
distances <- c(100, 200, 500, 1000, 2000, 5000, 10000, 20000, 30000, 40000, 50000, 80000, 100000, 120000, 150000)

# Read counts from the files
counts <- c(839, 809, 720, 617, 500, 410, 400, 370, 340, 308, 275, 225, 213, 194, 185)

# Create data frame
data <- data.frame(distance = distances, count = counts)

# Plot
library(ggplot2)
ggplot(data, aes(x = distance, y = count)) +
    geom_line() +
    geom_point() +
    scale_x_log10(breaks = distances, labels = distances) +
    labs(x = "Merge Distance (bp)", y = "Number of Merged Regions", 
         title = "Effect of Merge Distance on Number of Regions") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("../peaky/M_all_blacklisted/merge_distance_plot.png", width = 10, height = 6)
