mkdir peaky

cd peaky

# I feel as though most files are already in close enough shape to use peaky from chicago


salloc --nodes=1 --time=480 --partition=compute --mem=160G
srun --jobid=5698754 --pty /bin/bash

#the fragment file should have the shape like the following

chrom chromStart chromEnd  ID
1     566754   573939 142
1     753009   760251 196

#the restriction map file from chicago looks like the following

head comp_bio/capture_c_adipocytes/chicago/trial_run_chicago/restriction_fragments.rmap

chr1    12414   12460   0
chr1    12464   12685   1
chr1    12689   12828   2
chr1    12832   13314   3
chr1    13318   13419   4
chr1    13423   13565   5

# so I will modify the chicago restriction map file to look like the fragment file

awk '{gsub(/^chr/, "", $1); print $1, $2, $3, $4}' comp_bio/capture_c_adipocytes/chicago/trial_run_chicago/restriction_fragments.rmap > restriction_fragments.bed


# as for the interactions file from peaky's structure:

baitID preyID N
53559  52096 1
53559  52107 1
53559  52112 1


# and what I have for chicago is:

head comp_bio/capture_c_adipocytes/chicago/trial_run_chicago/trial/trial_updated.chinput

##        samplename=comp_bio/capture_c_adipocytes/chicago/trial_run_chicago/trial        bamname=DAAP6matureadipocytes_M1      baitmapfile=comp_bio/capture_c_adipocytes/chicago/baited_fragments.bed_4col.txt      digestfile=comp_bio/capture_c_adipocytes/chicago/genome.digest.bed
baitID  otherEndID      N       otherEndLen     distSign
886743  60679   1       49      60085
886743  146072  1       58      145432
886743  158481  1       102     157892
886743  318332  1       42      317923
886743  412161  1       108     411937


# so I will modify the chicago interactions file to look like the peaky interactions file
# remove the first row of comments
# and headers

#tail -n +3 comp_bio/capture_c_adipocytes/chicago/trial_run_chicago/trial/trial_updated.chinput  |  awk '{print $1, $2, $3}'  > interactions.tsv


tail -n +2 comp_bio/capture_c_adipocytes/main_capture_analysis/DAAP6mADs_M1_full.chinput  |  awk '{print $1, $2, $3}'  > DAAP6mADs_M1_interactions.tsv

# let's try and run peaky



install.packages("devtools")
library(devtools)
install_github("cqgd/R2BGLiMS")
install_github("cqgd/pky")

mamba create -n peaky python=3.9
mamba activate peaky
mamba install r-devtools
mamba install libssl-dev
mamba install libcurl4-openssl-dev
mamba install r-curl



### run peaky.R

# all baits
LC_ALL=C sort -k1,1 -k2,2n plot_all_mpcc.bedGraph > sorted_plot_all_mpcc.bedGraph
bedGraphToBigWig sorted_plot_all_mpcc.bedGraph ../ref_genomes/hg38.chrom.sizes all_mpcc.bw

# sig baits only
LC_ALL=C sort -k1,1 -k2,2n plot_sig_mpcc.bedGraph > sorted_plot_sig_mpcc.bedGraph
bedGraphToBigWig sorted_plot_sig_mpcc.bedGraph ../ref_genomes/hg38.chrom.sizes sig_mpcc.bw


### -------------------------------- ###
### TSS and gene body intersections with peaky mpcc
### -------------------------------- ###

cat ../peaky/plot_sig_mpcc.bedGraph | awk '{print $1 "\t" $2 "\t" $3}' > sig_mpcc.bed

bedtools intersect -a sig_mpcc.bed -b gene_tss_coordinates_5kbp.bed > sig_mpcc_gene_tss_5kbp.bed

bedtools intersect -a sig_mpcc.bed -b hg38_all_genes.bed -wo > sig_mpcc_gene_bodies.bed

bedtools intersect -a sig_mpcc.bed -b hg38_all_genes_extended.bed > sig_mpcc_gene_bodies_5kbp.bed
