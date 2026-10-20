module avail
# no homer


# make new conda env
cd main_capture_analysis
#conda create -n homercapture 

conda activate homercapture

#conda install bioconda::homer


### Finding Enriched Motifs in Genomic Regions (findMotifsGenome.pl)

findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]


# extract from master_fragment_key.bed columns 4,5,6 and 7, then empty column, then a column of only 0s or only 1s

cut -f4-7 master_fragment_key.bed | awk '{print $0 "\t\t0"}' > cre_regions_zero_stranded.bed
cut -f4-7 master_fragment_key.bed | awk '{print $0 "\t\t1"}' > cre_regions_ones_stranded.bed


findMotifsGenome.pl cre_regions_zero_stranded.bed comp_bio/capture_c_adipocytes/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna homer_output/ -size 200 -mask

# other strand
findMotifsGenome.pl cre_regions_ones_stranded.bed comp_bio/capture_c_adipocytes/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna homer_output_ones/ -size 200 -mask