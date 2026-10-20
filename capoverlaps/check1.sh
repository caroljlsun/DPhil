

### Setup
cd comp_bio/capture_c_adipocytes/capoverlaps

# reusing another env with the newish bedtools

mamba activate fullcap

### RERUN RERUN RERUN RERUN RERUN ###

### filtered and moved over to my working directory
zcat whradjbmi.giant-ukbb.meta-analysis.females.23May2018_harmonized.tsv.gz | awk -F'\t' '$11 < 0.00000005'| gzip >  comp_bio/capture_c_adipocytes/capoverlaps/whradjbmi_combined_filtered.gz 

### unzipped
gunzip -k whradjbmi_combined_filtered.gz

### cut out columns I think are interesting/useful? and convert to a bed-ish file

cat whradjbmi_combined_filtered | awk -F'\t' -v OFS='\t' '{print "chr"$4, $5, $5+1, $3, $9,  $1, $11}' >  whradjbmi_combined_extra_filtered.bed

### liftover ?
genetic_department/projects/TargetLookup/liftOver whradjbmi_combined_extra_filtered.bed genetic_department/projects/TargetLookup/hg19ToHg38.over.chain.gz whradjbmi_combined_hg38_extra_filtered.bed unmapped_whradjbmi_combined_hg38_extra_filtered.bed 

# the unmapped file is empty, and the wc -l matches between the old and new file

### get the updated rsids? I guess?
# using this: https://genome.sph.umich.edu/wiki/LiftOver#Lift_dbSNP_rs_numbers
# getting the rsids into a format that the python script needs it to be
zcat whradjbmi_combined_filtered.gz| awk -F'\t' '{print $3}' | sed 's/rs//g' | awk -F : '{print $1}' > whradjbmi_combined_filtered_rs.txt

# running the script
python LiftRsNumber.py whradjbmi_combined_filtered_rs.txt > whradjbmi_combined_filtered_lifted_rs.txt

# replace the old column of rsids with the new

awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[NR]=$1;next}{$4=a[FNR]}1' full_lifted_rsid.txt whradjbmi_combined_hg38_extra_filtered.bed > whradjbmi_combined_hg38_extra_filtered_lifted.bed

### Intersect with epigenetics

#just moving them over - it probably isn't best practice but the files seem relatively small?
cp genetic_department/share/GATL_GMGZ/CRE_AdiposeTissue/H3K27ac-ChromatinStates_AdiposeTissue_GMGZ.bed H3K27ac_peaks.bed

cp genetic_department/share/GATL_GMGZ/GSE178794_adiposeTissue_consensusPeaks-from11samples_unionPeaksInAtLeast3samples_b38.bed atac.bed

# bedtools intersect twice to keep information from the files, just in case
bedtools intersect -a whradjbmi_combined_hg38_extra_filtered_lifted.bed -b H3K27ac_peaks.bed -wa| bedtools intersect -b H3K27ac_peaks.bed -a - -wb > whradjbmi_combined_hg38_extra_H3K27ac.bed

wc -l 
4037

bedtools intersect -a whradjbmi_combined_hg38_extra_H3K27ac.bed -b atac.bed -wa| bedtools intersect -b atac.bed -a - -wb > whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed
wc -l 
709


