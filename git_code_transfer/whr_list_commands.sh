################################
### Setup ###
################################

cd /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps

# reusing another env with the newish bedtools

mamba activate fullcap


################################
### filtered and moved over to my working directory ###
################################
cd /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/harmonized_data/harmonized_gwas
zcat whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz | wc -l #27374931
zcat whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz | awk -F'\t' '$11 < 0.00000005'| gzip >  /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/whradjbmi_combined_filtered.gz 

zcat /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/whradjbmi_combined_filtered.gz | wc -l #54362
################################
### unzipped ###
################################

cd /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps

gunzip -k whradjbmi_combined_filtered.gz

### cut out columns I think are interesting/useful? and convert to a bed-ish file

cat whradjbmi_combined_filtered | awk -F'\t' -v OFS='\t' '{print "chr"$4, $5, $5+1, $3, $9,  $1, $11}' >  whradjbmi_combined_extra_filtered.bed

################################
### liftover  ###
################################

/nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/liftOver whradjbmi_combined_extra_filtered.bed /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/hg19ToHg38.over.chain.gz whradjbmi_combined_hg38_extra_filtered.bed unmapped_whradjbmi_combined_hg38_extra_filtered.bed 

#Deleted in new
#chr20   45796954        45796955        rs3091475:T:C   0.0138  W       4.399e-12


################################
### get the updated rsids###
################################

# using this: https://genome.sph.umich.edu/wiki/LiftOver#Lift_dbSNP_rs_numbers
# getting the rsids into a format that the python script needs it to be
zcat whradjbmi_combined_filtered.gz| awk -F'\t' '{print $3}' | sed 's/rs//g' | awk -F : '{print $1}' > whradjbmi_combined_filtered_rs.txt

# running the script
python2 LiftRsNumber.py whradjbmi_combined_filtered_rs.txt > whradjbmi_combined_filtered_lifted_rs.txt

# get the beginning of the rsid
cat whradjbmi_combined_filtered_lifted_rs.txt | awk -F'\t' '{print "rs"$2}' > rs_appended.txt

# get the ends of the rsid

cat whradjbmi_combined_filtered | awk -F'\t' '{print $3}' | awk -F: -v OFS=':' '{print "",$2, $3}' > rsid_ends.txt

paste -d '' rs_appended.txt rsid_ends.txt > full_lifted_rsid.txt

# remove that one awkard missing rsid

grep -v 'rs3091475:T:C' full_lifted_rsid.txt > full_lifted_rsid_cleaned.txt

# replace the old column of rsids with the new

awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[NR]=$1;next}{$4=a[FNR]}1' full_lifted_rsid_cleaned.txt whradjbmi_combined_hg38_extra_filtered.bed > whradjbmi_combined_hg38_extra_filtered_lifted.bed

################################
### Intersect with epigenetics ###
################################

#just moving them over - it probably isn't best practice but the files seem relatively small?
cp /nfs_home/projects/departments/nnrco/genetic_department/share/GATL_GMGZ/CRE_AdiposeTissue/H3K27ac-ChromatinStates_AdiposeTissue_GMGZ.bed H3K27ac_peaks.bed

cp /nfs_home/projects/departments/nnrco/genetic_department/share/GATL_GMGZ/GSE178794_adiposeTissue_consensusPeaks-from11samples_unionPeaksInAtLeast3samples_b38.bed atac.bed

# bedtools intersect twice to keep information from the files, just in case
bedtools intersect -a whradjbmi_combined_hg38_extra_filtered_lifted.bed -b H3K27ac_peaks.bed -wa| bedtools intersect -b H3K27ac_peaks.bed -a - -wb > whradjbmi_combined_hg38_extra_H3K27ac.bed

wc -l 
5289

bedtools intersect -a whradjbmi_combined_hg38_extra_H3K27ac.bed -b atac.bed -wa| bedtools intersect -b atac.bed -a - -wb > whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed
wc -l 
966

################################
### make H3K27ac only and ATAC only interesected files ###
################################

bedtools intersect -a whradjbmi_combined_hg38_extra_filtered_lifted.bed -b H3K27ac_peaks.bed -wa| bedtools intersect -b H3K27ac_peaks.bed -a - -wb > whradjbmi_combined_hg38_extra_H3K27ac_only.bed
wc -l 
5289

bedtools intersect -a whradjbmi_combined_hg38_extra_filtered_lifted.bed -b atac.bed -wa| bedtools intersect -b atac.bed -a - -wb > whradjbmi_combined_hg38_extra_ATAC_only.bed
wc -l 
1250

# set aside for future use, if not enough regions for capture c


################################
### remove if within 5 kbp of a gene ###
################################

#full list of genes in hg38 is here /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/gencode.v26.GRCh38.genes_for_gene.protein_coding.txt
# cp /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/gencode.v26.GRCh38.genes_for_gene.protein_coding.txt hg38_all_genes.bed
# bedtools slop -i hg38_all_genes.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b 5000 > hg38_all_genes_extended.bed
# bedtools intersect -v -a whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed -b hg38_all_genes_extended.bed | wc -l
# bedtools intersect -v -a whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed -b hg38_all_genes_extended.bed > whradjbmi_combined_hg38_intergenic.bed
# ### merging snps within 50+/- kbp of each other
# bedtools merge -d 50000 -i whradjbmi_combined_hg38_intergenic.bed > whradjbmi_combined_hg38_intergenic_merged.bed
# # or?
# bedtools cluster -d 50000 -i whradjbmi_combined_hg38_intergenic.bed > whradjbmi_combined_hg38_intergenic_clustered.bed
# wc -l 168

################################
### never mind!! don't remove genic snps, rerun merge and clump ###
################################

# move a list of hg38 genes over 
cp /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/gencode.v26.GRCh38.genes_for_gene.protein_coding.txt hg38_all_genes.bed

# merge the snps within 50+/- kbp of each other
bedtools merge -d 50000 -i whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed > whradjbmi_combined_hg38_merged.bed
wc -l 370

# cluster the snps within 50+/- kbp of each other
# gives me the identity of the merged snps in each cluster
bedtools cluster -d 50000 -i whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed > whradjbmi_combined_hg38_clustered.bed
wc -l 966

################################
### assess complexity by looking at how many genes are close by ###
################################

# uni supervisor suggested I find the genes close by the clusters, not the snps

bedtools slop -i whradjbmi_combined_hg38_merged.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b 500000 > whradjbmi_combined_hg38_merged_500kbp.bed

# +/- 500kbp, so a million after expr 

bedtools intersect -wa -wb -a whradjbmi_combined_hg38_merged_500kbp.bed -b hg38_all_genes.bed > whradjbmi_combined_hg38_merged_500kbp_genes.bed
wc -l 6313


# check the length of the unique snps
cat whradjbmi_combined_hg38_merged_500kbp_genes.bed | awk '{print $1,$2,$3}' | uniq | wc -l
369

# I think one of them simply didn't intersect at all, 369 = 370 - 1
# need to find which one that is...
bedtools intersect -wa -wb -v -a whradjbmi_combined_hg38_merged_500kbp.bed -b hg38_all_genes.bed > whradjbmi_combined_hg38_merged_500kbp_nogenes.bed 
# empty 4th col, so I'll just add a column with NA
sed -i 's/$/\tNA/' whradjbmi_combined_hg38_merged_500kbp_nogenes.bed

# yep one of the snps didn't have intersections +/- 500kbp genes

# unique, process columns such that you get one column with all the intersected genes in one place

awk 'BEGIN{OFS="\t"; FS=OFS} {a[$1 OFS $2 OFS $3]=a[$1 OFS $2 OFS $3]","$10} END {for (i in a) {split(i, arr, OFS); if (a[i] ~ /^,/) {a[i] = substr(a[i], 2)}; print arr[1], arr[2], arr[3], a[i]}}' whradjbmi_combined_hg38_merged_500kbp_genes.bed > whradjbmi_combined_hg38_merged_500kbp_genes_concatenated.bed
wc -l 369

# append the missing snp/cluster

cat whradjbmi_combined_hg38_merged_500kbp_nogenes.bed >> whradjbmi_combined_hg38_merged_500kbp_genes_concatenated.bed
wc -l 370

################################
### Run capsequm on snps ###
################################
# need to make a bed file with only 4 columns: chr start stop unique_name

awk 'BEGIN{OFS="\t"; FS=OFS} {print $1, $2, $3, $1"_"$2"_"$3"_"$4"_cluster_"$17"_whradjbmi_combined_H3K27ac_ATAC" }' whradjbmi_combined_hg38_clustered.bed | head
awk 'BEGIN{OFS="\t"; FS=OFS} {print $1, $2, $3, $1"_"$2"_"$3"_"$4"_cluster_"$17"_whradjbmi_combined_H3K27ac_ATAC" }' whradjbmi_combined_hg38_clustered.bed > whradjbmi_combined_hg38_clustered_4col.bed

cd /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh whradjbmi_combined_hg38_extra_H3K27ac_ATAC_final /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/whradjbmi_combined_hg38_clustered_4col.bed
#/fsx/home_dirs/czjs/doctoral_work/repeat-annotation/RepeatMasker/RepeatMasker -species "Homo sapiens" -pa 4 /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/whradjbmi_combined_hg38_extra_H3K27ac_ATAC_3/oligo_seqs.fa

################################
### Making tidy snp file ###
################################

## add to the snp files whether or not an oligo could be designed for it

cd /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/whradjbmi_combined_hg38_extra_H3K27ac_ATAC_final

wc -l oligo_info.txt
1433

wc -l /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/whradjbmi_combined_hg38_clustered_4col.bed
966

# capsequm designs 2 probes per region, 2*966 = 1932 , 1932-1433 = 499, 499/2 =~ 250 regions either don't have any probe, or have half the probes etc.
# find which ones have probes and which don't, turn that to a binary yes no column

## structure of oligo_info.txt

# chr     start   stop    fragment_start  fragment_stop   side_of_fragment        sequence        total_number_of_alignments      density_score    repeat_length   repeat_class    GC%     associations
# chr1    3031529 3031599 3031529 3031673 L       GATCCCTTGGGGGGCCGGCCAGCCACAAGGACGCCCAGCAGGCGACTCCTGGATTCACCCCCTCCCCACC  1       1.0     0NA      71      chr1_3031614_3031615_rs10909864T:C_whradjbmi_combined_H3K27ac_ATAC,
# chr1    3031603 3031673 3031529 3031673 R       TCTCTGGACCTGAGGCCCGCAGTGCTCTGGCAGGGCCAGCAATGACGCTGACGCTGGGAGAGATGTGATC  1       1.0     0NA      63      chr1_3031614_3031615_rs10909864T:C_whradjbmi_combined_H3K27ac_ATAC,

## structure of whradjbmi_clustered

# chr1    3031614 3031615 rs10909864T:C   -0.0215 w       1.629e-10       chr1    3031036 3033036 peak_85 chr1    3031539 3032321 peak89  7       1
# chr1    9292645 9292646 rs114709597C:G  0.0226  w       6.55e-12        chr1    9289341 9298741 peak_228        chr1    9292558 9293469 peak197 6       2
# chr1    11805747        11805748        rs13306561A:G   0.017   w       4.194e-08       chr1    11801343        11807743        peak_318        chr1    11805513        11806521   peak264 11      3

## Make key
cd /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps

awk 'BEGIN{OFS="\t"; FS=OFS} {print $1"_"$2"_"$3"_"$4"_cluster_"$17"_whradjbmi_combined_H3K27ac_ATAC" }' whradjbmi_combined_hg38_clustered.bed > temp.txt
paste whradjbmi_combined_hg38_clustered.bed temp.txt > whradjbmi_combined_hg38_clustered_keyed.bed
rm temp.txt

## sort before join
# 18th column is the key
# awk '{print $18}' whradjbmi_combined_hg38_clustered_keyed.bed| head # quick check

sort -k18,18 whradjbmi_combined_hg38_clustered_keyed.bed > whradjbmi_combined_hg38_clustered_sorted.bed

# 13th column is the key
# awk '{print $13}' oligo_info.txt| head # quick check
cp /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/whradjbmi_combined_hg38_extra_H3K27ac_ATAC_final/oligo_info.txt whradjbmi_oligo_info.txt
sort -k13,13 whradjbmi_oligo_info.txt > whradjbmi_oligo_info_sorted.txt
sed 's/,$//' whradjbmi_oligo_info_sorted.txt > whradjbmi_oligo_info_sorted_no_comma.txt
sed -i '1d' whradjbmi_oligo_info_sorted_no_comma.txt

## join 

join -1 18 -2 13 -t $'\t' -a 1 -o 1.18 whradjbmi_combined_hg38_clustered_sorted.bed whradjbmi_oligo_info_sorted_no_comma.txt | uniq -c > joined_counts.txt

# count the number of regions with 1, 2, or none probes
awk '{if ($1 == 1) print $0, "one"; else if ($1 == 2) print $0, "both"; else print $0, "none";}' joined_counts.txt > final_output.txt
# this is overkill, just use the first column
grep "one" final_output.txt | wc -l
# 250, as it should be

paste whradjbmi_combined_hg38_clustered_sorted.bed <(awk -F' '  '{print $1}' final_output.txt) | head
paste whradjbmi_combined_hg38_clustered_sorted.bed <(awk -F' '  '{print $1}' final_output.txt) >  whradjbmi_combined_hg38_clustered_sorted_final.bed # check the sizes of this....

wc -l whradjbmi_combined_hg38_clustered_sorted_final.bed # 966, correct




################################
### Gene expression check ###
################################

# files here https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/tpms-by-tissue/gene_tpm_2017-06-05_v8_adipose_subcutaneous.gct.gz

# wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/tpms-by-tissue/gene_tpm_2017-06-05_v8_adipose_subcutaneous.gct.gz

# nvm, the .gct file is too weird
#https://www.proteinatlas.org/about/download
wget https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip

# processed in gene_expression.R

sort -k10,10 whradjbmi_combined_hg38_merged_500kbp_genes.bed > whradjbmi_combined_hg38_merged_500kbp_genes_sorted.bed
sort -k1,1 expressed_genes.txt > expressed_genes_sorted.txt


join -a1 -e 'not_expressed' -o '0,2.1' -1 10 -2 1 whradjbmi_combined_hg38_merged_500kbp_genes_sorted.bed expressed_genes_sorted.txt | grep not_expressed
join -a1 -e 'not_expressed'  -t $'\t' -o '1.1,1.2,1.3,1.4,1.5,1.6,0,2.1' -1 10 -2 1 whradjbmi_combined_hg38_merged_500kbp_genes_sorted.bed expressed_genes_sorted.txt > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed.bed


wc -l whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed.bed # 6313
head whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed.bed

chr3 137687690 138687691 chr3 138123718 138132387 A4GNT not_expressed
chr3 137881411 138889511 chr3 138123718 138132387 A4GNT not_expressed
chr12 124338449 125352411 chr12 125065379 125143320 AACS AACS


head whradjbmi_combined_hg38_clustered.bed
chr1    3031614 3031615 rs10909864T:C   -0.0215 w       1.629e-10       chr1    3031036 3033036 peak_85 chr1    3031539 3032321 peak89  7       1
chr1    9292645 9292646 rs114709597C:G  0.0226  w       6.55e-12        chr1    9289341 9298741 peak_228        chr1    9292558 9293469 peak197 6       2
chr1    11805747        11805748        rs13306561A:G   0.017   w       4.194e-08       chr1    11801343        11807743        peak_318        chr1    11805513        11806521        peak264 11      3
chr1    11806126        11806127        rs13306560C:T   -0.0292 w       3.164e-08       chr1    11801343        11807743        peak_318        chr1    11805513        11806521        peak264 11      3
wc -l 706

bedtools intersect -b whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed.bed -a whradjbmi_combined_hg38_clustered_sorted_final.bed -u| wc -l

awk 'BEGIN{OFS="\t"; FS=OFS} {a[$1 OFS $2 OFS $3]=a[$1 OFS $2 OFS $3]","$8} END {for (i in a) {split(i, arr, OFS); if (a[i] ~ /^,/) {a[i] = substr(a[i], 2)}; print arr[1], arr[2], arr[3], a[i]}}' whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed.bed > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed_concatenated.bed

bedtools intersect -wb -u -a whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed_concatenated.bed -b whradjbmi_combined_hg38_merged_500kbp_genes_concatenated.bed | bedtools intersect -a - -b whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed_concatenated.bed | head
bedtools intersect -wo -a whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed_concatenated.bed -b whradjbmi_combined_hg38_merged_500kbp_genes_concatenated.bed | awk '!seen[$4]++' | wc -l




################################
### gene expression check with DAAP6 ###
################################

## formatting
# replace commas with tabs
tr ',' '\t' < DAAP6_RNAseq_results_APs_v_mADs.csv > DAAP6_RNAseq_results_APs_v_mADs.tsv
# remove header
sed -i '1d' DAAP6_RNAseq_results_APs_v_mADs.tsv


## sort before join
# sort -k10,10 whradjbmi_combined_hg38_merged_500kbp_genes.bed > whradjbmi_combined_hg38_merged_500kbp_genes_sorted.bed # already sorted
sort -k1,1 DAAP6_RNAseq_results_APs_v_mADs.tsv > DAAP6_RNAseq_results_APs_v_mADs_sorted.txt


join -a1 -e 'not_expressed' -o '0,2.1' -1 10 -2 1 whradjbmi_combined_hg38_merged_500kbp_genes_sorted.bed DAAP6_RNAseq_results_APs_v_mADs_sorted.txt | grep not_expressed
join -a1 -e 'not_expressed'  -t $'\t' -o '1.1,1.2,1.3,1.4,1.5,1.6,0,2.1,2.6' -1 10 -2 1 whradjbmi_combined_hg38_merged_500kbp_genes_sorted.bed DAAP6_RNAseq_results_APs_v_mADs_sorted.txt > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed



wc -l whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed # 6313
head whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed

awk 'BEGIN{OFS="\t"; FS=OFS} {a[$1 OFS $2 OFS $3]=a[$1 OFS $2 OFS $3]","$8} END {for (i in a) {split(i, arr, OFS); if (a[i] ~ /^,/) {a[i] = substr(a[i], 2)}; print arr[1], arr[2], arr[3], a[i]}}' whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated.bed


## check consistency? between gtex and daap6?

whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed.bed
whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed

# Sort the files based on the 7th field
sort -k7,7 whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed.bed > whradjbmi_combined_GTEX_sorted.bed
sort -k7,7 whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed > whradjbmi_combined_DAAP6_sorted.bed

# Join the files based on the 7th field and check if the 8th fields are identical

join -j 7 whradjbmi_combined_GTEX_sorted.bed whradjbmi_combined_DAAP6_sorted.bed | sort -u -k1,1 | wc -l # 3222 are either shared or not expressed

join -j 7 whradjbmi_combined_GTEX_sorted.bed whradjbmi_combined_DAAP6_sorted.bed | awk '$8 == $16' |  sort -u -k1,1  | wc -l # 840 out of 4188 are not expressed 


################################
### Making tidy cluster file ###
################################

# make a bed file with unique names and the number of snps in the cluster, and the genes close by 500 kbp, and whether they are expressed or not
# in daap6 and in gtex, + how many probes got designed etc


## All genes, regardless of gene expression:
# recover the original region sizes to encounter less pain later on:

bedtools slop -i whradjbmi_combined_hg38_merged_500kbp_genes_concatenated.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b -500000 > whradjbmi_combined_hg38_merged_500kbp_genes_concatenated_OG.bed

## GTEX first
# recover the original region sizes to encounter less pain later on:
bedtools slop -i whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed_concatenated.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b -500000 > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed_concatenated_OG.bed
# change into gene expression column into one that only has the expressed genes, not including the "not expressed"

awk 'BEGIN{OFS="\t"; FS=OFS} {gsub(/,not_expressed/, "", $4); print $1, $2, $3, $4}' whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed_concatenated_OG.bed > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed_concatenated_OG_cleaned.bed # wc -l 369

# count the numbers of expressed gtex genes
awk -F'\t' '{print gsub(/,/, "", $4)}' whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed_concatenated_OG_cleaned.bed > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed_concatenated_gene_count.bed # wc -l 369

## DAAP6 in house next
# recover the original region sizes to encounter less pain later on:
bedtools slop -i whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b -500000 > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated_OG.bed

# change into gene expression column into one that only has the expressed genes, not including the "not expressed"

awk 'BEGIN{OFS="\t"; FS=OFS} {gsub(/,not_expressed/, "", $4); print $1, $2, $3, $4}' whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated_OG.bed > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated_OG_cleaned.bed # wc -l 369

## Biggest chain of intersects you've made in your life

bedtools intersect -a whradjbmi_combined_hg38_clustered_sorted_final.bed -b whradjbmi_combined_hg38_merged_500kbp_genes_concatenated_OG.bed -wb | bedtools intersect -a - -b  whradjbmi_combined_hg38_merged_500kbp_genes_sorted_expressed_concatenated_OG_cleaned.bed  -wb | bedtools intersect -a - -b  whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated_OG_cleaned.bed  -wb  > whradjbmi_long.bed

chr10   102600761       102600762       rs4917976:C:T   -0.0098 W       1.411e-08       chr10   102595043       102605443       peak_5527       chr10   102600699       102601318       peak6436        4       224     chr10_102600761_102600762_rs4917976:C:T_cluster_224_whradjbmi_combined_H3K27ac_ATAC        2       chr10   102600761       102600762       LDB1,ELOVL3,FBXL15,CUEDC2,C10orf95,CYP17A1,PPRC1,NOLC1,NFKB2,PSD,MFSD13A,TRIM8,ARL3,SFXN2,BORCS7,PITX3,GBF1,ACTR1A,SUFU,AS3MT,CNNM2,NT5C2,WBP1L    chr10   102600761       102600762       ACTR1A,ARL3,BORCS7,CUEDC2,GBF1,LDB1,NFKB2,NOLC1,NT5C2,PPRC1,SUFU,TRIM8,WBP1L    chr10   102600761       102600762       ACTR1A,ARL3,BORCS7,CNNM2,CUEDC2,CYP17A1,ELOVL3,FBXL15,GBF1,LDB1,MFSD13A,NFKB2,NOLC1,NT5C2,PITX3,PPRC1,PSD,SFXN2,SUFU,TRIM8,WBP1L


## Pick out the fields of interest + filter by unique cluster + print $0

# interesting fields:


awk 'BEGIN{OFS="\t"; FS=OFS} {print $20, $21, $22, $17, $19, $23, $27, $31 } ' whradjbmi_long.bed > whradjbmi_long_filtered.bed

head 

chr10   102775388       102782989       226     2       FBXL15,CUEDC2,C10orf95,CYP17A1,RPEL1,NFKB2,PSD,MFSD13A,TRIM8,ARL3,SFXN2,BORCS7,GBF1,ACTR1A,SUFU,AS3MT,CNNM2,NT5C2,INA,WBP1L     ACTR1A,ARL3,BORCS7,CUEDC2,GBF1,NFKB2,NT5C2,SUFU,TRIM8,WBP1L        ACTR1A,ARL3,BORCS7,CNNM2,CUEDC2,CYP17A1,FBXL15,GBF1,INA,MFSD13A,NFKB2,NT5C2,PSD,SFXN2,SUFU,TRIM8,WBP1L
chr10   102775388       102782989       226     1       FBXL15,CUEDC2,C10orf95,CYP17A1,RPEL1,NFKB2,PSD,MFSD13A,TRIM8,ARL3,SFXN2,BORCS7,GBF1,ACTR1A,SUFU,AS3MT,CNNM2,NT5C2,INA,WBP1L     ACTR1A,ARL3,BORCS7,CUEDC2,GBF1,NFKB2,NT5C2,SUFU,TRIM8,WBP1L        ACTR1A,ARL3,BORCS7,CNNM2,CUEDC2,CYP17A1,FBXL15,GBF1,INA,MFSD13A,NFKB2,NT5C2,PSD,SFXN2,SUFU,TRIM8,WBP1L
chr10   112973123       113053289       227     2       HABP2,VTI1A,TCF7L2      not_expressed,TCF7L2,VTI1A      HABP2,TCF7L2,VTI1A
chr10   112973123       113053289       227     1       HABP2,VTI1A,TCF7L2      not_expressed,TCF7L2,VTI1A      HABP2,TCF7L2,VTI1A
chr10   112973123       113053289       227     2       HABP2,VTI1A,TCF7L2      not_expressed,TCF7L2,VTI1A      HABP2,TCF7L2,VTI1A
## first reduce the line count by counting/collpasing by the number of probes designed...

awk 'BEGIN{OFS="\t"; FS=OFS} {a[$1 OFS $2 OFS $3]=a[$1 OFS $2 OFS $3]","$5} END {for (i in a) {split(i, arr, OFS); if (a[i] ~ /^,/) {a[i] = substr(a[i], 2)}; print arr[1], arr[2], arr[3], a[i]}}' whradjbmi_long_filtered.bed > whradjbmi_long_filtered_collapsed.bed #| wc -l # 367

awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, gsub(/2/, "", $4)}' whradjbmi_long_filtered_collapsed.bed > whradjbmi_long_probe_2_count.bed
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, gsub(/1/, "", $4)}' whradjbmi_long_filtered_collapsed.bed > whradjbmi_long_probe_1_count.bed

## collapse 6th (gtex) and 7th (daap6) fields correctly

awk 'BEGIN{OFS="\t"; FS=OFS} {!seen[$1 OFS $2 OFS $3 OFS $7]++} END {for (key in seen) {print key}}' whradjbmi_long_filtered.bed > whradjbmi_long_filtered_collapsed_gtex_gene_list.bed
awk 'BEGIN{FS=OFS="\t"} {gsub(/not_expressed,/, "", $4); gsub(/not_expressed/, "", $4); print}' whradjbmi_long_filtered_collapsed_gtex_gene_list.bed > whradjbmi_long_filtered_collapsed_gtex_gene_list_cleaned.bed
# need another one that has non empty 4th columns
awk 'BEGIN{FS=OFS="\t"} {gsub(/not_expressed,/, "NA,", $4); gsub(/not_expressed/, "NA", $4); print}' whradjbmi_long_filtered_collapsed_gtex_gene_list.bed > whradjbmi_long_filtered_collapsed_gtex_gene_list_cleaned_nonempty.bed 


awk 'BEGIN{OFS="\t"; FS=OFS} {!seen[$1 OFS $2 OFS $3 OFS $8]++} END {for (key in seen) {print key}}' whradjbmi_long_filtered.bed > whradjbmi_long_filtered_collapsed_DAAP6_gene_list.bed
awk 'BEGIN{FS=OFS="\t"} {gsub(/not_expressed,/, "", $4); gsub(/not_expressed/, "", $4); print}' whradjbmi_long_filtered_collapsed_DAAP6_gene_list.bed > whradjbmi_long_filtered_collapsed_DAAP6_gene_list_cleaned.bed
awk 'BEGIN{FS=OFS="\t"} {gsub(/not_expressed,/, "NA,", $4); gsub(/not_expressed/, "NA", $4); print}' whradjbmi_long_filtered_collapsed_DAAP6_gene_list.bed > whradjbmi_long_filtered_collapsed_DAAP6_gene_list_cleaned_nonempty.bed 


# then count the number of expressed genes according to gtex

awk -F'\t' 'BEGIN{OFS="\t"} {n=split($4, a, ","); print $1, $2, $3, n}' whradjbmi_long_filtered_collapsed_gtex_gene_list_cleaned.bed > whradjbmi_long_gtex_collapsed_count.bed

# then count the number of expressed genes according to in house daap6 data

awk -F'\t' 'BEGIN{OFS="\t"} {n=split($4, a, ","); print $1, $2, $3, n}' whradjbmi_long_filtered_collapsed_DAAP6_gene_list_cleaned.bed > whradjbmi_long_DAAP6_collapsed_count.bed

# then count the number of snps that map onto the region

awk 'BEGIN{OFS="\t"; FS=OFS} {a[$1 OFS $2 OFS $3]=a[$1 OFS $2 OFS $3]","$4} END {for (i in a) {split(i, arr, OFS); if (a[i] ~ /^,/) {a[i] = substr(a[i], 2)}; print arr[1], arr[2], arr[3], a[i]}}' whradjbmi_long_filtered.bed > whradjbmi_long_filtered_snps_collapsed.bed #| wc -l # 437

awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, gsub(/,/, "", $4)+1}' whradjbmi_long_filtered_snps_collapsed.bed > whradjbmi_long_filtered_collapsed_snp_count.bed


## create the main collapsed file that I'll paste everything to:

awk 'BEGIN{OFS="\t"; FS=OFS} {!seen[$1 OFS $2 OFS $3 OFS $6]++} END {for (key in seen) {print key}}' whradjbmi_long_filtered.bed > whradjbmi_long_filtered_main.bed
awk -F'\t' 'BEGIN{OFS="\t"} {n=split($4, a, ","); print $1, $2, $3, n}' whradjbmi_long_filtered_main.bed > whradjbmi_long_filtered_main_count.bed


# modified files can be intersected instead

bedtools intersect -a whradjbmi_long_filtered_main.bed -b whradjbmi_long_filtered_main_count.bed -wb | bedtools intersect -a - -b whradjbmi_long_filtered_collapsed_snp_count.bed -wb | bedtools intersect -a - -b whradjbmi_long_probe_1_count.bed -wb | bedtools intersect -a - -b whradjbmi_long_probe_2_count.bed -wb | bedtools intersect -a - -b whradjbmi_long_filtered_collapsed_gtex_gene_list_cleaned_nonempty.bed -wb | bedtools intersect -a - -b whradjbmi_long_gtex_collapsed_count.bed -wb | bedtools intersect -a - -b whradjbmi_long_filtered_collapsed_DAAP6_gene_list_cleaned_nonempty.bed -wb | bedtools intersect -a - -b whradjbmi_long_DAAP6_collapsed_count.bed -wb |  awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $8, $12, $16, $20, $24, $28, $32, $36 }' > whradjbmi_combined_hg38_clustered_sorted_final_tidy_OG_regions.bed 

# sort the file I guess that

bedtools sort -i whradjbmi_combined_hg38_clustered_sorted_final_tidy_OG_regions.bed > whradjbmi_combined_hg38_clustered_sorted_final_tidy_OG_regions_pretty.bed

cat whradjbmi_combined_hg38_clustered_sorted_final_tidy_OG_regions_pretty.bed > whradjbmi_combined_hg38_clustered_sorted_final_tidy_OG_regions_pretty.tsv
#awk 'BEGIN{OFS="\t"; FS=OFS} {print NF}' whradjbmi_combined_hg38_clustered_sorted_final_tidy_OG_regions.bed | sort | uniq -c

## next steps

# number of regions within 20kbp of each other

# merge the snps within 50+/- kbp of each other
bedtools merge -d 200000 -i whradjbmi_combined_hg38_clustered_sorted_final_tidy_OG_regions_pretty.bed > whradjbmi_combined_hg38_regions_close.bed
wc -l 243 
# add -c and -o options 


# number of genes close by

## filter t2d_oligo_info.txt density score < 40

################################
### gene set pathways etc check
################################


cat HALLMARK_ADIPOGENESIS.v2023.2.Hs.grp KEGG* REACTOME* WP_* | wc -l #746

cat HALLMARK_ADIPOGENESIS.v2023.2.Hs.grp KEGG* REACTOME* WP_* | grep -v '_'| wc -l #25 

cat HALLMARK_ADIPOGENESIS.v2023.2.Hs.grp KEGG* REACTOME* WP_* |  grep -v '_' > gene_sets.txt


sort -u gene_sets.txt | wc -l # 522
sort -u gene_sets.txt > gene_sets_sorted.txt

# binary yes/no

## sort before join
# sort -k10,10 whradjbmi_combined_hg38_merged_500kbp_genes.bed > whradjbmi_combined_hg38_merged_500kbp_genes_sorted.bed # already sorted

join -a1 -e 'not_in_geneset' -o '0,2.1' -1 10 -2 1 whradjbmi_combined_hg38_merged_500kbp_genes_sorted.bed gene_sets_sorted.txt | grep -v not_in_geneset | head


join -a1 -e 'not_in_geneset'  -t $'\t' -o '1.1,1.2,1.3,1.4,1.5,1.6,0,2.1' -1 10 -2 1 whradjbmi_combined_hg38_merged_500kbp_genes_sorted.bed gene_sets_sorted.txt > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_genesets.bed



wc -l whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_genesets.bed # 6313
head whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_genesets.bed

awk 'BEGIN{OFS="\t"} {key = $1 OFS $2 OFS $3; if (key in a) a[key] = a[key] "," $8; else a[key] = $8} END {for (key in a) {split(key, arr, OFS); print arr[1], arr[2], arr[3], a[key]}}' whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_genesets.bed > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_genesets_concatenated.bed

bedtools slop -i whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_genesets_concatenated.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b -500000 > whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_genesets_concatenated_OG.bed

awk 'BEGIN{OFS="\t"; FS=OFS} {gsub(/,not_in_geneset/, "", $4); print $1, $2, $3, $4}' whradjbmi_combined_hg38_merged_500kbp_genes_sorted_DAAP6_genesets_concatenated_OG.bed > whradjbmi_combined_hg38_merged_500kbp_genesets_OG_cleaned.bed # wc -l 441


awk 'BEGIN{FS=OFS="\t"} {gsub(/not_in_geneset,/, "", $4); gsub(/not_in_geneset/, "", $4); print}' whradjbmi_combined_hg38_merged_500kbp_genesets_OG_cleaned.bed > whradjbmi_combined_hg38_pathways.bed
awk 'BEGIN{FS=OFS="\t"} {gsub(/not_in_geneset,/, "NA,", $4); gsub(/not_in_geneset/, "NA", $4); print}' whradjbmi_combined_hg38_merged_500kbp_genesets_OG_cleaned.bed > whradjbmi_combined_hg38_pathways_nonempty.bed 

#count

awk -F'\t' 'BEGIN{OFS="\t"} {n=split($4, a, ","); print $1, $2, $3, n}' whradjbmi_combined_hg38_pathways.bed > whradjbmi_long_pathways_collapsed_count.bed


bedtools intersect -a whradjbmi_combined_hg38_clustered_sorted_final_tidy_OG_regions.bed -b whradjbmi_long_pathways_collapsed_count.bed -wb |  awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $16}' > whradjbmi_combined_hg38_tidy_genesets.bed

cat whradjbmi_combined_hg38_tidy_genesets.bed > whradjbmi_combined_hg38_tidy_genesets.tsv



################################
### eQTLS ###
################################

# files here /novo/projects/departments/nnrco/genetic_department/share/for_GATL/collabs/Carol
# move some over?

cp /novo/projects/departments/nnrco/genetic_department/share/for_GATL/collabs/Carol/Adipose_Subcutaneous_snps.hg19.bed Adipose_Subcutaneous_snps.hg19.bed 
cp /novo/projects/departments/nnrco/genetic_department/share/for_GATL/collabs/Carol/AdipoExpress_snps.hg19.bed AdipoExpress_snps.hg19.bed
cp /novo/projects/departments/nnrco/genetic_department/share/for_GATL/collabs/Carol/FUSION_ge_adipose_naive_snps.hg19.bed FUSION_ge_adipose_naive_snps.hg19.bed
cp /novo/projects/departments/nnrco/genetic_department/share/for_GATL/collabs/Carol/TwinsUK_ge_fat_snps.hg19.bed TwinsUK_ge_fat_snps.hg19.bed



# all in hg19 so need to convert to hg38

/nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/liftOver AdipoExpress_snps.hg19.bed /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/hg19ToHg38.over.chain.gz AdipoExpress_snps_lifted.bed  unmapped_AdipoExpress_snps.bed
/nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/liftOver Adipose_Subcutaneous_snps.hg19.bed /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/hg19ToHg38.over.chain.gz Adipose_Subcutaneous_snps_lifted.bed  unmapped_Adipose_Subcutaneous_snps.bed
/nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/liftOver FUSION_ge_adipose_naive_snps.hg19.bed /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/hg19ToHg38.over.chain.gz FUSION_ge_adipose_naive_snps_lifted.bed  unmapped_FUSION_ge_adipose_naive_snps.bed
/nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/liftOver TwinsUK_ge_fat_snps.hg19.bed /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/hg19ToHg38.over.chain.gz TwinsUK_ge_fat_snps_lifted.bed  unmapped_TwinsUK_ge_fat_snps.bed


# 
bedtools intersect -a Adipose_Subcutaneous_snps_lifted.bed -b whradjbmi_combined_hg38_clustered.bed > Adipose_Subcutaneous_snps_lifted_intersected_whradjbmi.bed # | head # wc -l #10 
bedtools intersect -a AdipoExpress_snps_lifted.bed -b whradjbmi_combined_hg38_clustered.bed > AdipoExpress_snps_lifted_intersected_whradjbmi.bed  # | head # wc -l #25
bedtools intersect -a FUSION_ge_adipose_naive_snps_lifted.bed -b whradjbmi_combined_hg38_clustered.bed > FUSION_ge_adipose_naive_snps_lifted_intersected_whradjbmi.bed  #| head # wc -l #11
bedtools intersect -a TwinsUK_ge_fat_snps_lifted.bed -b whradjbmi_combined_hg38_clustered.bed > TwinsUK_ge_fat_snps_lifted_intersected_whradjbmi.bed # | head # wc -l #16

################################
### TSSs ###
################################


wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz

/nfs_home/users/czjs/miniconda3/envs/fullcap/bin/python TSS_generator.py

head gene_tss_coordinates.bed 

# expand by 1 and 2 kbp

bedtools slop -i gene_tss_coordinates.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b -500 > gene_tss_coordinates_1kbp.bed
bedtools slop -i gene_tss_coordinates.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b -1000 > gene_tss_coordinates_2kbp.bed

# intersect -v

bedtools intersect -v -b gene_tss_coordinates_1kbp.bed -a whradjbmi_combined_hg38_clustered_sorted_final.bed >  whradjbmi_combined_hg38_1kbp_TSS.bed #wc -l #680
cat whradjbmi_combined_hg38_1kbp_TSS.bed > whradjbmi_combined_hg38_1kbp_TSS.tsv

bedtools intersect -v -b gene_tss_coordinates_2kbp.bed -a whradjbmi_combined_hg38_clustered_sorted_final.bed > whradjbmi_combined_hg38_2kbp_TSS.bed #wc -l #573
cat whradjbmi_combined_hg38_2kbp_TSS.bed > whradjbmi_combined_hg38_2kbp_TSS.tsv



################################
### Finemapping overlap ###
################################

# file is here: /novo/projects/departments/nnrco/genetic_department/share/for_LWMR/tde_hepatocytes_v1/nngene_internal.ukbb_EUR_LGLM_2023.whradjbmi_sexcombined/cs_annotated.tsv

cp /novo/projects/departments/nnrco/genetic_department/share/for_LWMR/tde_hepatocytes_v1/nngene_internal.ukbb_EUR_LGLM_2023.whradjbmi_sexcombined/cs_annotated.tsv whradjbmi_combined_hg38_finemapping.bed

#sed -i '1d' whradjbmi_combined_hg38_finemapping.bed # remove 1st row, once. 


### -------------------------------- ###
### dealing with problematic oligos
### -------------------------------- ###

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/main_attempt_120824/25829_1/oligo_info.txt whr_oligos_round1.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/main_attempt_120824/25829_1/leftright_coords_round1.bed whr_oligos_round1_redesign.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/main_attempt_120824/25829_2/oligo_info.txt whr_oligos_round2.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/main_attempt_120824/25829_2/leftright_coords_round2.bed whr_oligos_round2_redesign.bed

# processed in processing_oligos.R, output is whr_pick_new_snps.bed
# further processed in processing_snps_whradjbmi.R


### building off of whradjbmi_combined_hg38_selected_oligos.bed

mkdir manual_left_right 

cp whradjbmi_combined_hg38_selected_oligos.bed manual_left_right/oligo_info.txt

cd manual_left_right

autogen_probes_left_right_modded.R 

# manually ran through the R script because it's a pain

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh whradjbmi_round2 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right/leftright_coords_round1.bed
squeue --user=czjs --iterate=30_seconds

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/whradjbmi_round2/oligo_info.txt whradjbmi_round2_oligo_info.txt

autogen_probes_left_right_modded.R

processing_snps_whradjbmi.R

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh whradjbmi_round3 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right/whradjbmi_combined_hg38_selected_snps_round2_capready.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/whradjbmi_round3/oligo_info.txt whradjbmi_reselected_oligo_info.txt

autogen_probes_left_right_modded.R

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh whradjbmi_reselected_round2 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right/leftright_coords_reselected_round1.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/whradjbmi_reselected_round2/oligo_info.txt whradjbmi_reselected_oligo_info_round2.txt



### -------------------------------- ###
### filtering out SNPs that are exonic
### -------------------------------- ###


bedtools intersect -v -a whradjbmi_combined_hg38_clustered_sorted_final.bed -b exon_coordinates.bed | wc -l #654
bedtools intersect -v -a whradjbmi_combined_hg38_2kbp_TSS.bed -b exon_coordinates.bed | wc -l #484 from 573


bedtools intersect -v -a whradjbmi_combined_hg38_2kbp_TSS.bed -b exon_coordinates.bed > whradjbmi_combined_hg38_clustered_sorted_final_no_exons_no_TSS.bed
# number of regions left: 233

wc -l whradjbmi_combined_hg38_clustered_sorted_final.bed #966, 482 snps removed

# save the intersection and have a look on ucsc?

bedtools intersect -a whradjbmi_combined_hg38_clustered_sorted_final.bed -b exon_coordinates.bed > whradjbmi_combined_hg38_clustered_sorted_final_exonic_check.bed

# rerun capsequm with intronic no tss snps

awk 'BEGIN{OFS="\t"; FS=OFS} {print $1, $2, $3, $1"_"$2"_"$3"_"$4"_cluster_"$17"_whradjbmi_combined_H3K27ac_ATAC" }' whradjbmi_combined_hg38_clustered_sorted_final_no_exons_no_TSS.bed > whradjbmi_combined_hg38_clustered_sorted_final_no_exons_no_TSS_4_col.bed


sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh whr_no_exons_no_TSS /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/whradjbmi_combined_hg38_clustered_sorted_final_no_exons_no_TSS_4_col.bed
# and then you need to filter to get the selected oligos...use processing_snps

### REPEAT/ REDO

### building off of whradjbmi_combined_hg38_selected_oligos.bed

mkdir manual_left_right_whr

cd manual_left_right_whr

autogen_probes_left_right_modded.R whradjbmi_combined_hg38_selected_oligos.bed

# manually ran through the R script because it's a pain

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh whradjbmi_round2 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_whr/leftright_coords_round1.bed
squeue --user=czjs --iterate=30_seconds

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/whradjbmi_round2/oligo_info.txt whradjbmi_round2_oligo_info.txt

autogen_probes_left_right_modded.R whradjbmi_round2_oligo_info.txt

processing_snps_whradjbmi.R

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh whradjbmi_round3 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_whr/whradjbmi_combined_hg38_selected_snps_round2_capready.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/whradjbmi_round3/oligo_info.txt whradjbmi_reselected_oligo_info.txt

autogen_probes_left_right_modded.R whradjbmi_reselected_oligo_info.txt

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh whradjbmi_reselected_round2 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_whr/leftright_coords_reselected_round1.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/whradjbmi_reselected_round2/oligo_info.txt whradjbmi_reselected_oligo_info_round2.txt

autogen_probes_left_right_modded.R whradjbmi_reselected_oligo_info_round2.txt

/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/processing_oligos_whr.R

