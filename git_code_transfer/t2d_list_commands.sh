
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
zcat T2D.Suzuki.Nature.EUR_hg19_harmonized.tsv.gz | wc -l #19303704
zcat T2D.Suzuki.Nature.EUR_hg19_harmonized.tsv.gz | awk -F'\t' '$11 < 0.00000005'| gzip >  /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/T2D_EUR_filtered.gz 
zcat /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/T2D_EUR_filtered.gz | wc -l # 64948
################################
### unzipped ###
################################

cd /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps

gunzip -k T2D_EUR_filtered.gz

### cut out columns I think are interesting/useful? and convert to a bed-ish file

cat T2D_EUR_filtered | awk -F'\t' -v OFS='\t' '{print "chr"$4, $5, $5+1, $3, $9,  $1, $11}' >  T2D_EUR_extra_filtered.bed

################################
### liftover  ###
################################

/nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/liftOver T2D_EUR_extra_filtered.bed /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/hg19ToHg38.over.chain.gz T2D_EUR_hg38_extra_filtered.bed unmapped_T2D_EUR_hg38_extra_filtered.bed 

# #Deleted in new
# chr6    51097410        51097411        rs9463672       -0.0236 T       1.118e-09
# #Deleted in new
# chr9    136142203       136142204       rs514659        -0.0425 T       7.174e-11
# #Deleted in new
# chr9    136143120       136143121       rs613534        -0.0389 T       1.508e-09
# #Deleted in new
# chr9    136143121       136143122       rs543968        -0.0389 T       1.546e-09
# #Deleted in new
# chr9    136144284       136144285       rs597988        0.0421  T       1.074e-10

################################
### get the updated rsids###
################################

# using this: https://genome.sph.umich.edu/wiki/LiftOver#Lift_dbSNP_rs_numbers
# getting the rsids into a format that the python script needs it to be
zcat T2D_EUR_filtered.gz| awk -F'\t' '{print $3}' | sed 's/rs//g' | awk -F : '{print $1}' > T2D_EUR_filtered_rs.txt

# running the script
python2 LiftRsNumber.py T2D_EUR_filtered_rs.txt > T2D_EUR_filtered_lifted_rs.txt

# get the beginning of the rsid
cat T2D_EUR_filtered_lifted_rs.txt | awk -F'\t' '{print "rs"$2}' > rs_appended.txt

# remove missing rsids

awk 'NR % 2 == 0 {print $4}' unmapped_T2D_EUR_hg38_extra_filtered.bed > rsids_to_remove.txt

grep -v -f rsids_to_remove.txt rs_appended.txt > full_lifted_rsid_cleaned.txt

# replace the old column of rsids with the new

awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[NR]=$1;next}{$4=a[FNR]}1' full_lifted_rsid_cleaned.txt T2D_EUR_hg38_extra_filtered.bed > T2D_EUR_hg38_extra_filtered_lifted.bed

################################
### Intersect with epigenetics ###
################################

#just moving them over - it probably isn't best practice but the files seem relatively small?
cp /nfs_home/projects/departments/nnrco/genetic_department/share/GATL_GMGZ/CRE_AdiposeTissue/H3K27ac-ChromatinStates_AdiposeTissue_GMGZ.bed H3K27ac_peaks.bed

cp /nfs_home/projects/departments/nnrco/genetic_department/share/GATL_GMGZ/GSE178794_adiposeTissue_consensusPeaks-from11samples_unionPeaksInAtLeast3samples_b38.bed atac.bed

# bedtools intersect twice to keep information from the files, just in case
bedtools intersect -a T2D_EUR_hg38_extra_filtered_lifted.bed -b H3K27ac_peaks.bed -wa| bedtools intersect -b H3K27ac_peaks.bed -a - -wb > T2D_EUR_hg38_extra_H3K27ac.bed

wc -l 
5346

bedtools intersect -a T2D_EUR_hg38_extra_H3K27ac.bed -b atac.bed -wa| bedtools intersect -b atac.bed -a - -wb > T2D_EUR_hg38_extra_H3K27ac_ATAC.bed
wc -l 
1008

################################
### make H3K27ac only and ATAC only interesected files ###
################################

bedtools intersect -a T2D_EUR_hg38_extra_filtered_lifted.bed -b H3K27ac_peaks.bed -wa| bedtools intersect -b H3K27ac_peaks.bed -a - -wb > T2D_EUR_hg38_extra_H3K27ac_only.bed
wc -l 
5346

bedtools intersect -a T2D_EUR_hg38_extra_filtered_lifted.bed -b atac.bed -wa| bedtools intersect -b atac.bed -a - -wb > T2D_EUR_hg38_extra_ATAC_only.bed
wc -l 
1300

# set aside for future use, if not enough regions for capture c


################################
### merge and clump ###
################################

# move a list of hg38 genes over 
cp /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/gencode.v26.GRCh38.genes_for_gene.protein_coding.txt hg38_all_genes.bed

# merge the snps within 50+/- kbp of each other
bedtools merge -d 50000 -i T2D_EUR_hg38_extra_H3K27ac_ATAC.bed > T2D_EUR_hg38_merged.bed
wc -l 442

# cluster the snps within 50+/- kbp of each other
# gives me the identity of the merged snps in each cluster
bedtools cluster -d 50000 -i T2D_EUR_hg38_extra_H3K27ac_ATAC.bed > T2D_EUR_hg38_clustered.bed
wc -l 1008

################################
### assess complexity by looking at how many genes are close by ###
################################

# uni supervisor suggested I find the genes close by the clusters, not the snps

bedtools slop -i T2D_EUR_hg38_merged.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b 500000 > T2D_EUR_hg38_merged_500kbp.bed

# +/- 500kbp, so a million after expr 

bedtools intersect -wa -wb -a T2D_EUR_hg38_merged_500kbp.bed -b hg38_all_genes.bed > T2D_EUR_hg38_merged_500kbp_genes.bed
wc -l 6866


# check the length of the unique snps
cat T2D_EUR_hg38_merged_500kbp_genes.bed | awk '{print $1,$2,$3}' | uniq | wc -l
441

# I think one of them simply didn't intersect at all, 441 = 442 - 1
# need to find which one that is...
bedtools intersect -wa -wb -v -a T2D_EUR_hg38_merged_500kbp.bed -b hg38_all_genes.bed > T2D_EUR_hg38_merged_500kbp_nogenes.bed 
# empty 4th col, so I'll just add a column with NA
sed -i 's/$/\tNA/' T2D_EUR_hg38_merged_500kbp_nogenes.bed

# yep one of the snps didn't have intersections +/- 500kbp genes

# unique, process columns such that you get one column with all the intersected genes in one place

awk 'BEGIN{OFS="\t"; FS=OFS} {a[$1 OFS $2 OFS $3]=a[$1 OFS $2 OFS $3]","$10} END {for (i in a) {split(i, arr, OFS); if (a[i] ~ /^,/) {a[i] = substr(a[i], 2)}; print arr[1], arr[2], arr[3], a[i]}}' T2D_EUR_hg38_merged_500kbp_genes.bed > T2D_EUR_hg38_merged_500kbp_genes_concatenated.bed
wc -l 441

# append the missing snp/cluster

cat T2D_EUR_hg38_merged_500kbp_nogenes.bed >> T2D_EUR_hg38_merged_500kbp_genes_concatenated.bed
wc -l 442 # perfect

################################
### Run capsequm on snps ###
################################
# need to make a bed file with only 4 columns: chr start stop unique_name

awk 'BEGIN{OFS="\t"; FS=OFS} {print $1, $2, $3, $1"_"$2"_"$3"_"$4"_cluster_"$17"_T2D_EUR_H3K27ac_ATAC" }' T2D_EUR_hg38_clustered.bed | head
awk 'BEGIN{OFS="\t"; FS=OFS} {print $1, $2, $3, $1"_"$2"_"$3"_"$4"_cluster_"$17"_T2D_EUR_H3K27ac_ATAC" }' T2D_EUR_hg38_clustered.bed > T2D_EUR_hg38_clustered_4col.bed

cd /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh T2D_EUR_hg38_extra_H3K27ac_ATAC_final /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/T2D_EUR_hg38_clustered_4col.bed
#/fsx/home_dirs/czjs/doctoral_work/repeat-annotation/RepeatMasker/RepeatMasker -species "Homo sapiens" -pa 4 /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/T2D_EUR_hg38_extra_H3K27ac_ATAC_3/oligo_seqs.fa


################################
### Making tidy snp file ###
################################

## add to the snp files whether or not an oligo could be designed for it

cd /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/T2D_EUR_hg38_extra_H3K27ac_ATAC_final

wc -l oligo_info.txt
1515

wc -l /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/T2D_EUR_hg38_clustered_4col.bed
1008

# capsequm designs 2 probes per region, 2*1008 = 2016 , 2016-1515 = 501, 501/2 =~ 251 regions either don't have any probe, or have half the probes etc.
# find which ones have probes and which don't, turn that to a binary yes no column


## Make key
cd /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps

awk 'BEGIN{OFS="\t"; FS=OFS} {print $1"_"$2"_"$3"_"$4"_cluster_"$17"_T2D_EUR_H3K27ac_ATAC" }' T2D_EUR_hg38_clustered.bed > temp.txt
paste T2D_EUR_hg38_clustered.bed temp.txt > T2D_EUR_hg38_clustered_keyed.bed
rm temp.txt

## sort before join
# 18th column is the key
# awk '{print $18}' T2D_EUR_hg38_clustered_keyed.bed| head # quick check

sort -k18,18 T2D_EUR_hg38_clustered_keyed.bed > T2D_EUR_hg38_clustered_sorted.bed

# 13th column is the key
# awk '{print $13}' oligo_info.txt| head # quick check
cp /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/T2D_EUR_hg38_extra_H3K27ac_ATAC_final/oligo_info.txt t2d_oligo_info.txt
sort -k13,13 t2d_oligo_info.txt > t2d_oligo_info_sorted.txt
sed 's/,$//' t2d_oligo_info_sorted.txt > t2d_oligo_info_sorted_no_comma.txt
sed -i '1d' t2d_oligo_info_sorted_no_comma.txt

## join 

join -1 18 -2 13 -t $'\t' -a 1 -o 1.18 T2D_EUR_hg38_clustered_sorted.bed t2d_oligo_info_sorted_no_comma.txt | uniq -c > joined_counts.txt

# count the number of regions with 1, 2, or none probes
awk '{if ($1 == 1) print $0, "one"; else if ($1 == 2) print $0, "both"; else print $0, "none";}' joined_counts.txt > final_output.txt
# this is overkill, just use the first column
grep "one" final_output.txt | wc -l
# 251

paste T2D_EUR_hg38_clustered_sorted.bed <(awk -F' '  '{print $1}' final_output.txt) | head
paste T2D_EUR_hg38_clustered_sorted.bed <(awk -F' '  '{print $1}' final_output.txt) >  T2D_EUR_hg38_clustered_sorted_final.bed

wc -l T2D_EUR_hg38_clustered_sorted_final.bed # 1008, correct



################################
### Gene expression check ###
################################

# files here https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/tpms-by-tissue/gene_tpm_2017-06-05_v8_adipose_subcutaneous.gct.gz

# wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/tpms-by-tissue/gene_tpm_2017-06-05_v8_adipose_subcutaneous.gct.gz

# nvm, the .gct file is too weird
#https://www.proteinatlas.org/about/download
wget https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip

# processed in gene_expression.R

sort -k10,10 T2D_EUR_hg38_merged_500kbp_genes.bed > T2D_EUR_hg38_merged_500kbp_genes_sorted.bed
sort -k1,1 expressed_genes.txt > expressed_genes_sorted.txt


join -a1 -e 'not_expressed' -o '0,2.1' -1 10 -2 1 T2D_EUR_hg38_merged_500kbp_genes_sorted.bed expressed_genes_sorted.txt | grep not_expressed
join -a1 -e 'not_expressed'  -t $'\t' -o '1.1,1.2,1.3,1.4,1.5,1.6,0,2.1' -1 10 -2 1 T2D_EUR_hg38_merged_500kbp_genes_sorted.bed expressed_genes_sorted.txt > T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed.bed


wc -l T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed.bed # 6866
head T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed.bed

chr1    32317982        33317983        chr1    33306766        33321098        A3GALT2 not_expressed
chr3    137848480       138889511       chr3    138123718       138132387       A4GNT   not_expressed
chr3    138109174       139109175       chr3    138123718       138132387       A4GNT   not_expressed

head T2D_EUR_hg38_clustered.bed
chr1    2213349 2213350 rs262695        -0.0256 T       4.893e-10       chr1    2211961 2215361 peak_66 chr1    2212445 2213351 peak61  9       1
chr1    2232129 2232130 rs263533        -0.0219 T       9.208e-09       chr1    2226361 2244961 peak_67 chr1    2231708 2232229 peak64  6       1
chr1    2234479 2234480 rs12061341      0.0422  T       2.09e-09        chr1    2226361 2244961 peak_67 chr1    2234119 2234783 peak65  10      1
wc -l 1008

# bedtools intersect -b T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed.bed -a T2D_EUR_hg38_clustered_sorted_final.bed -u| head

awk 'BEGIN{OFS="\t"; FS=OFS} {a[$1 OFS $2 OFS $3]=a[$1 OFS $2 OFS $3]","$8} END {for (i in a) {split(i, arr, OFS); if (a[i] ~ /^,/) {a[i] = substr(a[i], 2)}; print arr[1], arr[2], arr[3], a[i]}}' T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed.bed > T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed_concatenated.bed



################################
### check overlap between whr and t2d regions ###
################################


# on the snp level

bedtools intersect -a T2D_EUR_hg38_clustered_sorted_final.bed -b T2D_EUR_hg38_clustered_sorted_final.bed | wc -l # 143

bedtools intersect -b T2D_EUR_hg38_clustered_sorted_final.bed -a T2D_EUR_hg38_clustered_sorted_final.bed > T2D_EUR_t2d_intersect_snp.bed
# swapped order because need the rs base change information

# on the region/cluster level? 


bedtools intersect -a T2D_EUR_hg38_merged.bed -b T2D_EUR_hg38_merged.bed | wc -l # 84
bedtools intersect -a T2D_EUR_hg38_merged.bed -b T2D_EUR_hg38_merged.bed > T2D_EUR_t2d_intersect_regions.bed

# the 2 above should have perfect overlap right
bedtools intersect -a T2D_EUR_t2d_intersect_regions.bed -b T2D_EUR_t2d_intersect_snp.bed | wc -l # 143
# perfect

# check the beta values of the snps in the intersected regions
# all in beta_effect.r 
# but the hist shows a fairly symmetrical curve centered on 0, +/- 0.15

################################
### check if coding snps are synonymous ###
################################

# snpnexus requires the file formatting to be like this

Type	Id	Position	Alelle1	Allele2	Strand
Chromosome	1	942451	T	C	1
Chromosome	3	9810376	-	GAT	1
Chromosome	7	25226951	TA	GTT	1


# let's try with T2D_EUR_t2d_intersect_snp first

awk '{gsub(/[^0-9]/, "", $1); split($4,a,":"); print "Chromosome", $1, $2, a[2], a[3], "1" }' T2D_EUR_t2d_intersect_snp.bed > T2D_EUR_t2d_intersect_snp_nexus_check.txt

################################
### file with only non-coding variants ###
################################

# remove if within 5 kbp of a gene 


bedtools slop -i hg38_all_genes.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b 5000 > hg38_all_genes_extended.bed


bedtools intersect -v -a T2D_EUR_t2d_intersect_snp.bed -b hg38_all_genes_extended.bed | wc -l #39

bedtools intersect -v -a T2D_EUR_t2d_intersect_regions.bed -b hg38_all_genes_extended.bed | wc -l #12

bedtools intersect -v -a T2D_EUR_hg38_merged.bed -b hg38_all_genes_extended.bed | wc -l #61

bedtools intersect -v -a T2D_EUR_hg38_merged.bed -b hg38_all_genes_extended.bed | wc -l #60


################################
### include intronic only, intronic + non-synonymous too ###
################################

# after uploading the nexus check file to: https://www.snp-nexus.org/v4/
# download the file under genomic coordinates and overlapping, and the second tab overlapped or nearest genes

awk 'BEGIN{OFS="\t"; FS=OFS} !seen[$1]++' near_gens_T2D_EUR_t2d_intersect_snp_nexus_check_3.txt | wc -l # 144, meaning 143 minus header, which is right

# unique
awk 'BEGIN{OFS="\t"; FS=OFS} !seen[$1]++' near_gens_T2D_EUR_t2d_intersect_snp_nexus_check_3.txt > T2D_EUR_t2d_intersect_nexus_unique.txt
#awk ' BEGIN{OFS="\t"; FS=OFS} {print $2, $3, $6}' T2D_EUR_t2d_intersect_nexus_unique.txt | head

# counts
wc -l T2D_EUR_t2d_intersect_nexus_unique.txt # 144
awk 'BEGIN{OFS="\t"; FS=OFS} $6 ~ /intronic/ {print}' T2D_EUR_t2d_intersect_nexus_unique.txt | wc -l # 90
awk 'BEGIN{OFS="\t"; FS=OFS} $6 ~ / syn/ {print}' T2D_EUR_t2d_intersect_nexus_unique.txt | wc -l # 1, the space is important so as to not capture the nonsyn
awk 'BEGIN{OFS="\t"; FS=OFS} $6 ~ /nonsyn/ {print}' T2D_EUR_t2d_intersect_nexus_unique.txt | wc -l # 1

### repeat for the other files e.g. whr only
awk '{gsub(/[^0-9]/, "", $1); split($4,a,":"); print "Chromosome", $1, $2, a[2], a[3], "1" }' T2D_EUR_hg38_clustered_sorted_final.bed > T2D_EUR_snp_nexus_check.txt

# unique
awk 'BEGIN{OFS="\t"; FS=OFS} !seen[$1]++' near_gens_T2D_EUR_snp_nexus_check.txt > T2D_EUR_intersect_nexus_unique.txt

# counts 

wc -l T2D_EUR_intersect_nexus_unique.txt # 966
awk 'BEGIN{OFS="\t"; FS=OFS} $6 ~ /intronic/ {print}' T2D_EUR_intersect_nexus_unique.txt | wc -l # 602
awk 'BEGIN{OFS="\t"; FS=OFS} $6 ~ / syn/ {print}' T2D_EUR_intersect_nexus_unique.txt | wc -l # 14, the space is important so as to not capture the nonsyn
awk 'BEGIN{OFS="\t"; FS=OFS} $6 ~ /nonsyn/ {print}' T2D_EUR_intersect_nexus_unique.txt | wc -l # 32


### t2d only
# awkardly paste back parts of the rsid from og main file: T2D_EUR_filtered

# Sort the files based on the columns to be joined
sort -k3,3 T2D_EUR_filtered > T2D_EUR_filtered_sorted.txt
sort -k4,4 T2D_EUR_hg38_clustered_sorted_final.bed > T2D_EUR_hg38_clustered_sorted_for_rsid.txt

# Join the files and extract the desired fields
join -1 3 -2 4 T2D_EUR_filtered_sorted.txt T2D_EUR_hg38_clustered_sorted_for_rsid.txt | wc -l # 1008, good

join -1 3 -2 4 T2D_EUR_filtered_sorted.txt T2D_EUR_hg38_clustered_sorted_for_rsid.txt | awk '{print "Chromosome", $4, $5, $6, $7 ,"1"}' > t2d_snp_nexus_check.txt


# unique
awk 'BEGIN{OFS="\t"; FS=OFS} !seen[$1]++' near_gens_t2d_snp_nexus_check.txt > t2d_intersect_nexus_unique.txt

# counts
wc -l t2d_intersect_nexus_unique.txt # 988, something is not quite right with this number... 
# snpnexus seems to only clock 987 variants from the list I gave it

awk 'BEGIN{OFS="\t"; FS=OFS} $6 ~ /intronic/ {print}' t2d_intersect_nexus_unique.txt | wc -l # 644
awk 'BEGIN{OFS="\t"; FS=OFS} $6 ~ / syn/ {print}' t2d_intersect_nexus_unique.txt | wc -l # 5, the space is important so as to not capture the nonsyn
awk 'BEGIN{OFS="\t"; FS=OFS} $6 ~ /nonsyn/ {print}' t2d_intersect_nexus_unique.txt | wc -l # 12 

################################
### gene expression check with DAAP6 ###
################################

## formatting
# replace commas with tabs
tr ',' '\t' < DAAP6_RNAseq_results_APs_v_mADs.csv > DAAP6_RNAseq_results_APs_v_mADs.tsv
# remove header
sed -i '1d' DAAP6_RNAseq_results_APs_v_mADs.tsv


## sort before join
# sort -k10,10 T2D_EUR_hg38_merged_500kbp_genes.bed > T2D_EUR_hg38_merged_500kbp_genes_sorted.bed # already sorted
sort -k1,1 DAAP6_RNAseq_results_APs_v_mADs.tsv > DAAP6_RNAseq_results_APs_v_mADs_sorted.txt


join -a1 -e 'not_expressed' -o '0,2.1' -1 10 -2 1 T2D_EUR_hg38_merged_500kbp_genes_sorted.bed DAAP6_RNAseq_results_APs_v_mADs_sorted.txt | grep not_expressed
join -a1 -e 'not_expressed'  -t $'\t' -o '1.1,1.2,1.3,1.4,1.5,1.6,0,2.1,2.6' -1 10 -2 1 T2D_EUR_hg38_merged_500kbp_genes_sorted.bed DAAP6_RNAseq_results_APs_v_mADs_sorted.txt > T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed



wc -l T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed # 6866
head T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed

awk 'BEGIN{OFS="\t"; FS=OFS} {a[$1 OFS $2 OFS $3]=a[$1 OFS $2 OFS $3]","$8} END {for (i in a) {split(i, arr, OFS); if (a[i] ~ /^,/) {a[i] = substr(a[i], 2)}; print arr[1], arr[2], arr[3], a[i]}}' T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed > T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated.bed


## check consistency? between gtex and daap6?

T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed.bed
T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed

# Sort the files based on the 7th field
sort -k7,7 T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed.bed > T2D_GTEX_sorted.bed
sort -k7,7 T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed.bed > T2D_DAAP6_sorted.bed

# Join the files based on the 7th field and check if the 8th fields are identical

join -j 7 T2D_GTEX_sorted.bed T2D_DAAP6_sorted.bed | sort -u -k1,1 | wc -l # 4188 are either shared or not expressed

join -j 7 T2D_GTEX_sorted.bed T2D_DAAP6_sorted.bed | awk '$8 == $16' |  sort -u -k1,1  | wc -l # 1142 out of 4188 are not expressed 



################################
### Making tidy cluster file ###
################################

# make a bed file with unique names and the number of snps in the cluster, and the genes close by 500 kbp, and whether they are expressed or not
# in daap6 and in gtex, + how many probes got designed etc


## All genes, regardless of gene expression:
# recover the original region sizes to encounter less pain later on:

bedtools slop -i T2D_EUR_hg38_merged_500kbp_genes_concatenated.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b -500000 > T2D_EUR_hg38_merged_500kbp_genes_concatenated_OG.bed

## GTEX first
# recover the original region sizes to encounter less pain later on:
bedtools slop -i T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed_concatenated.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b -500000 > T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed_concatenated_OG.bed
# change into gene expression column into one that only has the expressed genes, not including the "not expressed"

awk 'BEGIN{OFS="\t"; FS=OFS} {gsub(/,not_expressed/, "", $4); print $1, $2, $3, $4}' T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed_concatenated_OG.bed > T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed_concatenated_OG_cleaned.bed # wc -l 441

# count the numbers of expressed gtex genes
awk -F'\t' '{print $4}' T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed_concatenated_cleaned.bed | tr ',' '\n' | sort | uniq | wc -l
awk -F'\t' '{print gsub(/,/, "", $4)}' T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed_concatenated_cleaned.bed > T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed_concatenated_gene_count.bed # wc -l 441

## DAAP6 in house next
# recover the original region sizes to encounter less pain later on:
bedtools slop -i T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b -500000 > T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated_OG.bed

# change into gene expression column into one that only has the expressed genes, not including the "not expressed"

awk 'BEGIN{OFS="\t"; FS=OFS} {gsub(/,not_expressed/, "", $4); print $1, $2, $3, $4}' T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated_OG.bed > T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated_OG_cleaned.bed # wc -l 441

## Biggest chain of intersects you've made in your life

bedtools intersect -a T2D_EUR_hg38_clustered_sorted_final.bed -b T2D_EUR_hg38_merged_500kbp_genes_concatenated_OG.bed -wb | bedtools intersect -a - -b  T2D_EUR_hg38_merged_500kbp_genes_sorted_expressed_concatenated_OG_cleaned.bed  -wb | bedtools intersect -a - -b  T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_expressed_concatenated_OG_cleaned.bed  -wb  > T2D_long.bed

chr10   102775388       102775389       rs2486757       -0.0225 T       1.76e-08        chr10   102772843       102788443       peak_5546       chr10   102775220       102775508       peak6448        4       56      chr10_102775388_102775389_rs2486757_cluster_56_T2D_EUR_H3K27ac_ATAC  2 chr10   102775388       102782989       FBXL15,CUEDC2,C10orf95,CYP17A1,RPEL1,NFKB2,PSD,MFSD13A,TRIM8,ARL3,SFXN2,BORCS7,GBF1,ACTR1A,SUFU,AS3MT,CNNM2,NT5C2,INA,WBP1L     chr10   102775388       102782989       ACTR1A,ARL3,BORCS7,CUEDC2,GBF1,NFKB2,NT5C2,SUFU,TRIM8,WBP1L     chr10   102775388     102782989       ACTR1A,ARL3,BORCS7,CNNM2,CUEDC2,CYP17A1,FBXL15,GBF1,INA,MFSD13A,NFKB2,NT5C2,PSD,SFXN2,SUFU,TRIM8,WBP1L



## Pick out the fields of interest + filter by unique cluster + print $0

# interesting fields:


awk 'BEGIN{OFS="\t"; FS=OFS} {print $20, $21, $22, $17, $19, $23, $27, $31 } ' T2D_long.bed > T2D_long_filtered.bed

head 

chr10   102775388       102782989       56      2       FBXL15,CUEDC2,C10orf95,CYP17A1,RPEL1,NFKB2,PSD,MFSD13A,TRIM8,ARL3,SFXN2,BORCS7,GBF1,ACTR1A,SUFU,AS3MT,CNNM2,NT5C2,INA,WBP1L     ACTR1A,ARL3,BORCS7,CUEDC2,GBF1,NFKB2,NT5C2,SUFU,TRIM8,WBP1L     ACTR1A,ARL3,BORCS7,CNNM2,CUEDC2,CYP17A1,FBXL15,GBF1,INA,MFSD13A,NFKB2,NT5C2,PSD,SFXN2,SUFU,TRIM8,WBP1L
chr10   102775388       102782989       56      1       FBXL15,CUEDC2,C10orf95,CYP17A1,RPEL1,NFKB2,PSD,MFSD13A,TRIM8,ARL3,SFXN2,BORCS7,GBF1,ACTR1A,SUFU,AS3MT,CNNM2,NT5C2,INA,WBP1L     ACTR1A,ARL3,BORCS7,CUEDC2,GBF1,NFKB2,NT5C2,SUFU,TRIM8,WBP1L     ACTR1A,ARL3,BORCS7,CNNM2,CUEDC2,CYP17A1,FBXL15,GBF1,INA,MFSD13A,NFKB2,NT5C2,PSD,SFXN2,SUFU,TRIM8,WBP1L
chr10   110864659       110864660       57      2       ADRA2A,DUSP5,SMC3,BBIP1,SHOC2,RBM20,PDCD4       not_expressed   ADRA2A,BBIP1,DUSP5,PDCD4,SHOC2,SMC3
chr10   110918899       110918900       58      2       ADRA2A,DUSP5,SMC3,BBIP1,SHOC2,RBM20,PDCD4       not_expressed   ADRA2A,BBIP1,DUSP5,PDCD4,SHOC2,SMC3
chr10   112973123       113074092       59      2       HABP2,VTI1A,TCF7L2      not_expressed,TCF7L2,VTI1A      HABP2,TCF7L2,VTI1A

## first reduce the line count by counting/collpasing by the number of probes designed...

awk 'BEGIN{OFS="\t"; FS=OFS} {a[$1 OFS $2 OFS $3]=a[$1 OFS $2 OFS $3]","$5} END {for (i in a) {split(i, arr, OFS); if (a[i] ~ /^,/) {a[i] = substr(a[i], 2)}; print arr[1], arr[2], arr[3], a[i]}}' T2D_long_filtered.bed > T2D_long_filtered_collapsed.bed #| wc -l # 437

awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, gsub(/2/, "", $4)}' T2D_long_filtered_collapsed.bed > T2D_long_probe_2_count.bed
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, gsub(/1/, "", $4)}' T2D_long_filtered_collapsed.bed > T2D_long_probe_1_count.bed

## collapse 6th (gtex) and 7th (daap6) fields correctly

awk 'BEGIN{OFS="\t"; FS=OFS} {!seen[$1 OFS $2 OFS $3 OFS $7]++} END {for (key in seen) {print key}}' T2D_long_filtered.bed > T2D_long_filtered_collapsed_gtex_gene_list.bed
awk 'BEGIN{FS=OFS="\t"} {gsub(/not_expressed,/, "", $4); gsub(/not_expressed/, "", $4); print}' T2D_long_filtered_collapsed_gtex_gene_list.bed > T2D_long_filtered_collapsed_gtex_gene_list_cleaned.bed
# need another one that has non empty 4th columns
awk 'BEGIN{FS=OFS="\t"} {gsub(/not_expressed,/, "NA,", $4); gsub(/not_expressed/, "NA", $4); print}' T2D_long_filtered_collapsed_gtex_gene_list.bed > T2D_long_filtered_collapsed_gtex_gene_list_cleaned_nonempty.bed 


awk 'BEGIN{OFS="\t"; FS=OFS} {!seen[$1 OFS $2 OFS $3 OFS $8]++} END {for (key in seen) {print key}}' T2D_long_filtered.bed > T2D_long_filtered_collapsed_DAAP6_gene_list.bed
awk 'BEGIN{FS=OFS="\t"} {gsub(/not_expressed,/, "", $4); gsub(/not_expressed/, "", $4); print}' T2D_long_filtered_collapsed_DAAP6_gene_list.bed > T2D_long_filtered_collapsed_DAAP6_gene_list_cleaned.bed
awk 'BEGIN{FS=OFS="\t"} {gsub(/not_expressed,/, "NA,", $4); gsub(/not_expressed/, "NA", $4); print}' T2D_long_filtered_collapsed_DAAP6_gene_list.bed > T2D_long_filtered_collapsed_DAAP6_gene_list_cleaned_nonempty.bed 


# then count the number of expressed genes according to gtex

awk -F'\t' 'BEGIN{OFS="\t"} {n=split($4, a, ","); print $1, $2, $3, n}' T2D_long_filtered_collapsed_gtex_gene_list_cleaned.bed > T2D_long_gtex_collapsed_count.bed

# then count the number of expressed genes according to in house daap6 data

awk -F'\t' 'BEGIN{OFS="\t"} {n=split($4, a, ","); print $1, $2, $3, n}' T2D_long_filtered_collapsed_DAAP6_gene_list_cleaned.bed > T2D_long_DAAP6_collapsed_count.bed

# then count the number of snps that map onto the region

awk 'BEGIN{OFS="\t"; FS=OFS} {a[$1 OFS $2 OFS $3]=a[$1 OFS $2 OFS $3]","$4} END {for (i in a) {split(i, arr, OFS); if (a[i] ~ /^,/) {a[i] = substr(a[i], 2)}; print arr[1], arr[2], arr[3], a[i]}}' T2D_long_filtered.bed > T2D_long_filtered_snps_collapsed.bed #| wc -l # 437

awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, gsub(/,/, "", $4)+1}' T2D_long_filtered_snps_collapsed.bed > T2D_long_filtered_collapsed_snp_count.bed


## create the main collapsed file that I'll paste everything to:

awk 'BEGIN{OFS="\t"; FS=OFS} {!seen[$1 OFS $2 OFS $3 OFS $6]++} END {for (key in seen) {print key}}' T2D_long_filtered.bed > T2D_long_filtered_main.bed
awk -F'\t' 'BEGIN{OFS="\t"} {n=split($4, a, ","); print $1, $2, $3, n}' T2D_long_filtered_main.bed > T2D_long_filtered_main_count.bed



# start pasting I guess

# paste -d '\t' glory_collapsed_main.bed glory_collapsed_snp_count.bed glory_sure_probe_1_count.bed glory_sure_probe_2_count.bed glory_collapsed_gtex_gene_list.bed glory_gtex_collapsed_count.bed glory_collapsed_DAAP6_gene_list.bed glory_DAAP6_collapsed_count.bed > T2D_EUR_hg38_clustered_sorted_final_tidy_OG_regions.bed

# modified files can be intersected instead

bedtools intersect -a T2D_long_filtered_main.bed -b T2D_long_filtered_main_count.bed -wb | bedtools intersect -a - -b T2D_long_filtered_collapsed_snp_count.bed -wb | bedtools intersect -a - -b T2D_long_probe_1_count.bed -wb | bedtools intersect -a - -b T2D_long_probe_2_count.bed -wb | bedtools intersect -a - -b T2D_long_filtered_collapsed_gtex_gene_list_cleaned_nonempty.bed -wb | bedtools intersect -a - -b T2D_long_gtex_collapsed_count.bed -wb | bedtools intersect -a - -b T2D_long_filtered_collapsed_DAAP6_gene_list_cleaned_nonempty.bed -wb | bedtools intersect -a - -b T2D_long_DAAP6_collapsed_count.bed -wb |  awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $8, $12, $16, $20, $24, $28, $32, $36 }' > T2D_EUR_hg38_clustered_sorted_final_tidy_OG_regions.bed 

# sort the file I guess that

bedtools sort -i T2D_EUR_hg38_clustered_sorted_final_tidy_OG_regions.bed > T2D_EUR_hg38_clustered_sorted_final_tidy_OG_regions_pretty.bed

cat T2D_EUR_hg38_clustered_sorted_final_tidy_OG_regions_pretty.bed > T2D_EUR_hg38_clustered_sorted_final_tidy_OG_regions_pretty.tsv
#awk 'BEGIN{OFS="\t"; FS=OFS} {print NF}' T2D_EUR_hg38_clustered_sorted_final_tidy_OG_regions.bed | sort | uniq -c

## next steps

# number of regions within 20kbp of each other

# merge the snps within 50+/- kbp of each other
bedtools merge -d 200000 -i T2D_EUR_hg38_clustered_sorted_final_tidy_OG_regions_pretty.bed > T2D_EUR_hg38_regions_close.bed
wc -l 325 
# add -c and -o options 


# number of genes close by

## filter t2d_oligo_info.txt density score < 40


################################
### eQTLS ###
################################

# files here /novo/projects/departments/nnrco/genetic_department/share/for_GATL/collabs/Carol
# move some over?

# cp /novo/projects/departments/nnrco/genetic_department/share/for_GATL/collabs/Carol/Adipose_Subcutaneous_snps.hg19.bed Adipose_Subcutaneous_snps.hg19.bed 
# cp /novo/projects/departments/nnrco/genetic_department/share/for_GATL/collabs/Carol/AdipoExpress_snps.hg19.bed AdipoExpress_snps.hg19.bed
# cp /novo/projects/departments/nnrco/genetic_department/share/for_GATL/collabs/Carol/FUSION_ge_adipose_naive_snps.hg19.bed FUSION_ge_adipose_naive_snps.hg19.bed
# cp /novo/projects/departments/nnrco/genetic_department/share/for_GATL/collabs/Carol/TwinsUK_ge_fat_snps.hg19.bed TwinsUK_ge_fat_snps.hg19.bed



# all in hg19 so need to convert to hg38

# /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/liftOver AdipoExpress_snps.hg19.bed /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/hg19ToHg38.over.chain.gz AdipoExpress_snps_lifted.bed  unmapped_AdipoExpress_snps.bed
# /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/liftOver Adipose_Subcutaneous_snps.hg19.bed /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/hg19ToHg38.over.chain.gz Adipose_Subcutaneous_snps_lifted.bed  unmapped_Adipose_Subcutaneous_snps.bed
# /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/liftOver FUSION_ge_adipose_naive_snps.hg19.bed /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/hg19ToHg38.over.chain.gz FUSION_ge_adipose_naive_snps_lifted.bed  unmapped_FUSION_ge_adipose_naive_snps.bed
# /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/liftOver TwinsUK_ge_fat_snps.hg19.bed /nfs_home/projects/departments/nnrco/genetic_department/projects/TargetLookup/hg19ToHg38.over.chain.gz TwinsUK_ge_fat_snps_lifted.bed  unmapped_TwinsUK_ge_fat_snps.bed


# 
bedtools intersect -a Adipose_Subcutaneous_snps_lifted.bed -b T2D_EUR_hg38_clustered.bed > Adipose_Subcutaneous_snps_lifted_intersected_T2D_EUR.bed # | head # wc -l #12 
bedtools intersect -a AdipoExpress_snps_lifted.bed -b T2D_EUR_hg38_clustered.bed > AdipoExpress_snps_lifted_intersected_T2D_EUR.bed  # | head # wc -l #21
bedtools intersect -a FUSION_ge_adipose_naive_snps_lifted.bed -b T2D_EUR_hg38_clustered.bed > FUSION_ge_adipose_naive_snps_lifted_intersected_T2D_EUR.bed  #| head # wc -l #7
bedtools intersect -a TwinsUK_ge_fat_snps_lifted.bed -b T2D_EUR_hg38_clustered.bed > TwinsUK_ge_fat_snps_lifted_intersected_T2D_EUR.bed # | head # wc -l #1


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
# sort -k10,10 T2D_EUR_hg38_merged_500kbp_genes.bed > T2D_EUR_hg38_merged_500kbp_genes_sorted.bed # already sorted


join -a1 -e 'not_in_geneset' -o '0,2.1' -1 10 -2 1 T2D_EUR_hg38_merged_500kbp_genes_sorted.bed gene_sets_sorted.txt | grep -v not_in_geneset | head


join -a1 -e 'not_in_geneset'  -t $'\t' -o '1.1,1.2,1.3,1.4,1.5,1.6,0,2.1' -1 10 -2 1 T2D_EUR_hg38_merged_500kbp_genes_sorted.bed gene_sets_sorted.txt > T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_genesets.bed



wc -l T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_genesets.bed # 6866
head T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_genesets.bed

awk 'BEGIN{OFS="\t"} {key = $1 OFS $2 OFS $3; if (key in a) a[key] = a[key] "," $8; else a[key] = $8} END {for (key in a) {split(key, arr, OFS); print arr[1], arr[2], arr[3], a[key]}}' T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_genesets.bed > T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_genesets_concatenated.bed

bedtools slop -i T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_genesets_concatenated.bed -g /fsx/home_dirs/czjs/doctoral_work/2023/guides/sizes.hg.38.fai -b -500000 > T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_genesets_concatenated_OG.bed

awk 'BEGIN{OFS="\t"; FS=OFS} {gsub(/,not_in_geneset/, "", $4); print $1, $2, $3, $4}' T2D_EUR_hg38_merged_500kbp_genes_sorted_DAAP6_genesets_concatenated_OG.bed > T2D_EUR_hg38_merged_500kbp_genesets_OG_cleaned.bed # wc -l 441


awk 'BEGIN{FS=OFS="\t"} {gsub(/not_in_geneset,/, "", $4); gsub(/not_in_geneset/, "", $4); print}' T2D_EUR_hg38_merged_500kbp_genesets_OG_cleaned.bed > T2D_EUR_hg38_pathways.bed
awk 'BEGIN{FS=OFS="\t"} {gsub(/not_in_geneset,/, "NA,", $4); gsub(/not_in_geneset/, "NA", $4); print}' T2D_EUR_hg38_merged_500kbp_genesets_OG_cleaned.bed > T2D_EUR_hg38_pathways_nonempty.bed 

#count

awk -F'\t' 'BEGIN{OFS="\t"} {n=split($4, a, ","); print $1, $2, $3, n}' T2D_EUR_hg38_pathways.bed > T2D_long_pathways_collapsed_count.bed


bedtools intersect -a T2D_EUR_hg38_clustered_sorted_final_tidy_OG_regions.bed -b T2D_long_pathways_collapsed_count.bed -wb |  awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $16}' > T2D_EUR_hg38_tidy_genesets.bed

cat T2D_EUR_hg38_tidy_genesets.bed > T2D_EUR_hg38_tidy_genesets.tsv


################################
### Finemapping overlap ###
################################

# file is here: /novo/projects/departments/nnrco/genetic_department/share/for_LWMR/fm_precomputation_v1/diagram.2024_Suzuki.T2D.EUR/cs_annotated.tsv

cp /novo/projects/departments/nnrco/genetic_department/share/for_LWMR/fm_precomputation_v1/diagram.2024_Suzuki.T2D.EUR/cs_annotated.tsv T2D_EUR_Suzuki_finemapping.bed
# sed -i '1d' T2D_EUR_Suzuki_finemapping.bed # remove 1st row, once. 

# match/ join to this file T2D_EUR_hg38_clustered_sorted_final.bed

# Sort files by the first column if not already sorted
sort -k1,1 T2D_EUR_Suzuki_finemapping.bed > T2D_EUR_Suzuki_finemapping_sorted.bed
sort -k1,1 T2D_EUR_hg38_clustered_sorted_final.bed > T2D_EUR_hg38_clustered_sorted_final_sorted.bed

# Join the files on the first column
join -1 1 -2 1 -t $'\t' T2D_EUR_Suzuki_finemapping_sorted.bed T2D_EUR_hg38_clustered_sorted_final_sorted.bed > joined_output.bed


# nngene_internal.ukbb_EUR_LGLM_2023.T2D_EUR_sexcombined
# TG2HDLcRatio_UKBB_EUR_Oliveri_2024

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

bedtools intersect -v -b gene_tss_coordinates_1kbp.bed -a T2D_EUR_hg38_clustered_sorted_final.bed >  T2D_EUR_hg38_1kbp_TSS.bed #wc -l #685
cat T2D_EUR_hg38_1kbp_TSS.bed > T2D_EUR_hg38_1kbp_TSS.tsv

bedtools intersect -v -b gene_tss_coordinates_2kbp.bed -a T2D_EUR_hg38_clustered_sorted_final.bed > T2D_EUR_hg38_2kbp_TSS.bed #wc -l #572
cat T2D_EUR_hg38_2kbp_TSS.bed > T2D_EUR_hg38_2kbp_TSS.tsv

### -------------------------------- ###
### dealing with problematic oligos
### -------------------------------- ###

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/main_attempt_120824/15419_1/oligo_info.txt t2d_oligos_round1.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/main_attempt_120824/15419_1/leftright_coords_round1.bed t2d_oligos_round1_redesign.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/main_attempt_120824/15419_2/oligo_info.txt t2d_oligos_round2.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/main_attempt_120824/15419_2/leftright_coords_round2.bed t2d_oligos_round2_redesign.bed

# whichever one is 8.30pm onwards\
\
### building off of t2d_hg38_selected_oligos.bed

mkdir manual_left_right_t2d

cd manual_left_right_t2d

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/T2D_EUR_hg38_selected_oligos.tsv t2d_round1_oligo_info.txt

autogen_probes_left_right_modded.R

# manually ran through the R script because it's a pain

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh t2d_round2 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_t2d/leftright_coords_round1.bed
squeue --user=czjs --iterate=30_seconds

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/t2d_round2/oligo_info.txt t2d_round2_oligo_info.txt

autogen_probes_left_right_modded.R

processing_snps_T2D.R

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh t2d_reselected_round2 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_t2d/t2d_hg38_selected_snps_round2_capready.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/t2d_reselected_round2/oligo_info.txt t2d_reselected_oligo_info_round2.txt

autogen_probes_left_right_modded.R

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh t2d_reselected_round3 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_t2d/leftright_coords_reselected_round1.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/t2d_reselected_round3/oligo_info.txt t2d_reselected_oligo_info_round3.txt

autogen_probes_left_right_modded.R

#### overlaps between whr and t2d #####

bedtools intersect -a final_clean_nondup_oligos_whr.bed -b final_clean_nondup_oligos_t2d.bed -wb > whr_t2d_overlap.bed

bedtools intersect -a whr_no_overlap_rows.bed -b t2d_no_overlap_rows.bed -wb


### -------------------------------- ###
### filtering out SNPs that are exonic
### -------------------------------- ###


bedtools intersect -v -a T2D_EUR_hg38_clustered_sorted_final.bed -b exon_coordinates.bed | wc -l #707
bedtools intersect -v -a T2D_EUR_hg38_2kbp_TSS.bed -b exon_coordinates.bed | wc -l #452 from 572


bedtools intersect -v -a T2D_EUR_hg38_2kbp_TSS.bed -b exon_coordinates.bed > T2D_EUR_hg38_clustered_sorted_final_no_exons_no_TSS.bed

wc -l T2D_EUR_hg38_clustered_sorted_final.bed #1008, 301 snps removed
wc -l T2D_EUR_hg38_clustered_sorted_final_no_exons_no_TSS.bed #452
# save the intersection and have a look on ucsc?

bedtools intersect -a T2D_EUR_hg38_clustered_sorted_final.bed -b exon_coordinates.bed > T2D_EUR_hg38_clustered_sorted_final_exonic_check.bed
# number of regions left: 236
# rerun capsequm with intronic no tss snps

awk 'BEGIN{OFS="\t"; FS=OFS} {print $1, $2, $3, $1"_"$2"_"$3"_"$4"_cluster_"$17"_T2D_EUR_H3K27ac_ATAC" }' T2D_EUR_hg38_clustered_sorted_final_no_exons_no_TSS.bed > T2D_EUR_hg38_clustered_sorted_final_no_exons_no_TSS_4_col.bed


sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh t2d_no_exons_no_TSS /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/T2D_EUR_hg38_clustered_sorted_final_no_exons_no_TSS_4_col.bed
# and then you need to filter to get the selected oligos...use processing_snps



### REPEAT/ REDO

### building off of T2D_EUR_hg38_selected_oligos.bed

mkdir manual_left_right_t2d

cd manual_left_right_t2d

autogen_probes_left_right_modded.R T2D_EUR_hg38_selected_oligos.bed

# manually ran through the R script because it's a pain

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh T2D_EUR_round2 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_t2d/leftright_coords_round1.bed
squeue --user=czjs --iterate=30_seconds

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/T2D_EUR_round2/oligo_info.txt t2d_round2_oligo_info.txt

autogen_probes_left_right_modded.R t2d_round2_oligo_info.txt

processing_snps_T2D_EUR.R

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh T2D_EUR_round3 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_t2d/t2d_hg38_selected_snps_round2_capready.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/T2D_EUR_round3/oligo_info.txt t2d_reselected_oligo_info.txt

autogen_probes_left_right_modded.R t2d_reselected_oligo_info.txt

sbatch /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh t2d_reselected_round2 /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/manual_left_right_t2d/leftright_coords_reselected_round1.bed

cp /novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/t2d_reselected_round2/oligo_info.txt t2d_reselected_oligo_info_round2.txt

autogen_probes_left_right_modded.R t2d_reselected_oligo_info_round2.txt

/novo/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/processing_oligos_t2d.R




#sed -i '1d' final_clean_nondup_oligos_whr.bed 

#sed -i '1d' final_clean_nondup_oligos_t2d.bed



bedtools intersect -a final_clean_nondup_oligos_whr.bed -b final_clean_nondup_oligos_t2d.bed -wb > whr_t2d_overlap.bed

bedtools intersect -a whr_no_overlap_rows.bed -b t2d_no_overlap_rows.bed -wb # good, nothing
