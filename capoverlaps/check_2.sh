
### Setup
cd comp_bio/capture_c_adipocytes/capoverlaps

# reusing another env with the newish bedtools

mamba activate fullcap

### make H3K27ac only and ATAC only interesected files, just in case:

bedtools intersect -a whradjbmi_combined_hg38_extra_filtered_lifted.bed -b H3K27ac_peaks.bed -wa| bedtools intersect -b H3K27ac_peaks.bed -a - -wb > whradjbmi_combined_hg38_extra_H3K27ac_only.bed
wc -l 
4037

bedtools intersect -a whradjbmi_combined_hg38_extra_filtered_lifted.bed -b atac.bed -wa| bedtools intersect -b atac.bed -a - -wb > whradjbmi_combined_hg38_extra_ATAC_only.bed
wc -l 
908

# set aside for future use, if not enough regions for capture c



### remove if within 5 kbp of a gene 

#full list of genes in hg38 is here genetic_department/projects/TargetLookup/gencode.v26.GRCh38.genes_for_gene.protein_coding.txt
# cp genetic_department/projects/TargetLookup/gencode.v26.GRCh38.genes_for_gene.protein_coding.txt hg38_all_genes.bed
# bedtools slop -i hg38_all_genes.bed -g doctoral_work/2023/guides/sizes.hg.38.fai -b 5000 > hg38_all_genes_extended.bed
# bedtools intersect -v -a whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed -b hg38_all_genes_extended.bed | wc -l
# bedtools intersect -v -a whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed -b hg38_all_genes_extended.bed > whradjbmi_combined_hg38_intergenic.bed
# ### merging snps within 50+/- kbp of each other
# bedtools merge -d 50000 -i whradjbmi_combined_hg38_intergenic.bed > whradjbmi_combined_hg38_intergenic_merged.bed
# # or?
# bedtools cluster -d 50000 -i whradjbmi_combined_hg38_intergenic.bed > whradjbmi_combined_hg38_intergenic_clustered.bed
# wc -l 168

### never mind!! don't remove genic snps, rerun merge and clump

# move a list of hg38 genes over 
cp genetic_department/projects/TargetLookup/gencode.v26.GRCh38.genes_for_gene.protein_coding.txt hg38_all_genes.bed

# merge the snps within 50+/- kbp of each other
bedtools merge -d 50000 -i whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed > whradjbmi_combined_hg38_merged.bed
wc -l 268

# cluster the snps within 50+/- kbp of each other
# gives me the identity of the merged snps in each cluster
bedtools cluster -d 50000 -i whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed > whradjbmi_combined_hg38_clustered.bed
wc -l 706


### assess complexity by looking at how many genes are close by 
# uni supervisor suggested I find the genes close by the clusters, not the snps

bedtools slop -i whradjbmi_combined_hg38_merged.bed -g doctoral_work/2023/guides/sizes.hg.38.fai -b 500000 > whradjbmi_combined_hg38_merged_500kbp.bed

# +/- 500kbp, so a million after expr 

bedtools intersect -wa -wb -a whradjbmi_combined_hg38_merged_500kbp.bed -b hg38_all_genes.bed > whradjbmi_combined_hg38_merged_500kbp_genes.bed
wc -l 4572


# check the length of the unique snps
cat whradjbmi_combined_hg38_merged_500kbp_genes.bed | awk '{print $1,$2,$3}' | uniq | wc -l
267

# I think one of them simply didn't intersect at all, 267 = 268 - 1
# need to find which one that is...
bedtools intersect -wa -wb -v -a whradjbmi_combined_hg38_merged_500kbp.bed -b hg38_all_genes.bed > whradjbmi_combined_hg38_merged_500kbp_nogenes.bed 


# yep one of the snps didn't have intersections +/- 500kbp genes

# unique, process columns such that you get one column with all the intersected genes in one place

awk 'BEGIN{OFS="\t"; FS=OFS} {a[$1 OFS $2 OFS $3]=a[$1 OFS $2 OFS $3]","$10} END {for (i in a) {split(i, arr, OFS); if (a[i] ~ /^,/) {a[i] = substr(a[i], 2)}; print arr[1], arr[2], arr[3], a[i]}}' whradjbmi_combined_hg38_merged_500kbp_genes.bed > whradjbmi_combined_hg38_merged_500kbp_genes_concatenated.bed
wc -l 267

# append the missing snp/cluster

cat whradjbmi_combined_hg38_merged_500kbp_nogenes.bed >> whradjbmi_combined_hg38_merged_500kbp_genes_concatenated.bed

### Run capsequm on snps

# need to make a bed file with only 4 columns: chr start stop unique_name

awk 'BEGIN{OFS="\t"; FS=OFS} {print $1, $2, $3, $1"_"$2"_"$3"_"$4"_whradjbmi_combined_H3K27ac_ATAC" }' whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed | head
awk 'BEGIN{OFS="\t"; FS=OFS} {print $1, $2, $3, $1"_"$2"_"$3"_"$4"_whradjbmi_combined_H3K27ac_ATAC" }' whradjbmi_combined_hg38_extra_H3K27ac_ATAC.bed > whradjbmi_combined_hg38_extra_H3K27ac_ATAC_4col.bed

cd comp_bio/capture_c_adipocytes/capsequm2/oligo
sbatch comp_bio/capture_c_adipocytes/capsequm2/guides/run_capsequm.sh

