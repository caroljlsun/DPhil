# copy of commands I used

#awk '{print $2+$5-500}' kristian_control_guide_regions_and_guides.txt # = the motif start, 500 bp because that's what kristian padded with
#awk '{print $2+$6-500 }' kristian_control_guide_regions_and_guides.txt # the mostif stop

awk '{print $1 "\t" $2+$5-500 "\t" $2+$6-500}' kristian_control_guide_regions_and_guides.txt > guides.bed

sed '1d' guides.bed > guides_no_head.bed

sort -k1,1V -k2,2n -k3,3n guides_no_head.bed > guides_no_head_sorted.bed

bedtools annotate -both -i contigs_no_head_sorted.bed -files guides_no_head_sorted.bed 

# check the differences of motifs start and stop


#awk '{print $6-$5}' kristian_control_guide_regions_and_guides.txt

conda deactivate  
conda remove --name capsequm-env --all
# conda create --name capsequm-env python>=3.4 pysam numpy pandas biopython repeatmasker blat star # nope doesn't work
cd ../doctoral_work/conda_stuff

conda env create -f capsequm.yml
# blah
cd ~/doctoral_work/2023/oligo
#python design.py Capture -f ~/hg19/Sequence/genome.fa -g hg19 -b viewpoints.bed -o 50 -e NlaIII -s ~/hg19/STAR/

python design.py Capture -f /nfs_home/users/czjs/doctoral_work/2023/guides/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz -g hg38 -b /nfs_home/users/czjs/doctoral_work/2023/guides/single_bp_regions.bed 



python design.py Capture -f /nfs_home/users/czjs/doctoral_work/2023/guides/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -g hg38 -b /nfs_home/users/czjs/doctoral_work/2023/guides/single_bp_regions.bed -s .conda/pkgs/star-2.7.10b-h9ee0642_0/bin

# error???

sort -k1,1V -k2,2n -k3,3n corrected_replaced_single_bp_regions.bed > sorted_corrected_replaced_single_bp_regions.bed
