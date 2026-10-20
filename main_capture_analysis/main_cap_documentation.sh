cd comp_bio/capture_c_adipocytes/main_capture_analysis


wget https://objectstorage.uk-london-1.oraclecloud.com/p/vnPapcM--KLQxT1GuZYE1ym9w0ivtKsWSSd3E2n-kd1juD_qdu6ZGXagYTslb0VS/n/cnyr09qj8zbo/b/england-data/o/out/CP2020081400001/X204SC25010669-Z01-F001/checkSize.xls
wget https://objectstorage.uk-london-1.oraclecloud.com/p/vnPapcM--KLQxT1GuZYE1ym9w0ivtKsWSSd3E2n-kd1juD_qdu6ZGXagYTslb0VS/n/cnyr09qj8zbo/b/england-data/o/out/CP2020081400001/X204SC25010669-Z01-F001/MD5.txt
wget https://objectstorage.uk-london-1.oraclecloud.com/p/vnPapcM--KLQxT1GuZYE1ym9w0ivtKsWSSd3E2n-kd1juD_qdu6ZGXagYTslb0VS/n/cnyr09qj8zbo/b/england-data/o/out/CP2020081400001/X204SC25010669-Z01-F001/Readme.html

wget https://objectstorage.uk-london-1.oraclecloud.com/p/vnPapcM--KLQxT1GuZYE1ym9w0ivtKsWSSd3E2n-kd1juD_qdu6ZGXagYTslb0VS/n/cnyr09qj8zbo/b/england-data/o/out/CP2020081400001/X204SC25010669-Z01-F001/X204SC25010669-Z01-F001_01.tar
wget https://objectstorage.uk-london-1.oraclecloud.com/p/vnPapcM--KLQxT1GuZYE1ym9w0ivtKsWSSd3E2n-kd1juD_qdu6ZGXagYTslb0VS/n/cnyr09qj8zbo/b/england-data/o/out/CP2020081400001/X204SC25010669-Z01-F001/X204SC25010669-Z01-F001_02.tar

md5sum X204SC25010669-Z01-F001_01.tar 
f2a18d1578d547ea5d5c44cf6c40dbaa  X204SC25010669-Z01-F001_01.tar # match!
md5sum X204SC25010669-Z01-F001_02.tar 
a9117951d902d6cd7dbcdb9f1be41659  X204SC25010669-Z01-F001_02.tar # match!
# md5sums from novogene
f2a18d1578d547ea5d5c44cf6c40dbaa  X204SC25010669-Z01-F001_01.tar
a9117951d902d6cd7dbcdb9f1be41659  X204SC25010669-Z01-F001_02.tar


# Extract the tar files
# tar -xvf X204SC25010669-Z01-F001_01.tar # gave me undetermined and P3
# tar -xvf X204SC25010669-Z01-F001_02.tar # gave me P1, M2, M1, M3, P2

mamba activate fullcap

# config files etc...

cp comp_bio/capture_c_adipocytes/pilot_capture_pt2_analysis/capcruncher_config.yml comp_bio/capture_c_adipocytes/main_capture_analysis/capcruncher_config.yml

# modify accordingly and rename files to

#DAAP6mADs_M1_R1 # etc

# mature adipocytes
mv M1_EKDL250000811-1A_22WLVCLT3_L3_1.fq.gz DAAP6mADs_M1_R1.fq.gz
mv M1_EKDL250000811-1A_22WLVCLT3_L3_2.fq.gz DAAP6mADs_M1_R2.fq.gz

mv M2_EKDL250000811-1A_22WLVCLT3_L3_1.fq.gz DAAP6mADs_M2_R1.fq.gz
mv M2_EKDL250000811-1A_22WLVCLT3_L3_2.fq.gz DAAP6mADs_M2_R2.fq.gz

mv M3_EKDL250000811-1A_22WLVCLT3_L3_1.fq.gz DAAP6mADs_M3_R1.fq.gz
mv M3_EKDL250000811-1A_22WLVCLT3_L3_2.fq.gz DAAP6mADs_M3_R2.fq.gz

# preadipocytes

mv P1_EKDL250000811-1A_22WLVCLT3_L3_1.fq.gz DAAP6pADs_P1_R1.fq.gz
mv P1_EKDL250000811-1A_22WLVCLT3_L3_2.fq.gz DAAP6pADs_P1_R2.fq.gz

mv P2_EKDL250000811-1A_22WLVCLT3_L3_1.fq.gz DAAP6pADs_P2_R1.fq.gz
mv P2_EKDL250000811-1A_22WLVCLT3_L3_2.fq.gz DAAP6pADs_P2_R2.fq.gz

mv P3_EKDL250000811-1A_22WLVCLT3_L3_1.fq.gz DAAP6pADs_P3_R1.fq.gz
mv P3_EKDL250000811-1A_22WLVCLT3_L3_2.fq.gz DAAP6pADs_P3_R2.fq.gz

# lol need to change the fq into fastq
mv DAAP6mADs_M1_R1.fq.gz DAAP6mADs_M1_R1.fastq.gz
mv DAAP6mADs_M1_R2.fq.gz DAAP6mADs_M1_R2.fastq.gz

mv DAAP6mADs_M2_R1.fq.gz DAAP6mADs_M2_R1.fastq.gz
mv DAAP6mADs_M2_R2.fq.gz DAAP6mADs_M2_R2.fastq.gz

mv DAAP6mADs_M3_R1.fq.gz DAAP6mADs_M3_R1.fastq.gz
mv DAAP6mADs_M3_R2.fq.gz DAAP6mADs_M3_R2.fastq.gz

mv DAAP6pADs_P1_R1.fq.gz DAAP6pADs_P1_R1.fastq.gz
mv DAAP6pADs_P1_R2.fq.gz DAAP6pADs_P1_R2.fastq.gz

mv DAAP6pADs_P2_R1.fq.gz DAAP6pADs_P2_R1.fastq.gz
mv DAAP6pADs_P2_R2.fq.gz DAAP6pADs_P2_R2.fastq.gz

mv DAAP6pADs_P3_R1.fq.gz DAAP6pADs_P3_R1.fastq.gz
mv DAAP6pADs_P3_R2.fq.gz DAAP6pADs_P3_R2.fastq.gz

# undetermined
mkdir -p undetermined
mv Undetermined_* undetermined/

# bring over and process the restriction fragment bed file


cp comp_bio/capture_c_adipocytes/capoverlaps/final_merged_oligos_non_overlapping.bed final_merged_oligos_non_overlapping.bed

# processing in R because it's easier...
# processing in processing_documentation.R to get unique_restriction_fragments.bed

capcruncher pipeline -c all --profile capc --cores 32  -j --rerun-incomplete --notemp --resources "mem_mb=118000" --keep-going
capcruncher pipeline --unlock

squeue --user=czjs --iterate=30_seconds

mamba install -c bioconda ucsc-bedtobigbed

# may need to edit checkpoint deduplication rule to have longer times to finish bc its massive

# also the pipeline stops each time I close my laptop rahhhhh

# tmux time again


tmux new -s capcruncher
mamba activate fullcap
capcruncher pipeline --unlock -c all --profile capc --cores 32  -j --rerun-incomplete --notemp --resources "mem_mb=118000" --keep-going

capcruncher pipeline -c all --profile capc --cores 32  -j --rerun-incomplete --notemp --resources "mem_mb=118000" --keep-going --keep-incomplete
tmux attach -t capcruncher

# unlock error
# need to go through each smk file and unlock them before
snakemake --snakefile /nfs_home/users/czjs/miniconda3/envs/fullcap/lib/python3.10/site-packages/capcruncher/pipeline/workflow/rules/align.smk --profile capc --configfile capcruncher_config.yml --unlock

Error executing rule count on cluster (jobid: 3165, external: 6122077, jobscript: comp_bio/capture_c_adipocytes/main_capture_analysis/.snakemake/tmp.wl8lpyjo/snakejob.count.3165.sh). For error details see the cluster log and the log files of the involved rule(s).
Removing output files of failed job count since they might be corrupted:
capcruncher_output/interim/pileups/counts_by_restriction_fragment/DAAP6mADs_M3.hdf5

# the log files for the rule count OOM'd out even with 16 gbs...
# increase memory and time

# need me a dpnii restriction map file
cp comp_bio/capture_c_adipocytes/main_capture_analysis/capcruncher_output/resources/restriction_fragments/genome.digest.bed genome.digest.bed 

# zip and send over

tar -czf mADs_M1_norm.tar.gz comp_bio/capture_c_adipocytes/main_capture_analysis/capcruncher_output/results/DAAP6mADs_M1/bigwigs/norm
tar -czf mADs_M2_norm.tar.gz comp_bio/capture_c_adipocytes/main_capture_analysis/capcruncher_output/results/DAAP6mADs_M2/bigwigs/norm
tar -czf mADs_M3_norm.tar.gz comp_bio/capture_c_adipocytes/main_capture_analysis/capcruncher_output/results/DAAP6mADs_M3/bigwigs/norm


tar -czf pADs_P1_norm.tar.gz comp_bio/capture_c_adipocytes/main_capture_analysis/capcruncher_output/results/DAAP6pADs_P1/bigwigs/norm
tar -czf pADs_P2_norm.tar.gz comp_bio/capture_c_adipocytes/main_capture_analysis/capcruncher_output/results/DAAP6pADs_P2/bigwigs/norm
tar -czf pADs_P3_norm.tar.gz comp_bio/capture_c_adipocytes/main_capture_analysis/capcruncher_output/results/DAAP6pADs_P3/bigwigs/norm


# Backup into NNEDL

cp comp_bio/capture_c_adipocytes/main_capture_analysis/X204SC25010669-Z01-F001_01.tar /nfs_home/projects/departments/compbio/sysbio/Projects/czjs_sequencing/X204SC25010669-Z01-F001_01.tar

cp comp_bio/capture_c_adipocytes/main_capture_analysis/X204SC25010669-Z01-F001_02.tar /nfs_home/projects/departments/compbio/sysbio/Projects/czjs_sequencing/X204SC25010669-Z01-F001_02.tar




### run peaky.R

cp comp_bio/capture_c_adipocytes/peaky/M1/plot_all_mpcc.bedGraph M1_plot_all_mpcc.bedGraph 
cp comp_bio/capture_c_adipocytes/peaky/M1/plot_sig_mpcc.bedGraph M1_plot_sig_mpcc.bedGraph


### -------------------------------- ###
### TSS and gene body intersections with peaky mpcc
### -------------------------------- ###

cat M1_plot_sig_mpcc.bedGraph | awk '{print $1 "\t" $2 "\t" $3}' > M1_sig_mpcc.bed

bedtools intersect -a M1_sig_mpcc.bed -b ../capoverlaps/gene_tss_coordinates_5kbp.bed > M1_sig_mpcc_gene_tss_5kbp.bed

bedtools intersect -a M1_sig_mpcc.bed -b ../capoverlaps/hg38_all_genes.bed -wo > M1_sig_mpcc_gene_bodies.bed

bedtools intersect -a M1_sig_mpcc.bed -b ../capoverlaps/hg38_all_genes_extended.bed > M1_sig_mpcc_gene_bodies_5kbp.bed


# quarto isn't working for some reason
# load a different environment and try again?

#conda create --name cloned_fullcap --clone fullcap
mamba activate cloned_fullcap
mamba uninstall bokeh --force 
pip install bokeh 
pip install --force-reinstall -v "bokeh==3.2.1"
mamba install panel 


#quarto         render         capcruncher_output/results/capcruncher_report.qmd         --to html         --execute         -P fastq_deduplication_path:$(realpath capcruncher_output/interim/statistics/deduplication/data/)         -P fastq_trimming_path:$(realpath capcruncher_output/interim/statistics/trimming/trimming.json)         -P fastq_flash_path:$(realpath capcruncher_output/interim/statistics/flash/flash.json)         -P fastq_digestion_path:$(realpath capcruncher_output/interim/statistics/digestion/data/)         -P reporter_filtering_path:$(realpath capcruncher_output/interim/statistics/filtering/data/)         -P reporter_deduplication_path:$(realpath capcruncher_output/interim/statistics/deduplication_final/data/)         -P reporter_cis_trans_path:$(realpath capcruncher_output/interim/statistics/cis_and_trans_reporters/data/)         --log capcruncher_output/logs/make_report.log         2> capcruncher_output/logs/make_report.log.err;
capcruncher pipeline -c all --profile capc --cores 32 --forcerun make_report --rerun-incomplete --notemp --resources "mem_mb=118000" --keep-going

IncompleteFilesException:
The files below seem to be incomplete. If you are sure that certain files are not incomplete, mark them as complete with

    snakemake --cleanup-metadata <filenames>

To re-generate the files rerun your command with the --rerun-incomplete flag.
Incomplete files:
capcruncher_output/results/DAAP6mADs_M1/bigwigs/raw/DAAP6mADs_M1_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6mADs_M1/bigwigs/raw/DAAP6mADs_M1_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6mADs_M1/bigwigs/norm/DAAP6mADs_M1_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6mADs_M1/bigwigs/norm/DAAP6mADs_M1_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6mADs_M2/bigwigs/raw/DAAP6mADs_M2_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6mADs_M2/bigwigs/raw/DAAP6mADs_M2_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6mADs_M2/bigwigs/norm/DAAP6mADs_M2_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6mADs_M2/bigwigs/norm/DAAP6mADs_M2_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6mADs_M3/bigwigs/raw/DAAP6mADs_M3_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6mADs_M3/bigwigs/raw/DAAP6mADs_M3_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6mADs_M3/bigwigs/norm/DAAP6mADs_M3_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6mADs_M3/bigwigs/norm/DAAP6mADs_M3_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P1/bigwigs/raw/DAAP6pADs_P1_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P1/bigwigs/raw/DAAP6pADs_P1_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P1/bigwigs/norm/DAAP6pADs_P1_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P1/bigwigs/norm/DAAP6pADs_P1_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P2/bigwigs/raw/DAAP6pADs_P2_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P2/bigwigs/raw/DAAP6pADs_P2_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P2/bigwigs/norm/DAAP6pADs_P2_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P2/bigwigs/norm/DAAP6pADs_P2_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P3/bigwigs/raw/DAAP6pADs_P3_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P3/bigwigs/raw/DAAP6pADs_P3_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P3/bigwigs/norm/DAAP6pADs_P3_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig
capcruncher_output/results/DAAP6pADs_P3/bigwigs/norm/DAAP6pADs_P3_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig
capcruncher_output/results/comparisons/counts_per_viewpoint/norm/chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.tsv
capcruncher_output/results/comparisons/counts_per_viewpoint/norm/chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.tsv

# seems like cluster 137 and 291 did not work very well 
#capcruncher pipeline -c all --profile capc --cores 32 --cleanup-metadata capcruncher_output/results/DAAP6mADs_M1/bigwigs/raw/DAAP6mADs_M1_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6mADs_M1/bigwigs/raw/DAAP6mADs_M1_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6mADs_M1/bigwigs/norm/DAAP6mADs_M1_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6mADs_M1/bigwigs/norm/DAAP6mADs_M1_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6mADs_M2/bigwigs/raw/DAAP6mADs_M2_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6mADs_M2/bigwigs/raw/DAAP6mADs_M2_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6mADs_M2/bigwigs/norm/DAAP6mADs_M2_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6mADs_M2/bigwigs/norm/DAAP6mADs_M2_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6mADs_M3/bigwigs/raw/DAAP6mADs_M3_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6mADs_M3/bigwigs/raw/DAAP6mADs_M3_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6mADs_M3/bigwigs/norm/DAAP6mADs_M3_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6mADs_M3/bigwigs/norm/DAAP6mADs_M3_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P1/bigwigs/raw/DAAP6pADs_P1_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P1/bigwigs/raw/DAAP6pADs_P1_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P1/bigwigs/norm/DAAP6pADs_P1_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P1/bigwigs/norm/DAAP6pADs_P1_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P2/bigwigs/raw/DAAP6pADs_P2_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P2/bigwigs/raw/DAAP6pADs_P2_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P2/bigwigs/norm/DAAP6pADs_P2_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P2/bigwigs/norm/DAAP6pADs_P2_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P3/bigwigs/raw/DAAP6pADs_P3_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P3/bigwigs/raw/DAAP6pADs_P3_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P3/bigwigs/norm/DAAP6pADs_P3_chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.bigWig capcruncher_output/results/DAAP6pADs_P3/bigwigs/norm/DAAP6pADs_P3_chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.bigWig capcruncher_output/results/comparisons/counts_per_viewpoint/norm/chr15_41791645_41791646_rs1197705:G:A_cluster_291_whradjbmi_combined_H3K27ac_ATAC.tsv capcruncher_output/results/comparisons/counts_per_viewpoint/norm/chr15_41791645_41791646_rs1197705_cluster_137_T2D_EUR_H3K27ac_ATAC.tsv 
# doesnt worrkk

# mamba install panel
#cd comp_bio/capture_c_adipocytes/main_capture_analysis && /nfs_home/users/czjs/miniconda3/envs/cloned_fullcap/bin/python3.10 -m snakemake --snakefile '/novo/users/czjs/miniconda3/envs/cloned_fullcap/lib/python3.10/site-packages/capcruncher/pipeline/workflow/Snakefile' --target-jobs 'make_report:' --cores 'all' --attempt 3 --force-use-threads  --resources 'mem_mb=1000' 'mem_mib=954' 'disk_mb=1000' 'disk_mib=954' --wait-for-files-file 'comp_bio/capture_c_adipocytes/main_capture_analysis/.snakemake/tmp.7m6mgkj2/snakejob.make_report.50.sh.waitforfilesfile.txt' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'params' 'mtime' 'code' 'software-env' 'input' --skip-script-cleanup  --use-conda  --conda-frontend 'mamba' --conda-base-path '/nfs_home/users/czjs/miniconda3' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --printshellcmds  --latency-wait 240 --scheduler 'greedy' --scheduler-solver-path '/nfs_home/users/czjs/miniconda3/envs/cloned_fullcap/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && exit 0 || exit 1

#worked!!! sort of...

salloc --nodes=1 --time=480 --partition=highmem --mem=72G
srun --jobid=6200800 --pty /bin/bash

mamba install plotly
mamba install panel
mamba install seaborn pathlib
mamba install json

cd comp_bio/capture_c_adipocytes/main_capture_analysis && /nfs_home/users/czjs/miniconda3/envs/cloned_fullcap/bin/python3.10 -m snakemake --snakefile '/novo/users/czjs/miniconda3/envs/cloned_fullcap/lib/python3.10/site-packages/capcruncher/pipeline/workflow/Snakefile' --target-jobs 'make_report:' --attempt 3 --resources 'mem_mb=16000' 'mem_mib=15258' 'disk_mb=16000' 'disk_mib=15258' --wait-for-files-file 'comp_bio/capture_c_adipocytes/main_capture_analysis/.snakemake/tmp.7m6mgkj2/snakejob.make_report.50.sh.waitforfilesfile.txt' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'params' 'mtime' 'code' 'software-env' 'input' --skip-script-cleanup  --use-conda  --conda-frontend 'mamba' --conda-base-path '/nfs_home/users/czjs/miniconda3' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --printshellcmds  --latency-wait 240 --scheduler 'greedy' --scheduler-solver-path '/nfs_home/users/czjs/miniconda3/envs/cloned_fullcap/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && exit 0 || exit 1


ERROR: pip's dependency resolver does not currently take into account all the packages that are installed. This behaviour is the source of the following dependency conflicts.
capcruncher 0.3.12 requires numpy<=1.26.4, but you have numpy 2.2.3 which is incompatible.
capcruncher 0.3.12 requires pandas<=2.1.2, but you have pandas 2.2.3 which is incompatible.
ibis-framework 8.0.0 requires numpy<2,>=1, but you have numpy 2.2.3 which is incompatible.
matplotlib 3.8.3 requires numpy<2,>=1.21, but you have numpy 2.2.3 which is incompatible.
scikit-learn 1.4.1.post1 requires numpy<2.0,>=1.19.5, but you have numpy 2.2.3 which is incompatible.
scipy 1.12.0 requires numpy<1.29.0,>=1.22.4, but you have numpy 2.2.3 which is incompatible.
statsmodels 0.14.1 requires numpy<2,>=1.18, but you have numpy 2.2.3 which is incompatible.

mamba install numpy=1.26.4 --no-deps
mamba install pandas=2.1.2 --no-deps
mamba install bokeh=3.2.1 --no-deps