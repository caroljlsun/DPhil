#!/bin/bash
#SBATCH --job-name=capsequm_seek
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --ntasks-per-node=1
#SBATCH --time 0-06:00
#SBATCH --mem=65000MB
#SBATCH --output=/nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/outputs/-%A-%a.out
#SBATCH --error=/nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/outputs/%A-%a.err
#SBATCH --mail-user=czjs@novonordisk.com
#SBATCH --mail-type=ALL

#---------------------------------------------
# Procedure: Script to run capsequm, an R script to process, capsequm again etc
# Project: Capture C
# User: CZJS
# Date: Sep-2023
#---------------------------------------------

### Set up


module load anaconda3/2021.05

#needed to make interactive conda env
eval "$(conda shell.bash hook)"

conda activate /nfs_home/users/czjs/miniconda3/envs/capsequm2

# make directory to store all the output

DIR_PATH="/nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/""main_attempt_""$(date "+%d%m%y")""/""$RANDOM"
mkdir $DIR_PATH -p

cd $DIR_PATH

### Counters
ROUND=1
GUIDES_LEFT=1
PATH_TO_BED_FILE="/nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/T2D_EUR_hg38_clustered_4col.bed"
ROUND_PATH="$DIR_PATH""_""$ROUND"

mkdir $ROUND_PATH -p
# python /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/design.py Capture -f /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -g hg38 -b /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capoverlaps/whradjbmi_combined_hg38_extra_H3K27ac_ATAC_4col.bed -s /data/nnrco/rnaseq_genome/Homo_sapiens.GRCh38/index_STAR 


### Loop

while [ $GUIDES_LEFT -gt 0 ]; #while the number of guides to design > 0
do
   
    echo "Round is "$ROUND "Guides left are" $GUIDES_LEFT
    
    cd $ROUND_PATH
    
    # python /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/design.py Capture -f /nfs_home/users/czjs/doctoral_work/2023/guides/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -g hg38 -b /nfs_home/users/czjs/doctoral_work/2023/guides/left_right_single_coords.bed -s /data/nnrco/rnaseq_genome/Homo_sapiens.GRCh38/index_STAR 
    
    # change so that it's iterative
    

    
    python /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/design.py Capture -f /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -g hg38 -b $PATH_TO_BED_FILE -s /data/nnrco/rnaseq_genome/Homo_sapiens.GRCh38/index_STAR 


     # Call plink
    module load R/4.2.0
    
    Rscript --vanilla /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/guides/autogen_probes_left_right_modded.R $ROUND_PATH $ROUND
    
    
     # read number of guide regions to design
    GUIDE_FILE_PATH="$ROUND_PATH""/""leftright_coords_""$ROUND"".txt"
    
    GUIDES_LEFT=$(cat $GUIDE_FILE_PATH)
    
    echo "Guides left" $GUIDES_LEFT
    
    echo "End of round" $ROUND
    
    PATH_TO_BED_FILE="$ROUND_PATH""/leftright_coords_round""$ROUND"".bed" #changes for the next round, offset by 1!!!

    echo $PATH_TO_BED_FILE
    
    let ROUND++ #add round counter
    
    ROUND_PATH="$DIR_PATH""_""$ROUND"
    
    mkdir $ROUND_PATH -p
    
    
    if [ $ROUND -gt 2 ]
    then
        break
    fi
done