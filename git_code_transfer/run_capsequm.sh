#!/bin/bash
#SBATCH --job-name=capsequm2
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
# Procedure: Script to run one specific line from capsequm
# Project: Capture C
# User: CZJS
# Date: Sep-2023
#---------------------------------------------

# Add missing import statement
import os

module load anaconda3/2021.05

#needed to make interactive conda env
eval "$(conda shell.bash hook)"

conda activate /nfs_home/users/czjs/miniconda3/envs/capsequm2

# arguments

output_folder=$1
file_to_process=$2
 
# make directory to store all the output

DIR_PATH="/nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/""$output_folder"
echo "output path is" $DIR_PATH

mkdir $DIR_PATH -p


cd $DIR_PATH

python /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/capsequm2/oligo/design.py Capture -f /nfs_home/projects/departments/nnrco/comp_bio/capture_c_adipocytes/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -g hg38 -b $file_to_process -s /data/nnrco/rnaseq_genome/Homo_sapiens.GRCh38/index_STAR 
