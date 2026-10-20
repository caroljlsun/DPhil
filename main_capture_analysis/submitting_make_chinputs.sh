#!/bin/bash
#SBATCH --job-name=chinputs
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --ntasks-per-node=1
#SBATCH --time 0-06:00
#SBATCH --mem=65000MB
#SBATCH --output=comp_bio/capture_c_adipocytes/main_capture_analysis/my_logs/-%A-%a.out
#SBATCH --error=comp_bio/capture_c_adipocytes/main_capture_analysis/my_logs/%A-%a.err
#SBATCH --mail-user=czjs@novonordisk.com
#SBATCH --mail-type=ALL

#---------------------------------------------
# Procedure: Script to make chinputs etc.
# Project: Capture C
# User: CZJS
# Date: Feb-2025
#---------------------------------------------
# follow the structure of DAAP6mADs_M1 etc.

PASSAGE=$1 

module load R/4.2.0
    
Rscript --vanilla comp_bio/capture_c_adipocytes/main_capture_analysis/making_chinputs_n_peakys.R $PASSAGE