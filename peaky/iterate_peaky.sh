#!/bin/bash
#SBATCH --job-name=peaky
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --cpus-per-task=64
#SBATCH --ntasks-per-node=1
#SBATCH --time 0-12:00
#SBATCH --mem=65000MB
#SBATCH --output=comp_bio/capture_c_adipocytes/peaky/logs/-%A-%a.out
#SBATCH --error=comp_bio/capture_c_adipocytes/peaky/logs/%A-%a.err
#SBATCH --mail-user=czjs@novonordisk.com
#SBATCH --mail-type=ALL



#---------------------------------------------
# Procedure: Script to reiterate peaky
# Project: Peaky
# User: CZJS
# Date: Jan-2025
#---------------------------------------------

# run like this: sbatch iterate_peaky.sh testing 250000 10000 5

module load singularity
module load rocker/rstudio/4.4.0/


rstudio-shell Rscript $1 $2 $3 $4