#!/bin/bash
#SBATCH --job-name=snakemake_job
#SBATCH --output=snakemake_job.out
#SBATCH --error=snakemake_job.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=118G


cd comp_bio/capture_c_adipocytes/main_capture_analysis

mamba activate fullcap

capcruncher pipeline -c all --profile capc --cores 32  -j --rerun-incomplete --notemp --resources "mem_mb=118000" --keep-going
