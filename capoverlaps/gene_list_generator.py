# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz

import pyranges as pr 
import pandas as pd

# df = pr.read_gtf("gencode.v39.annotation.nochr.gtf").df

# Read the GTF file
df = pr.read_gtf("gencode.v39.annotation.gtf.gz").df

# Filter for genes and select relevant columns
gene_coordinates = df.query('Feature=="gene" and (gene_type=="protein_coding" or gene_type=="processed_pseudogene" or gene_type=="unprocessed_pseudogene")')[["Chromosome", "Start", "End", "Strand", "gene_name", "gene_type"]].drop_duplicates()

# Rename columns to match BED format
gene_coordinates = gene_coordinates.rename(columns={"Start": "start", "End": "stop"})

# Reorder columns according to BED format: Chromosome, Start, Stop, Name
gene_coordinates_bed = gene_coordinates[["Chromosome", "start", "stop", "gene_name"]]

# Save to BED file
gene_coordinates_bed.to_csv("gene_body_coordinates.bed", sep="\t", header=False, index=False)