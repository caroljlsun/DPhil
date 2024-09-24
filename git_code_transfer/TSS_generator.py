# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz

import pyranges as pr 
import pandas as pd

# df = pr.read_gtf("gencode.v39.annotation.nochr.gtf").df

df = pr.read_gtf("gencode.v39.annotation.gtf.gz").df

gene_coordinates = df.query('Feature=="gene" and (gene_type=="protein_coding" or gene_type=="processed_pseudogene" or gene_type=="unprocessed_pseudogene")')[["Chromosome","Start","End", "Strand", "gene_name", "gene_type"]].drop_duplicates()

gene_coordinates["TSS"] = [x if s=="+" else y for x,y,s in zip(gene_coordinates.Start, gene_coordinates.End, gene_coordinates.Strand)]

gene_coordinates = gene_coordinates[["Chromosome", "TSS", "gene_name","gene_type"]].drop_duplicates()

# Assuming TSS as both start and end for simplicity; adjust if needed
gene_coordinates = gene_coordinates[["Chromosome", "TSS", "gene_name", "gene_type"]]
gene_coordinates["stop"] = gene_coordinates["TSS"]+1 
gene_coordinates = gene_coordinates.rename(columns={"TSS": "start"})


# Reorder according to minimum BED format: Chromosome, Start, Stop, Name
gene_coordinates_bed = gene_coordinates[["Chromosome", "start", "stop", "gene_name"]]

# Save to BED filegene_coordinates = gene_coordinates.rename(columns={"TSS": "Start"})
gene_coordinates_bed.to_csv("gene_tss_coordinates.bed", sep="\t", header=False, index=False)