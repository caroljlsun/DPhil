# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz

import pyranges as pr 
import pandas as pd

# df = pr.read_gtf("gencode.v39.annotation.nochr.gtf").df

df = pr.read_gtf("gencode.v39.annotation.gtf.gz").df

# Filter for exons
exon_coordinates = df.query('Feature=="exon"')[["Chromosome", "Start", "End", "Strand", "gene_name", "gene_type", "exon_number"]].drop_duplicates()


#print(exon_coordinates.head())

# Reorder according to minimum BED format: Chromosome, Start, Stop, Name
#exon_coordinates_bed = exon_coordinates[["Chromosome", "Start", "End", "gene_name", "exon_number"]]

# Save to BED file
exon_coordinates.to_csv("exon_coordinates.bed", sep="\t", header=False, index=False)