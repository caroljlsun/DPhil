#--------------------------------------------- #
# Procedure: Gene expression GTEX data
# Project: Capture-C
# User: CZJS
# Date: MAY-2024
#--------------------------------------------- #


### ----------------------------------------- ###
### Libraries ----
### ----------------------------------------- ###

easypackages::libraries( "generics",
                         "ggplot2",
                         "plyr",
                         "lubridate",
                         "Rmisc",
                         "tidyr",
                         "readxl",
                         "dplyr",
                         "data.table",
                         "gridExtra",
                         "tidyverse")


### ----------------------------------------- ###
### Set up ----
### ----------------------------------------- ###

setwd("comp_bio/capture_c_adipocytes/capoverlaps/")


### ----------------------------------------- ###
### Import and wrangle ----
### ----------------------------------------- ###

# Read the .bed file
bed_data <- read.table("whradjbmi_combined_hg38_merged_500kbp_genes.bed", header = FALSE, stringsAsFactors = FALSE)

# Get gene names and ens names

bed_genes <- bed_data$V10
bed_ens <- bed_data$V8

# Read the tsv file
gtex_data <- read.table("rna_tissue_gtex.tsv", skip = 1, sep = "\t", stringsAsFactors = FALSE, fill = T, header = F)


# Subset by adipose tissue

gtex_adipose <- gtex_data[gtex_data$V3 == "adipose tissue",]
gtex_adipose_genes <- gtex_adipose$V2

# Exploring data

head(gtex_adipose)

hist(gtex_adipose$V4) # what on earth is with the axis 
max(gtex_adipose$V4) # 41K??? 
which(gtex_adipose$V4 == 41921.1)
# 661399 ENSG00000228253 MT-ATP8 adipose tissue 41921.1 51142.1 48792.8

# let's limit the hist a bit 

ggplot(data = gtex_adipose, aes(x = V4)) + 
  geom_histogram(bins = 100) + 
  theme_minimal()+
  xlim(0, 20)+
  ylim(0, 3000)

# Filter by expression

gtex_adipose_filtered <- gtex_adipose[gtex_adipose$V4 > 5,] # decidedly arbitrary


# Check which genes from the .bed file are present in the .tsv file


expressed_genes <- bed_genes[bed_genes %in% gtex_adipose_filtered$V2]

# Print the expressed genes
print(expressed_genes)

# check that the standard ones exist:

match("PPARG", expressed_genes)
#--------------------------------------------- #
# Procedure: Summarise grs
# Project: Obesity ML
# User: CZJS
# Date: JUL-2022
#--------------------------------------------- #


### ----------------------------------------- ###
### Libraries ----
### ----------------------------------------- ###

easypackages::libraries( "generics",
                         "ggplot2",
                         "plyr",
                         "lubridate",
                         "Rmisc",
                         "tidyr",
                         "readxl",
                         "dplyr",
                         "data.table",
                         "gridExtra",
                         "tidyverse")


### ----------------------------------------- ###
### Set up ----
### ----------------------------------------- ###

setwd("comp_bio/capture_c_adipocytes/capoverlaps/")


### ----------------------------------------- ###
### Import and wrangle ----
### ----------------------------------------- ###

# Read the .bed file
bed_data <- read.table("whradjbmi_combined_hg38_merged_500kbp_genes.bed", header = FALSE, stringsAsFactors = FALSE)

# Get gene names and ens names

bed_genes <- bed_data$V10
bed_ens <- bed_data$V8

# Read the tsv file
gtex_data <- read.table("rna_tissue_gtex.tsv", skip = 1, sep = "\t", stringsAsFactors = FALSE, fill = T, header = F)


# Subset by adipose tissue

gtex_adipose <- gtex_data[gtex_data$V3 == "adipose tissue",]
gtex_adipose_genes <- gtex_adipose$V2

# Exploring data

head(gtex_adipose)

hist(gtex_adipose$V4) # what on earth is with the axis 
max(gtex_adipose$V4) # 41K??? 
which(gtex_adipose$V4 == 41921.1)
# 661399 ENSG00000228253 MT-ATP8 adipose tissue 41921.1 51142.1 48792.8

# let's limit the hist a bit 

ggplot(data = gtex_adipose, aes(x = V4)) + 
  geom_histogram(bins = 100) + 
  theme_minimal()+
  xlim(0, 20)+
  ylim(0, 3000)

# Filter by expression

gtex_adipose_filtered <- gtex_adipose[gtex_adipose$V4 > 5,] # decidedly arbitrary

# check that the standard ones exist:

match("LEP", gtex_adipose_filtered)

which(gtex_adipose$V2 == "LEP")
gtex_adipose[13406,]

### debug later ###
# Check which genes from the .bed file are present in the .tsv file

expressed_genes <- bed_genes[bed_genes %in% gtex_adipose_filtered$V2]

# Check the expressed genes

print(unique(expressed_genes))


write.table(unique(expressed_genes), "expressed_genes.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = F)
