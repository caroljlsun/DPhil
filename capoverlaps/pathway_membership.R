library(tidyverse)
library(data.table)

setwd("comp_bio/capture_c_adipocytes/capoverlaps/")

# List of files to read
files <- list.files(pattern = "HALLMARK_ADIPOGENESIS.v2023.2.Hs.grp|KEGG.*|REACTOME.*|WP_.*")

whr_list <- fread("whradjbmi_combined_hg38_tidy_genesets.tsv")
t2d_list <- fread("T2D_EUR_hg38_tidy_genesets.tsv")

# Read and concatenate files
data_list <- lapply(files, function(file) fread(file, skip = 1))

combined_data <- bind_rows(data_list)

long_data <- combined_data %>%
    pivot_longer(cols = starts_with("#"), names_to = "Pathway", values_to = "Gene") %>% 
    drop_na()

clean_long_data <- long_data %>%
    mutate(Pathway = str_extract(Pathway, "[^/]+$"))

# Count the occurrences of each pathway
pathway_counts <- clean_long_data %>%
    count(Pathway, sort = TRUE)

# Plot the data as a bar chart
ggplot(pathway_counts, aes(x = reorder(Pathway, n), y = n, fill = Pathway)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip coordinates for better readability
    labs(title = "Gene Set Pathway Counts", x = "Pathway", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")  # Remove legend if not needed



# Split the genes column in the SNP file into individual genes
snp_data_long <- whr_list %>%
    separate_rows(V4, sep = ",") %>% 
    select(c(V1, V2, V3, V4))

# Intersect the genes with genesets
intersected_data <- snp_data_long %>%
    inner_join(clean_long_data, by = c("V4" = "Gene"))

# Aggregate the results to show which genesets each SNP may belong to
# aggregated_data <- intersected_data %>%
#     group_by(SNP) %>%
#     summarize(Genesets = paste(unique(Pathway), collapse = ","))

# Count the occurrences of each pathway
intersected_pathway_counts <- intersected_data %>%
    count(Pathway, sort = TRUE)

# Plot the data as a bar chart
ggplot(intersected_pathway_counts, aes(x = reorder(Pathway, n), y = n, fill = Pathway)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip coordinates for better readability
    labs(title = "WHR Gene Set Pathway Counts", x = "Pathway", y = "Count") +
    theme_minimal() +
    ylim(0,100)+
    theme(legend.position = "none")  # Remove legend if not needed

