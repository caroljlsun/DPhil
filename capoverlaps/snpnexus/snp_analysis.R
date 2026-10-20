library(tidyverse)
library(data.table)

setwd("comp_bio/capture_c_adipocytes/capoverlaps/snpnexus/")


foo <- fread("t2d_pervariant.tsv")
buzz <- fread("pathway_t2d.txt")


foo_cleaned <- foo %>% 
    select(id, UCSC) %>% 
    mutate(UCSC = str_remove(UCSC, "^CONSEQ="), 
           annotation = str_extract(UCSC, "^[^(]+")) %>% 
    mutate(clean_annotation = str_replace(annotation, "non-coding intronic", "intronic"))

ggplot(foo_cleaned, aes(x = clean_annotation)) +
    geom_bar() +
    labs(title = "Bar Plot of T2D snp annotations", x = "Annotation", y = "Count") 

## counted df

bar <- as.data.table(table(foo_cleaned$clean_annotation))

bar <- bar %>%
    group_by(V1) %>%
    summarise(N = sum(N)) %>%
    mutate(percentage = N / sum(N) * 100)

library(plotly)
library(viridis)

colors <- inferno(length(unique(bar$V1)))

p <- plot_ly(bar, labels = ~V1, values = ~N, type = 'pie',textposition = 'outside',textinfo = 'label+percent', marker = list(colors = colors), textfont = list(size = 18)) %>% 
    layout(legend = list(font = list(size = 18)))
p
top10_buzz <- buzz %>%
    arrange(`p-Value`) %>%  # Sort by p-value in ascending order
    slice_head(n = 10)      # Select the top 10 rows

# Plot the data
ggplot(top10_buzz, aes(x = reorder(Description, -`p-Value`), y = `p-Value`, fill = `Parent(s)`)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip coordinates for better readability
    labs(title = "Plot of Most Significant T2D GO Annotations", x = "GO", y = "p-value") +
    theme_minimal() 
# foo_cleaned <- foo %>%
#     select(id, UCSC) %>%
#     mutate(
#         UCSC = str_remove(UCSC, "^CONSEQ="),  # Remove the CONSEQ= prefix
#         annotation = str_extract_all(UCSC, "[^,]+(?=\\()")) %>% 
#     select(-UCSC)
# 
# df <- foo_cleaned %>%
#     mutate(
#         class = case_when(
#             sapply(annotation, function(x) all(x == "intronic")) ~ "non-coding",
#             sapply(annotation, function(x) any(x %in% c("exonic", "coding"))) ~ "coding",
#             sapply(annotation, function(x) all(x == "")) ~ NA,
#             TRUE ~ "non-coding"
#         )
#     )




foo <- fread("whr_pervariant.tsv")
buzz <- fread("pathway_whr.txt")


foo_cleaned <- foo %>% 
    select(id, UCSC) %>% 
    mutate(UCSC = str_remove(UCSC, "^CONSEQ="), 
           annotation = str_extract(UCSC, "^[^(]+")) %>% 
    mutate(clean_annotation = str_replace(annotation, "non-coding intronic", "intronic"))

ggplot(foo_cleaned, aes(x = clean_annotation)) +
    geom_bar() +
    labs(title = "Bar Plot of WHR snp annotations", x = "Annotation", y = "Count") 

top10_buzz <- buzz %>%
    arrange(`p-Value`) %>%  # Sort by p-value in ascending order
    slice_head(n = 10)      # Select the top 10 rows

# Plot the data
ggplot(top10_buzz, aes(x = reorder(Description, -`p-Value`), y = `p-Value`, fill = `Parent(s)`)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip coordinates for better readability
    labs(title = "Plot of Most Significant WHR GO Annotations", x = "GO", y = "p-value") +
    theme_minimal() 
