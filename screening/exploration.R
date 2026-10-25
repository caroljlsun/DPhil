#---------------------------------------------
# Procedure: Script to analyse harmony output
# Project: Screening
# User: CZJS
# Date: Jan-2025
#---------------------------------------------


# Setup -------------------------------------------------------------------


library(tidyverse)
library(data.table)
library(FactoMineR)
library(factoextra)
library(Hmisc)
library(bit64)
library(NbClust)
library(pheatmap)
library(umap)



setwd("comp_bio/capture_c_adipocytes/screening/")

# Import data ---------------------------------------------------

pilot_squared <- fread("Evaluation_2_PlateResults.txt", skip = 8) %>% 
    mutate(rowcol = paste0(Row, Column)) %>% 
    select(-"Time [s]", -Compound, -Concentration, -"Cell Type", -V21, -"Cell Count", -Timepoint, -"Number of Analyzed Fields")


gene_key <- fread("GENE_KEY.csv")%>% 
    mutate(rowcol = paste0(Row, Column)) %>% 
    select(-Row, -Column)

# Wrangle -----------------------------------------------------------------

pilot_squared <- merge(pilot_squared, gene_key, by = "rowcol")

genes_of_interest <- pilot_squared %>% 
    filter(str_detect(Gene, "NTC|PPARG|IRS1|PAPOLA|POLA1|PTEN"))
# Heatmap -----------------------------------------------------------------

# convert i64 type to numeric
genes_of_interest$`Nuclei - Intensity Cell Alexa 568 Sum - Sum per Well` <- genes_of_interest$`Nuclei - Intensity Cell Alexa 568 Sum - Sum per Well`/1

num_goi <- subset(genes_of_interest, select = -c(Gene, rowcol, Row, Column))

norm_goi <- scale(num_goi)
                  
goi_pca <- princomp(norm_goi)

summary(goi_pca)

fviz_eig(goi_pca, addlabels = T)

fviz_pca_var(goi_pca, col.var = "black")
fviz_pca_ind(goi_pca)

set.seed(100)

goi_PCA <- PCA(norm_goi)
goi_hcpc<-HCPC(goi_PCA, nb.clust=4,consol=FALSE,graph=FALSE)

plot.HCPC(goi_hcpc, choice = "tree", rect = T)

clustered_f <- goi_hcpc$data.clust

foo <- clustered_f%>% 
    group_by(clust) %>% 
    summarise(across( everything(), ~mean(.x, na.rm = T) ),  ) 

pheatmap(as.matrix(t(foo[,-1])),
         scale = "row",
         labels_col = foo$clust,
         angle_col = 0,
         cluster_cols = F)


qux <- subset(genes_of_interest, select = -c(rowcol, Row, Column))


row.names(qux) <- make.names(qux$Gene, T) 
qux_matrix <- as.matrix(qux[,-12], rownames = make.names(qux$Gene, T) )
norm_qux <- scale(qux_matrix)


pheatmap(norm_qux, angle_col = 45)


# jaccards similarity
# merge ntcs into one
# merge gene duplicates

# normalised and rownamed

genes_of_interest$`Nuclei - Intensity Cell Alexa 568 Sum - Sum per Well` <- genes_of_interest$`Nuclei - Intensity Cell Alexa 568 Sum - Sum per Well`/1

genes_of_interest <- as.data.frame(genes_of_interest)

rownames(genes_of_interest) <- make.names(genes_of_interest$Gene, unique = T)

num_goi <- genes_of_interest %>% select(-c(Gene, rowcol, Row, Column))

goi_pca <- PCA(num_goi, scale = T)

summary(goi_pca)

fviz_eig(goi_pca, addlabels = T)

fviz_pca_var(goi_pca, col.var = "black")
fviz_pca_ind(goi_pca)


# normalised, rownames, 2 way median polish

library(platetools)

# polish then select...

df <- data.frame(well = num_to_well(1:96),
                 vals = rnorm(96))

raw_map(data = df$vals,
        well = df$well,
        plate = 96)

b_score(data = df$vals,
        well = df$well,
        plate = 96)

b_map(data = df$vals,
      well = df$well,
      plate = 96)

pilot_squared$well <- paste0(LETTERS[pilot_squared$Row], sprintf("%02d", pilot_squared$Column))

raw_map(data = pilot_squared$`Lipid_droplets - Intensity Lipid_droplets Alexa 568 Mean - Mean per Well`, 
        well = pilot_squared$well, 
        plate = 96)


b_map(data = pilot_squared$`Lipid_droplets - Intensity Lipid_droplets Alexa 568 Mean - Mean per Well`, 
        well = pilot_squared$well, 
        plate = 96)


b_scored <- b_score(data = pilot_squared$`Lipid_droplets - Intensity Lipid_droplets Alexa 568 Mean - Mean per Well`, 
                    well = pilot_squared$well, 
                    plate = 96)


# Define the columns you want to apply the b_score function to
columns_to_score <- c(colnames(pilot_squared[,4:14]))

# Apply the b_score function to each column using lapply
b_scored_list <- lapply(columns_to_score, function(column) {
    b_score(data = pilot_squared[[column]], 
            well = pilot_squared$well, 
            plate = 96)
})

# Combine the results back into the dataframe
names(b_scored_list) <- columns_to_score
b_scored_df <- as.data.frame(b_scored_list) %>% 
    select(-matches("\\.well$"))

# Add the well column back to the dataframe
b_scored_df$well <- pilot_squared$well
b_scored_df$Gene <- pilot_squared$Gene

# Display the first few rows to check the result
head(b_scored_df)

b_scored_goi <- b_scored_df  %>% 
    filter(str_detect(Gene, "NTC|PPARG|IRS1|PAPOLA|POLA1|PTEN"))

rownames(b_scored_goi) <- make.names(b_scored_goi$Gene, unique = T)

num_goi <- b_scored_goi %>% select(-c(Gene, well))

goi_pca <- PCA(num_goi, scale = T)

summary(goi_pca)

fviz_eig(goi_pca, addlabels = T)

fviz_pca_var(goi_pca,axes = c(2,3), col.var = "black")
fviz_pca_ind(goi_pca, axes = c(1,2), col.ind = "cos2") +
    scale_color_gradient2(low="blue", mid="white", 
                          high="red", midpoint=0.50)+
    theme_minimal()

# take first 4 dimensions and then umap that


# normalised, rownamed and merged by gene duplicates, ntcs, using median


# Filter, group by Gene, and calculate the median for each group
median_df <- pilot_squared %>%
    filter(grepl("NTC|PPARG|IRS1|PAPOLA|POLA1|PTEN", Gene)) %>%
    mutate(Gene = gsub("NTC[0-9]+", "NTC", Gene)) %>%
    group_by(Gene) %>%
    summarize_all(median, na.rm = TRUE) %>% 
    as.data.frame()

row.names(median_df) <- make.names(median_df$Gene, T) 
med_matrix <- as.matrix(median_df[,-c(1,2,3,4, 16, 17)])
norm_med <- scale(med_matrix)


pheatmap(norm_med[,-9])



# UMAP --------------------------------------------------------------------

library(plotly)

qux_umap <- umap(norm_qux)

qux_layout <- qux_umap$layout
qux_layout <- data.frame(qux_layout)
qux_plot_data <- cbind(qux_layout, rownames(norm_qux))



fig <- plot_ly(qux_plot_data, x = ~X1, y = ~X2, color = ~rownames(norm_qux), type = 'scatter', mode = 'markers')%>%  
    layout(
        plot_bgcolor = "#e5ecf6",
        legend=list(title=list(text='features')), 
        xaxis = list( 
            title = "0"),  
        yaxis = list( 
            title = "1")) 


fig