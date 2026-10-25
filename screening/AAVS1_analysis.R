library(tidyverse)
library(data.table)
library(FactoMineR)
library(factoextra)
library(Hmisc)
library(bit64)
library(NbClust)
library(pheatmap)
library(umap)
library(platetools)
library(plotly)
library(corrplot)
library(ggcorrplot)
library(heatmaply)


### Analyse the nuclei


plate_1 <- fread("AAVS1_nuclei_count.txt", skip = 8) %>% 
    mutate(well = paste0(LETTERS[Row], sprintf("%02d", Column)))

#removing the 1st and last column because it is empty

plate_1 <- plate_1[!plate_1$Column %in% c(1,12)]


raw_map(data = plate_1$`Nuclei - Number of Objects`, 
        well = plate_1$well, 
        plate = 96)



