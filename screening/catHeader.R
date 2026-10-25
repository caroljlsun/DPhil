# -------------------------------------------------------- #
# Procedure: Header function for RMDs
# Project: Obesity ML
# User: CZJS
# Date: MAR-2022
# -------------------------------------------------------- #

# made because I want many plots under tabsets
# because having a million plots isn't great
# but data exploration is all about that so


# yeah this is a function that helps you print in rmds a bunch 
# of tabsets

# intake is a list of plots you have made



easypackages::libraries("tidyverse",
                        "data.table"
)



# --------------------------------- #
# Function ----
# --------------------------------- #

# taken from stackexchange at : https://stackoverflow.com/questions/53444333/dynamic-tabsets-with-multiple-plots-r-markdown

# specify the number of levels you want, it's the previous header + 1

catHeader <- function(text = "", level = 3) {
    
    cat(paste0("\n\n", # Makes 2 new lines
               paste(rep("#", level),collapse = ""), #Pastes #'s "level" times
               " ", text, "\n")) # Names the new section on the same line, then
    #inserts a new line.
} 



# --------------------------------- #
# Example usage ----
# --------------------------------- #

# paste the following into an RMD file 

# ```{ r example}
# 
# source("comp_bio/ml_obesity/project_data/obesity_Imaging/sandbox/testouts/Imaging/functions/custom_distributions_1.R")
# 
# load("comp_bio/ml_obesity/project_data/obesity_Imaging/sandbox/testouts/Imaging/all_t2.RData")
# 
# 
# 
# to_plots <- c("BMI_2", "Sex", "Total_BMD_.bone_mineral_density.", "Total_tissue_fat_percentage", "WHR", "BMI_3", "Sex", "Total_BMD_.bone_mineral_density..1", "Total_tissue_fat_percentage.1", "WHR")[1:5]
# 
# plotcompile <- list()
# 
# for (i in 2:length(to_plots)) {
#   
#   name <- paste0(to_plots[i])
#   plotcompile[[name]] <- custom_dist(dataset = all_t2,
#                                      plotx = "BMI_2",
#                                      ploty = to_plots[i])
#   
# }
# 
# for (j in seq_along(plotcompile)) {
#   
#   tmp <- plotcompile[[j]]
#   catHeader(to_plots[j+1], 4)
#   print(tmp)
#   
# }
# 
# ```