#---------------------------------------------
# Procedure: Script to run SPADE
# Project: Screening
# User: CZJS
# Date: Feb-2025
#---------------------------------------------


# Setup -------------------------------------------------------------------
install.packages('biclust')
install.packages('data.table')
install.packages('diptest')
install.packages('evtree')
install.packages('ggdendro')
install.packages("ggfortify")
install.packages('ggplot2')
install.packages('ggrepel')
install.packages('ggRandomForests')
install.packages('gridExtra')
install.packages('gtable')
install.packages('gtools')
install.packages('igraph')
install.packages('MASS')
install.packages('packcircles')
install.packages('plyr')
install.packages("randomForestSRC")
install.packages('reshape2')

install.packages("BiocManager")
BiocManager::install()
BiocManager::install("flowCore")

#source("http://bioconductor.org/biocLite.R") #for older R versions
#biocLite(suppressUpdates = TRUE) #for older R versions
#biocLite("flowCore", suppressUpdates = TRUE) #for older R versions

install.packages('devtools')
library("devtools")
install_github('tchitchek-lab/SPADEVizR')

library("SPADEVizR")

# probably need to ask hpc people for help in installing the packages...
