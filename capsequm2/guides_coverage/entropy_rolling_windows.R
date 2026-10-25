library("tidyverse")
library("readr")
#library("Hmisc")

path_current <- "~/doctoral_work/2023/guides/guides_coverage/"
setwd(dir = path_current)

guides_raw <- read_tsv("guides_no_head_sorted.bed", col_names = F)
contigs_raw <- read_tsv("contigs_no_head_sorted.bed", col_names = F)
guides_in_context <- read_tsv("guides_in_context.txt", col_names = T)

### correct guides 
# magic number of subtracting 500 from the motifs - bc the coordinates were
# given incorrectly, all off by 500bp

guides_context_corrected <- guides_in_context
guides_context_corrected$motif_start <- guides_in_context$motif_start - 500
guides_context_corrected$motif_stop <- guides_in_context$motif_stop - 500

base_contig <- matrix(0, nrow = 500, ncol = 1)

test_tibble <- as_tibble(base_contig, .name_repair = "universal")

base_contig[guides_context_corrected$motif_start[1]:guides_context_corrected$motif_stop[1]] <- 1

### time for some loops
# convert into lapply functions later


list_of_matrices <- c()
list_of_contigs <- c()
temp_rows <- 0

for (i in seq(1, dim(contigs_raw)[1])){
    ### make matrices of appropriate size per contig
    
    temp_rows <- contigs_raw$X3[i]-contigs_raw$X2[i]
    
    list_of_matrices[i] <- tibble(matrix(0, nrow = temp_rows, ncol = 1))
    
    ### match contigs and guides, then assign 1
    
    for (j in seq(1, dim(guides_context_corrected)[1])){
        
        if(contigs_raw$X4[i] == guides_context_corrected$contig[j] &
           guides_context_corrected$motif_start[j] > 0 &
           guides_context_corrected$motif_start[j] < temp_rows &
           guides_context_corrected$motif_stop[j] < temp_rows){
            
            list_of_matrices[[i]][guides_context_corrected$motif_start[j]:guides_context_corrected$motif_stop[j]] <- list_of_matrices[[i]][guides_context_corrected$motif_start[j]:guides_context_corrected$motif_stop[j]] + 1
            
            
        }else if(contigs_raw$X4[i] == guides_context_corrected$contig[j] &
                 guides_context_corrected$motif_start[j] > 0 &
                 guides_context_corrected$motif_start[j] < temp_rows &
                 guides_context_corrected$motif_stop[j] >= temp_rows){

            list_of_matrices[[i]][guides_context_corrected$motif_start[j]:temp_rows] <- list_of_matrices[[i]][guides_context_corrected$motif_start[j]:temp_rows] + 1

        }
        
    }
    
   # list_of_matrices[i] <- as.data.frame(tibble(list_of_matrices[i]))
}

### ok then plot? coverage? 

#plot(y = list_of_matrices[[2]], x = seq(1:length(list_of_matrices[[2]])))


library(runner)
library(entropy)

sample_matrix <- list_of_matrices[[5]]


sample_runner <- runner(sample_matrix, 
       k = 10,
       f = entropy)

plot(y = sample_runner, x = seq(1:length(sample_runner)))

sample_dataframe <- as.data.frame(sample_runner)

#ggplot(sample_dataframe, aes(y = sample_runner, x = 1:dim(sample_dataframe)[1])) +
#    geom_line()


# reiterate across matrices

list_of_entropy <- c()

list_of_entropy <-  lapply(list_of_matrices, function(x) runner(x, k = 10, f = entropy))
plot(y = list_of_entropy[[2]], x = seq(1:length(list_of_entropy[[2]])))


### work out distances of guides within contig?

test_distances <- guides_context_corrected %>% 
    arrange(contig, desc(motif_start)) %>% 
    group_by(contig) %>% 
    mutate(differences = motif_start-lead(motif_start))

mean_distances <- test_distances %>%
    group_by(contig) %>% 
    summarise(mean = mean(differences, na.rm = T)) 


ggplot(mean_distances, aes(x = contig, y = mean))+
    geom_point()+
    geom_abline(intercept = mean(mean_distances$mean), slope = 0)+
    theme_minimal()
