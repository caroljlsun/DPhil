

setwd("comp_bio/capture_c_adipocytes/main_capture_analysis")

library(UMI4Cats)


# Make UMI-4C object including grouping by condition


my_col_data <- data.frame(sampleID = c("pADs_P1_43123", "mADs_M1_43123", "pADs_P2_43123", "mADs_M2_43123"),
                          condition = c("pre", "mature", "pre", "mature"),
                          replicate = c("1", "1", "2", "2"),
                          viewpoint = c("43123", "43123", "43123", "43123"),
                          file = c("umi_counts/pADs_P1_filled_sorted.tsv", "umi_counts/mADs_M1_filled_sorted.tsv", "umi_counts/pADs_P2_filled_sorted.tsv", "umi_counts/mADs_M2_filled_sorted.tsv"))

mu_umi <- makeUMI4C(
    colData = my_col_data,
    viewpoint_name = "43123",
    grouping = "condition",
    bait_expansion = 1e6,
    bait_exclusion = 15000
)

plotUMI4C(mu_umi, grouping="condition",TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, dgram_function = "quotient")

## 5) Get significant interactions 
# Generate windows
my_win_frags <- makeWindowFragments(mu_umi, n_frags=15, sliding = 0.8)

# Call interactions
my_gr <-  callInteractions(umi4c = mu_umi, 
                           design = ~condition, 
                           query_regions = my_win_frags, 
                           padj_threshold = 0.1, 
                           zscore_threshold=1)


# Plot interactions
my_all <- plotInteractionsUMI4C(mu_umi, my_gr,TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, grouping = NULL, significant=FALSE, xlim=c(15.5e6, 16.5e6)) # Plot all regions


my_sign <- plotInteractionsUMI4C(mu_umi, my_gr,TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, grouping = NULL, significant=TRUE, xlim=c(15.5e6, 16.5e6)) # Plot only significantly interacting regions
cowplot::plot_grid(my_all, my_sign, ncol=2, labels=c("All", "Significant"))

# Obtain unique significant interactions
my_inter <- getSignInteractions(my_gr)

## 6) Differential testing ----------------------
# Fisher test
umi_fisher <- fisherUMI4C(mu_umi, query_regions = my_inter,
                          filter_low = 20, 
                          grouping="condition")
plotUMI4C(umi_fisher, ylim = c(0, 50), grouping="condition")

# DESeq2 Wald Test
umi_wald <- waldUMI4C(umi4c=mu_umi,
                      query_regions = my_inter,
                      design=~condition)

plotUMI4C(umi_wald, ylim = c(0, 10),TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, grouping="condition")
