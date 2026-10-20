


library(UMI4Cats)

## 0) Download example data -------------------------------
hg19_example_data <- downloadUMI4CexampleData(out_dir = "hg19_example_data")

## 1) Generate Digested genome ----------------------------
# The selected RE in this case is DpnII (|GATC), so the cut_pos is 0, and the res_enz "GATC".
hg19_dpnii <- digestGenome(
    cut_pos = 0,
    res_enz = "GATC",
    name_RE = "DpnII",
    ref_gen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    out_path = file.path("hg19_digested_genome/")
)



## 2) Process UMI-4C fastq files --------------------------
raw_dir <- file.path(path, "CIITA", "fastq")

contactsUMI4C(
    fastq_dir = raw_dir,
    wk_dir = file.path(path, "CIITA"),
    bait_seq = "GGACAAGCTCCCTGCAACTCA",
    bait_pad = "GGACTTGCA",
    res_enz = "GATC",
    cut_pos = 0,
    digested_genome = hg19_dpnii,
    bowtie_index = file.path(path, "ref_genome", "ucsc.hg19.chr16"),
    ref_gen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    threads = 5
)

## 3) Get filtering and alignment stats -------------------
statsUMI4C(wk_dir = file.path(path, "CIITA"))

## 4) Analyze the results ---------------------------------
# Load sample processed file paths
files <- list.files(file.path(path, "CIITA", "count"),
                    pattern = "*_counts.tsv.gz",
                    full.names = TRUE
)

# Create colData including all relevant information
colData <- data.frame(
    sampleID = gsub("_counts.tsv.gz", "", basename(files)),
    file = files,
    stringsAsFactors = FALSE
)

library(tidyr)
colData <- colData %>%
    separate(sampleID,
             into = c("condition", "replicate", "viewpoint"),
             remove = FALSE
    )

# Make UMI-4C object including grouping by condition
umi <- makeUMI4C(
    colData = colData,
    viewpoint_name = "CIITA",
    grouping = "condition",
    bait_expansion = 2e6
)

# Plot replicates
plotUMI4C(umi, grouping=NULL)

## 5) Get significant interactions 
# Generate windows
win_frags <- makeWindowFragments(umi, n_frags=8)

# Call interactions
gr <-  callInteractions(umi4c = umi, 
                        design = ~condition, 
                        query_regions = win_frags, 
                        padj_threshold = 0.01, 
                        zscore_threshold=2)

# Plot interactions
all <- plotInteractionsUMI4C(umi, gr, grouping = NULL, significant=FALSE, xlim=c(10.75e6, 11.1e6)) # Plot all regions
sign <- plotInteractionsUMI4C(umi, gr, grouping = NULL, significant=TRUE, xlim=c(10.75e6, 11.1e6)) # Plot only significantly interacting regions
cowplot::plot_grid(all, sign, ncol=2, labels=c("All", "Significant"))

# Obtain unique significant interactions
inter <- getSignInteractions(gr)

## 6) Differential testing ----------------------
# Fisher test
umi_fisher <- fisherUMI4C(umi, query_regions = inter,
                          filter_low = 20, 
                          grouping="condition")
plotUMI4C(umi_fisher, ylim = c(0, 10), grouping="condition")

# DESeq2 Wald Test
umi_wald <- waldUMI4C(umi4c=umi,
                      query_regions = inter,
                      design=~condition)

plotUMI4C(umi_wald, ylim = c(0, 10), grouping="condition")
