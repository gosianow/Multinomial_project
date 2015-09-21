# BioC devel

# Run DM_0.2.3 analysis 

# Created 18 Sep 2015


setwd("/home/gosia/multinomial_project/simulations_sqtl3_hsapiens_noDE_noNull")


##############################################
### DM analysis 
##############################################

library(DM)
library(ggplot2)
library(GenomicRanges)


results_dir <- "dm_0_2_3/"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)


# simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, sep = "\t", as.is = TRUE)

# save(simulation_details, file = "3_truth/simulation_details.Rdata")


load("3_truth/simulation_details.Rdata")


counts <- simulation_details[, grep("isoformCount", colnames(simulation_details))]

ns <- length(grep("isoformCount", colnames(simulation_details)))

colnames(counts) <- paste0("s", 1:ns)

gene_id <- simulation_details$gene_id
feature_id <- simulation_details$transcript_id
sample_id <- colnames(counts)
group <- rep(c("G1", "G2"), times = c(20, 60))

### Create dmDSdata object
data <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)

plotData(data, out_dir = results_dir)

### Filtering
data <- dmFilter(data, min_samps_gene_expr = 60, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.1, max_features = Inf)


plotData(data, out_dir = paste0(results_dir, "filtered_"))




### Estimate dispersion
data <- dmDispersion(data, disp_interval = c(0, 1e+03), verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

# save(data, file = paste0(results_dir, "dmDSdispersion.Rdata"))
# load(paste0(results_dir, "dmDSdispersion.Rdata"))

plotDispersion(data, out_dir = results_dir)

### Fit full model
data <- dmFit(data, BPPARAM = BiocParallel::MulticoreParam(workers = 5))


### Fit null model & LR test
data <- dmLRT(data, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

plotLRT(data, out_dir = results_dir)

results <- results(data)

write.table(results, paste0(results_dir, "results.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(data, file = paste0(results_dir, "dmDSLRT.Rdata"))

### Plot fits for interesting genes
oo <- order(results$pvalue, decreasing = FALSE)
gene_id <- results[oo, "gene_id"][1:5]

plotFit(data, gene_id = gene_id, out_dir = results_dir)



##############################################
### Check filtering and truth
##############################################


load("3_truth/simulation_details.Rdata")


truth2 <- simulation_details[, c("gene_id", "transcript_id", "gene_ds_status", "transcript_ds_status")]

ds_info <- truth2[truth2$transcript_ds_status == 1, c(1, 2)]
colnames(ds_info) <- c("gene_id", "feature_id")


info <- counts(data)[, c(1,2)]

table(info$feature_id %in% ds_info$feature_id)
table(unique(info$gene_id) %in% unique(ds_info$gene_id))


DM:::dmDS_plotDataInfo(info, ds_info, out_dir = results_dir)






































