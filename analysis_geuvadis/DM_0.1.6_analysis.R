# BioC 3.1

# Created 8 July 2015
# Analyse data from DM_0.1.5_filtering


##############################################################################################################

setwd("/home/Shared/data/seq/geuvadis/")

library(DM)

library(BiocParallel)
library(limma)
source("/home/gosia/R/multinomial_project/package_devel/0_my_printHead.R")
library(edgeR)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)




# setwd("/Users/gosia/geuvadis/")
# 
# library(DM)
# library(BiocParallel)
# library(limma)
# source("/Users/gosia/R/multinomial_project/package_devel/0_my_printHead.R")
# library(edgeR)
# 
# library(ggplot2)
# library(reshape2)
# library(gridExtra)
# library(RColorBrewer)
# 
# Rfiles <- list.files("/Users/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
# for(i in Rfiles) source(i)


##########################################################################################
# Run the DMsQTL pipeline
##########################################################################################


######### run on DM_0_1_5 data 
out.dir <- "dm_0_1_6_analysis/dm_0_1_5_data/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

BPPARAM <- MulticoreParam(workers = 2, stop.on.error = TRUE, log = TRUE)


load(paste0("dm_0_1_5_data/dgeSQTL.RData"))

object.size(dgeSQTL)


data <- new("dmSQTLdata", counts = dgeSQTL$counts, genotypes = dgeSQTL$genotypes, samples = data.frame(sample_id = dgeSQTL$samples$sample_names))




xx <- do.call(rbind, data@counts)
object.size(data@counts)
object.size(xx)


yy <- do.call(rbind, data@genotypes)
object.size(data@genotypes)
object.size(yy)




### subset 100 genes
geneList <- names(dgeSQTL$counts)[1:10]

data <- new("dmSQTLdata", counts = dgeSQTL$counts[geneList], genotypes = dgeSQTL$genotypes[geneList], samples = data.frame(sample_id = dgeSQTL$samples$sample_names))

BPPARAM <- MulticoreParam(workers = 3)





## Input files: transcript expression, gene location and genotype information
data_dir <- "data/"


### read counts 
counts_path= paste0(data_dir, "expression/trExpCount_CEU.tsv")
counts_raw <- read.table(counts_path, header = TRUE, as.is = TRUE)

gene_id_counts <- counts_raw$geneId
feature_id_counts <- counts_raw$trId
counts <- as.matrix(counts_raw[, -c(1:2)])


### read genotypes
chr = "19"
genotypes_raw <- read.table(paste0(data_dir, "genotypes/snps_CEU_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)

snp_ranges <- GenomicRanges::GRanges(genotypes_raw$chr, IRanges::IRanges(genotypes_raw$start, genotypes_raw$end))
names(snp_ranges) <- genotypes_raw$snpId

snp_id_genotypes <- genotypes_raw$snpId

genotypes <- as.matrix(genotypes_raw[, -c(1:4)])


### read ranges
library(rtracklayer)
genes_path = paste0(data_dir, "annotation/genes_noChr.bed")
gene_ranges = import(genes_path)
names(gene_ranges) <- mcols(gene_ranges)$name


sample_id <- colnames(genotypes)

window <- 5e3




data <- dmSQTLdataFromRanges(counts, gene_id_counts, feature_id_counts, gene_ranges, genotypes, snp_id_genotypes, snp_ranges, sample_id, window = 5e3)

dmSQTL_plotData(data, out_dir = "./")


data_filter <- dmSQTL_filter(data, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.1, max_features = Inf, minor_allel_freq = 0.05, BPPARAM = MulticoreParam(workers = 5))


dmSQTL_plotData(data_filter, out_dir = "./filtering_")


library(pryr)
data_org <- data

BPPARAM <- MulticoreParam(workers = 3, stop.on.error = TRUE)

# tracemem(data)
data <- dm_estimateMeanExpression(data_org, BPPARAM = BPPARAM)







######### commonDispersion
common_dispersion <- dispersion <- dmSQTL_estimateCommonDispersion(data, disp_adjust = TRUE, disp_interval = c(0, 1e+5), disp_tol = 1e-01, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)


common_dispersion <- 6

dispersion <- dmSQTL_estimateTagwiseDispersion(data, disp_adjust = TRUE, disp_mode = c("optimize", "optim", "constrOptim", "grid")[4], disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = common_dispersion, disp_init_weirMoM = FALSE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = c("none", "common", "trended")[1], disp_prior_df = 10, disp_span = 0.3, prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)


dmSQTL_plotDispersion(dispersion, mean_expression, common_dispersion = common_dispersion, out_dir = "./")




fit <- dmSQTL_fit(data, dispersion, prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)



table <- dmSQTL_test(fit, verbose=FALSE, BPPARAM = BPPARAM)


dm_plotTable(table, out_dir = "./")




gene_id <- "ENSG00000007341.12"
snp_id <- "snp_1_113108919"


Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)


dmSQTL_plotFit(data, gene_id, snp_id, fit = fit, table = table, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[1], order = TRUE, plot_full = ifelse(is.null(fit), FALSE, TRUE), plot_null = ifelse(is.null(fit), FALSE, TRUE), out_dir = "./")


data@counts[[gene_id]]




dmSQTL_plotFit(data, gene_id, snp_id, fit = fit, table = table, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[5], order = TRUE, plot_full = ifelse(is.null(fit), FALSE, TRUE), plot_null = ifelse(is.null(fit), FALSE, TRUE), out_dir = "./pres_5")



##########################################################################################
# check if snps from sQTLseekeR paper are significant 
##########################################################################################


snp <- "snp_19_41937095"
gene <- "ENSG00000105341.11"


snp <- "snp_5_96244549"
gene <- "ENSG00000164308.12"













