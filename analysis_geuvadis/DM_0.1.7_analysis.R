# BioC 3.1

# Created 31 July 2015
# Analyse data from DM_0.1.5_filtering


##############################################################################################################

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


setwd("/home/Shared/data/seq/geuvadis/")

library(GenomicRanges)
library(BiocParallel)
library(edgeR)


library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

library(rtracklayer)

library(DM)


Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)


##########################################################################################
# Run the DMsQTL pipeline
##########################################################################################


######### run on DM_0_1_7 data 
out_dir <- "dm_0_1_7_analysis/dm_0_1_7_data/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)



## Input files: transcript expression, gene location and genotype information
data_dir <- "data/"


### read counts 
counts_path <- paste0(data_dir, "expression/trExpCount_CEU.tsv")
counts_raw <- read.table(counts_path, header = TRUE, as.is = TRUE)

gene_id_counts <- counts_raw$geneId
feature_id_counts <- counts_raw$trId
counts <- as.matrix(counts_raw[, -c(1:2)])


### read genotypes
chr = "19"
genotypes_raw <- read.table(paste0(data_dir, "genotypes/snps_CEU_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(chr, nrow(genotypes_raw)), IRanges::IRanges(genotypes_raw$start, genotypes_raw$end))
names(snp_ranges) <- genotypes_raw$snpId


snp_id_genotypes <- genotypes_raw$snpId

genotypes <- as.matrix(genotypes_raw[, -c(1:4)])


### read ranges

genes_path = paste0(data_dir, "annotation/genes_noChr.bed")
gene_ranges = rtracklayer::import(genes_path)
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name


sample_id <- colnames(genotypes)

window <- 5e3






data_org <- dmSQTLdataFromRanges(counts, gene_id_counts, feature_id_counts, gene_ranges, genotypes, snp_id_genotypes, snp_ranges, sample_id, window = 5e3)


# data <- dmSQTLdata(counts, gene_id_counts, feature_id_counts, genotypes, gene_id_genotypes, snp_id_genotypes, sample_id)


# x <- data; min_samps_gene_expr = 70; min_gene_expr = 1; min_samps_feature_prop = 5; min_feature_prop = 0.1; max_features = Inf; minor_allel_freq = 0.2; BPPARAM = MulticoreParam(workers = 10)

data <- dmSQTLfilter(data_org, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.1, max_features = Inf, minor_allel_freq = 0.05, BPPARAM = MulticoreParam(workers = 10))


# data@genotypes <- MatrixList(lapply(data@genotypes, function(g){ g[1, , drop = FALSE]}))


dmSQTLplotData(data, out_dir = out_dir)


save(data, file = paste0(out_dir, "data.RData"))

load(paste0(out_dir, "data.RData"))

# x <- data; mean_expression = TRUE; common_dispersion = TRUE; tagwise_dispersion = TRUE; disp_adjust = TRUE; disp_mode = c("optimize", "optim", "constrOptim", "grid")[4]; disp_interval = c(0, 1e+5); disp_tol = 1e-08; disp_init = 100; disp_init_weirMoM = TRUE; disp_grid_length = 21; disp_grid_range = c(-10, 10); disp_moderation = c("none", "common", "trended")[1]; disp_prior_df = 10; disp_span = 0.3; prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2]; prop_tol = 1e-12; verbose = FALSE; BPPARAM = MulticoreParam(workers = 10)


data_disp <- dmSQTLdispersion(data, mean_expression = TRUE, common_dispersion = TRUE, tagwise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = c("optimize", "optim", "constrOptim", "grid")[4], disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = c("none", "common", "trended")[1], disp_prior_df = 10, disp_span = 0.3, prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers = 10))




save(data_disp, file = paste0(out_dir, "data_disp.RData"))


dmSQTLplotDispersion(data_disp, out_dir = out_dir)


# x <- data_disp; dispersion = "tagwise_dispersion"; prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2]; prop_tol = 1e-12; verbose = FALSE; BPPARAM = MulticoreParam(workers=10)


data_fit <- dmSQTLfit(data_disp, dispersion = "tagwise_dispersion", prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers = 10))



snp_id <- "snp_19_55578437"
gene_id <- "ENSG00000131037.8"


# x <- data_fit; plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[2]; order = TRUE; plot_full = TRUE; plot_nunll = TRUE

dmSQTLplotFit(data_fit, gene_id, snp_id, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[2], order = TRUE, plot_full = TRUE, plot_nunll = TRUE, out_dir = out_dir)


snp_id <- "snp_19_58785739"



data_test <- dmSQTLtest(data_fit, BPPARAM = MulticoreParam(workers=10))


dmSQTLplotTest(data_test, out_dir = out_dir)


snp_id <- "snp_19_55578437"
gene_id <- "ENSG00000131037.8"

dmSQTLplotFit(data_test, gene_id, snp_id, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[2], order = TRUE, plot_full = TRUE, plot_nunll = TRUE, out_dir = out_dir)


##########################################################################################
# check if snps from sQTLseekeR paper are significant 
##########################################################################################


snp <- "snp_19_41937095"
gene <- "ENSG00000105341.11"


snp <- "snp_5_96244549"
gene <- "ENSG00000164308.12"













