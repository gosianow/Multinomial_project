# BioC Devel

# Created 10 Sep 2015

##############################################################################################################

setwd("/home/Shared/data/seq/geuvadis/")


library(DM)
library(ggplot2)

out_dir <- "dm_0_2_3_analysis/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


########################################################
# sqtl data for chr1
########################################################

## Input files: transcript expression, gene location and genotype information
data_dir <- "data/"


### read counts 
counts_path <- paste0(data_dir, "expression/trExpCount_CEU.tsv")
counts_raw <- read.table(counts_path, header = TRUE, as.is = TRUE)

gene_id <- counts_raw$geneId
feature_id <- counts_raw$trId
counts <- counts_raw[, -c(1:2)]


### read ranges

genes_path = paste0(data_dir, "annotation/genes_noChr.bed")
gene_ranges = rtracklayer::import(genes_path)
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name


### read genotypes
chr <- "19"
genotypes_raw <- read.table(paste0(data_dir, "genotypes/snps_CEU_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(chr, nrow(genotypes_raw)), IRanges::IRanges(genotypes_raw$start, genotypes_raw$end))
names(snp_ranges) <- genotypes_raw$snpId
snp_id <- genotypes_raw$snpId

genotypes <- genotypes_raw[, -c(1:4)]


sample_id <- colnames(genotypes)

window <- 5e3


### DM SQTL analysis

d <- dmSQTLdataFromRanges(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window)


all(names(d@counts) == names(d@genotypes))


d <- dmFilter(d, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.05, max_features = Inf, minor_allel_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotData(d, paste0(out_dir, "chr",chr, "_"))

save(d, file = paste0(out_dir, "chr",chr, "_d.Rdata"))


d <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotDispersion(d, paste0(out_dir, "chr",chr, "_"))


d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


d <- dmLRT(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotLRT(d, out_dir = paste0(out_dir, "chr",chr, "_"))


save(d, file = paste0(out_dir, "chr",chr, "_d.Rdata"))



chr <- "19"
load(paste0(out_dir, "chr",chr, "_d.Rdata"))



# plotFit(d, gene_id, snp_id, out_dir = paste0(out_dir, "chr",chr, "_"), plot_type = "boxplot1")























