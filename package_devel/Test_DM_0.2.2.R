# BioC Devel

# Created 26 Aug 2015


# library(devtools); library(GenomicRanges); library(BiocParallel); library(edgeR)

# Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE); for(i in Rfiles) source(i)



########################################################
# data from Mark
########################################################


setwd("/home/Shared/tmp/gosia/")

library(DM)

load("d.Rdata")

load("ds.Rdata")



load("first_hurdle.Rdata")


# o <- order(lookup$gene_id[m], lookup$transcript_id[m]) 

# d <- dmDSdata(counts = as.matrix(countsw[o,7:12]), gene_id_counts = lookup$gene_id[m][o], feature_id_counts = rownames(countsw)[o], sample_id = colnames(countsw)[7:12], group = as.factor(df$condition)[7:12])

# system.time(d <- dmDSfilter(d, max_features = 5))


d <- dmDSdata(counts = as.matrix(countsw), gene_id_counts = lookup$gene_id[m], feature_id_counts = rownames(countsw), sample_id = colnames(countsw), group = as.factor(df$condition))


d <- dmFilter(d)


# save(d, file = "d.Rdata")

plotData(d, out_dir = "./")


plot(d, out_dir = "./")


ds <- dmDispersion(d, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


# save(ds, file = "ds.Rdata")

plotDispersion(ds, out_dir = "./")



df <- dmFit(ds, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotFit(df, gene_id = "ENSG00000004059", plot_type = "barplot", out_dir = "./")



pairwise_comparison <- matrix(c(1, 2, 2, 3), nrow = 2, byrow = TRUE)

dt <- dmLRT(df, pairwise_comparison = pairwise_comparison, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotLRT(dt, out_dir = "./")

head(results(dt))


########################################################
# sqtl data
########################################################

# library(devtools); library(GenomicRanges); library(BiocParallel); library(edgeR); library(ggplot2); library(reshape2)

# Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE); for(i in Rfiles) source(i)

#######################################

library(DM)

setwd("/home/gosia/R/multinomial_project/package_devel/DM")


## Input files: transcript expression, gene location and genotype information
data_dir <- "data-raw/geuvadis/data/"


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


d <- dmSQTLdataFromRanges(counts, gene_id_counts, feature_id_counts, gene_ranges, genotypes, snp_id_genotypes, snp_ranges, sample_id, window = 5e3)


all(names(d@counts) == names(d@genotypes))


d <- dmFilter(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotData(d, "~/")


# save(d, file = "data-raw/d.Rdata")



ds <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


# save(ds, file = "data-raw/ds.Rdata")


plotDispersion(ds, "~/")


df <- dmFit(ds, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


snp_id <- "snp_19_54704760"
gene_id <- "ENSG00000170889.9"

plotFit(df, gene_id, snp_id, out_dir = "~/")



dt <- dmLRT(df, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotLRT(dt, out_dir = "~/")


plotFit(dt, gene_id, snp_id, out_dir = "~/", plot_type = "boxplot2")






























