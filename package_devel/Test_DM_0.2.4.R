# BioC Devel

# Created 21 Sep 2015
# - use unique snps per gene

##############################################################################################################


######################################################################
# SQTL - geuvadis data for chr1
######################################################################


setwd("/home/Shared/data/seq/geuvadis/")


library(DM)
library(ggplot2)

out_dir <- "dm_0_2_4_test/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

chr <- "19"


### Source all R files in DM package
Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
for(i in 1:length(Rfiles)) source(Rfiles[i])



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

genotypes_raw <- read.table(paste0(data_dir, "genotypes/snps_CEU_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(chr, nrow(genotypes_raw)), IRanges::IRanges(genotypes_raw$start, genotypes_raw$end))
names(snp_ranges) <- genotypes_raw$snpId
snp_id <- genotypes_raw$snpId

genotypes <- genotypes_raw[, -c(1:4)]


sample_id <- colnames(genotypes)

window <- 5e3

BPPARAM = BiocParallel::MulticoreParam(workers = 5)


### DM SQTL analysis

d <- dmSQTLdataFromRanges(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

plotData(d, paste0(out_dir, "chr",chr, "_"))



d <- dmFilter(d, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.05, max_features = Inf, minor_allel_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 5))


plotData(d, paste0(out_dir, "chr",chr, "_filtered_"))


save(d, file = paste0(out_dir, "chr",chr, "_d.Rdata"))
load(paste0(out_dir, "chr",chr, "_d.Rdata"))

BPPARAM = BiocParallel::MulticoreParam(workers = 10)


ds <- dmDispersion(d, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotDispersion(ds, paste0(out_dir, "chr",chr, "_"))


df <- dmFit(ds, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


save(df, file = paste0(out_dir, "chr",chr, "_df.Rdata"))

load(paste0(out_dir, "chr",chr, "_df.Rdata"))

gene_id <- "ENSG00000131037.8"
snp_id <- "snp_19_55603388"

plotFit(df, gene_id, snp_id, out_dir = paste0(out_dir, "chr",chr, "_"), plot_type = "boxplot1")




dt <- dmLRT(df, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotLRT(dt, out_dir = paste0(out_dir, "chr",chr, "_"))



plotFit(dt, gene_id, snp_id, out_dir = paste0(out_dir, "chr",chr, "_"), plot_type = "boxplot2")


res <- results(dt)

res <- res[order(res$pvalue, decreasing = FALSE), ]



plotFit(dt, gene_id = res[1:5, "gene_id"], snp_id = res[1:5, "snp_id"], out_dir = paste0(out_dir, "chr",chr, "_"), plot_type = "boxplot1")


table(res$df)


gene_id = res[1, "gene_id"]
snp_id = res[1, "snp_id"]
block_id <- res[1, "block_id"]


df@fit_full[[gene_id]][[block_id]]


















######## DM_0.2.3 plots of unique snps

d <- dmSQTLdataFromRanges(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window)

plotData(d, paste0(out_dir, "chr",chr, "_"))

nr_snps <- width(d@genotypes)

nr_snps_unique <- sapply(1:length(d@genotypes), function(i){
  
  nrow(unique(d@genotypes[[i]]))
  
})


sum(nr_snps)
sum(nr_snps_unique)



df <- data.frame(nr_snps = log10(nr_snps), nr_snps_unique = log10(nr_snps_unique), fraction = nr_snps_unique/nr_snps)

ggp <- ggplot(data = df, aes(x = nr_snps, y = nr_snps_unique)) + 
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed", size =  0.5)


ggp2 <- ggplot(data = df, aes(x = nr_snps, y = fraction)) + 
  geom_point(alpha = 0.5, size = 1.5) +
  
  
  pdf(paste0(out_dir, "chr",chr, "_", "nr_snps_unique.pdf"))

print(ggp)
print(ggp2)

dev.off()


### filtering

d <- dmFilter(d, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.05, max_features = Inf, minor_allel_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotData(d, paste0(out_dir, "chr",chr, "_filtered_"))


nr_snps <- width(d@genotypes)

nr_snps_unique <- sapply(1:length(d@genotypes), function(i){
  
  nrow(unique(d@genotypes[[i]]))
  
})


sum(nr_snps)
sum(nr_snps_unique)



df <- data.frame(nr_snps = log10(nr_snps), nr_snps_unique = log10(nr_snps_unique), fraction = nr_snps_unique/nr_snps)

ggp <- ggplot(data = df, aes(x = nr_snps, y = nr_snps_unique)) + 
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed", size =  0.5)


ggp2 <- ggplot(data = df, aes(x = nr_snps, y = fraction)) + 
  geom_point(alpha = 0.5, size = 1.5) +
  
  
  pdf(paste0(out_dir, "chr",chr, "_filtered_", "nr_snps_unique.pdf"))

print(ggp)
print(ggp2)

dev.off()







########################################################
# DS - data from Mark
########################################################


setwd("/home/Shared/tmp/gosia/")

library(DM)

load("first_hurdle.Rdata")



# ### Source all R files in DM package
# Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
# for(i in 1:length(Rfiles)) source(Rfiles[i])



counts = countsw
gene_id = lookup$gene_id[m]
feature_id = rownames(countsw)
sample_id = colnames(countsw)
group = as.factor(df$condition)


d <- dmDSdata(counts = as.matrix(countsw), gene_id = lookup$gene_id[m], feature_id = rownames(countsw), sample_id = colnames(countsw), group = as.factor(df$condition))


d <- dmFilter(d)


d <- d[1:2000, ]

# save(d, file = "d.Rdata")

plotData(d, out_dir = "./")


### dispersion

ds_none <- dmDispersion(d, verbose = TRUE, disp_moderation = "none", BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotDispersion(ds_none, out_dir = "./none_")




ds_trended <- dmDispersion(d, verbose = TRUE, disp_moderation = "trended", BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotDispersion(ds_trended, out_dir = "./trended_")


### fit and LRT

ds <- ds_none


df <- dmFit(ds, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

dt <- dmLRT(df, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


out_dir = "./"
pvalues <- res$pvalue



plotLRT(dt, out_dir = "./")


res <- results(dt)

res <- res[order(res$pvalue, decreasing = FALSE), ]

head(res)


plotFit(dt, gene_id = head(res)$gene_id, out_dir = "./")






















