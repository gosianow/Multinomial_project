# BioC Devel

# Created 19 Aug 2015


# library(devtools); library(GenomicRanges); library(BiocParallel); library(edgeR)

# Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE); for(i in Rfiles) source(i)


########################################################
# simulation data
########################################################

setwd("/home/gosia/R/multinomial_project/package_devel/DM")


# create metadata file
metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata



htseq <- read.table("data-raw/simulations_sim5_drosophila_noDE_noNull/2_counts/dexseq_nomerge/htseq_counts.txt", header = TRUE, as.is = TRUE)
htseq <- htseq[!grepl(pattern = "_", htseq$group_id), ]


counts <- as.matrix(htseq[,-1])
group_id <- htseq[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id_counts <- group_split[, 1]
feature_id_counts <- group_split[, 2]
sample_id = metadata$sample_id
group = metadata$group


counts_spl <- split.data.frame(counts, gene_id_counts)

x <- counts_spl
 
object <- MatrixList(x)


object.size(counts_spl)

object.size(object)


########################################################
# data from Mark 
########################################################


setwd("/home/Shared/tmp/gosia/")

load("first_hurdle.Rdata")

mixx <- sample.int(length(lookup$gene_id))

countsw <- countsw[lookup$transcript_id, ]
counts <- countsw[mixx, ]
gene_id_counts <- lookup$gene_id[mixx]
feature_id_counts <- lookup$transcript_id[mixx]
sample_id = df$sample
group = df$condition


d <- dmDSdata(counts, gene_id_counts, feature_id_counts, sample_id, group)


system.time(d <- dmDSfilter(d))


dd <- dmDSdispersion(x = d, mean_expression = TRUE, common_dispersion = TRUE, tagwise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 10, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))




########################################################
# data from Mark - subset 
########################################################


setwd("/home/Shared/tmp/gosia/")

library(DM)

load("first_hurdle.Rdata")

o <- order(lookup$gene_id[m], lookup$transcript_id[m]) 


d <- dmDSdata(counts = as.matrix(countsw[o,7:12]), gene_id_counts = lookup$gene_id[m][o], 
             feature_id_counts = rownames(countsw)[o], sample_id = colnames(countsw)[7:12], 
             group = as.factor(df$condition)[7:12])



system.time(d <- dmDSfilter(d, max_features = 5))



ds <- dmDSdispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

dmDSplotDispersion(ds, out_dir = "./")



df <- dmDSfit(ds, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


dmDSplotFit(df, gene_id = "ENSG00000203685", plot_type = "barplot", out_dir = "./")


# save(df, file = "df.Rdata")

dt <- dmDStest(df)


head(dt@table)


########################################################
# sqtl data
########################################################

library(devtools); library(GenomicRanges); library(BiocParallel); library(edgeR); library(ggplot2); library(reshape2)

Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE); for(i in Rfiles) source(i)

#######################################

library(DM)

setwd("/home/gosia/R/multinomial_project/package_devel/DM")

load("data-raw/data_sqtl.Rdata")





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

library(rtracklayer)

genes_path = paste0(data_dir, "annotation/genes_noChr.bed")
gene_ranges = rtracklayer::import(genes_path)
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name


sample_id <- colnames(genotypes)

window <- 5e3




d <- dmSQTLdataFromRanges(counts, gene_id_counts, feature_id_counts, gene_ranges, genotypes, snp_id_genotypes, snp_ranges, sample_id, window = 5e3)



all(names(d@counts) == names(d@genotypes))




d <- dmSQTLfilter(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

dmSQTLplotData(d, "~/")


# save(d, file = "data-raw/data_sqtl.Rdata")


d <- dmSQTLdispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


dmSQTLplotDispersion(d, "~/")


df <- dmSQTLfit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


snp_id <- "snp_19_54704760"
gene_id <- "ENSG00000170889.9"

dmSQTLplotFit(df, gene_id, snp_id, "~/")



dt <- dmSQTLtest(df, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


dmSQTLplotTest(dt, out_dir = "~/")


dmSQTLplotFit(dt, gene_id, snp_id, out_dir = "~/")


######## data from package

### counts
 head(dataSQTL_counts)
 counts <- as.matrix(dataSQTL_counts[, -1])

 group_id <- dataSQTL_counts[,1]
 group_split <- limma::strsplit2(group_id, ":")
 gene_id_counts <- group_split[, 1]
 feature_id_counts <- group_split[, 2]

 ### gene_ranges
 dataSQTL_gene_ranges
 gene_ranges <- dataSQTL_gene_ranges
 names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name


 ### genotypes
 head(dataSQTL_genotypes)
 genotypes <- as.matrix(dataSQTL_genotypes[, -(1:4)])

 snp_id_genotypes <- dataSQTL_genotypes$snp_id

 snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(dataSQTL_genotypes$chr),
 IRanges::IRanges(dataSQTL_genotypes$start, dataSQTL_genotypes$end))
 names(snp_ranges) <- dataSQTL_genotypes$snp_id

 all(colnames(counts) == colnames(genotypes))

 sample_id <- colnames(counts)

 ### create dmSQTLdata object
 data <- dmSQTLdataFromRanges(counts, gene_id_counts, feature_id_counts,
 gene_ranges, genotypes, snp_id_genotypes, snp_ranges, sample_id,
 window = 5e3)



 dmSQTLplotData(data, "~/")


data <- dmSQTLfilter(data, BPPARAM = BiocParallel::MulticoreParam(workers = 5))


dmSQTLplotData(data, "~/")

## Not run:

### This part is quite time consuming thus if possible, increase the number of cores.
data <- dmSQTLdispersion(data, BPPARAM = BiocParallel::MulticoreParam(workers = 5))
## End(Not run)



dmSQTLplotDispersion(data, "~/")

d <- dmSQTLfit(data, BPPARAM = BiocParallel::MulticoreParam(workers = 5))


snp_id <- "snp_19_12796435"
gene_id <- "ENSG00000132004.7"

dmSQTLplotFit(d, gene_id, snp_id)



d <- dmSQTLtest(d)


dmSQTLplotTest(d, "~/")






























