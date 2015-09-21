# BioC Devel

# Created 10 Sep 2015
# Updated 21 Sep 2015
# - test DS analysis: data from Mark 

##############################################################################################################


########################################################
# data from Mark
########################################################


setwd("/home/Shared/tmp/gosia/")

library(DM)

load("first_hurdle.Rdata")





# o <- order(lookup$gene_id[m], lookup$transcript_id[m]) 

# d <- dmDSdata(counts = as.matrix(countsw[o,7:12]), gene_id_counts = lookup$gene_id[m][o], feature_id_counts = rownames(countsw)[o], sample_id = colnames(countsw)[7:12], group = as.factor(df$condition)[7:12])

# system.time(d <- dmDSfilter(d, max_features = 5))


counts = countsw
gene_id = lookup$gene_id[m]
feature_id = rownames(countsw)
sample_id = colnames(countsw)
group = as.factor(df$condition)


d <- dmDSdata(counts = as.matrix(countsw), gene_id = lookup$gene_id[m], feature_id = rownames(countsw), sample_id = colnames(countsw), group = as.factor(df$condition))


d <- dmFilter(d)


# save(d, file = "d.Rdata")

plotData(d, out_dir = "./")


d <- d[, 1:6]



ds_none <- dmDispersion(d, verbose = TRUE, disp_moderation = "none", BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotDispersion(ds, out_dir = "./none_")


ds_trended <- dmDispersion(d, verbose = TRUE, disp_moderation = "trended", BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotDispersion(ds, out_dir = "./trended_")



# save(ds, file = "ds.Rdata")

plotDispersion(ds, out_dir = "./")



df <- dmFit(ds, BPPARAM = BiocParallel::MulticoreParam(workers = 10))



plotFit(df, gene_id = "ENSG00000004059", plot_type = "barplot", out_dir = "./")



pairwise_comparison <- matrix(c(1, 2, 2, 3), nrow = 2, byrow = TRUE)

dt <- dmLRT(df, pairwise_comparison = pairwise_comparison, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotLRT(dt, out_dir = "./")

head(results(dt))

















######################################################################
# sqtl data for chr19
######################################################################


setwd("/home/Shared/data/seq/geuvadis/")


library(DM)
library(ggplot2)

out_dir <- "dm_0_2_3_test/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)



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


########## plot the number of unique snps

chr <- "19"
load(paste0(out_dir, "chr",chr, "_d.Rdata"))


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


####### find blocks of equivalent snps 

i = 2

dim(d@genotypes[[i]])
dim(unique(d@genotypes[[i]]))


matching_snps <- match(data.frame(t(d@genotypes[[i]])), data.frame(t(d@genotypes[[i]])))

d@genotypes[[i]][matching_snps == 12, ]





########################################################
# Check the dispersion estimates when using null model
########################################################


d <- dmSQTLdataFromRanges(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window)

d <- dmFilter(d, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.05, max_features = Inf, minor_allel_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotData(d, paste0(out_dir, "chr",chr, "_"))


dn <- d

genotypes <- dn@genotypes

### keep only one SNP per gene
inds <- 1:length(genotypes)
names(inds) <- names(genotypes)

genotypesn <- DM:::MatrixList(lapply(inds, function(g){ 
# g = 1
snp <- genotypes[[g]][1, , drop = FALSE]
snp[,] <- 1
return(snp)
 }))


dn@genotypes <- genotypesn


dn <- dmDispersion(dn, common_dispersion = TRUE, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotDispersion(dn, paste0(out_dir, "chr",chr, "_dn_"))


save(dn, file = paste0(out_dir, "chr",chr, "_dn.Rdata"))




chr <- "19"
load(paste0(out_dir, "chr",chr, "_d.Rdata"))


x <- d

genewise_dispersion = x@genewise_dispersion; mean_expression = x@mean_expression; nr_features = width(x@counts); common_dispersion = x@common_dispersion

w <- sapply(genewise_dispersion, length)

mean_expression <- rep(mean_expression, w)
nr_features <- rep(nr_features, w)

df <- data.frame(mean_expression = log10(mean_expression + 1), dispersion = log10(unlist(genewise_dispersion)), nr_features = nr_features)


wn <- sapply(dn@genewise_dispersion, length)

dfn <- data.frame(mean_expression = log10(rep(dn@mean_expression, wn) + 1), dispersion = log10(unlist(dn@genewise_dispersion)), nr_features = rep(width(dn@counts), wn))


df_quant <- quantile(na.omit(df$nr_features), probs = 0.95)

ggp2 <- ggplot(df, aes(x = mean_expression, y = dispersion, colour = nr_features )) +
theme_bw() +
xlab("Log10(mean expression)") +
ylab("Log10(dispersion)") +
geom_point(size = 1.5, alpha = 0.5, na.rm = TRUE) +
theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") +
guides(colour = guide_colorbar(barwidth = 20, barheight = 0.5)) +
scale_colour_gradient(limits = c(2, df_quant), breaks = seq(2, df_quant, 1), low = "royalblue2", high="red2", name = "Number of features", na.value = "red2") +
geom_point(data = dfn, aes(x = mean_expression, y = dispersion), colour = "black", size = 1)


if(length(common_dispersion)){
  ggp2 <- ggp2 + geom_hline(yintercept = log10(common_dispersion), colour = "black", linetype = "dashed", size =  0.5)
}


pdf(paste0(out_dir, "chr",chr, "_", "dispersion_vs_mean_2.pdf"))

print(ggp2)

dev.off()



########################################################
# Check the fittings using null model dispersion
########################################################


dndisp <- unlist(dn@genewise_dispersion)
dnmean <- dn@mean_expression

which(dndisp > 150)

dim(d@genotypes[["ENSG00000105193.3"]])

dim(d@genotypes[["ENSG00000105379.4"]])


### one with max disp
which.max(unlist(d@genewise_dispersion))

table(d@genotypes[["ENSG00000105379.4"]]["snp_19_51852143", ])

plotFit(d, gene_id = "ENSG00000105379.4", snp_id = "snp_19_51852143", out_dir = paste0(out_dir, "chr",chr, "_"), plot_type = "boxplot1")


### one with min disp
which.min(unlist(d@genewise_dispersion[["ENSG00000105379.4"]]))


table(d@genotypes[["ENSG00000105379.4"]]["snp_19_51857093", ])

plotFit(d, gene_id = "ENSG00000105379.4", snp_id = "snp_19_51857093", out_dir = paste0(out_dir, "chr",chr, "_"), plot_type = "boxplot1")






### plot the error

full_dispersion = unlist(d@genewise_dispersion)
w <- sapply(d@genewise_dispersion, length)

dndisp <- unlist(dn@genewise_dispersion)
names(dndisp) <- limma::strsplit2(names(dndisp), ".snp_")[, 1]

setdiff(names(dndisp), names(w))


null_dispersion = rep(dndisp[names(w)], w)

error <- full_dispersion - null_dispersion

df <- data.frame(null_dispersion = null_dispersion, error = error)


ggp <- ggplot(data = df, aes(x = null_dispersion, y = error)) +
geom_point(alpha = 0.5, size = 1.5) +
geom_hline(yintercept = 0, colour = "red", linetype = "dashed", size =  0.5)

pdf(paste0(out_dir, "chr",chr, "_", "disp_error.pdf"))

print(ggp)

dev.off()






d2 <- d

d2@common_dispersion <- dn@common_dispersion

w <- sapply(d@genewise_dispersion, length)
d2@genewise_dispersion <- relist(rep(unlist(dn@genewise_dispersion), w), d@genewise_dispersion)


plotDispersion(d2, paste0(out_dir, "chr",chr, "_d2_"))


d2 <- dmFit(d2, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


d2 <- dmLRT(d2, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotLRT(d2, out_dir = paste0(out_dir, "chr",chr, "_d2_"))


save(d2, file = paste0(out_dir, "chr",chr, "_d2.Rdata"))



chr <- "19"
load(paste0(out_dir, "chr",chr, "_d2.Rdata"))




resd <- results(d)
resd2 <- results(d2)

df <- data.frame(pvalues_full = (resd$pvalue), pvalues_null = (resd2$pvalue))

ggp <- ggplot(data = df, aes(x = pvalues_full, y = pvalues_null)) + 
geom_point(alpha = 0.5, size = 1.5) +
geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed", size =  0.5)


pdf(paste0(out_dir, "chr",chr, "_", "diff_pvalues.pdf"))

print(ggp)

dev.off()






########################################################
# Keep unique snps
########################################################


d <- dmSQTLdataFromRanges(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window)

d <- dmFilter(d, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.05, max_features = Inf, minor_allel_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

d3 <- d

### keep only unique SNPs per gene
inds <- 1:length(d@genotypes)
names(inds) <- names(d@genotypes)

genotypes3 <- DM:::MatrixList(lapply(inds, function(g){ 
# g = 1
unique(d@genotypes[[g]])
 }))


d3@genotypes <- genotypes3

plotData(d3, paste0(out_dir, "chr",chr, "_d3_"))

save(d3, file = paste0(out_dir, "chr",chr, "_d3.Rdata"))


d3 <- dmDispersion(d3, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotDispersion(d3, paste0(out_dir, "chr",chr, "_d3_"))


d3 <- dmFit(d3, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


d3 <- dmLRT(d3, BPPARAM = BiocParallel::MulticoreParam(workers = 10))


plotLRT(d3, out_dir = paste0(out_dir, "chr",chr, "_d3_"))

save(d3, file = paste0(out_dir, "chr",chr, "_d3.Rdata"))






















