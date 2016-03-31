
# R32dev

library(devtools)

load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")

##############################################################################
# GEUVADIS 
##############################################################################

library(DRIMSeq)
library(ggplot2)
library(limma)
library(GenomicRanges)

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/Shared/data/seq/geuvadis'
workers=10
population='CEU'
chr='22'


setwd(rwd)

out_dir <- "drimseq_0_3_3_analysis_f/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


data_dir <- "data/"

########################################################
# sqtl analysis per chromosome
########################################################

out_name <- paste0(out_dir, population, "_chr",chr, "_")


load(paste0(out_name, "d.Rdata"))


out_dir <- "drimseq_0_3_3_analysis_f/test_deviance/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


out_name <- paste0(out_dir, population, "_chr", chr, "_")




## LR test
dlr <- dmTest(d, test = "lr", BPPARAM = BiocParallel::MulticoreParam(workers = workers))

reslr <- results(dlr)

plotTest(dlr, out_dir = paste0(out_name, "lr_"))


## Fql test
dfql <- dmTest(d, test = "fql", BPPARAM = BiocParallel::MulticoreParam(workers = workers))

resfql <- results(dfql)

plotTest(dfql, out_dir = paste0(out_name, "fql_"))


## F test
df <- dmTest(d, test = "f", BPPARAM = BiocParallel::MulticoreParam(workers = workers))

resf <- results(df)

plotTest(df, out_dir = paste0(out_name, "f_"))






resm <- merge(reslr[, c("gene_id", "snp_id", "pvalue")], resf[, c("gene_id", "snp_id", "pvalue")], by = c("gene_id", "snp_id"), suffixes = c(".lr",".f"))

png(paste0(out_name, "scatter_pvalues_f.png"))

plot(resm[, "pvalue.lr"], resm[, "pvalue.f"], xlim = c(0, 1), ylim = c(0, 1), xlab = "LR", ylab = "F")
abline(a = 0, b = 1, col = "blue")

dev.off()




resm <- merge(reslr[, c("gene_id", "snp_id", "pvalue")], resfql[, c("gene_id", "snp_id", "pvalue")], by = c("gene_id", "snp_id"), suffixes = c(".lr",".fql"))

png(paste0(out_name, "scatter_pvalues_fql.png"))

plot(resm[, "pvalue.lr"], resm[, "pvalue.fql"], xlim = c(0, 1), ylim = c(0, 1), xlab = "LR", ylab = "Fql")
abline(a = 0, b = 1, col = "blue")

dev.off()


table(resfql$f < 0)

neg_genes <- resfql[resfql$f < 0, c("gene_id", "snp_id")]
neg_genes <- neg_genes[complete.cases(neg_genes), ]
dim(neg_genes)


genes <- neg_genes[!duplicated(neg_genes[, "gene_id"]), ]
dim(genes)



gene <- genes[1, "gene_id"]
snp <- genes[1, "snp_id"]



for(i in 1:10){
  
  plotFit(dfql, gene_id = genes[i, "gene_id"], snp_id = genes[i, "snp_id"], out_dir = paste0(out_name, "dfql", "_neg_dev_no_zeros_", i, "_"))
  
}


############################################################################


res <- results(d)

table(res$f < 0)

res <- res[order(res$f, decreasing = FALSE), ]

gene <- "ENSG00000131100.8"

d@fit_full[[gene]]
d@fit_full[[gene]]@metadata


d@fit_null[[gene]]@metadata


block <- "block_6"

d@fit_null[[gene]]@metadata[block, ]

genotypes <- d@genotypes[[gene]][block, ]


pi <- d@fit_null[[gene]][[block]][, "null"]
gamma0 <- d@genewise_dispersion[[gene]][block]
y <- d@counts[[gene]][, !is.na(genotypes)]


DRIMSeq:::dm_devG(pi = pi[-length(pi)], gamma0 = gamma0, y = y)

DRIMSeq:::dm_likG(pi = pi[-length(pi)], gamma0 = gamma0, y = y)

DRIMSeq:::dm_lik(pi = pi[-length(pi)], gamma0 = gamma0, y = y)



res[res$gene_id == gene & res$block_id == block, ]

snp <- "snp_22_18071444"


plotFit(d, gene_id = gene, snp_id = snp, plot_type = "boxplot1", out_dir = out_name)





# ### DRIMSeq SQTL analysis
# 
# d <- dmSQTLdataFromRanges(counts = counts_raw[, -c(1:2)], gene_id = counts_raw$geneId, feature_id = counts_raw$trId, gene_ranges = gene_ranges, genotypes = genotypes_raw[, -c(1:4)], snp_id = genotypes_raw$snpId, snp_ranges = snp_ranges, sample_id = colnames(genotypes_raw[, -c(1:4)]), window = 5e3, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# 
# d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5, min_samps_feature_prop = 0, minor_allele_freq = 5, min_gene_expr = 10, min_feature_expr = 10, min_feature_prop = 0, max_features = Inf, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# 
# plotData(d, out_dir = out_name)
# 
# 
# d <- dmDispersion(d, verbose = TRUE,  speed = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# 
# plotDispersion(d, out_dir = out_name)
# 
# 
# d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# 
# ## LR test
# d <- dmTest(d, test = "lr", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# plotTest(d, out_dir = paste0(out_name, "lr_"))
# 
# res <- results(d)
# 
# write.table(res, file = paste0(out_name, "results_lr.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# 
# res <- res[order(res$pvalue, decreasing = FALSE), ]
# 
# plotFit(d, gene_id = res[1, "gene_id"], snp_id = res[1, "snp_id"], out_dir = paste0(out_name, "lr_"))
# 
# 
# 
# ## F test
# d <- dmTest(d, test = "f", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# plotTest(d, out_dir = paste0(out_name, "f_"))
# 
# res <- results(d)
# 
# write.table(res, file = paste0(out_name, "results_f.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# 
# res <- res[order(res$pvalue, decreasing = FALSE), ]
# 
# plotFit(d, gene_id = res[1, "gene_id"], snp_id = res[1, "snp_id"], out_dir = paste0(out_name, "f_"))
# 
# save(d, file = paste0(out_name, "d.Rdata"))







