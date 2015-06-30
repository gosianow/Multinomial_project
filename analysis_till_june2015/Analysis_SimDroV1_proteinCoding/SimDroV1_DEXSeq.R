# BioC 3.0
# Crated 3 Feb 2015
# Midified

# Run DEXSeq on bitseq counts

#####################################################################

setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding")

# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


out.dir <- "DEXSeq_1.10.8/bitseq_1.10.0/"
dir.create(out.dir, showWarnings=F, recursive=T)


library("DEXSeq")
library("edgeR")
library("BiocParallel")


nCores <- 10
BPPARAM = MulticoreParam(workers=nCores)
countfiles <-  paste0("BitSeq_1.10.0/BitSeq_counts", metadata$SampleName1, ".txt")


dxr   <- DEXSeqDataSetFromHTSeq(countfiles = countfiles, sampleData = metadata , design = ~ sample + exon + condition:exon)


dxr  <- estimateSizeFactors(dxr)
dxr  <- estimateDispersions(dxr  ,BPPARAM=BPPARAM )
dxr   <- testForDEU(dxr ,BPPARAM=BPPARAM)


DEXSeq_res_exon <-  DEXSeqResults(dxr)

DEXSeq_res_exon <- as.data.frame(DEXSeq_res_exon[, 1:7])
DEXSeq_res_exon <- data.frame(groupID = rownames(DEXSeq_res_exon),DEXSeq_res_exon)
DEXSeq_res_exon$geneID <- strsplit2(as.character(DEXSeq_res_exon$groupID),":")[,1] 


DEXSeq_res_gene <- perGeneQValue(DEXSeq_res_exon, p = "pvalue")

DEXSeq_res_gene <- as.data.frame(DEXSeq_res_gene)
DEXSeq_res_gene <- data.frame(rownames(DEXSeq_res_gene),DEXSeq_res_gene)
colnames(DEXSeq_res_gene) <- c("Gene","padjust")


save(dxr , file=paste0(out.dir,"dexseq_1.10.8_bitseq.RData"))


write.table(DEXSeq_res_exon, file=paste0(out.dir,"dexseq_1.10.8_exon_bitseq_results.txt") , row.names=FALSE, sep="\t", quote = FALSE)
write.table(DEXSeq_res_gene, file=paste0( out.dir,"dexseq_1.10.8_gene_bitseq_results.txt") ,row.names=FALSE,sep="\t", quote = FALSE)

dim(DEXSeq_res_gene)
sum(DEXSeq_res_gene$padjust < 0.05)


dim(DEXSeq_res_exon[!is.na(DEXSeq_res_exon$padj),])





















