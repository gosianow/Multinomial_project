# BioC 3.0
# Crated 17 June 2015
# Midified

# Run DEXSeq on kallisto and bitseq counts

#####################################################################

setwd("/home/gosia/Multinomial_project/Simulations_drosophila_Sim5_noDE_noNull")

# create metadata file
metadata <- data.frame(sampleName1 = paste0("sample_",1:6), sampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


out.dir <- "DEXSeq_1.10.8/"
dir.create(out.dir, showWarnings=F, recursive=T)


library("DEXSeq")
library("edgeR")
library("BiocParallel")


nCores <- 10
BPPARAM = MulticoreParam(workers = nCores)
countfiles <-  paste0("2_counts/kallisto/kallisto",1:6, ".txt")


dxr   <- DEXSeqDataSetFromHTSeq(countfiles = countfiles, sampleData = metadata , design = ~ sample + exon + condition:exon)


dxr  <- estimateSizeFactors(dxr)
dxr  <- estimateDispersions(dxr, BPPARAM=BPPARAM )
dxr   <- testForDEU(dxr, BPPARAM=BPPARAM)


DEXSeq_res_exon <-  DEXSeqResults(dxr)

DEXSeq_res_gene <- perGeneQValue(DEXSeq_res_exon, p = "pvalue")


DEXSeq_res_exon <- as.data.frame(DEXSeq_res_exon[, c("groupID", "featureID", "pvalue", "padj")])
colnames(DEXSeq_res_exon) <- c("gene_id", "feature_id", "PValue", "adjPValue")
DEXSeq_res_exon$feature_id <- substring(DEXSeq_res_exon$feature_id, 2)

DEXSeq_res_exon$group_id <- paste0(DEXSeq_res_exon$gene_id, ":", DEXSeq_res_exon$feature_id)


DEXSeq_res_gene <- as.data.frame(DEXSeq_res_gene)
DEXSeq_res_gene <- data.frame(rownames(DEXSeq_res_gene),DEXSeq_res_gene)
colnames(DEXSeq_res_gene) <- c("gene_id", "adjPValue")


save(dxr , file=paste0(out.dir,"dexseq_kallisto.RData"))


write.table(DEXSeq_res_exon, file=paste0(out.dir,"dexseq_kallisto_exon.txt") , row.names=FALSE, sep="\t", quote = FALSE)
write.table(DEXSeq_res_gene, file=paste0( out.dir,"dexseq_kallisto.txt") ,row.names=FALSE,sep="\t", quote = FALSE)






















