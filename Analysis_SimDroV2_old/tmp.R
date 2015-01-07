# source("/home/gosia/R/R_Multinomial_project/Analysis_SimDroV2/tmp.R")


#######################################################
# 
# Created 04 Sep 2014 

# comprare the DM version 4 commonDispersion and adjustement performance with other methods and DM dirmult

# Update 04 Sep 2014:



#######################################################
# BioC 2.14

setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2")

# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata



##########################################################################
# run new DM (common dispersion, and adjustement)  on FC
##########################################################################

library(edgeR)
library(parallel)
library(dirmult)
library(limma)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version4/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v4.R"))



##############################################################################################################
# run new DM (common dispersion, and adjustement)  on FC
##############################################################################################################

out.dir <- "DM_v4/fc/"
dir.create(out.dir, showWarnings=F, recursive=T)


#### load FC

fc <- read.table("featureCounts/featureCounts.txt")

counts <- fc[,2:7]
colnames(counts) <- metadata$SampleName1
gene.id <- strsplit2(fc[,1], ":")[,1]
ete.id <- fc[,1]

dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))

name1 <- "fc"
dge <- dgeOrg
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)





###############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM4adj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]


## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-03, mcCores=20, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

dgeDM <- dmTest(dgeDM, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###############################################################
name1 <- "fc_g0_s3_keep0s_subsetInf_DM4adj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 3
dge <- dge[keep,]


## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-03, mcCores=20, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

dgeDM <- dmTest(dgeDM, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))





##############################################################################################################
# run new DM (common dispersion, and adjustement)  on htseq counts
##############################################################################################################

out.dir <- "DM_v4/htseq/"
dir.create(out.dir, showWarnings=F, recursive=T)


#### load htseq counts

library(DEXSeq)

ecs <- DEXSeqDataSetFromHTSeq(countfiles = paste0("htseqCounts/htseq_samp_", metadata$SampleName1, ".txt"), sampleData= data.frame(row.names = metadata$SampleName, condition = metadata$condition), design = ~sample + exon + condition:exon)

counts <- featureCounts(ecs)

colnames(counts) <- metadata$SampleName1
gene.id <- strsplit2(rownames(counts), ":")[,1]
ete.id <- rownames(counts)

dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))

name1 <- "htseq"
dge <- dgeOrg
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)







################################################################
name1 <- "htseq_g0_s3_keep0s_subsetInf_DM4adj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 3
dge <- dge[keep,]


## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-03, mcCores=20, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

dgeDM <- dmTest(dgeDM, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




################################################################
name1 <- "htseq_g0_s4_keep0s_subsetInf_DM4adj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]


## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-03, mcCores=20, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

dgeDM <- dmTest(dgeDM, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))









































