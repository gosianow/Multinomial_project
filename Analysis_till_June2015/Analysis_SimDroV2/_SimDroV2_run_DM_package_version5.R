#######################################################
# 
# Created 22 Oct 2014 

# run DM version 5 (SQTL) on new Simulations_drosophila_V2

# Update 29 Oct 2014:
# in dmFit dispersion = "commonDispersion"


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

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v5.R"))





##############################################################################################################
# run DM version 5 (common dispersion + adjustement) on FC
##############################################################################################################

out.dir <- "DM_v5/fc/"
dir.create(out.dir, showWarnings=F, recursive=T)


#### load FC

fc <- read.table("featureCounts/all/featureCounts.txt")

counts <- fc[,2:7]
colnames(counts) <- metadata$SampleName1
gene_id <- strsplit2(fc[,1], ":")[,1]
ete_id <- fc[,1]

dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id))

name1 <- "fc"
dge <- dgeOrg
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)




###fc_g0_s4_keep0s_subsetInf_DM5############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)

mcCores <- 20

## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###fc_g0_s4_keep0s_subsetInf_DM5adj############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5adj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


nlevels(dge$genes$gene_id)

mcCores <- 20

## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))



###fc_g0_s4_keep0s_subsetInf_DM5adjM2############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5adjM2"

dgeDM$commonDispersion <- dgeDM$commonDispersion/2

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###fc_g0_s4_keep0s_subsetInf_DM5adjM4############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5adjM4"

dgeDM$commonDispersion <- dgeDM$commonDispersion/2

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###fc_g0_s4_keep0s_subsetInf_DM5adjM8############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5adjM8"

dgeDM$commonDispersion <- dgeDM$commonDispersion/2

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))





###fc_g0_s4_keep0s_subsetInf_DM5TG############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TG"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)

mcCores <- 20

## run DM pipeline
dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=TRUE, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion="tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###fc_g0_s4_keep0s_subsetInf_DM5TGadj############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TGadj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


mcCores <- 20

## run DM pipeline
dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion="tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))





###fc_g0_s4_keep0s_subsetInf_DM5TGo############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TGo"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)

mcCores <- 30

## run DM pipeline
dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE, modeDisp="optim", initDisp = 100)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=TRUE, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion="tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###fc_g0_s4_keep0s_subsetInf_DM5TGoadj############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TGoadj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


mcCores <- 30

## run DM pipeline
dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE, modeDisp="optim", initDisp = 100)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion="tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))



###fc_g0_s4_keep0s_subsetInf_DM5TGco############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TGco"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)

mcCores <- 20

## run DM pipeline
dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE, modeDisp="constrOptim", initDisp = 100)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=TRUE, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion="tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###fc_g0_s4_keep0s_subsetInf_DM5TGcoadj############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TGcoadj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


mcCores <- 20

## run DM pipeline
dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE, modeDisp="constrOptim", initDisp = 100)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion="tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))



##############################################################################################################
# run DM version 5 (common dispersion + adjustement)  on htseq counts
##############################################################################################################

out.dir <- "DM_v5/htseq/"
dir.create(out.dir, showWarnings=F, recursive=T)


#### load htseq counts

library(DEXSeq)

ecs <- DEXSeqDataSetFromHTSeq(countfiles = paste0("htseqCounts/htseq_samp_", metadata$SampleName1, ".txt"), sampleData= data.frame(row.names = metadata$SampleName, condition = metadata$condition), design = ~sample + exon + condition:exon)

counts <- featureCounts(ecs)

colnames(counts) <- metadata$SampleName1
gene_id <- strsplit2(rownames(counts), ":")[,1]
ete_id <- rownames(counts)

dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id))

name1 <- "htseq"
dge <- dgeOrg
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)





###htseq_g0_s4_keep0s_subsetInf_DM5#############################################################
name1 <- "htseq_g0_s4_keep0s_subsetInf_DM5"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

mcCores <- 20


## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###htseq_g0_s4_keep0s_subsetInf_DM5adj#############################################################
name1 <- "htseq_g0_s4_keep0s_subsetInf_DM5adj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

mcCores <- 20


## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))






##############################################################################################################
# run DM version 5 (common dispersion + adjustement) on FC by Gordon
##############################################################################################################

out.dir <- "DM_v5/fcG/"
dir.create(out.dir, showWarnings=F, recursive=T)


#### load FC

fc <- read.table("featureCounts_Gordon/all/featureCounts.txt")

counts <- fc[,2:7]
colnames(counts) <- metadata$SampleName1
gene_id <- strsplit2(fc[,1], ":")[,1]
ete_id <- fc[,1]

dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id))

name1 <- "fc"
dge <- dgeOrg
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)




###fcG_g0_s4_keep0s_subsetInf_DM5############################################################
name1 <- "fcG_g0_s4_keep0s_subsetInf_DM5"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)

mcCores <- 10

## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###fcG_g0_s4_keep0s_subsetInf_DM5adj############################################################
name1 <- "fcG_g0_s4_keep0s_subsetInf_DM5adj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


nlevels(dge$genes$gene_id)

mcCores <- 10

## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))








##############################################################################################################
# run DM version 5 (common dispersion + adjustement) on RSEM RSEM_exp_count
##############################################################################################################

out.dir <- "DM_v5/rsem/"
dir.create(out.dir, showWarnings=F, recursive=T)


#### load RSEM

rsem <- read.table("RSEM/RSEM_exp_count.xls", header = TRUE)

counts <- round(rsem[,4:9])
colnames(counts) <- metadata$SampleName1
gene_id <- rsem[,"gene_id"]
ete_id <- rsem[,"transcript_id"]

dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id))

name1 <- "rsem"
dge <- dgeOrg
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)




###rsem_g0_s4_keep0s_subsetInf_DM5############################################################
name1 <- "rsem_g0_s4_keep0s_subsetInf_DM5"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)

mcCores <- 20

## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###rsem_g0_s4_keep0s_subsetInf_DM5adj############################################################
name1 <- "rsem_g0_s4_keep0s_subsetInf_DM5adj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


nlevels(dge$genes$gene_id)

mcCores <- 20

## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




###rsem_g0_s4_keep0s_subsetInf_DM5adjM2############################################################
name1 <- "rsem_g0_s4_keep0s_subsetInf_DM5adjM2"


dgeDM$commonDispersion <- dgeDM$commonDispersion/2

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))






























