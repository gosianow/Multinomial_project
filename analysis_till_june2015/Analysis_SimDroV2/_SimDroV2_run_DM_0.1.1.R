##############################################################################

# BioC 3.0
# Created 28 Nov 2014:

# Run DM_0.1.1


##############################################################################


setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2")

# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata



out.dir <- "DM_0.1.1/"
dir.create(out.dir, showWarnings=F, recursive=T)


library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)

#### function to simulate data from DM
source("/home/gosia/R/R_Multinomial_project/Analysis_SimDM/simulate_from_DM.R")





##############################################################################################################
# run DM on FC 
##############################################################################################################

out.dir <- "DM_0.1.1/fc/"
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




### fc_DM_gridCommon ############################################################
name1 <- "fc_DM_gridCommon"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)

mcCores <- 20



## run DM pipeline
# dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

dgeDM <- dge
dgeDM$commonDispersion <- 5792.15274446394

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)


dgeDM <- dmEstimateTagwiseDisp(dgeDM, group = NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 15, gridRange = c(-7, 7), trend = c("none", "commonDispersion", "trendedDispersion")[2], priorDf = 4, span = 0.3, mcCores = mcCores, verbose = FALSE)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)


dgeDM <- dmFit(dgeDM, group=NULL, dispersion= "tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)


dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))


plot(log10(dgeDM$meanExpr), log10(dgeDM$tagwiseDispersion))




### fc_DM_gridTrend ############################################################
name1 <- "fc_DM_gridTrend"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)

mcCores <- 20



## run DM pipeline
# dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

dgeDM <- dge
dgeDM$commonDispersion <- 5792.15274446394

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)


dgeDM <- dmEstimateTagwiseDisp(dgeDM, group = NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 15, gridRange = c(-7, 7), trend = c("none", "commonDispersion", "trendedDispersion")[3], priorDf = 4, span = 0.3, mcCores = mcCores, verbose = FALSE)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)


dgeDM <- dmFit(dgeDM, group=NULL, dispersion= "tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)


dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))


plot(log10(dgeDM$meanExpr), log10(dgeDM$tagwiseDispersion))







