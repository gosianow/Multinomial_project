######################################################
# BioC 2.14

# Created 04 Nov 2014 
# Analysis of Kim_adenocarcinoma data

# Update 19 Nov 2014:
# plots of TREND dispersion vs. mean 
#######################################################



setwd("/home/Shared/data/seq/Kim_adenocarcinoma/")


metadata <- read.table("metadata/Malgorzata Nowicka2014-11-04GSE37764.csv", stringsAsFactors=F, sep=",", header=T) 

metadataOrg <- metadata <- metadata[metadata$X == "RNA-seq",]



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



#### load FC

fc <- read.table("featureCounts/fc.txt", header = T)

colnames(fc) <- gsub("_s.bam","",gsub("X.home.Shared.data.seq.Kim_adenocarcinoma.bam_insilicodb.", "",colnames(fc)))






# Null_normal1
model <- "Null_normal1"

out.dir <- paste0("DM_v5/fc/diff_out_", model)
dir.create(out.dir, showWarnings=F, recursive=T)

metadata <-  metadataOrg[metadataOrg$Tissue.Type == "normal", ]
metadata$condition = c(rep("C1", 3), rep("C2", 3))


counts <- fc[, metadata[, "ids"]]
gene_id <- strsplit2(fc[,1], ":")[,1]
ete_id <- fc[,1]


dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id))

name1 <- "fc"
dge <- dgeOrg
# write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)




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




###fc_g0_s4_keep0s_subsetInf_DM5TGoWadj############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TGoWadj"

mcCores <- 20

## run DM pipeline
dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE, modeDisp="optim", initDisp = 10, initWeirMoM = TRUE)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)
write.table(dgeDM$initDispersion, paste0(out.dir, "/",name1,"_initDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion="tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))









########### Null_tumor1 ###########
model <- "Null_tumor1"

out.dir <- paste0("DM_v5/fc/diff_out_", model)
dir.create(out.dir, showWarnings=F, recursive=T)

metadata <-  metadataOrg[metadataOrg$Tissue.Type == "tumor", ]
metadata$condition = c(rep("C1", 3), rep("C2", 3))


counts <- fc[, metadata[, "ids"]]
gene_id <- strsplit2(fc[,1], ":")[,1]
ete_id <- fc[,1]


dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id))

name1 <- "fc"
dge <- dgeOrg
# write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)




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




###fc_g0_s4_keep0s_subsetInf_DM5TGoWadj############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TGoWadj"

mcCores <- 10

## run DM pipeline
dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE, modeDisp="optim", initDisp = 10, initWeirMoM = TRUE)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)
write.table(dgeDM$initDispersion, paste0(out.dir, "/",name1,"_initDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion="tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))







#######################################################
### calculate the mean gene expression
#######################################################

out.dir <- "PLOTS_DM_v5_TREND_dispVSmean/"
dir.create(out.dir, recursive = T, showWarnings = FALSE)


# model <- "diff_out_Null_normal1"
model <- "diff_out_Null_tumor1"

dir.create(paste0(out.dir,"/", model))

load(paste0("DM_v5/fc/", model ,"/fc_g0_s4_keep0s_subsetInf_DM5adj_dgeDM.RData"))


meanExpr <- sapply(dgeDM$counts, function(g){ mean(colSums(g)) } )

meanExpr <- data.frame(gene_id = names(meanExpr), meanExpr = meanExpr)

head(meanExpr)

table <- meanExpr

#######################################################
# plot dispersion vs mean
#######################################################



### load common dispersions
cDisp <- read.table(paste0("DM_v5/fc/",model,"/fc_g0_s4_keep0s_subsetInf_DM5adj_commonDispersion.txt"))


files <- list.files(path = paste0("DM_v5/fc/", model), pattern = "_results.xls" )
files <- files[grepl(pattern = "TG", files)]

TGmethods <- gsub(pattern = "_results.xls", replacement = "" , files)


for( i in 1:length(TGmethods)){
  # i = 1
  
  tDisp <- read.table(paste0("DM_v5/fc/", model ,"/",TGmethods[i],"_tagwiseDispersion.txt"))
  tName <- paste0(TGmethods[i],"_tagwiseDispersion")
  colnames(tDisp) <- c("gene_id", tName)
  
  
  table <- unique(merge(table, tDisp, by = "gene_id", all.x=TRUE))
  
  pdf(paste0(out.dir, "/", model, "/TREMD_mean_vs_gamma-",TGmethods[i],".pdf"))
  
  smoothScatter(log10(table$meanExpr), log10(table[,tName]), xlab="log10 mean gene expression", ylab="log10 gamma +", nrpoints = Inf, colramp=colorRampPalette(c("white", "grey40")), pch = 19, cex=0.6)
  abline(h = log10(cDisp), col = "red")

  dev.off()
  
}














































