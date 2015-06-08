##############################################################################################################
# created 14 Nov 2014
# BioC 3.0
# test DM package for standard differential splicing pipeline
# Test modeDip = grid, priorDf, gridLength

##############################################################################################################


#######################################################
# create dge object from FC counts from Drosophila
#######################################################

setwd("/home/gosia/Multinomial_project/DM_package_devel/")


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


library(limma)

fc <- read.table("featureCounts/featureCounts.txt")

counts <- fc[,2:7]
colnames(counts) <- metadata$SampleName1
gene_id <- strsplit2(fc[,1], ":")[,1]
ete_id <- fc[,1]

library(edgeR)
dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id))

colnames(dgeOrg$counts) <- dgeOrg$samples$group
rownames(dgeOrg$counts) <- dgeOrg$genes$ete_id


######### filtering
name1 <- "SimDroV2_fc"
dge <- dgeOrg
keep <- rowSums(cpm(dge)>0) >= 4
dge <- dge[keep,]



# dge <- dge[1:2000, ]

dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id))

dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR



##############################################################################################################
# run DM pipeline for aternative splicing
##############################################################################################################

out.dir <- "DM_0.1.1_ModeratedDispersion/"
dir.create(out.dir, showWarnings=F, recursive=T)

# save.image(paste0(out.dir, "testWorkspace.RData"))



#######################################################
# test DM_0.1.1
#######################################################

### If all above is fine, load this workspace 
setwd("/home/gosia/Multinomial_project/DM_package_devel/")
out.dir <- "DM_0.1.1_ModeratedDispersion/"
load(paste0(out.dir, "testWorkspace.RData"))


library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")


### Source all R files in DM package
Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in 1:length(Rfiles)) source(Rfiles[i])



######################## Test modeDip = grid 


dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, subset=Inf, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=20, verbose=FALSE)
### Common dispersion for full data set = 3741.896


gridLength = 15
gridRange = c(-7, 7)



splinePts <- seq(from = gridRange[1], to = gridRange[2], length = gridLength)
splinePts
splineDisp <- dgeDM$commonDispersion * 2^splinePts
splineDisp



dgeDMnone <- dgeDM <- dmEstimateTagwiseDisp(dgeDM, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp=c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = FALSE, gridLength=gridLength, gridRange=gridRange, trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf=priorDf, span=span, mcCores=20, verbose=FALSE)


pdf(paste0(out.dir, "DispVsMean_grid_none.pdf"))

plot(log10(dgeDM$meanExpr+1), log10(dgeDM$tagwiseDispersion), pch=".", ylim=c(1,6))
abline(h=log10(dgeDM$commonDisp), col="deeppink", lwd = 2, lty=2)

plot(log10(dgeDM$meanExpr+1), dgeDM$tagwiseDispersion, pch=".")
abline(h=dgeDM$commonDisp, col="deeppink", lwd = 2, lty=2)

dev.off()



priorDf = 4 ### priorN <- priorDf/(nlibs - ngroups)
span = 0.3



dgeDM <- dmEstimateTagwiseDisp(dgeDM, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp=c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = FALSE, gridLength=gridLength, gridRange=gridRange, trend = c("none", "commonDispersion", "trendedDispersion")[2], priorDf = 10, span = 0.3, mcCores = 20, verbose=FALSE, plot=FALSE)


pdf(paste0(out.dir, "DispVsMean_grid_commonDispersion.pdf"))

plot(log10(dgeDM$meanExpr+1), log10(dgeDM$tagwiseDispersion), pch=".", ylim=c(1,6))
abline(h=log10(dgeDM$commonDisp), col="deeppink", lwd = 2, lty=2)

plot(log10(dgeDM$meanExpr+1), dgeDM$tagwiseDispersion, pch=".")
abline(h=dgeDM$commonDisp, col="deeppink", lwd = 2, lty=2)

dev.off()





#### if plot == TRUE

plotGenes <- names(dgeDM$meanExpr)[which(dgeDM$meanExpr > 10 & dgeDM$meanExpr < 100 & dgeDMnone$tagwiseDispersion > 1e5)]
subset <- 10


pdf(paste0(out.dir, "/logLiks.pdf"))

for(i in plotGenes[seq(subset)]){
  
  
  
  ylim <- c(min(c(dgeDM$plotLoglik[plotGenes[i],], dgeDM$plotLoglik0[plotGenes[i],], dgeDM$plotModeration[plotGenes[i],])) ,max(c(dgeDM$plotLoglik[plotGenes[i],], dgeDM$plotLoglik0[plotGenes[i],], dgeDM$plotModeration[plotGenes[i],])))
  
#   plot(log10(dgeDM$plotSplineDisp), log10(dgeDM$plotLoglik[plotGenes[i],] - ylim[1] + 1), type="l", lwd=3, col = "black", ylim = log10(ylim - ylim[1] +1) )
#   lines( log10(dgeDM$plotSplineDisp), log10(dgeDM$plotLoglik0[plotGenes[i],]- ylim[1] + 1), type="l", lwd=3, col = "black", lty=2)
#   lines( log10(dgeDM$plotSplineDisp), log10(dgeDM$plotModeration[plotGenes[i],]- ylim[1] + 1), type="l", lwd=3, col = "grey", lty=2)
#   abline(v = log10(dgeDM$plotOutDisp[plotGenes[i]]), col = "pink")
#   abline(v = log10(dgeDM$commonDispersion), col="grey", lty=2)
  
r <- 1:15

  plot(dgeDM$plotSplineDisp[r], dgeDM$plotLoglik[plotGenes[i],r], type="l", lwd=3, col = "black", ylim=ylim )
  lines( dgeDM$plotSplineDisp[r], dgeDM$plotLoglik0[plotGenes[i],r], type="l", lwd=3, col = "black", lty=2)
  lines( dgeDM$plotSplineDisp[r], dgeDM$plotModeration[plotGenes[i],r], type="l", lwd=3, col = "grey", lty=2)
  abline(v = dgeDM$plotOutDisp[plotGenes[i]], col = "pink")
  abline(v = dgeDM$commonDispersion, col="grey", lty=2)
abline( v = dgeDNnone$tagwiseDispersion[] )
  
plot(dgeDM$plotSplineDisp[r], dgeDM$plotLoglik[plotGenes[i],r], type="l", lwd=3, col = "black", )
abline(v = dgeDM$plotOutDisp[plotGenes[i]], col = "pink")
abline(v = dgeDM$commonDispersion, col="grey", lty=2)

plot( dgeDM$plotSplineDisp[r], dgeDM$plotLoglik0[plotGenes[i],r], type="l", lwd=3, col = "black", lty=2)
abline(v = dgeDM$plotOutDisp[plotGenes[i]], col = "pink")
abline(v = dgeDM$commonDispersion, col="grey", lty=2)

plot( dgeDM$plotSplineDisp[r], dgeDM$plotModeration[plotGenes[i],r], type="l", lwd=3, col = "grey", lty=2)
abline(v = dgeDM$plotOutDisp[plotGenes[i]], col = "pink")
abline(v = dgeDM$commonDispersion, col="grey", lty=2)
}

dev.off()







dgeDM <- dmEstimateTagwiseDisp(dgeDM, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp=c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = FALSE, gridLength=gridLength, gridRange=gridRange, trend = c("none", "commonDispersion", "trendedDispersion")[3], priorDf=priorDf, span=span, mcCores=20, verbose=FALSE)


pdf(paste0(out.dir, "DispVsMean_grid_trendedDispersion.pdf"))

plot(log10(dgeDM$meanExpr+1), log10(dgeDM$tagwiseDispersion), pch=".", ylim=c(1,6))
abline(h=log10(dgeDM$commonDisp), col="deeppink", lwd = 2, lty=2)

plot(log10(dgeDM$meanExpr+1), dgeDM$tagwiseDispersion, pch=".")
abline(h=dgeDM$commonDisp, col="deeppink", lwd = 2, lty=2)

dev.off()




################ Test other


dgeDM <- dmEstimateTagwiseDisp(dgeDM, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp=c("optimize", "optim", "constrOptim", "grid")[3], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength=11, gridRange=c(-3, 3), trend = c("none", "commonDispersion", "trendedDispersion")[2], priorDf=10, span=0.3, mcCores=20, verbose=FALSE)


pdf(paste0(out.dir, "DispVsMean_constrOptim_Weir.pdf"))

plot(log10(dgeDM$meanExpr+1), log10(dgeDM$tagwiseDispersion), pch=".")
abline(h=log10(dgeDM$commonDisp), col="deeppink", lwd = 2, lty=2)

plot(log10(dgeDM$meanExpr+1), dgeDM$tagwiseDispersion, pch=".")
abline(h=dgeDM$commonDisp, col="deeppink", lwd = 2, lty=2)


dev.off()




dgeDM <- dmEstimateTagwiseDisp(dgeDM, group=NULL, adjust = FALSE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp=c("optimize", "optim", "constrOptim", "grid")[1], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength=11, gridRange=c(-3, 3), trend = c("none", "commonDispersion", "trendedDispersion")[2], priorDf=10, span=0.3, mcCores=20, verbose=FALSE)


pdf(paste0(out.dir, "DispVsMean_optimize5.pdf"))

plot(log10(dgeDM$meanExpr+1), log10(dgeDM$tagwiseDispersion), pch=".")
abline(h=log10(dgeDM$commonDisp), col="deeppink", lwd = 2, lty=2)

plot(log10(dgeDM$meanExpr+1), dgeDM$tagwiseDispersion, pch=".")
abline(h=dgeDM$commonDisp, col="deeppink", lwd = 2, lty=2)


dev.off()



#######################################################
# test DM_1.0
#######################################################

### If all above is fine, load this workspace 
setwd("/home/gosia/Multinomial_project/DM_package_devel/")
out.dir <- "DM_output/"
load(paste0(out.dir, "testWorkspace.RData"))


library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

dgeDM <- dmFit(dge, group=NULL, dispersion=3000, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)


### Source all R files in DM package
Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in 1:length(Rfiles)) source(Rfiles[i])



debug(dmLogLikkm1)

dgeDM <- dmFit(dge, group=NULL, dispersion=3000, mode="FisherScoring", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 1) 
# Error in names(pi) <- rownames(y) :
#   'names' attribute [8] must be the same length as the vector [4]


undebug(dmLogLikkm1)


#######################################################


































































