# BioC 3.1

# Created 2 June 2015
# Analyse data from DM_0.1.5_filtering with DM_0.1.5 (New dgeSQTL object and use of BiocParallel)


##############################################################################################################

setwd("/home/Shared/data/seq/geuvadis/")

library(DM)
library(BiocParallel)

library(limma)
source("/home/gosia/R/multinomial_project/package_devel/0_my_printHead.R")

library(edgeR)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)


Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)




##########################################################################################
# Run the DMsQTL pipeline
##########################################################################################


######### run on DM_0_1_5 data 
out.dir <- "dm_0_1_5_analysis/Results_Data_DM_0_1_5_TagwiseDisp_gridNone/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

BPPARAM <- MulticoreParam(workers = 10)


load(paste0("DM_0_1_5_Data/dgeSQTL.RData"))




# ### subset 100 genes
# dgeSQTL_org <- dgeSQTL
# geneList <- names(dgeSQTL_org$counts)[1:1000]
# dgeSQTL$counts <- dgeSQTL_org$counts[geneList]
# dgeSQTL$genotypes <- dgeSQTL_org$genotypes[geneList]
# 
# dgeSQTL$commonDispersion <- 4

# dgeSQTL <- dgeSQTL_org




######### commonDispersion
cat(paste0("Estimate common dispersion \n"))

dgeSQTL <- dmSQTLEstimateCommonDisp(dgeSQTL, adjustDisp = TRUE, intervalDisp = c(0, 1e+3), tolDisp = 1e-01, modeProp = "constrOptim2G", tolProp = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)

save(dgeSQTL, file = paste0(out.dir, "dgeSQTL.RData"))

commonDispersion <- dgeSQTL$commonDispersion
save(commonDispersion, file = paste0(out.dir, "commonDispersion.RData"))


######### tagwiseDispersion
cat(paste0("Estimate tagwise dispersion \n"))

load(paste0(out.dir, "commonDispersion.RData"))
dgeSQTL$commonDispersion <- commonDispersion

dgeSQTL <- dmSQTLEstimateTagwiseDisp(dgeSQTL, adjustDisp = TRUE, modeDisp = c("optimize", "optim", "constrOptim", "grid")[4], intervalDisp = c(0, 1e+5), tolDisp = 1e-08,  initDisp = 10, initWeirMoMDisp = TRUE, gridLengthDisp = 15, gridRangeDisp = c(-7, 7), trendDisp = c("none", "commonDispersion", "trendedDispersion")[1], priorDfDisp = 10, spanDisp = 0.3, modeProp = c( "constrOptim2", "constrOptim2G", "FisherScoring")[2], tolProp = 1e-12, verbose = FALSE, plot = FALSE, BPPARAM = BPPARAM)


save(dgeSQTL, file = paste0(out.dir, "dgeSQTL.RData"))


######### DM fitting 
cat(paste0("Fitting DM \n"))

dgeSQTL <- dmSQTLFit(dgeSQTL, model = c("full", "null")[1], dispersion = c("commonDispersion", "tagwiseDispersion")[2], modeProp=c("constrOptim2", "constrOptim2G", "FisherScoring")[2], tolProp = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)

save(dgeSQTL, file = paste0(out.dir, "dgeSQTL.RData"))


table(unlist(lapply(dgeSQTL$fitFull, function(g){sum(sapply(g, is.null))})))



######### LR testing
cat(paste0("Testing \n"))

dgeSQTL <- dmSQTLTest(dgeSQTL, dispersion = c("commonDispersion", "tagwiseDispersion")[2], modeProp="constrOptim2G", tolProp = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)


save(dgeSQTL, file = paste0(out.dir, "dgeSQTL.RData"))

write.table(dgeSQTL$table, paste0(out.dir, "results.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




### plot histogram of p-values
pdf(paste0(out.dir, "Histogram_pValues.pdf"))
hist(dgeSQTL$table$pValue, breaks = 100, main = "", col = "hotpink", xlab = "p-values")
dev.off()


### plot dispersion versus mean
commonDispersion <- log10(dgeSQTL$commonDispersion)
tagwiseDispersion = log10(unlist(dgeSQTL$tagwiseDispersion))

meanExpr = log10(dgeSQTL$meanExpr[unlist(lapply(names(dgeSQTL$tagwiseDispersion), function(g){rep(g, length(dgeSQTL$tagwiseDispersion[[g]]))}))] +1 )

nrTrans <- dgeSQTL$table[, c("geneID", "snpID", "nrGroups", "df")]
rownames(nrTrans) <- paste0(nrTrans$geneID, ".", nrTrans$snpID)
nrTrans$nrTrans <- nrTrans$df/(nrTrans$nrGroups - 1) + 1

df <- data.frame(tagwiseDispersion, meanExpr, nrTrans = nrTrans[names(tagwiseDispersion), "nrTrans"])


ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = nrTrans)) +
theme_bw() +
xlab("Log10 meanExpr") +
ylab("Log10 tagwiseDispersion") +
geom_point(size = 1, alpha = 0.5) +
geom_hline(aes(yintercept =  commonDispersion), colour = "deeppink", linetype="dashed", size = 1)+
theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.title = element_text(size=14, face="bold"), legend.text = element_text(size = 14)) +
scale_colour_gradient(limits=c(1, 10), breaks = c(2, 4, 6, 8, 10), low = "blue", high="red", name = "nrTrans", na.value = "red")


pdf(paste0(out.dir, "DispersionVersusMean.pdf"), 7, 5)
print(ggp2)
dev.off()



























