# BioC 3.0

# Created 13 Jan 2014
# Modyfied 13 Jan 2014

# Test the DM sQTL pipeline


##############################################################################################################
# Test DM_0.1.1 (this is exactly DM version 5)
##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")


library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")


load("DMv5_sQTL_test/dgeSQTL.RData")
# load("DMv5_sQTL_analysis/DM_genotypes/dgeSQTL.RData")

### subset the dgeSQTL object 
keep.genes <- intersect(unique(dgeSQTL$genes$gene_id), unique(dgeSQTL$SNPs$gene_id))[1:20]
dgeSQTL$counts <- dgeSQTL$counts[keep.genes]
dgeSQTL$genes <- dgeSQTL$genes[dgeSQTL$genes$gene_id %in% names(dgeSQTL$counts), ]
dgeSQTL$SNPs <- dgeSQTL$SNPs[dgeSQTL$SNPs$gene_id %in% names(dgeSQTL$counts), ]
dgeSQTL$genotypes <- dgeSQTL$genotypes[dgeSQTL$SNPs$gene_id %in% names(dgeSQTL$counts), ]


out.dir <- "DM_0.1.1_sQTL_test/"
dir.create(out.dir, showWarnings=F, recursive=T)


dgeSQTL <- dmSQTLFit(dgeSQTL, model="full", dispersion=3, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 10)


dgeSQTL <- dmSQTLTest(dgeSQTL, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 10)


dgeSQTL <- dmSQTLEstimateCommonDisp(dgeSQTL, adjust = TRUE, subset=Inf, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=10, verbose=FALSE)


save.image(paste0(out.dir, "workspace.RData"))

dgeSQTL11 <- dgeSQTL

##############################################################################################################
# Test DM_0.1.2 - implement tagwise dispersion 
##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")


out.dir <- "DM_0.1.2_sQTL_test/"
dir.create(out.dir, showWarnings=F, recursive=T)


load("DMv5_sQTL_test/dgeSQTL.RData") ## Chr19 only
# load("DMv5_sQTL_analysis/DM_genotypes/dgeSQTL.RData") ## Full data set

### subset the dgeSQTL object 
keep.genes <- intersect(unique(dgeSQTL$genes$gene_id), unique(dgeSQTL$SNPs$gene_id))[1:20]
dgeSQTL$counts <- dgeSQTL$counts[keep.genes]
dgeSQTL$genes <- dgeSQTL$genes[dgeSQTL$genes$gene_id %in% names(dgeSQTL$counts), ]
dgeSQTL$genotypes <- dgeSQTL$genotypes[dgeSQTL$SNPs$gene_id %in% names(dgeSQTL$counts), ]
dgeSQTL$SNPs <- dgeSQTL$SNPs[dgeSQTL$SNPs$gene_id %in% names(dgeSQTL$counts), ]


### Source all R files in DM package
Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)


dgeSQTL$commonDispersion <- 1

# debug(dmOneGeneManyGroups)

dgeSQTL <- dmSQTLFit(dgeSQTL, model="full", dispersion=c("commonDispersion", "tagwiseDispersion")[1], mode=c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 15)

# undebug(dmOneGeneManyGroups)


dgeSQTL <- dmSQTLTest(dgeSQTL, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 15)


### check how many snps are called significant 
dim(dgeSQTL$table)
sum(dgeSQTL$table$FDR < 0.05, na.rm = TRUE)




dgeSQTL <- dmSQTLEstimateCommonDisp(dgeSQTL, adjust = TRUE, subset=Inf, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+4), tol = 1e-00, mcCores=10, verbose=FALSE)


dgeSQTL$commonDispersion


dgeSQTL <- dmSQTLEstimateTagwiseDisp(dgeSQTL, adjust = TRUE, mode = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 15, gridRange = c(-7, 7), trend = c("none", "commonDispersion", "trendedDispersion")[2], priorDf = 10, span = 0.3, mcCores = 20, verbose = FALSE, plot = FALSE)



dgeSQTL$tagwiseDispersion



































