
# BioC 14

# Created 07 Oct 2014
# Modyfied 28 Oct 2014


##############################################################################################################

# test for sQTL pipeline

##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

load("DM_sQTL_test/dgeSQTL.RData")


library(edgeR)
library(parallel)

library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v5.R"))

out.dir <- "DM_sQTL_test/"
dir.create(out.dir, showWarnings=F, recursive=T)




dgeSQTL <- dmSQTLFit(dgeSQTL, model="full", dispersion=3000, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 30)


dgeSQTL <- dmSQTLTest(dgeSQTL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 30)




dgeSQTL <- dmSQTLEstimateCommonDisp(dgeSQTL, adjust = FALSE, subset=Inf, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-03, mcCores=30, verbose=FALSE)




##############################################################################################################

# test for standard differential splicing pipeline

##############################################################################################################


#######################################################
# create dge object from FC counts from Drosophila
#######################################################

setwd("/home/gosia/Multinomial_project/DM_package/")


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


library(limma)

fc <- read.table("PLOTS4/featureCounts/featureCounts.txt")

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

# dge <- dge[1:5000, ]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id))

dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


# save(dge, file="/home/Shared/tmp/gosia/Example_DM/dge.RData")



#######################################################
# run DM pipeline for aternative splicing
#######################################################
## constrOptim does not work!!!

library(edgeR)
library(parallel)

library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")


source(paste0(Rdir, "dmFunctions_v5_back1.R"))

source(paste0(Rdir, "dmFunctions_v5_back2.R"))

source(paste0(Rdir, "dmFunctions_v5.R"))


out.dir <- "PLOTS5/"
dir.create(out.dir, showWarnings=F, recursive=T)


############### run 

length(unique(as.character(dge$genes$gene.id)))
length(unique(dge$genes$gene.id))






system.time(dgeDM <- dmFit(dge, group=NULL, dispersion=3000, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20))


system.time(dgeDM <- dmFit(dge, group=NULL, dispersion=3000, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20))




system.time(dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2", mcCores=30, verbose=FALSE))

system.time(dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2G", mcCores=30, verbose=FALSE))



system.time(dgeDMadj <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2", mcCores=30, verbose=FALSE))

system.time(dgeDMadj <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", mcCores=30, verbose=FALSE))







dgeDM <- dgeDMadj

dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)


dgeDM <- dmTest(dgeDM, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)


head(dgeDM$table)


save(dgeDM, file = paste0(out.dir, "/", name1, "_dgeDM",".Rdata"))
write.table(dgeDM$table, paste0(out.dir, "/", name1, "_table",".xls"), sep = "\t", quote = F, row.names = F)


























