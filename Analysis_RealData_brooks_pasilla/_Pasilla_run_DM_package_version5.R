#######################################################
# source("/home/gosia/R/R_Multinomial_project/Analysis_RealData_brooks_pasilla/_Pasilla_run_DM_package_version5.R")

# Created 04 Nov 2014 

# comprare the DM version 4 commonDispersion and adjustement performance with other methods and DM dirmult

# Update 05 Nov 2014:

# + run DM version 5 (SQTL) with all different versions of tagwise dispersion estimation
# move DM results to Shared directory

#######################################################
# BioC 2.14


setwd("/home/Shared/data/seq/brooks_pasilla")

# create metadata file
sri.org <- read.table("/home/Shared/data/seq/brooks_pasilla/download/SraRunInfo.csv", stringsAsFactors=F, sep=",", header=T)
keep <- grep(paste("GSM4611", 76:82, sep="", collapse="|"), sri.org$SampleName)
sri <- sri.org[keep,]

# manual trimming of reads to have equal lenght 
which.trim <- sri[sri$SampleName=="GSM461179", c("Run")]

# cmd <- paste0("java -jar /usr/local/software/Trimmomatic-0.30/trimmomatic-0.30.jar  SE -threads 20 download/", which.trim , ".fastq.gz download/", which.trim , ".trimed.fastq.gz CROP:40 \n")
# cat(cmd)
# for(i in 1:length(cmd))
#   system(cmd[i])

sri[sri$Run %in% which.trim, "avgLength"] <- 40
sri[sri$Run %in% which.trim, "Run"] <- paste0(sri[sri$Run %in% which.trim, "Run"], ".trimed")

sri$LibraryName = gsub("S2_DRSC_","",sri$LibraryName) # trim label
metadata = unique(sri[,c("LibraryName","LibraryLayout", "SampleName", "avgLength" )])

for(i in seq_len(nrow(metadata))) {
  rw = (sri$LibraryName==metadata$LibraryName[i])
  if(metadata$LibraryLayout[i]=="PAIRED") {
    metadata$fastq1[i] = paste0("download/",sri$Run[rw],"_1.fastq.gz",collapse=",")
    metadata$fastq2[i] = paste0("download/",sri$Run[rw],"_2.fastq.gz",collapse=",")
    metadata$ReadLength[i] <- metadata$avgLength[i] / 2
    metadata$SRR[i] <- paste0(sri$Run[rw], collapse=",")
  } else {
    metadata$fastq1[i] = paste0("download/",sri$Run[rw],".fastq.gz",collapse=",")
    metadata$fastq2[i] = ""
    metadata$ReadLength[i] <- metadata$avgLength[i]
    metadata$SRR[i] <- paste0(sri$Run[rw], collapse=",")
  }
}

metadata$condition = "CTL"
metadata$condition[grep("RNAi",metadata$LibraryName)] = "KD"
metadata$shortname = paste( seq_len(nrow(metadata)), substr(metadata$condition,1,2),  substr(metadata$LibraryLayout,1,2), metadata$ReadLength, sep=".")

metadata$color[metadata$condition == "CTL"] <- "chartreuse3"
metadata$color[metadata$condition == "KD"] <- "darkorchid3"

metadata.org <- metadata




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

fc <- read.table("featureCounts/fc.txt", header = T)

colnames(fc) <- gsub("_s.bam","",gsub("sam.", "",colnames(fc)))

counts <- fc[, metadata$SampleName]
gene_id <- strsplit2(fc[,1], ":")[,1]
ete_id <- fc[,1]


dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id))

name1 <- "fc"
dge <- dgeOrg
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)




# ###fc_g0_s4_keep0s_subsetInf_DM5adj############################################################
# name1 <- "fc_g0_s4_keep0s_subsetInf_DM5adj"
# dge <- dgeOrg
# keep <- rowSums(cpm(dge) > 0) >= 5
# dge <- dge[keep,]
# dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
# dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
# dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR
# 
# 
# nlevels(dge$genes$gene_id)
# 
# mcCores <- 20
# 
# ## run DM pipeline
# dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)
# 
# write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)
# 
# dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
# 
# dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
# 
# write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
# save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))
# 
# 
# 
# 
# 
# ###fc_g0_s4_keep0s_subsetInf_DM5TGadj############################################################
# name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TGadj"
# dge <- dgeOrg
# keep <- rowSums(cpm(dge) > 0) >= 4
# dge <- dge[keep,]
# dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
# dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
# dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR
# 
# 
# mcCores <- 20
# 
# ## run DM pipeline
# dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)
# 
# write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)
# 
# dgeDM <- dmFit(dgeDM, group=NULL, dispersion="tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
# 
# dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
# 
# write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
# save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))





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




###fc_g0_s4_keep0s_subsetInf_DM5TGo10adj############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TGo10adj"

mcCores <- 30

## run DM pipeline
dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE, modeDisp="optim", initDisp = 10)

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




###fc_g0_s4_keep0s_subsetInf_DM5TGco10adj############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TGco10adj"

mcCores <- 20

## run DM pipeline
dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE, modeDisp="constrOptim", initDisp = 10)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion="tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))



###fc_g0_s4_keep0s_subsetInf_DM5TGcoWadj############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM5TGcoWadj"

mcCores <- 20

## run DM pipeline
dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE, modeDisp="constrOptim", initDisp = 10, initWeirMoM = TRUE)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)
write.table(dgeDM$initDispersion, paste0(out.dir, "/",name1,"_initDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)


dgeDM <- dmFit(dgeDM, group=NULL, dispersion="tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))













































