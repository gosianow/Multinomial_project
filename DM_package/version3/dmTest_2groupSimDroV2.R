# BioC 14
# Created 31 Aug 2014

# Updated 03 Sep 2014: + plot profile likelihood for htseq counts and for fc counts


setwd("/home/gosia/Multinomial_project/DM_package/")


#######################################################


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


#######################################################
# create dge object from FC counts from Drosophila
#######################################################


library(limma)

fc <- read.table("PLOTS3/featureCounts/featureCounts.txt")

counts <- fc[,2:7]
colnames(counts) <- metadata$SampleName1
gene.id <- strsplit2(fc[,1], ":")[,1]
ete.id <- fc[,1]

library(edgeR)
dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))


######### filtering
name1 <- "SimDroV2_fc"
dge <- dgeOrg
keep <- rowSums(cpm(dge)>0) >= 4
dge <- dge[keep,]


length(unique(dge$genes$gene.id))



#######################################################
# create dge object from HTSEQ counts from Drosophila
#######################################################

library(DEXSeq)
library("edgeR")

ecs <- DEXSeqDataSetFromHTSeq(countfiles = paste0("PLOTS3/htseqCounts/htseq_samp_", metadata$SampleName1, ".txt"), sampleData= data.frame(row.names = metadata$SampleName, condition = metadata$condition), design = ~sample + exon + condition:exon)

counts <- featureCounts(ecs)

colnames(counts) <- metadata$SampleName1
gene.id <- strsplit2(rownames(counts), ":")[,1]
ete.id <- rownames(counts)

dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))


#########  filtering
name1 <- "SimDroV2_htseq"
dge <- dgeOrg
keep <- rowSums(cpm(dge)>0) >= 4
dge <- dge[keep,]

length(unique(dge$genes$gene.id))



#######################################################
# run DM pipeline
#######################################################
## constrOptim does not work!!!

library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version3/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions.R"))

out.dir <- "PLOTS3/2groupSimDroV2/"
dir.create(out.dir, showWarnings=F, recursive=T)


############### run 

length(unique(as.character(dge$genes$gene.id)))
length(unique(dge$genes$gene.id))


dgeDM <- dmEstimateCommonDisp(dge, group=dge$samples$group, subset=1000, adjust = TRUE, mode = "constrOptim2", mcCores=30, interval = c(0, 1e+5), tol = 1e-05, verbose=FALSE)

dgeDM$commonDispersion
log(dgeDM$commonDispersion)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", mcCores=30, verbose=FALSE)

dgeDM <- dmTest(dgeDM, mode="constrOptim2", mcCores=30, verbose=FALSE)
  
head(dgeDM$table)
  

save(dgeDM, file = paste0(out.dir, "/", name1, "_dgeDM",".Rdata"))
write.table(dgeDM$table, paste0(out.dir, "/", name1, "_table",".xls"), sep = "\t", quote = F, row.names = F)




#######################################################
# plot profile likelihood for common dispersion
#######################################################
library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version3/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions.R"))

out.dir <- "PLOTS3/2groupSimDroV2/"
dir.create(out.dir, showWarnings=F, recursive=T)


################### run

adjust = FALSE
mode = "constrOptim2"
mcCores = 30
verbose = FALSE


gX <- seq(0, 1e+8, by = 100)
gX <- gX[-1]
ploglikY <- rep(0, length(gX))

group <- dge$samples$group
colnames(dge$counts) <- dge$samples$group
rownames(dge$counts) <- dge$genes$ete.id

group <- as.factor(group)
ngroups <- nlevels(group)

## only for 1 group
# for(gr in 1){
  gr=1
  grIndx <- which(group==(levels(group)[gr]))
  
  y <- split(as.data.frame(dge$counts[, grIndx, drop=FALSE]), as.character(dge$genes$gene.id))

for(i in 1:length(ploglikY))
  ploglikY[i] <- dmAdjustedProfileLik(gamma0=gX[i], y, adjust=adjust, mode=mode, mcCores=mcCores, common = TRUE, verbose=verbose) 

write.table(data.frame(gX=gX, ploglikY=ploglikY), paste0(out.dir, "/", name1,"_gr",gr, "_profileLogLik",".txt"), quote = F, sep = "\t", row.names = F)

pdf(paste0(out.dir, "/", name1,"_gr",gr, "_profileLogLik",".pdf"))
plot(gX, ploglikY, type="l", col="deeppink", lwd=4, xlab="gamma +", ylab="Profile loglikelihood")
dev.off()



## adjusted for htseq counts
pdf(paste0(out.dir, "/", name1,"_gr",gr, "_profileLogLik",".pdf"))
plot(gX, ploglikY, type="l", col="deeppink", lwd=4, xlab="gamma +", ylab="Profile loglikelihood")
plot(gX, -log(-ploglikY), type="l", col="deeppink", lwd=4, xlab="gamma +", ylab="- log -Profile loglikelihood")
plot(gX[ploglikY != 0], ploglikY[ploglikY != 0], type="l", col="deeppink", lwd=4, xlab="gamma +", ylab="Profile loglikelihood")
plot(gX[gX > 1e+5 & ploglikY != 0], ploglikY[gX > 1e+5 & ploglikY != 0], type="l", col="deeppink", lwd=4, xlab="gamma +", ylab="Profile loglikelihood")
dev.off()

# }















