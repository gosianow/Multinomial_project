


# BioC 14
# Created 03 Sep 2014: test dmEstimateCommonDisp in dmFunctions_v4.R
# Analysis of drosophila simulation data
# Updated 03 Sep 2014:

# + pooled dispersion 



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

setwd("/home/gosia/Multinomial_project/DM_package/")


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


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




######### check the difference in commonDispersion for different subset size
name1 <- "htseq_g0_s3_keep0s_subsetInf_DM"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 3
dge <- dge[keep,]



library("edgeR")
dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))


######### standard filtering
name1 <- "SimDroV2"
dge <- dgeOrg
keep <- rowSums(cpm(dge)>0) >= 3
dge <- dge[keep,]
# dge$counts[ dge$counts == 0 ] <- 1

length(unique(dge$genes$gene.id))

out.dir <- "PLOTS3/2groupSimDroV2/"
dir.create(out.dir, showWarnings=F, recursive=T)



#######################################################
# run DM pipeline
#######################################################
## constrOptim does not work!!!

library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version4/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v4.R"))

out.dir <- "PLOTS4/2groupSimDroV2/"
dir.create(out.dir, showWarnings=F, recursive=T)


############### run 

length(unique(as.character(dge$genes$gene.id)))
length(unique(dge$genes$gene.id))


# dge <- dge[1:100, ]
# dge$genes$gene.id <- as.factor(as.character(dge$genes$gene.id))


# dgeDM <- dmFit(dge, group=NULL, dispersion=3000, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)


dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-05, mcCores=20, verbose=FALSE)


dgeDMadj <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-05, mcCores=1, verbose=FALSE)

dgeDM <- dgeDMadj

dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)


dgeDM <- dmTest(dgeDM, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

  
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


# gX <- seq(0, 1e+3, by = 10)
gX <- seq(1e+3, 1e+5, by = 20)
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

pdf(paste0(out.dir, "/", name1,"_gr",gr, "_profileLogLik",".pdf"))
plot(gX, ploglikY, type="l", col="deeppink", lwd=4, xlab="gamma +", ylab="Profile loglikelihood")
dev.off()

# }


gX1 <- gX
ploglikY1 <- ploglikY





#######################################################
# test dmAdj for Inf values
#######################################################



dge = dgeFit

ynames <- names(y)

adjtmp <- mclapply(seq(length(y)), function(g){  
  
  a <- adjCROneGeneManyGroups(y = y[[g]], group, gamma0 = gamma0, piH = dge$fit[[g]]$piH)  
  
  if(!is.null(a))
    names(a) <- ynames[g]
  
  return(a)
  
}, mc.cores=mcCores)

at <- unlist(adjtmp)

atInf <- at[at == Inf]

gInf <- names(atInf)

yOrg <- y

y <- yOrg[["FBgn0000054"]]

yOrg[gInf]


g <- "FBgn0000054"
y = yOrg[[g]]; gamma0 = gamma0; piH = dge$fit[[g]]$piH


# SOLUTION:
# a <- adjCROneGeneGroup(y = y[, grIndx, drop = FALSE], gamma0, piH = piH[, lgroups[gr]])


#######################################################
# 
#######################################################























