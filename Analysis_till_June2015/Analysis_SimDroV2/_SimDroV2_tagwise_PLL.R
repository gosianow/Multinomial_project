

# BioC 2.14

# Created 28 Oct 2014
# Modyfied 28 Oct 2014



#######################################################
# create dge object from FC counts from Drosophila
#######################################################

setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2")


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


library(limma)

fc <- read.table("/home/gosia/Multinomial_project/DM_package/PLOTS4/featureCounts/featureCounts.txt")

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

dge <- dge[1:1000, ]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id))

dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR




#######################################################
# run DM pipeline for aternative splicing
#######################################################


library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")

source(paste0(Rdir, "dmFunctions_v5.R"))



out.dir <- "PLOTS_DM_v5_TagwisePLL/"
dir.create(out.dir, showWarnings=F, recursive=T)


#######################################################

group=NULL; adjust = TRUE; mode = "constrOptim2G"; epsilon = 1e-05; maxIte = 1000; interval = c(0, 1e+5); tol = 1e-00; mcCores=20; verbose=FALSE; modeDisp=c("optimize", "optim", "constrOptim")[2]; initDisp = 100


y <- dge$counts
genes <- names(y)

if(is.null(group)) group <- dge$samples$group
group <- as.factor(group)
ngroups <- nlevels(group)
lgroups <- levels(group)

igroups <- list()
for(gr in 1:ngroups){
  # gr=2
  igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
  
}

g <- 2
y[[g]]
parGamma0 <- seq(100, 800, by = 10)


g <- 1
y[[g]]
# parGamma0 <- seq(100, 100000, by = 100)
# parGamma0 <- seq(100, 10000, by = 10)
parGamma0 <- seq(100, 2000, by = 10)


g <- 3
y[[g]]
parGamma0 <- seq(50, 500, by = 10)


g <- 4
y[[g]]
parGamma0 <- seq(100, 100000, by = 100)


g <- 6
y[[g]]
parGamma0 <- seq(100, 100000, by = 100)




tgPLL <- rep(0, length(parGamma0))

for(i in 1:length(tgPLL)){
  
  tgPLL[i] <- dmAdjustedProfileLikTG(gamma0 = parGamma0[i], y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)
  
  
}



pdf(paste0(out.dir,"tgPLL_g", g, ".pdf"))

plot(parGamma0, tgPLL, type = "l" , col="deeppink", lwd=4, xlab="gamma +", ylab="Profile loglikelihood")

plot(parGamma0[-1], diff(tgPLL) > 0, col = "deeppink2", lwd=4, xlab="gamma +", ylab="Monotonicity of profile loglikelihood")

dev.off()


# ploted genes:
# - optimize()
# FBgn0000014 FBgn0000015 FBgn0000017 FBgn0000018 FBgn0000022 FBgn0000024
# 99987.5632    401.5831    118.1388  99999.3633          NA  99999.5388 

# - optim()
# FBgn0000014 FBgn0000015 FBgn0000017 FBgn0000018 FBgn0000022 FBgn0000024
# 167.6487    401.5470    100.0003  33028.3749          NA  30592.7970

# - constrOptim()
# FBgn0000014 FBgn0000015 FBgn0000017 FBgn0000018 FBgn0000022 FBgn0000024
# 1212.3080    403.1098    117.5000   6012.1467          NA  13171.6480







