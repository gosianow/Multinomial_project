##############################################################################################################
# created 17 Nov 2014
# BioC 3.0
# test the parallel solutions on DM_0.1.1

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



#dge <- dge[1:20000, ]

dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id))

dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR



##############################################################################################################
# run DM pipeline for aternative splicing
##############################################################################################################

out.dir <- "DM_output/"
dir.create(out.dir, showWarnings=F, recursive=T)

# save.image(paste0(out.dir, "testWorkspace.RData"))



#######################################################
# test DM_0.1.1
#######################################################

### If all above is fine, load this workspace 
setwd("/home/gosia/Multinomial_project/DM_package_devel/")
out.dir <- "DM_output/"
load(paste0(out.dir, "testWorkspace.RData"))

### load DM
library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")


### Source all R files in DM package
Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in 1:length(Rfiles)) source(Rfiles[i])



######################## Test modeDip = grid 


# dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, subset=Inf, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=20, verbose=FALSE)

dgeDM <- dge
dgeDM$commonDispersion <- 1365.728



gridLength = 11
gridRange = c(-5, 5)



splinePts <- seq(from = gridRange[1], to = gridRange[2], length = gridLength)
splinePts
splineDisp <- dgeDM$commonDispersion * 2^splinePts
splineDisp


priorDf = 4 ### priorN <- priorDf/(nlibs - ngroups)
span = 0.3



time_mclapply <- system.time( dgeDM <- dmEstimateTagwiseDisp(dgeDM, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp=c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = FALSE, gridLength=gridLength, gridRange=gridRange, trend = c("none", "commonDispersion", "trendedDispersion")[2], priorDf=priorDf, span=span, mcCores=20, verbose=FALSE) )

time_mclapply


### function dmEstimateTagwiseDisp:
group=NULL; adjust = FALSE; mode = "constrOptim2G"; epsilon = 1e-05; maxIte = 1000; modeDisp=c("optimize", "optim", "constrOptim", "grid")[4]; interval = c(0, 1e+5); tol = 1e-00;  initDisp = 10; initWeirMoM = FALSE; gridLength=11; gridRange=c(-5, 5); trend = c("none", "commonDispersion", "trendedDispersion")[2]; priorDf=4; span=0.3; verbose=FALSE

mcCores=10


dge <- dgeDM

y <- dge$counts
genes <- names(y)
ngenes <- length(y)

if(is.null(group)) group <- dge$samples$group
group <- as.factor(group)
ngroups <- nlevels(group)
lgroups <- levels(group)
nlibs <- length(group)

igroups <- list()
for(gr in 1:ngroups){
  # gr=2
  igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
  
}






mcCores=5

### calculate mean expression of genes 
time <- system.time(meanExpr1 <- unlist(mclapply(seq(ngenes), function(g){
  # message("working")
  sum(y[[g]]) / nlibs 
  },  mc.cores=mcCores)) )
time 
head(meanExpr1)






time <- system.time(meanExpr <- bpvec(seq(ngenes), function(g){ 
  
  lg <- length(g)
  mexp <- numeric(lg)
  for(i in seq(lg))
    mexp[i] <- sum(y[[g[i]]]) / nlibs 
  
  return(mexp)
  
  }, AGGREGATE=c, BPPARAM = MulticoreParam(mcCores)) )

time
head(meanExpr)



all.equal(meanExpr1, meanExpr)




names(meanExpr) <- genes
dge$meanExpr <- meanExpr



### genrate spline dispersion
splinePts <- seq(from = gridRange[1], to = gridRange[2], length = gridLength)
splineDisp <- dge$commonDispersion * 2^splinePts




### calculate the likelihood for each gene at the spline dispersion points

mcCores=20

time <- system.time( loglik0L <- mclapply(seq(length(y)), function(g){
  # g = 1
  
  ll <- numeric(gridLength)
  
  for(i in seq(gridLength)){
    # i = 1 
    out <- dmAdjustedProfileLikTG(gamma0 = splineDisp[i], y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)
    if(is.null(out))
      return(NULL)
    
    ll[i] <- out
    
  }
  
  return(ll)
  
}, mc.cores = mcCores) )

time


############################### BiocParallel - bplapply

library(BiocParallel)

mcCores <- 20

time <- system.time( loglik0L <- bplapply(seq(length(y)), function(g){
  # g = 1
  
  ll <- numeric(gridLength)
  
  for(i in seq(gridLength)){
    # i = 1 
    out <- dmAdjustedProfileLikTG(gamma0 = splineDisp[i], y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)
    if(is.null(out))
      return(NULL)
    
    ll[i] <- out
    
  }
  
  return(ll)
  
}, BPPARAM = MulticoreParam(mcCores) ) )

time

names(loglik0L) <- genes

system.time( loglik0 <- do.call(rbind, loglik0L) )

loglik0.1 <- loglik0


############################### BiocParallel 

library(BiocParallel)

mcCores <- 20

time <- system.time( loglik0L <- bpvec(seq(length(y)), function(g){
  # g = 20:40

  lg <- length(g)
  ll <- matrix(0, nrow = lg, ncol = gridLength)
  
  for(j in seq(lg)){
    # j = 1
    
    for(i in seq(gridLength)){
      # i = 1 
      
      out <- dmAdjustedProfileLikTG(gamma0 = splineDisp[i], y = y[[g[j]]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)
      
      if(is.null(out))
        ll[j, i] <- NA
      else 
        ll[j, i] <- out

      
    }
 
  }
  
  return(ll)
  
}, AGGREGATE=rbind, BPPARAM = MulticoreParam(mcCores) ) )

time



### vectorizing
mcCores = 2

bpvec(1:10, function(v) {
  message("working") ## 10 tasks, 4 messages 
  sqrt(v) ## <- this function must be vectorized 
}, AGGREGATE=c, BPPARAM = MulticoreParam(mcCores))









































